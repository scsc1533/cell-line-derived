#!/usr/bin/env python3
"""Stable length-gene selection for the exact published-SVC stacking model.

This wrapper adds:
  1. optional reuse of a precomputed all-sample D_target table;
  2. Top-N direction-consistent length genes per OOF training fold;
  3. final stability selection by frequency across all OOF folds.

It reuses run_stacking_exact_SVC.py for the exact expression SVC and the rest
of the stacking workflow. Keep this file, run_stacking_exact_SVC.py,
stacking_expression_length.py and EML.py in the same directory.
"""

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import stacking_expression_length as core
import run_stacking_exact_SVC as exact


def preparse_stability_args():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--precomputed-d-target", default=None)
    p.add_argument("--length-top-n-htr8", type=int, default=20)
    p.add_argument("--length-top-n-k562", type=int, default=20)
    p.add_argument("--length-top-n-liver", type=int, default=5)
    p.add_argument("--length-stability-threshold", type=float, default=0.60)
    p.add_argument(
        "--length-ranking-direction",
        choices=["negative", "any"],
        default="negative",
        help="negative requires median D_target(GDM)-median D_target(Healthy)<0",
    )
    # Consume these here so the expected number of selection calls is known,
    # then add them back for the original parser used by exact.main().
    p.add_argument("--oof-folds", type=int, default=5)
    p.add_argument("--oof-repeats", type=int, default=5)
    known, remaining = p.parse_known_args()

    if min(
        known.length_top_n_htr8,
        known.length_top_n_k562,
        known.length_top_n_liver,
    ) < 1:
        p.error("All per-pool --length-top-n-* values must be >=1")
    if not (0 < known.length_stability_threshold <= 1):
        p.error("--length-stability-threshold must be in (0,1]")
    if known.oof_folds < 2 or known.oof_repeats < 1:
        p.error("OOF folds must be >=2 and repeats must be >=1")

    remaining.extend(["--oof-folds", str(known.oof_folds)])
    remaining.extend(["--oof-repeats", str(known.oof_repeats)])

    # The original parser marks raw length inputs as required. When a complete
    # precomputed table is supplied, inject harmless placeholders; the raw-data
    # functions are replaced below and these paths are never opened.
    if known.precomputed_d_target:
        required_placeholders = {
            "--gene-score": "PRECOMPUTED_D_TARGET",
            "--reference": "PRECOMPUTED_D_TARGET",
            "--length": "PRECOMPUTED_D_TARGET",
        }
        for flag, value in required_placeholders.items():
            if flag not in remaining:
                remaining.extend([flag, value])

    sys.argv = [sys.argv[0]] + remaining
    return known


STABILITY_ARGS = preparse_stability_args()
EXPECTED_OOF_SELECTIONS = (
    STABILITY_ARGS.oof_folds * STABILITY_ARGS.oof_repeats
)
REQUIRED_SELECTION_COUNT = int(math.ceil(
    EXPECTED_OOF_SELECTIONS * STABILITY_ARGS.length_stability_threshold
))

TOP_N_BY_COMPARISON = {
    "HTR_8_SVneo_vs_rest": STABILITY_ARGS.length_top_n_htr8,
    "K562_vs_rest": STABILITY_ARGS.length_top_n_k562,
    "HepG2_Hep3B2.1_7_vs_rest": STABILITY_ARGS.length_top_n_liver,
}


def load_precomputed_distance(path):
    x = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    required = {"comparison", "sample", "Gene", "D_target"}
    missing = required.difference(x.columns)
    if missing:
        raise ValueError(
            f"Precomputed D_target file is missing columns: {sorted(missing)}"
        )
    x = x[x["comparison"].isin(core.COMPARISON_FEATURE_STEM)].copy()
    x["sample"] = x["sample"].astype(str)
    x["Gene"] = x["Gene"].astype(str)
    x["D_target"] = pd.to_numeric(x["D_target"], errors="coerce")
    if "plasma_total_count" in x.columns:
        x["plasma_total_count"] = pd.to_numeric(
            x["plasma_total_count"], errors="coerce"
        )
    if x.empty:
        raise ValueError("No rows from the three candidate pools in D_target file")
    return x


PRECOMPUTED_DISTANCE = None
if STABILITY_ARGS.precomputed_d_target:
    PRECOMPUTED_DISTANCE = load_precomputed_distance(
        STABILITY_ARGS.precomputed_d_target
    )
    print(
        "Reusing precomputed D_target:",
        STABILITY_ARGS.precomputed_d_target,
        f"({len(PRECOMPUTED_DISTANCE)} rows)",
    )


# If requested, replace the expensive raw length calculation. exact.main()
# still writes a copy into its new output directory for reproducibility, but it
# does not reread or recalculate the raw length distributions.
if PRECOMPUTED_DISTANCE is not None:
    def prepare_candidate_pools_from_precomputed(_file, _fdr):
        return (
            PRECOMPUTED_DISTANCE[["comparison", "Gene"]]
            .drop_duplicates()
            .reset_index(drop=True)
        )

    def skip_raw_length_read(*_args, **_kwargs):
        return pd.DataFrame(), []

    def return_precomputed_distance(*_args, **_kwargs):
        return PRECOMPUTED_DISTANCE.copy()

    core.prepare_candidate_pools = prepare_candidate_pools_from_precomputed
    core.read_length_counts = skip_raw_length_read
    exact.compute_all_d_target = return_precomputed_distance


SELECTION_STATE = {
    "call_index": 0,
    "history": [],
    "stable_all": None,
    "stable_selected": None,
}


def calculate_gene_statistics(distance_dt, train_ids, y_by_sample, args):
    ids = set(map(str, train_ids))
    z = distance_dt[distance_dt["sample"].astype(str).isin(ids)].copy()
    if "plasma_total_count" in z.columns:
        # Supports raising the count threshold when reusing a table generated
        # with a lower threshold. A lower new threshold cannot restore rows
        # already removed from the precomputed table.
        z = z[z["plasma_total_count"] >= args.min_plasma_gene_count]
    z["y"] = z["sample"].astype(str).map(y_by_sample)

    rows = []
    for (comparison, gene), g in z.groupby(
        ["comparison", "Gene"], sort=False
    ):
        healthy = g.loc[g["y"] == 0, "D_target"].dropna().to_numpy()
        gdm = g.loc[g["y"] == 1, "D_target"].dropna().to_numpy()
        if (
            len(healthy) >= args.length_min_group_n
            and len(gdm) >= args.length_min_group_n
        ):
            try:
                p_value = mannwhitneyu(
                    healthy, gdm, alternative="two-sided"
                ).pvalue
            except ValueError:
                p_value = np.nan
        else:
            p_value = np.nan
        healthy_median = np.median(healthy) if len(healthy) else np.nan
        gdm_median = np.median(gdm) if len(gdm) else np.nan
        rows.append({
            "comparison": comparison,
            "Gene": gene,
            "n_healthy": len(healthy),
            "n_gdm": len(gdm),
            "healthy_median_D_target": healthy_median,
            "gdm_median_D_target": gdm_median,
            "gdm_minus_healthy_D_target": gdm_median - healthy_median,
            "wilcox_p_D_target": p_value,
        })
    return pd.DataFrame(rows)


def select_top_n_in_one_fold(stats, args):
    selected = []
    for comparison in core.COMPARISON_FEATURE_STEM:
        s = stats[
            (stats["comparison"] == comparison)
            & stats["wilcox_p_D_target"].notna()
            & (stats["n_healthy"] >= args.length_min_group_n)
            & (stats["n_gdm"] >= args.length_min_group_n)
        ].copy()
        if STABILITY_ARGS.length_ranking_direction == "negative":
            s = s[s["gdm_minus_healthy_D_target"] < 0].copy()
        s["absolute_median_effect"] = s[
            "gdm_minus_healthy_D_target"
        ].abs()
        s = s.sort_values(
            ["wilcox_p_D_target", "absolute_median_effect", "Gene"],
            ascending=[True, False, True],
        )
        top_n = TOP_N_BY_COMPARISON[comparison]
        keep = s.head(top_n).copy()
        if keep.empty:
            raise ValueError(
                f"No eligible length genes in OOF fold for {comparison}. "
                "Check count coverage, --length-min-group-n or direction rule."
            )
        keep["selection_rule"] = "OOF_direction_ranked_top_n"
        keep["configured_top_n_for_pool"] = top_n
        keep["rank_within_pool"] = np.arange(1, len(keep) + 1)
        selected.append(keep)
    return pd.concat(selected, ignore_index=True)


def build_stable_final_selection(full_training_stats):
    history = pd.concat(SELECTION_STATE["history"], ignore_index=True)
    counts = (
        history[["comparison", "Gene", "selection_fold_id"]]
        .drop_duplicates()
        .groupby(["comparison", "Gene"], as_index=False)
        .agg(selection_count=("selection_fold_id", "nunique"))
    )
    counts["total_training_folds"] = EXPECTED_OOF_SELECTIONS
    counts["selection_frequency"] = (
        counts["selection_count"] / counts["total_training_folds"]
    )
    counts["required_selection_count"] = REQUIRED_SELECTION_COUNT
    counts["stability_threshold"] = (
        STABILITY_ARGS.length_stability_threshold
    )
    counts["configured_top_n_for_pool"] = counts["comparison"].map(
        TOP_N_BY_COMPARISON
    )

    all_candidates = full_training_stats.merge(
        counts,
        on=["comparison", "Gene"],
        how="left",
    )
    all_candidates["selection_count"] = (
        all_candidates["selection_count"].fillna(0).astype(int)
    )
    all_candidates["total_training_folds"] = EXPECTED_OOF_SELECTIONS
    all_candidates["selection_frequency"] = (
        all_candidates["selection_count"] / EXPECTED_OOF_SELECTIONS
    )
    all_candidates["required_selection_count"] = REQUIRED_SELECTION_COUNT
    all_candidates["stability_threshold"] = (
        STABILITY_ARGS.length_stability_threshold
    )
    all_candidates["configured_top_n_for_pool"] = all_candidates[
        "comparison"
    ].map(TOP_N_BY_COMPARISON)
    all_candidates["stable_selected"] = (
        all_candidates["selection_count"] >= REQUIRED_SELECTION_COUNT
    )

    stable = all_candidates[all_candidates["stable_selected"]].copy()
    missing_pools = [
        comparison
        for comparison in core.COMPARISON_FEATURE_STEM
        if comparison not in set(stable["comparison"])
    ]
    if missing_pools:
        max_frequency = (
            all_candidates.groupby("comparison")["selection_frequency"]
            .max()
            .to_dict()
        )
        raise ValueError(
            "No gene reaches the requested stability threshold in: "
            f"{missing_pools}. Maximum observed frequencies: {max_frequency}. "
            "Lower --length-stability-threshold using training data only."
        )
    stable["selection_rule"] = "stable_across_OOF_training_folds"
    stable = stable.sort_values(
        ["comparison", "selection_count", "wilcox_p_D_target"],
        ascending=[True, False, True],
    )
    SELECTION_STATE["stable_all"] = all_candidates
    SELECTION_STATE["stable_selected"] = stable
    return stable, all_candidates


def stable_select_length_genes(distance_dt, train_ids, y_by_sample, args):
    """State-aware selector: first N calls are OOF; next call is final stable."""
    call_index = SELECTION_STATE["call_index"]
    stats = calculate_gene_statistics(
        distance_dt, train_ids, y_by_sample, args
    )
    if call_index < EXPECTED_OOF_SELECTIONS:
        selected = select_top_n_in_one_fold(stats, args)
        selected["selection_fold_id"] = call_index + 1
        SELECTION_STATE["history"].append(selected.copy())
        SELECTION_STATE["call_index"] += 1
        return selected, stats

    if call_index == EXPECTED_OOF_SELECTIONS:
        SELECTION_STATE["call_index"] += 1
        return build_stable_final_selection(stats)

    raise RuntimeError(
        "Length selector was called more times than expected. "
        f"Expected {EXPECTED_OOF_SELECTIONS} OOF calls plus one final call."
    )


core.select_length_genes = stable_select_length_genes


def find_output_arg():
    if "--output" not in sys.argv:
        return None
    i = sys.argv.index("--output")
    return sys.argv[i + 1] if i + 1 < len(sys.argv) else None


if __name__ == "__main__":
    exact.main()
    output_dir = find_output_arg()
    if output_dir:
        config = {
            "precomputed_d_target": STABILITY_ARGS.precomputed_d_target,
            "oof_folds": STABILITY_ARGS.oof_folds,
            "oof_repeats": STABILITY_ARGS.oof_repeats,
            "total_training_folds": EXPECTED_OOF_SELECTIONS,
            "top_n_per_pool_per_fold": TOP_N_BY_COMPARISON,
            "stability_threshold": (
                STABILITY_ARGS.length_stability_threshold
            ),
            "required_selection_count": REQUIRED_SELECTION_COUNT,
            "ranking_direction": STABILITY_ARGS.length_ranking_direction,
            "ranking_order": (
                "Wilcoxon P ascending, absolute median effect descending"
            ),
        }
        with open(
            os.path.join(output_dir, "00_length_stability_config.json"),
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(config, f, ensure_ascii=False, indent=2)
        if SELECTION_STATE["stable_all"] is not None:
            SELECTION_STATE["stable_all"].to_csv(
                os.path.join(
                    output_dir,
                    "05_length_stability_frequency_all_candidates.tsv",
                ),
                sep="\t",
                index=False,
            )
        if SELECTION_STATE["stable_selected"] is not None:
            SELECTION_STATE["stable_selected"].to_csv(
                os.path.join(
                    output_dir,
                    "05_final_stable_length_genes.tsv",
                ),
                sep="\t",
                index=False,
            )
