#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Parallel FIMO motif scanning for RBP motif databases.
#
# Main analysis:
#   CisBP-RNA human direct + inferred motifs
#
# Supplementary analysis:
#   ATtRACT human non-mutated motifs
#
# FASTA input:
#   FIMO-friendly nameOnly FASTA headers, e.g.
#   >ENSG00000119661|DNAL1|three_prime_UTR
#
# Scheduling:
#   All selected jobs, including ALL_BACKGROUND and feature jobs, are run in
#   parallel using MAX_JOBS. This is suitable after restricting the background
#   to the EMD-tested gene universe.
#
# Default mode:
#   RERUN_FAILED_ONLY=1
#   The script identifies failed/missing/incomplete jobs and reruns only those.
#   A successful job must have:
#     - fimo_status.tsv with status ok or ok_existing
#     - non-empty fimo.tsv
#
# Usage:
#   bash 11_run_FIMO_for_RBP_motif_databases_parallel.sh
#
# Optional:
#   MAX_JOBS=6 bash 11_run_FIMO_for_RBP_motif_databases_parallel.sh
#   RERUN_FAILED_ONLY=0 bash 11_run_FIMO_for_RBP_motif_databases_parallel.sh
# ============================================================

tool_dir="/opt/software/miniforge/envs/motif/bin"
fimo="${tool_dir}/fimo"

fasta_dir="./04_GENCODE_v43_UTR_CDS_FASTA_EMDtested_Hep3B_background_FIMO_nameOnly"
motif_root="./05_offline_motifs"
outdir="./06_FIMO_RBP_motif_scan_parallel_nameOnly_EMDtested_Hep3B_background"

cisbp_meme="${motif_root}/CisBP_RNA/CisBP_RNA_human_direct_and_inferred.meme"
attract_meme="${motif_root}/ATtRACT/ATtRACT_human_nonmutated_all_RBP.meme"

fimo_p_thresh="1e-4"
max_jobs="${MAX_JOBS:-6}"
rerun_failed_only="${RERUN_FAILED_ONLY:-1}"

mkdir -p "${outdir}"

if [[ ! -x "${fimo}" ]]; then
  echo "ERROR: fimo not found or not executable: ${fimo}" >&2
  exit 1
fi

if [[ ! -d "${fasta_dir}" ]]; then
  echo "ERROR: FASTA directory not found: ${fasta_dir}" >&2
  exit 1
fi

for motif_file in "${cisbp_meme}" "${attract_meme}"; do
  if [[ ! -s "${motif_file}" ]]; then
    echo "ERROR: motif MEME file not found or empty: ${motif_file}" >&2
    exit 1
  fi
done

if ! [[ "${max_jobs}" =~ ^[0-9]+$ ]] || [[ "${max_jobs}" -lt 1 ]]; then
  echo "ERROR: MAX_JOBS must be a positive integer. Current value: ${max_jobs}" >&2
  exit 1
fi

if ! [[ "${rerun_failed_only}" =~ ^[01]$ ]]; then
  echo "ERROR: RERUN_FAILED_ONLY must be 0 or 1. Current value: ${rerun_failed_only}" >&2
  exit 1
fi

manifest="${outdir}/fimo_scan_manifest.tsv"
joblist="${outdir}/fimo_joblist.tsv"
run_joblist="${outdir}/fimo_joblist_to_run.tsv"
logdir="${outdir}/logs"
mkdir -p "${logdir}"

sanitize_name() {
  echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g' | sed 's/_\+/_/g' | sed 's/^_\|_$//g'
}

add_job() {
  local database="$1"
  local analysis_role="$2"
  local motif_file="$3"
  local region_type="$4"
  local set_type="$5"
  local feature_group="$6"
  local fasta_file="$7"

  local feature_safe
  local fimo_dir
  local fimo_tsv

  feature_safe="$(sanitize_name "${feature_group}")"
  fimo_dir="${outdir}/${database}/${region_type}/${set_type}/${feature_safe}"
  fimo_tsv="${fimo_dir}/fimo.tsv"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
    "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" >> "${joblist}"
}

build_jobs_for_database() {
  local database="$1"
  local analysis_role="$2"
  local motif_file="$3"
  local region_type
  local background_fa
  local feature_fa
  local base
  local feature_group

  for region_type in three_prime_UTR five_prime_UTR CDS; do
    background_fa="${fasta_dir}/background_EMDtested_Hep3B2.1_7_vs_rest_${region_type}.fa"
    add_job "${database}" "${analysis_role}" "${motif_file}" "${region_type}" "background" "ALL_BACKGROUND" "${background_fa}"

    shopt -s nullglob
    feature_fastas=("${fasta_dir}"/*_feature_"${region_type}".fa)
    shopt -u nullglob

    if [[ "${#feature_fastas[@]}" -eq 0 ]]; then
      echo "WARNING: no feature FASTA files found for region: ${region_type}" >&2
    fi

    for feature_fa in "${feature_fastas[@]}"; do
      base="$(basename "${feature_fa}")"
      feature_group="${base%_feature_${region_type}.fa}"
      add_job "${database}" "${analysis_role}" "${motif_file}" "${region_type}" "feature" "${feature_group}" "${feature_fa}"
    done
  done
}

job_success_status() {
  local fimo_dir="$1"
  local fimo_tsv="$2"
  local status_file="${fimo_dir}/fimo_status.tsv"
  local status

  if [[ ! -s "${fimo_tsv}" || ! -s "${status_file}" ]]; then
    return 1
  fi

  status="$(awk -F '\t' 'NR == 1 {print $10}' "${status_file}")"
  [[ "${status}" == "ok" || "${status}" == "ok_existing" ]]
}

build_run_joblist() {
  printf "database\tanalysis_role\tregion_type\tset_type\tfeature_group\tfasta_file\tmotif_file\tfimo_dir\tfimo_tsv\n" > "${run_joblist}"

  while IFS=$'\t' read -r database analysis_role region_type set_type feature_group fasta_file motif_file fimo_dir fimo_tsv; do
    if [[ "${rerun_failed_only}" -eq 0 ]]; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
        "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" >> "${run_joblist}"
      continue
    fi

    if ! job_success_status "${fimo_dir}" "${fimo_tsv}"; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
        "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" >> "${run_joblist}"
    fi
  done < <(tail -n +2 "${joblist}")
}

run_one_job() {
  local database="$1"
  local analysis_role="$2"
  local region_type="$3"
  local set_type="$4"
  local feature_group="$5"
  local fasta_file="$6"
  local motif_file="$7"
  local fimo_dir="$8"
  local fimo_tsv="$9"
  local force_job_rerun="${10:-0}"

  local status_file="${fimo_dir}/fimo_status.tsv"
  local log_file="${logdir}/$(sanitize_name "${database}_${region_type}_${set_type}_${feature_group}").log"

  mkdir -p "${fimo_dir}"

  if [[ ! -s "${fasta_file}" ]]; then
    echo "Skip empty FASTA: ${fasta_file}" > "${log_file}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tskipped_empty_fasta\n" \
      "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
      "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" > "${status_file}"
    return 0
  fi

  if [[ "${force_job_rerun}" -eq 0 && -s "${fimo_tsv}" ]]; then
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tok_existing\n" \
      "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
      "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" > "${status_file}"
    return 0
  fi

  echo "Running FIMO: database=${database}, region=${region_type}, set=${set_type}, group=${feature_group}" > "${log_file}"

  if "${fimo}" \
    --oc "${fimo_dir}" \
    --verbosity 1 \
    --thresh "${fimo_p_thresh}" \
    --norc \
    --no-qvalue \
    "${motif_file}" \
    "${fasta_file}" >> "${log_file}" 2>&1; then

    if [[ -s "${fimo_tsv}" ]]; then
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tok\n" \
        "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
        "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" > "${status_file}"
    else
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmissing_fimo_tsv\n" \
        "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
        "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" > "${status_file}"
    fi
  else
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tfailed\n" \
      "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
      "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" > "${status_file}"
  fi
}

regenerate_manifest() {
  printf "database\tanalysis_role\tregion_type\tset_type\tfeature_group\tfasta_file\tmotif_file\tfimo_dir\tfimo_tsv\tstatus\n" > "${manifest}"

  while IFS=$'\t' read -r database analysis_role region_type set_type feature_group fasta_file motif_file fimo_dir fimo_tsv; do
    status_file="${fimo_dir}/fimo_status.tsv"
    if [[ -s "${status_file}" ]]; then
      cat "${status_file}" >> "${manifest}"
    else
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tnot_run\n" \
        "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" \
        "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" >> "${manifest}"
    fi
  done < <(tail -n +2 "${joblist}")
}

printf "database\tanalysis_role\tregion_type\tset_type\tfeature_group\tfasta_file\tmotif_file\tfimo_dir\tfimo_tsv\n" > "${joblist}"
build_jobs_for_database "CisBP_RNA" "main" "${cisbp_meme}"
build_jobs_for_database "ATtRACT" "supplementary" "${attract_meme}"
build_run_joblist

total_jobs="$(( $(wc -l < "${joblist}") - 1 ))"
jobs_to_run="$(( $(wc -l < "${run_joblist}") - 1 ))"
background_jobs_to_run="$(
  awk -F '\t' 'NR > 1 && $4 == "background" && $5 == "ALL_BACKGROUND" {n++} END {print n + 0}' "${run_joblist}"
)"
feature_jobs_to_run="$(
  awk -F '\t' 'NR > 1 && !($4 == "background" && $5 == "ALL_BACKGROUND") {n++} END {print n + 0}' "${run_joblist}"
)"

echo "Total FIMO jobs in design: ${total_jobs}"
echo "RERUN_FAILED_ONLY: ${rerun_failed_only}"
echo "Jobs selected to run: ${jobs_to_run}"
echo "ALL_BACKGROUND jobs selected: ${background_jobs_to_run}"
echo "Feature/non-background jobs selected: ${feature_jobs_to_run}"
echo "MAX_JOBS: ${max_jobs}"
echo "Output directory: ${outdir}"
echo "Selected job list: ${run_joblist}"

if [[ "${jobs_to_run}" -eq 0 ]]; then
  echo "No failed/missing/incomplete jobs detected. Regenerating manifest only."
  regenerate_manifest
else
  force_job_rerun="0"
  if [[ "${rerun_failed_only}" -eq 1 ]]; then
    force_job_rerun="1"
  fi

  echo "Running all selected jobs in parallel..."
  while IFS=$'\t' read -r database analysis_role region_type set_type feature_group fasta_file motif_file fimo_dir fimo_tsv; do
    run_one_job "${database}" "${analysis_role}" "${region_type}" "${set_type}" "${feature_group}" "${fasta_file}" "${motif_file}" "${fimo_dir}" "${fimo_tsv}" "${force_job_rerun}" &

    while [[ "$(jobs -rp | wc -l | tr -d ' ')" -ge "${max_jobs}" ]]; do
      sleep 2
    done
  done < <(tail -n +2 "${run_joblist}")

  wait || true
  regenerate_manifest
fi
failed_n="$(awk -F '\t' 'NR > 1 && $10 == "failed" {n++} END {print n + 0}' "${manifest}")"
missing_n="$(awk -F '\t' 'NR > 1 && $10 == "missing_fimo_tsv" {n++} END {print n + 0}' "${manifest}")"
not_run_n="$(awk -F '\t' 'NR > 1 && $10 == "not_run" {n++} END {print n + 0}' "${manifest}")"
ok_n="$(awk -F '\t' 'NR > 1 && ($10 == "ok" || $10 == "ok_existing") {n++} END {print n + 0}' "${manifest}")"

echo "Done."
echo "Manifest: ${manifest}"
echo "OK jobs: ${ok_n}"
echo "Failed jobs: ${failed_n}"
echo "Missing fimo.tsv jobs: ${missing_n}"
echo "Not-run jobs: ${not_run_n}"

if [[ "${failed_n}" -gt 0 || "${missing_n}" -gt 0 || "${not_run_n}" -gt 0 ]]; then
  echo "WARNING: some FIMO jobs are still failed/missing/not_run. Check selected jobs and logs:" >&2
  echo "  ${run_joblist}" >&2
  echo "  ${logdir}" >&2
fi