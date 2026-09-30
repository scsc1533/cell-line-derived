#!/usr/bin/env bash
set -euo pipefail

# ======================== User-adjustable paths ========================
INPUT_DIR="./"
CLIN_FILE="${INPUT_DIR}/groupinfo2.txt"
EXPRESSION_FILE="${INPUT_DIR}/merged_matrix.txt"
EXPRESSION_GENE_FILE="${INPUT_DIR}/features.txt"

PRECOMPUTED_D_TARGET="/data/work/01_2603cell_culture/05_model/04_mix_model/GDM_model_stacking_seed58/03_all_samples_gene_level_D_target.tsv.gz"

OUTPUT_DIR="/data/work/01_2603cell_culture/05_model/04_mix_model/GDM_model_stacking_stable_seed58"

# ======================== Published expression model ========================
SEED=58
TEST_SIZE=0.4
EXPRESSION_DETECTION_FREQUENCY=0.3
EXPRESSION_SCALER="standard"

# ======================== Length stability selection ========================
# 5 folds x 5 repeats = 25 OOF training folds
OOF_FOLDS=5
OOF_REPEATS=5

# Per-fold Top N can be configured independently for each candidate pool.
LENGTH_TOP_N_HTR8=30
LENGTH_TOP_N_K562=30
LENGTH_TOP_N_LIVER=5

# 0.60 across 25 folds means selection_count >= ceiling(25*0.60) = 15.
LENGTH_STABILITY_THRESHOLD=0.6

# negative: require median D_target(GDM) < median D_target(Healthy)
# any: do not restrict the median-effect direction
LENGTH_RANKING_DIRECTION="negative"

# A sample-gene record contributes only when its total count is >=10.
MIN_PLASMA_GENE_COUNT=10

# A gene must have at least this many valid samples in each disease group
# within the current OOF training fold.
LENGTH_MIN_GROUP_N=5

# A sample needs at least this many selected genes to calculate one compressed
# HTR8/Liver/K562 feature.
MIN_GENES_PER_FEATURE=1

BOOTSTRAP_N=2000

# ======================== Run ========================
python 3_run_stacking_exact_SVC_stable.py \
  --clin "${CLIN_FILE}" \
  --expression "${EXPRESSION_FILE}" \
  --expression-gene-file "${EXPRESSION_GENE_FILE}" \
  --expression-gene-column feature \
  --output "${OUTPUT_DIR}" \
  --precomputed-d-target "${PRECOMPUTED_D_TARGET}" \
  --group-column group \
  --positive-label 1 \
  --split-label-column label \
  --eligible-label train_test \
  --freq-threshold "${EXPRESSION_DETECTION_FREQUENCY}" \
  --scaler "${EXPRESSION_SCALER}" \
  --test-size "${TEST_SIZE}" \
  --seed "${SEED}" \
  --models SVC \
  --oof-folds "${OOF_FOLDS}" \
  --oof-repeats "${OOF_REPEATS}" \
  --length-top-n-htr8 "${LENGTH_TOP_N_HTR8}" \
  --length-top-n-k562 "${LENGTH_TOP_N_K562}" \
  --length-top-n-liver "${LENGTH_TOP_N_LIVER}" \
  --length-stability-threshold "${LENGTH_STABILITY_THRESHOLD}" \
  --length-ranking-direction "${LENGTH_RANKING_DIRECTION}" \
  --min-plasma-gene-count "${MIN_PLASMA_GENE_COUNT}" \
  --length-min-group-n "${LENGTH_MIN_GROUP_N}" \
  --min-genes-per-feature "${MIN_GENES_PER_FEATURE}" \
  --bootstrap "${BOOTSTRAP_N}"
