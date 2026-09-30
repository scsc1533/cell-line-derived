suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Prepare EMD/Wasserstein-defined feature gene sets for RBP motif enrichment
#
# Input:
#   all_comparisons_gene_scores.tsv
#   GENCODE GTF for SYMBOL -> ENSEMBL mapping
#
# Candidate rule:
#   perm_fdr_wasserstein < fdr_cutoff
#
# Output columns:
#   feature_group
#   SYMBOL
#   ENSEMBL
#   gene_type
#   EMD_score
#   perm_p_wasserstein
#   perm_fdr_wasserstein
#   rank_wasserstein
# ============================================================

# -------------------------
# 1. Input and output paths
# -------------------------

emd_file <- "./01_length_distribution_one_vs_rest/all_comparisons_gene_scores.tsv"
gtf_file <- "./GRCh38_GENCODE_20231021/gtf_gff/gencode.v43.annotation.gtf"

outdir <- "./01_EMD_feature_gene_sets"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 2. Parameters
# -------------------------

fdr_cutoff <- 0.05

# -------------------------
# 3. Helper functions
# -------------------------

extract_gtf_attr <- function(attr, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  out <- str_match(attr, pattern)[, 2]
  out
}

make_feature_group <- function(comparison) {
  x <- sub("_vs_rest$", "", comparison)
  x <- gsub("^Hep3B2\\.1_7$", "Hep3B2.1-7", x)
  x <- gsub("^HTR_8_SVneo$", "HTR-8", x)
  x <- gsub("^HTR-8_SVneo$", "HTR-8", x)
  x <- gsub("^HTR-8/SVneo$", "HTR-8", x)

  if (grepl("HepG2.*Hep3B2|Hep3B2.*HepG2", x)) {
    return("liver-like")
  }
  if (grepl("^K562$", x)) {
    return("K562-like")
  }
  if (grepl("HTR", x)) {
    return("HTR-8-like")
  }

  paste0(x, "-like")
}

write_tsv <- function(x, file) {
  write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# -------------------------
# 4. Read EMD result table
# -------------------------

emd <- fread(emd_file, data.table = FALSE, check.names = FALSE)

required_emd_cols <- c(
  "comparison",
  "Gene",
  "wasserstein_specificity_score",
  "perm_p_wasserstein",
  "perm_fdr_wasserstein",
  "rank_wasserstein"
)

missing_emd_cols <- setdiff(required_emd_cols, colnames(emd))
if (length(missing_emd_cols) > 0) {
  stop("Missing columns in EMD result table: ", paste(missing_emd_cols, collapse = ", "))
}

emd$wasserstein_specificity_score <- as.numeric(emd$wasserstein_specificity_score)
emd$perm_p_wasserstein <- as.numeric(emd$perm_p_wasserstein)
emd$perm_fdr_wasserstein <- as.numeric(emd$perm_fdr_wasserstein)
emd$rank_wasserstein <- as.numeric(emd$rank_wasserstein)

# -------------------------
# 5. Read GENCODE GTF gene annotation
# -------------------------

gtf_gene <- fread(
  gtf_file,
  sep = "\t",
  header = FALSE,
  data.table = FALSE,
  quote = "",
  skip = "#",
  col.names = c(
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )
) %>%
  filter(feature == "gene")

gtf_gene$ENSEMBL <- sub("\\..*$", "", extract_gtf_attr(gtf_gene$attribute, "gene_id"))
gtf_gene$SYMBOL <- extract_gtf_attr(gtf_gene$attribute, "gene_name")
gtf_gene$gene_type <- extract_gtf_attr(gtf_gene$attribute, "gene_type")

gene_map <- gtf_gene %>%
  filter(!is.na(SYMBOL), SYMBOL != "", !is.na(ENSEMBL), ENSEMBL != "") %>%
  select(SYMBOL, ENSEMBL, gene_type) %>%
  distinct()

multi_symbol_map <- gene_map %>%
  count(SYMBOL, name = "ensembl_n") %>%
  filter(ensembl_n > 1) %>%
  arrange(desc(ensembl_n), SYMBOL)

write_tsv(
  gene_map,
  file.path(outdir, "GENCODE_v43_SYMBOL_to_ENSEMBL_gene_map.tsv")
)

write_tsv(
  multi_symbol_map,
  file.path(outdir, "GENCODE_v43_SYMBOL_with_multiple_ENSEMBL.tsv")
)

# -------------------------
# 6. Prepare simplified feature gene table
# -------------------------

feature_all <- emd %>%
  transmute(
    feature_group = vapply(comparison, make_feature_group, character(1)),
    SYMBOL = Gene,
    EMD_score = wasserstein_specificity_score,
    perm_p_wasserstein = perm_p_wasserstein,
    perm_fdr_wasserstein = perm_fdr_wasserstein,
    rank_wasserstein = rank_wasserstein
  ) %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != "",
    is.finite(EMD_score),
    is.finite(perm_fdr_wasserstein)
  )

feature_candidate <- feature_all %>%
  filter(perm_fdr_wasserstein < fdr_cutoff) %>%
  arrange(feature_group, perm_fdr_wasserstein, rank_wasserstein, desc(EMD_score))

feature_candidate_mapped <- feature_candidate %>%
  left_join(gene_map, by = "SYMBOL") %>%
  select(
    feature_group,
    SYMBOL,
    ENSEMBL,
    gene_type,
    EMD_score,
    perm_p_wasserstein,
    perm_fdr_wasserstein,
    rank_wasserstein
  ) %>%
  arrange(feature_group, perm_fdr_wasserstein, rank_wasserstein, desc(EMD_score), SYMBOL)

feature_unmapped <- feature_candidate_mapped %>%
  filter(is.na(ENSEMBL) | ENSEMBL == "") %>%
  distinct(feature_group, SYMBOL, EMD_score, perm_p_wasserstein, perm_fdr_wasserstein, rank_wasserstein)

feature_duplicate_symbol <- feature_candidate_mapped %>%
  filter(!is.na(ENSEMBL), ENSEMBL != "") %>%
  count(feature_group, SYMBOL, name = "mapped_row_n") %>%
  filter(mapped_row_n > 1) %>%
  arrange(feature_group, desc(mapped_row_n), SYMBOL)

write_tsv(
  feature_candidate_mapped,
  file.path(outdir, "EMD_feature_genes_wasserstein_FDR005_with_ENSEMBL.tsv")
)

write_tsv(
  feature_unmapped,
  file.path(outdir, "EMD_feature_genes_wasserstein_FDR005_unmapped_SYMBOL.tsv")
)

write_tsv(
  feature_duplicate_symbol,
  file.path(outdir, "EMD_feature_genes_wasserstein_FDR005_SYMBOL_multiple_ENSEMBL.tsv")
)

feature_group_count <- feature_candidate_mapped %>%
  distinct(feature_group, SYMBOL) %>%
  count(feature_group, name = "feature_gene_n") %>%
  arrange(feature_group)

mapped_count <- feature_candidate_mapped %>%
  mutate(mapped = !is.na(ENSEMBL) & ENSEMBL != "") %>%
  distinct(feature_group, SYMBOL, mapped) %>%
  count(feature_group, mapped, name = "gene_n") %>%
  arrange(feature_group, desc(mapped))

write_tsv(
  feature_group_count,
  file.path(outdir, "EMD_feature_gene_count_by_group.tsv")
)

write_tsv(
  mapped_count,
  file.path(outdir, "EMD_feature_gene_mapping_count_by_group.tsv")
)

summary_df <- data.frame(
  item = c(
    "Input EMD rows",
    "Candidate FDR cutoff",
    "Candidate SYMBOL x feature_group rows",
    "Mapped output rows",
    "Unmapped SYMBOL rows",
    "SYMBOLs with multiple ENSEMBL in GTF",
    "Feature groups"
  ),
  value = c(
    nrow(emd),
    fdr_cutoff,
    nrow(feature_candidate),
    nrow(feature_candidate_mapped),
    nrow(feature_unmapped),
    nrow(feature_duplicate_symbol),
    paste(sort(unique(feature_candidate_mapped$feature_group)), collapse = ";")
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_df,
  file.path(outdir, "EMD_feature_gene_prepare_summary.tsv")
)

cat("Input EMD rows:", nrow(emd), "\n")
cat("Candidate FDR cutoff:", fdr_cutoff, "\n")
cat("Candidate SYMBOL x feature_group rows:", nrow(feature_candidate), "\n")
cat("Mapped output rows:", nrow(feature_candidate_mapped), "\n")
cat("Unmapped SYMBOL rows:", nrow(feature_unmapped), "\n")
cat("Feature groups:", paste(sort(unique(feature_candidate_mapped$feature_group)), collapse = ", "), "\n")
cat("Done.\n")
cat("Output directory:\n", outdir, "\n")

