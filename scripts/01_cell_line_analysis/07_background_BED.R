suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Extract gene-level merged 3'UTR / 5'UTR / CDS BED from GENCODE v43 GTF
# and prepare EMD-tested background genes for motif enrichment.
#
# Region strategy:
#   gene-level merged regions across all transcripts of the same gene.
#
# Background gene definition:
#   genes from all_comparisons_gene_scores.tsv with comparison == Hep3B2.1_7_vs_rest.
#   This defines the common EMD-tested gene universe for motif enrichment.
#
# Main outputs:
#   1. all gene-level merged region BED files
#   2. EMD-tested common background gene BED files
#   3. feature-group-specific feature/background BED files
# ============================================================

# -------------------------
# 1. Input and output paths
# -------------------------

gtf_file <- "./GRCh38_GENCODE_20231021/gtf_gff/gencode.v43.annotation.gtf"
emd_gene_score_file <- "./all_comparisons_gene_scores.tsv"
feature_gene_file <- "./EMD_feature_genes_wasserstein_FDR005_with_ENSEMBL.tsv"

outdir <- "./03_GENCODE_v43_UTR_CDS_BED_EMDtested_Hep3B_background"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 2. Parameters
# -------------------------

background_comparison <- "Hep3B2.1_7_vs_rest"
background_label <- "EMDtested_Hep3B2.1_7_vs_rest"

region_levels <- c("three_prime_UTR", "five_prime_UTR", "CDS")

# -------------------------
# 3. Helper functions
# -------------------------

extract_gtf_attr <- function(attr, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  str_match(attr, pattern)[, 2]
}

clean_ensembl <- function(x) {
  sub("\\..*$", "", x)
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

reduce_interval_vectors <- function(start0, end) {
  keep <- is.finite(start0) & is.finite(end) & end > start0
  start0 <- start0[keep]
  end <- end[keep]

  if (length(start0) == 0) {
    return(data.frame(start0 = numeric(0), end = numeric(0)))
  }

  ord <- order(start0, end)
  start0 <- start0[ord]
  end <- end[ord]

  out_start <- numeric(0)
  out_end <- numeric(0)
  current_start <- start0[1]
  current_end <- end[1]

  if (length(start0) > 1) {
    for (i in 2:length(start0)) {
      if (start0[i] <= current_end) {
        current_end <- max(current_end, end[i])
      } else {
        out_start <- c(out_start, current_start)
        out_end <- c(out_end, current_end)
        current_start <- start0[i]
        current_end <- end[i]
      }
    }
  }

  out_start <- c(out_start, current_start)
  out_end <- c(out_end, current_end)

  data.frame(start0 = out_start, end = out_end)
}

reduce_region_dt <- function(dt) {
  dt <- as.data.table(dt)
  dt <- dt[is.finite(start0) & is.finite(end) & end > start0]

  reduced <- dt[
    ,
    {
      x <- reduce_interval_vectors(start0, end)
      .(start0 = x$start0, end = x$end)
    },
    by = .(seqname, strand, ENSEMBL, SYMBOL, gene_type, region_type)
  ]

  reduced[
    ,
    region_length := sum(end - start0),
    by = .(ENSEMBL, region_type)
  ]

  setorder(reduced, seqname, start0, end, strand, ENSEMBL, region_type)
  reduced[]
}

make_bed <- function(dt) {
  dt <- as.data.table(dt)
  bed <- dt[
    ,
    .(
      chrom = seqname,
      start = as.integer(start0),
      end = as.integer(end),
      name = paste(ENSEMBL, SYMBOL, region_type, sep = "|"),
      score = 0,
      strand = strand,
      ENSEMBL = ENSEMBL,
      SYMBOL = SYMBOL,
      gene_type = gene_type,
      region_type = region_type
    )
  ]
  setorder(bed, chrom, start, end, strand, ENSEMBL)
  bed[]
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

write_bed <- function(dt, file) {
  write.table(
    dt,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

# -------------------------
# 4. Read GTF
# -------------------------

gtf <- fread(
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
)

gtf <- gtf %>%
  filter(feature %in% c("gene", "transcript", "exon", "CDS", "UTR", "five_prime_UTR", "three_prime_UTR", "five_prime_utr", "three_prime_utr"))

gtf$gene_id_raw <- extract_gtf_attr(gtf$attribute, "gene_id")
gtf$transcript_id_raw <- extract_gtf_attr(gtf$attribute, "transcript_id")
gtf$ENSEMBL <- clean_ensembl(gtf$gene_id_raw)
gtf$transcript_id <- clean_ensembl(gtf$transcript_id_raw)
gtf$SYMBOL <- extract_gtf_attr(gtf$attribute, "gene_name")
gtf$gene_type <- extract_gtf_attr(gtf$attribute, "gene_type")

gene_map <- gtf %>%
  filter(feature == "gene") %>%
  select(ENSEMBL, SYMBOL, gene_type, seqname, start, end, strand) %>%
  filter(!is.na(ENSEMBL), ENSEMBL != "", !is.na(SYMBOL), SYMBOL != "") %>%
  distinct()

write_tsv(
  gene_map,
  file.path(outdir, "GENCODE_v43_gene_map.tsv")
)

# -------------------------
# 5. Extract transcript-level CDS and UTR
# -------------------------

cds_tx <- gtf %>%
  filter(feature == "CDS", !is.na(transcript_id), transcript_id != "") %>%
  transmute(
    seqname,
    start = as.integer(start),
    end = as.integer(end),
    strand,
    ENSEMBL,
    transcript_id,
    SYMBOL,
    gene_type,
    region_type = "CDS",
    start0 = as.integer(start) - 1L
  )

cds_boundary <- cds_tx %>%
  group_by(transcript_id) %>%
  summarise(
    cds_min = min(start, na.rm = TRUE),
    cds_max = max(end, na.rm = TRUE),
    .groups = "drop"
  )

direct_utr <- gtf %>%
  filter(feature %in% c("five_prime_UTR", "three_prime_UTR", "five_prime_utr", "three_prime_utr")) %>%
  transmute(
    seqname,
    start = as.integer(start),
    end = as.integer(end),
    strand,
    ENSEMBL,
    transcript_id,
    SYMBOL,
    gene_type,
    region_type = case_when(
      feature %in% c("five_prime_UTR", "five_prime_utr") ~ "five_prime_UTR",
      feature %in% c("three_prime_UTR", "three_prime_utr") ~ "three_prime_UTR",
      TRUE ~ NA_character_
    ),
    start0 = as.integer(start) - 1L
  ) %>%
  filter(!is.na(region_type))

generic_utr <- gtf %>%
  filter(feature == "UTR", !is.na(transcript_id), transcript_id != "") %>%
  left_join(cds_boundary, by = "transcript_id") %>%
  filter(is.finite(cds_min), is.finite(cds_max)) %>%
  mutate(
    region_type = case_when(
      strand == "+" & end < cds_min ~ "five_prime_UTR",
      strand == "+" & start > cds_max ~ "three_prime_UTR",
      strand == "-" & end < cds_min ~ "three_prime_UTR",
      strand == "-" & start > cds_max ~ "five_prime_UTR",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(region_type)) %>%
  transmute(
    seqname,
    start = as.integer(start),
    end = as.integer(end),
    strand,
    ENSEMBL,
    transcript_id,
    SYMBOL,
    gene_type,
    region_type,
    start0 = as.integer(start) - 1L
  )

utr_tx <- bind_rows(direct_utr, generic_utr) %>%
  distinct(seqname, start, end, strand, ENSEMBL, transcript_id, SYMBOL, gene_type, region_type, start0)

region_tx <- bind_rows(cds_tx, utr_tx) %>%
  filter(
    region_type %in% region_levels,
    !is.na(ENSEMBL),
    ENSEMBL != "",
    !is.na(SYMBOL),
    SYMBOL != "",
    end > start0
  )

region_merged <- reduce_region_dt(region_tx)

write_tsv(
  region_merged,
  file.path(outdir, "GENCODE_v43_gene_level_merged_3UTR_5UTR_CDS_regions.tsv")
)

for (region_i in region_levels) {
  bed_i <- make_bed(region_merged[region_type == region_i])
  write_bed(
    bed_i,
    file.path(outdir, paste0("GENCODE_v43_gene_level_merged_", region_i, ".bed"))
  )
}

# -------------------------
# 6. Define common background genes from EMD-tested gene universe
# -------------------------

emd_score <- fread(emd_gene_score_file, data.table = FALSE, check.names = FALSE)

required_emd_cols <- c("comparison", "Gene")
missing_emd_cols <- setdiff(required_emd_cols, colnames(emd_score))
if (length(missing_emd_cols) > 0) {
  stop("Missing columns in EMD gene score table: ", paste(missing_emd_cols, collapse = ", "))
}

background_source <- emd_score %>%
  filter(comparison == background_comparison) %>%
  mutate(SYMBOL = Gene)

optional_background_cols <- c(
  "n_valid_samples",
  "target_mean_count",
  "rest_mean_count",
  "perm_p_wasserstein",
  "perm_fdr_wasserstein",
  "rank_wasserstein"
)

for (col_i in optional_background_cols) {
  if (!col_i %in% colnames(background_source)) {
    background_source[[col_i]] <- NA
  }
}

background_source <- background_source %>%
  select(comparison, SYMBOL, all_of(optional_background_cols)) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(SYMBOL, .keep_all = TRUE)

if (nrow(background_source) == 0) {
  stop("No genes found for background_comparison: ", background_comparison)
}

background_genes <- background_source %>%
  left_join(gene_map %>% select(SYMBOL, ENSEMBL, gene_type), by = "SYMBOL") %>%
  filter(!is.na(ENSEMBL), ENSEMBL != "") %>%
  distinct(SYMBOL, ENSEMBL, gene_type, .keep_all = TRUE)

background_regions <- region_merged %>%
  semi_join(background_genes, by = c("ENSEMBL", "SYMBOL", "gene_type"))

write_tsv(
  background_source,
  file.path(outdir, paste0("background_genes_", background_label, "_source_SYMBOL.tsv"))
)

write_tsv(
  background_genes,
  file.path(outdir, paste0("background_genes_", background_label, "_mapped.tsv"))
)

write_tsv(
  background_regions,
  file.path(outdir, paste0("background_regions_", background_label, "_gene_level_merged.tsv"))
)

for (region_i in region_levels) {
  bed_i <- make_bed(background_regions[background_regions$region_type == region_i, ])
  write_bed(
    bed_i,
    file.path(outdir, paste0("background_", background_label, "_", region_i, ".bed"))
  )
}
# -------------------------
# 7. Prepare feature/background BED per EMD feature group
# -------------------------

if (file.exists(feature_gene_file)) {
  feature_gene <- fread(feature_gene_file, data.table = FALSE, check.names = FALSE)

  required_feature_cols <- c("feature_group", "SYMBOL", "ENSEMBL")
  missing_feature_cols <- setdiff(required_feature_cols, colnames(feature_gene))
  if (length(missing_feature_cols) > 0) {
    stop("Missing columns in feature gene table: ", paste(missing_feature_cols, collapse = ", "))
  }

  feature_gene <- feature_gene %>%
    filter(!is.na(ENSEMBL), ENSEMBL != "", !is.na(SYMBOL), SYMBOL != "") %>%
    distinct(feature_group, SYMBOL, ENSEMBL, .keep_all = TRUE)

  feature_gene_before_common_background_n <- nrow(feature_gene)

  feature_gene <- feature_gene %>%
    semi_join(background_genes %>% distinct(ENSEMBL, SYMBOL), by = c("ENSEMBL", "SYMBOL"))

  write_tsv(
    feature_gene,
    file.path(outdir, paste0("EMD_feature_genes_within_common_background_", background_label, ".tsv"))
  )

  feature_groups <- sort(unique(feature_gene$feature_group))

  for (group_i in feature_groups) {
    group_safe <- sanitize_filename(group_i)

    group_feature_gene <- feature_gene %>%
      filter(feature_group == group_i) %>%
      distinct(ENSEMBL, SYMBOL)

    group_feature_regions <- region_merged %>%
      semi_join(group_feature_gene, by = c("ENSEMBL", "SYMBOL"))

    group_background_regions <- background_regions %>%
      anti_join(group_feature_gene, by = c("ENSEMBL", "SYMBOL"))

    write_tsv(
      group_feature_regions,
      file.path(outdir, paste0(group_safe, "_feature_regions_gene_level_merged.tsv"))
    )

    write_tsv(
      group_background_regions,
      file.path(outdir, paste0(group_safe, "_background_regions_EMDtested_Hep3B2.1_7_vs_rest_gene_level_merged.tsv"))
    )

    for (region_i in region_levels) {
      feature_bed_i <- make_bed(group_feature_regions[group_feature_regions$region_type == region_i, ])
      background_bed_i <- make_bed(group_background_regions[group_background_regions$region_type == region_i, ])

      write_bed(
        feature_bed_i,
        file.path(outdir, paste0(group_safe, "_feature_", region_i, ".bed"))
      )

      write_bed(
        background_bed_i,
        file.path(outdir, paste0(group_safe, "_background_EMDtested_Hep3B2.1_7_vs_rest_", region_i, ".bed"))
      )
    }
  }

  feature_region_count <- feature_gene %>%
    left_join(
      region_merged %>%
        distinct(ENSEMBL, SYMBOL, region_type),
      by = c("ENSEMBL", "SYMBOL")
    ) %>%
    count(feature_group, region_type, name = "gene_region_n") %>%
    arrange(feature_group, region_type)

  write_tsv(
    feature_region_count,
    file.path(outdir, "EMD_feature_gene_region_count_by_group.tsv")
  )
} else {
  warning("Feature gene file does not exist. Feature-group-specific BED files were skipped: ", feature_gene_file)
}

# -------------------------
# 8. Summary
# -------------------------

region_count <- region_merged %>%
  distinct(ENSEMBL, SYMBOL, region_type) %>%
  count(region_type, name = "gene_n") %>%
  arrange(region_type)

background_region_count <- background_regions %>%
  distinct(ENSEMBL, SYMBOL, region_type) %>%
  count(region_type, name = "background_gene_n") %>%
  arrange(region_type)

write_tsv(
  region_count,
  file.path(outdir, "GENCODE_v43_region_gene_count.tsv")
)

write_tsv(
  background_region_count,
  file.path(outdir, "background_region_gene_count.tsv")
)

summary_df <- data.frame(
  item = c(
    "GTF rows retained",
    "GENCODE genes in map",
    "Transcript-level region rows",
    "Gene-level merged region intervals",
    "Background source file",
    "Background comparison",
    "Background source SYMBOLs",
    "Background mapped genes",
    "Background merged region intervals"
  ),
  value = c(
    nrow(gtf),
    nrow(gene_map),
    nrow(region_tx),
    nrow(region_merged),
    emd_gene_score_file,
    background_comparison,
    nrow(background_source),
    nrow(background_genes),
    nrow(background_regions)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_df,
  file.path(outdir, "extract_regions_and_background_summary.tsv")
)

cat("GTF rows retained:", nrow(gtf), "\n")
cat("GENCODE genes in map:", nrow(gene_map), "\n")
cat("Transcript-level region rows:", nrow(region_tx), "\n")
cat("Gene-level merged region intervals:", nrow(region_merged), "\n")
cat("Background comparison:", background_comparison, "\n")
cat("Background source SYMBOLs:", nrow(background_source), "\n")
cat("Background mapped genes:", nrow(background_genes), "\n")
cat("Background merged region intervals:", nrow(background_regions), "\n")
cat("Done.\n")
cat("Output directory:\n", outdir, "\n")
