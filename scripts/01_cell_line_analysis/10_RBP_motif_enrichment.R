suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Summarize FIMO hits to gene level and test RBP motif enrichment.
#
# Main statistical unit:
#   gene-level has_motif_hit
#
# Tests:
#   Fisher exact test for each database x feature_group x region_type x RBP
#   Motif-level Fisher tests are also exported as supplementary results.
#
# Feature/background gene sets:
#   Parsed from FASTA headers generated from gene-level merged BED intervals.
#   Background genes are scanned once per region and then feature genes are
#   removed per feature_group before enrichment testing.
# ============================================================

# -------------------------
# 1. Input and output paths
# -------------------------

manifest_file <- "./06_FIMO_RBP_motif_scan_parallel_nameOnly_EMDtested_Hep3B_background/fimo_scan_manifest.tsv"

candidate_rbp_file <- "./04_RBP_p009_candidate_heatmap_cell_samples/logFC2/RBP_candidates_from_effect_size_table_wilcox_p005_by_group.tsv"

cisbp_metadata_file <- "./05_offline_motifs/CisBP_RNA/CisBP_RNA_human_all_motif_metadata.tsv"
attract_metadata_file <- "./05_offline_motifs/ATtRACT/ATtRACT_human_nonmutated_motif_metadata.tsv"

outdir <- "./07_RBP_motif_enrichment_new_back"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 2. Parameters
# -------------------------

fimo_hit_p_cutoff <- 1e-4
primary_database <- "CisBP_RNA"
supplementary_database <- "ATtRACT"

# -------------------------
# 3. Helper functions
# -------------------------

normalize_symbol <- function(x) {
  toupper(trimws(as.character(x)))
}

map_target_to_feature_group <- function(x) {
  if (grepl("^K562$|K562", x, ignore.case = TRUE)) {
    return("K562-like")
  }
  if (grepl("HTR", x, ignore.case = TRUE)) {
    return("HTR-8-like")
  }
  if (grepl("Hep|liver", x, ignore.case = TRUE)) {
    return("liver-like")
  }
  paste0(gsub("[^A-Za-z0-9]+", "_", x), "-like")
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

parse_sequence_name <- function(sequence_name) {
  x <- as.character(sequence_name)
  x <- sub("^>", "", x)
  x <- sub("\\s.*$", "", x)
  x <- sub("::.*$", "", x)

  parts <- strsplit(x, "\\|")
  out <- lapply(parts, function(p) {
    p <- c(p, rep(NA_character_, max(0, 3 - length(p))))
    data.frame(
      ENSEMBL = p[1],
      SYMBOL = p[2],
      region_type_from_name = p[3],
      stringsAsFactors = FALSE
    )
  })

  bind_rows(out)
}

read_fasta_gene_table <- function(fasta_file) {
  if (!file.exists(fasta_file) || file.info(fasta_file)$size == 0) {
    return(data.frame(ENSEMBL = character(0), SYMBOL = character(0)))
  }

  headers <- readLines(fasta_file, warn = FALSE)
  headers <- headers[grepl("^>", headers)]

  if (length(headers) == 0) {
    return(data.frame(ENSEMBL = character(0), SYMBOL = character(0)))
  }

  parse_sequence_name(headers) %>%
    filter(!is.na(ENSEMBL), ENSEMBL != "", !is.na(SYMBOL), SYMBOL != "") %>%
    distinct(ENSEMBL, SYMBOL)
}

read_fimo_hits <- function(fimo_file, database, region_type, set_type, feature_group, motif_key_map) {
  empty <- data.frame(
    database = character(0),
    region_type = character(0),
    set_type = character(0),
    feature_group = character(0),
    motif_id = character(0),
    RBP_SYMBOL = character(0),
    ENSEMBL = character(0),
    SYMBOL = character(0),
    best_pvalue = numeric(0),
    hit_count = integer(0),
    stringsAsFactors = FALSE
  )

  if (!file.exists(fimo_file) || file.info(fimo_file)$size == 0) {
    return(empty)
  }

  fimo <- tryCatch(
    fread(fimo_file, data.table = FALSE, check.names = FALSE, fill = TRUE),
    error = function(e) NULL
  )

  if (is.null(fimo) || nrow(fimo) == 0 || !"motif_id" %in% colnames(fimo)) {
    return(empty)
  }

  p_col <- intersect(c("p-value", "p.value", "pvalue"), colnames(fimo))[1]
  if (is.na(p_col) || !"sequence_name" %in% colnames(fimo)) {
    return(empty)
  }

  fimo$pvalue_numeric <- suppressWarnings(as.numeric(fimo[[p_col]]))

  fimo <- fimo %>%
    filter(
      !is.na(motif_id),
      motif_id != "",
      !grepl("^#", motif_id),
      !is.na(sequence_name),
      sequence_name != "",
      is.finite(pvalue_numeric),
      pvalue_numeric <= fimo_hit_p_cutoff
    )

  if (nrow(fimo) == 0) {
    return(empty)
  }

  seq_info <- parse_sequence_name(fimo$sequence_name)

  fimo_gene <- bind_cols(
    fimo %>% select(motif_id, pvalue_numeric),
    seq_info
  ) %>%
    filter(!is.na(ENSEMBL), ENSEMBL != "", !is.na(SYMBOL), SYMBOL != "") %>%
    left_join(motif_key_map, by = c("motif_id" = "motif_key"))

  if (nrow(fimo_gene) == 0) {
    return(empty)
  }

  fimo_gene %>%
    filter(!is.na(RBP_SYMBOL), RBP_SYMBOL != "") %>%
    group_by(database = database, region_type = region_type, set_type = set_type, feature_group = feature_group, motif_id, RBP_SYMBOL, ENSEMBL, SYMBOL) %>%
    summarise(
      best_pvalue = min(pvalue_numeric, na.rm = TRUE),
      hit_count = n(),
      .groups = "drop"
    )
}

make_fisher_result <- function(a, feature_n, c, background_n) {
  b <- feature_n - a
  d <- background_n - c

  if (feature_n <= 0 || background_n <= 0 || a < 0 || b < 0 || c < 0 || d < 0) {
    return(list(odds_ratio = NA_real_, p_value = NA_real_))
  }

  mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  test <- tryCatch(
    fisher.test(mat, alternative = "greater"),
    error = function(e) NULL
  )

  if (is.null(test)) {
    return(list(odds_ratio = NA_real_, p_value = NA_real_))
  }

  list(
    odds_ratio = unname(test$estimate),
    p_value = test$p.value
  )
}

format_gene_list <- function(x) {
  x <- sort(unique(x[!is.na(x) & x != ""]))
  paste(x, collapse = ";")
}

make_motif_key_map <- function(motif_meta) {
  key_cols <- intersect(
    c("MEME_motif_id", "Motif_ID", "Matrix_id", "Matrix_ID", "matrix_id"),
    colnames(motif_meta)
  )

  if (length(key_cols) == 0) {
    stop("No usable motif ID columns found in motif metadata.")
  }

  key_list <- lapply(key_cols, function(key_col) {
    motif_meta %>%
      transmute(
        database = database,
        motif_key = as.character(.data[[key_col]]),
        motif_key_source = key_col,
        MEME_motif_id = MEME_motif_id,
        RBP_SYMBOL = RBP_SYMBOL,
        RBP_SYMBOL_norm = RBP_SYMBOL_norm
      )
  })

  bind_rows(key_list) %>%
    filter(!is.na(motif_key), motif_key != "") %>%
    distinct(database, motif_key, MEME_motif_id, RBP_SYMBOL_norm, .keep_all = TRUE)
}

diagnose_fimo_file <- function(fimo_file, database, region_type, set_type, feature_group, motif_key_map) {
  out <- data.frame(
    database = database,
    region_type = region_type,
    set_type = set_type,
    feature_group = feature_group,
    fimo_tsv = fimo_file,
    file_exists = file.exists(fimo_file),
    file_size = ifelse(file.exists(fimo_file), file.info(fimo_file)$size, NA_real_),
    raw_rows = 0L,
    p_pass_rows = 0L,
    parsed_gene_rows = 0L,
    motif_id_matched_rows = 0L,
    unique_fimo_motif_ids = 0L,
    unique_matched_rbps = 0L,
    example_fimo_motif_ids = "",
    stringsAsFactors = FALSE
  )

  if (!file.exists(fimo_file) || file.info(fimo_file)$size == 0) {
    return(out)
  }

  fimo <- tryCatch(
    fread(fimo_file, data.table = FALSE, check.names = FALSE, fill = TRUE),
    error = function(e) NULL
  )

  if (is.null(fimo) || nrow(fimo) == 0 || !"motif_id" %in% colnames(fimo)) {
    return(out)
  }

  p_col <- intersect(c("p-value", "p.value", "pvalue"), colnames(fimo))[1]
  if (is.na(p_col) || !"sequence_name" %in% colnames(fimo)) {
    out$raw_rows <- nrow(fimo)
    return(out)
  }

  fimo$pvalue_numeric <- suppressWarnings(as.numeric(fimo[[p_col]]))
  fimo2 <- fimo %>%
    filter(
      !is.na(motif_id),
      motif_id != "",
      !grepl("^#", motif_id),
      !is.na(sequence_name),
      sequence_name != "",
      is.finite(pvalue_numeric),
      pvalue_numeric <= fimo_hit_p_cutoff
    )

  seq_info <- parse_sequence_name(fimo2$sequence_name)
  joined <- bind_cols(
    fimo2 %>% select(motif_id),
    seq_info
  ) %>%
    filter(!is.na(ENSEMBL), ENSEMBL != "", !is.na(SYMBOL), SYMBOL != "") %>%
    left_join(motif_key_map, by = c("motif_id" = "motif_key"))

  out$raw_rows <- nrow(fimo)
  out$p_pass_rows <- nrow(fimo2)
  out$parsed_gene_rows <- nrow(joined)
  out$motif_id_matched_rows <- sum(!is.na(joined$RBP_SYMBOL))
  out$unique_fimo_motif_ids <- n_distinct(fimo2$motif_id)
  out$unique_matched_rbps <- n_distinct(joined$RBP_SYMBOL[!is.na(joined$RBP_SYMBOL)])
  out$example_fimo_motif_ids <- paste(head(sort(unique(as.character(fimo2$motif_id))), 10), collapse = ";")

  out
}

# -------------------------
# 4. Read metadata
# -------------------------

manifest <- fread(manifest_file, data.table = FALSE, check.names = FALSE)
candidate <- fread(candidate_rbp_file, data.table = FALSE, check.names = FALSE)
cisbp_meta <- fread(cisbp_metadata_file, data.table = FALSE, check.names = FALSE)
attract_meta <- fread(attract_metadata_file, data.table = FALSE, check.names = FALSE)

required_manifest_cols <- c("database", "region_type", "set_type", "feature_group", "fasta_file", "fimo_tsv", "status")
required_candidate_cols <- c("SYMBOL", "target_group")
required_motif_cols <- c("MEME_motif_id", "RBP_SYMBOL")

missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))
missing_candidate_cols <- setdiff(required_candidate_cols, colnames(candidate))
missing_cisbp_cols <- setdiff(required_motif_cols, colnames(cisbp_meta))
missing_attract_cols <- setdiff(required_motif_cols, colnames(attract_meta))

if (length(missing_manifest_cols) > 0) {
  stop("Missing columns in FIMO manifest: ", paste(missing_manifest_cols, collapse = ", "))
}
if (length(missing_candidate_cols) > 0) {
  stop("Missing columns in candidate RBP table: ", paste(missing_candidate_cols, collapse = ", "))
}
if (length(missing_cisbp_cols) > 0) {
  stop("Missing columns in CisBP-RNA motif metadata: ", paste(missing_cisbp_cols, collapse = ", "))
}
if (length(missing_attract_cols) > 0) {
  stop("Missing columns in ATtRACT motif metadata: ", paste(missing_attract_cols, collapse = ", "))
}

motif_meta <- bind_rows(
  cisbp_meta %>%
    mutate(database = primary_database),
  attract_meta %>%
    mutate(database = supplementary_database)
) %>%
  mutate(
    RBP_SYMBOL = as.character(RBP_SYMBOL),
    RBP_SYMBOL_norm = normalize_symbol(RBP_SYMBOL)
  ) %>%
  distinct(database, MEME_motif_id, .keep_all = TRUE)

motif_key_map <- make_motif_key_map(motif_meta)

candidate_rbp <- candidate %>%
  mutate(
    RBP_SYMBOL = SYMBOL,
    RBP_SYMBOL_norm = normalize_symbol(SYMBOL),
    feature_group = vapply(target_group, map_target_to_feature_group, character(1))
  ) %>%
  distinct(target_group, feature_group, RBP_SYMBOL, RBP_SYMBOL_norm)

candidate_motif_coverage <- merge(
  candidate_rbp,
  data.frame(database = unique(motif_meta$database), stringsAsFactors = FALSE),
  by = NULL
) %>%
  left_join(
    motif_meta %>%
      distinct(database, RBP_SYMBOL_norm) %>%
      mutate(has_database_motif = TRUE),
    by = c("database", "RBP_SYMBOL_norm")
  ) %>%
  mutate(has_database_motif = ifelse(is.na(has_database_motif), FALSE, has_database_motif)) %>%
  select(database, target_group, feature_group, RBP_SYMBOL, has_database_motif) %>%
  arrange(database, feature_group, RBP_SYMBOL)

write_tsv(
  candidate_motif_coverage,
  file.path(outdir, "candidate_RBP_motif_database_coverage.tsv")
)

# -------------------------
# 5. Read FASTA gene sets
# -------------------------

manifest_ok <- manifest %>%
  filter(status %in% c("ok", "ok_existing"))

fasta_gene_sets <- manifest_ok %>%
  distinct(region_type, set_type, feature_group, fasta_file) %>%
  rowwise() %>%
  mutate(gene_table = list(read_fasta_gene_table(fasta_file))) %>%
  ungroup()

gene_set_long <- bind_rows(lapply(seq_len(nrow(fasta_gene_sets)), function(i) {
  genes_i <- fasta_gene_sets$gene_table[[i]]
  if (nrow(genes_i) == 0) {
    return(NULL)
  }
  genes_i %>%
    mutate(
      region_type = fasta_gene_sets$region_type[i],
      set_type = fasta_gene_sets$set_type[i],
      feature_group = fasta_gene_sets$feature_group[i]
    )
}))

write_tsv(
  gene_set_long,
  file.path(outdir, "FASTA_gene_sets_parsed_from_headers.tsv")
)

# -------------------------
# 6. Read and summarize FIMO hits
# -------------------------

fimo_diagnostics <- bind_rows(lapply(seq_len(nrow(manifest_ok)), function(i) {
  diagnose_fimo_file(
    fimo_file = manifest_ok$fimo_tsv[i],
    database = manifest_ok$database[i],
    region_type = manifest_ok$region_type[i],
    set_type = manifest_ok$set_type[i],
    feature_group = manifest_ok$feature_group[i],
    motif_key_map = motif_key_map %>% filter(database == manifest_ok$database[i])
  )
}))

write_tsv(
  fimo_diagnostics,
  file.path(outdir, "FIMO_reading_diagnostics.tsv")
)

hit_list <- vector("list", nrow(manifest_ok))

for (i in seq_len(nrow(manifest_ok))) {
  db_i <- manifest_ok$database[i]
  motif_meta_i <- motif_meta %>%
    filter(database == db_i)
  motif_key_map_i <- motif_key_map %>%
    filter(database == db_i)

  message(
    "Reading FIMO hits: ",
    db_i, " ",
    manifest_ok$region_type[i], " ",
    manifest_ok$set_type[i], " ",
    manifest_ok$feature_group[i]
  )

  hit_list[[i]] <- read_fimo_hits(
    fimo_file = manifest_ok$fimo_tsv[i],
    database = db_i,
    region_type = manifest_ok$region_type[i],
    set_type = manifest_ok$set_type[i],
    feature_group = manifest_ok$feature_group[i],
    motif_key_map = motif_key_map_i
  )
}

motif_gene_hits <- bind_rows(hit_list) %>%
  mutate(RBP_SYMBOL_norm = normalize_symbol(RBP_SYMBOL))

candidate_hits <- motif_gene_hits %>%
  inner_join(
    candidate_rbp %>% select(feature_group, RBP_SYMBOL_norm) %>% distinct(),
    by = c("feature_group", "RBP_SYMBOL_norm")
  )

write_tsv(
  candidate_hits,
  file.path(outdir, "candidate_RBP_motif_gene_level_hits.tsv")
)

# -------------------------
# 7. Fisher enrichment tests
# -------------------------

feature_sets <- gene_set_long %>%
  filter(set_type == "feature") %>%
  distinct(feature_group, region_type, ENSEMBL, SYMBOL)

background_sets <- gene_set_long %>%
  filter(set_type == "background", feature_group == "ALL_BACKGROUND") %>%
  distinct(region_type, ENSEMBL, SYMBOL)

test_grid <- expand.grid(
  database = sort(unique(motif_meta$database)),
  feature_group = sort(unique(feature_sets$feature_group)),
  region_type = c("three_prime_UTR", "five_prime_UTR", "CDS"),
  stringsAsFactors = FALSE
)

rbp_result_list <- list()
motif_result_list <- list()

for (i in seq_len(nrow(test_grid))) {
  db_i <- test_grid$database[i]
  group_i <- test_grid$feature_group[i]
  region_i <- test_grid$region_type[i]

  feature_genes <- feature_sets %>%
    filter(feature_group == group_i, region_type == region_i) %>%
    distinct(ENSEMBL, SYMBOL)

  background_genes <- background_sets %>%
    filter(region_type == region_i) %>%
    anti_join(feature_genes, by = c("ENSEMBL", "SYMBOL")) %>%
    distinct(ENSEMBL, SYMBOL)

  candidate_rbp_i <- candidate_rbp %>%
    filter(feature_group == group_i) %>%
    distinct(RBP_SYMBOL, RBP_SYMBOL_norm)

  motif_meta_i <- motif_meta %>%
    filter(database == db_i) %>%
    inner_join(candidate_rbp_i, by = "RBP_SYMBOL_norm")

  if (nrow(feature_genes) == 0 || nrow(background_genes) == 0 || nrow(motif_meta_i) == 0) {
    next
  }

  feature_hits_i <- motif_gene_hits %>%
    filter(
      database == db_i,
      region_type == region_i,
      set_type == "feature",
      feature_group == group_i
    )

  background_hits_i <- motif_gene_hits %>%
    filter(
      database == db_i,
      region_type == region_i,
      set_type == "background",
      feature_group == "ALL_BACKGROUND"
    ) %>%
    anti_join(feature_genes, by = c("ENSEMBL", "SYMBOL"))

  rbps_to_test <- motif_meta_i %>%
    distinct(RBP_SYMBOL = RBP_SYMBOL.x, RBP_SYMBOL_norm)

  rbp_rows <- lapply(seq_len(nrow(rbps_to_test)), function(j) {
    rbp_norm <- rbps_to_test$RBP_SYMBOL_norm[j]
    rbp_symbol <- rbps_to_test$RBP_SYMBOL[j]

    feature_hit_genes <- feature_hits_i %>%
      filter(RBP_SYMBOL_norm == rbp_norm) %>%
      distinct(ENSEMBL, SYMBOL)

    background_hit_genes <- background_hits_i %>%
      filter(RBP_SYMBOL_norm == rbp_norm) %>%
      distinct(ENSEMBL, SYMBOL)

    a <- nrow(feature_hit_genes)
    c <- nrow(background_hit_genes)
    fisher <- make_fisher_result(a, nrow(feature_genes), c, nrow(background_genes))

    motif_ids <- motif_meta_i %>%
      filter(RBP_SYMBOL_norm == rbp_norm) %>%
      pull(MEME_motif_id) %>%
      unique()

    data.frame(
      database = db_i,
      feature_group = group_i,
      region_type = region_i,
      RBP_SYMBOL = rbp_symbol,
      motif_n = length(motif_ids),
      feature_gene_n = nrow(feature_genes),
      background_gene_n = nrow(background_genes),
      feature_hit_gene_n = a,
      background_hit_gene_n = c,
      odds_ratio = fisher$odds_ratio,
      p_value = fisher$p_value,
      hit_feature_genes = format_gene_list(feature_hit_genes$SYMBOL),
      motif_ids = paste(sort(motif_ids), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  rbp_result_list[[length(rbp_result_list) + 1L]] <- bind_rows(rbp_rows)

  motifs_to_test <- motif_meta_i %>%
    distinct(MEME_motif_id, RBP_SYMBOL = RBP_SYMBOL.x, RBP_SYMBOL_norm)

  motif_rows <- lapply(seq_len(nrow(motifs_to_test)), function(j) {
    motif_id_i <- motifs_to_test$MEME_motif_id[j]
    rbp_symbol <- motifs_to_test$RBP_SYMBOL[j]

    feature_hit_genes <- feature_hits_i %>%
      filter(.data$motif_id == .env$motif_id_i) %>%
      distinct(ENSEMBL, SYMBOL)

    background_hit_genes <- background_hits_i %>%
      filter(.data$motif_id == .env$motif_id_i) %>%
      distinct(ENSEMBL, SYMBOL)

    a <- nrow(feature_hit_genes)
    c <- nrow(background_hit_genes)
    fisher <- make_fisher_result(a, nrow(feature_genes), c, nrow(background_genes))

    data.frame(
      database = db_i,
      feature_group = group_i,
      region_type = region_i,
      RBP_SYMBOL = rbp_symbol,
      motif_id = motif_id_i,
      feature_gene_n = nrow(feature_genes),
      background_gene_n = nrow(background_genes),
      feature_hit_gene_n = a,
      background_hit_gene_n = c,
      odds_ratio = fisher$odds_ratio,
      p_value = fisher$p_value,
      hit_feature_genes = format_gene_list(feature_hit_genes$SYMBOL),
      stringsAsFactors = FALSE
    )
  })

  motif_result_list[[length(motif_result_list) + 1L]] <- bind_rows(motif_rows)
}

rbp_enrichment <- bind_rows(rbp_result_list) %>%
  group_by(database, feature_group, region_type) %>%
  mutate(FDR_by_database_group_region = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(database) %>%
  mutate(FDR_by_database_global = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(database, feature_group, region_type, FDR_by_database_group_region, desc(odds_ratio), desc(feature_hit_gene_n))

motif_enrichment <- bind_rows(motif_result_list) %>%
  group_by(database, feature_group, region_type) %>%
  mutate(FDR_by_database_group_region = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  group_by(database) %>%
  mutate(FDR_by_database_global = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(database, feature_group, region_type, FDR_by_database_group_region, desc(odds_ratio), desc(feature_hit_gene_n))

# -------------------------
# 8. Write outputs
# -------------------------

write_tsv(
  rbp_enrichment,
  file.path(outdir, "RBP_level_motif_enrichment_Fisher_results.tsv")
)

write_tsv(
  motif_enrichment,
  file.path(outdir, "motif_level_enrichment_Fisher_results.tsv")
)

write_tsv(
  rbp_enrichment %>% filter(database == primary_database),
  file.path(outdir, "MAIN_CisBP_RNA_RBP_level_motif_enrichment.tsv")
)

write_tsv(
  rbp_enrichment %>% filter(database == supplementary_database),
  file.path(outdir, "SUPPLEMENT_ATtRACT_RBP_level_motif_enrichment.tsv")
)

summary_df <- data.frame(
  item = c(
    "FIMO manifest rows",
    "Candidate target_group x RBP rows",
    "Candidate RBP motif coverage rows",
    "Parsed FASTA gene set rows",
    "Candidate motif gene-level hit rows",
    "RBP-level enrichment rows",
    "Motif-level enrichment rows",
    "FIMO diagnostic p-pass rows",
    "FIMO diagnostic motif-id matched rows",
    "FIMO hit p-value cutoff"
  ),
  value = c(
    nrow(manifest),
    nrow(candidate_rbp),
    nrow(candidate_motif_coverage),
    nrow(gene_set_long),
    nrow(candidate_hits),
    nrow(rbp_enrichment),
    nrow(motif_enrichment),
    sum(fimo_diagnostics$p_pass_rows, na.rm = TRUE),
    sum(fimo_diagnostics$motif_id_matched_rows, na.rm = TRUE),
    fimo_hit_p_cutoff
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  summary_df,
  file.path(outdir, "RBP_motif_enrichment_summary.tsv")
)

cat("FIMO manifest rows:", nrow(manifest), "\n")
cat("Candidate target_group x RBP rows:", nrow(candidate_rbp), "\n")
cat("Parsed FASTA gene set rows:", nrow(gene_set_long), "\n")
cat("Candidate motif gene-level hit rows:", nrow(candidate_hits), "\n")
cat("RBP-level enrichment rows:", nrow(rbp_enrichment), "\n")
cat("Motif-level enrichment rows:", nrow(motif_enrichment), "\n")
cat("Done.\n")
cat("Output directory:\n", outdir, "\n")


