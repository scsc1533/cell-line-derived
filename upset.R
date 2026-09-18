library(data.table)
library(tidyverse)
library(ComplexUpset)

base_dir <- "/data/work/01_2603cell_culture/06_xinyue/01_mlRNA_cfRNA_vs_cell_debris"

wilcox_file <- "/data/work/01_2603cell_culture/03_analysis/01_Component/02_secreted/02_tpm_Component_comparison/all_celltypes_cfRNA_vs_debris_cell_Wilcoxon_results.txt"

outdir <- file.path(base_dir, "result_overlap")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- c(
  "Hep3B2.1-7",
  "HepG2",
  "K562",
  "HTR-8_SVneo",
  "HEK293T"
)


############################
# 1. Wilcoxon initial filter
############################

wilcox <- fread(wilcox_file, sep = "\t")

wilcox[
  Cell_type == "HTR-8/SVneo",
  Cell_type := "HTR-8_SVneo"
]

wilcox_filter <- wilcox[
  abs(log2FC) > 1 &
    p_adjusted < 0.05 &
    (
      direction == "up_in_cfRNA" & mean_cfRNA > 1 |
      direction == "up_in_debris_cell" & mean_debris_cell > 1
    )
]

fwrite(
  wilcox_filter,
  file.path(outdir, "01_wilcoxon_filtered.tsv"),
  sep = "\t"
)

cat("Wilcoxon direction count:\n")
print(table(wilcox_filter$direction))


############################
# 2. Method-supported genes
############################

get_support_gene <- function(cell, deg_dir) {
  
  f <- file.path(
    base_dir,
    cell,
    paste0(
      "gene_overlap_detailed_",
      deg_dir,
      "_logFC1_adj0.05.csv"
    )
  )
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(character(0))
  }
  
  x <- fread(f)
  
  method_cols <- c("DESeq2", "edgeR", "limma_voom")
  
  x[, (method_cols) := lapply(.SD, as.logical), .SDcols = method_cols]
  
  x[, method_n := rowSums(.SD, na.rm = TRUE), .SDcols = method_cols]
  
  unique(x[method_n >= 2, Gene])
}


support <- list()

for (cell in cells) {
  support[[paste(cell, "up", sep = "_")]] <- get_support_gene(cell, "up")
  support[[paste(cell, "down", sep = "_")]] <- get_support_gene(cell, "down")
}


############################
# 3. Overlap + upset
############################

overlap_plot <- function(df, dir_filter, prefix) {
  
  x <- df[direction == dir_filter]
  
  cat("\n", prefix, "\n")
  cat("direction:", dir_filter, "\n")
  cat("rows:", nrow(x), "\n")
  cat("genes:", length(unique(x$gene)), "\n")
  
  if (nrow(x) == 0) {
    warning("No genes for ", dir_filter)
    return(NULL)
  }
  
  sets <- lapply(
    cells,
    function(cell) unique(x[Cell_type == cell, gene])
  )
  
  names(sets) <- cells
  
  print(sapply(sets, length))
  
  genes <- unique(unlist(sets))
  
  long <- rbindlist(
    lapply(
      genes,
      function(g) {
        data.table(
          gene = g,
          cell_type = cells,
          in_set = sapply(sets, function(s) g %in% s)
        )
      }
    )
  )
  
  fwrite(
    long,
    file.path(outdir, paste0(prefix, "_overlap_long.tsv")),
    sep = "\t"
  )
  
  mat <- long %>%
    pivot_wider(
      names_from = cell_type,
      values_from = in_set,
      values_fill = FALSE
    )
  
  pdf(
    file.path(outdir, paste0(prefix, "_upset.pdf")),
    width = 15,
    height = 9
  )
  
  print(
    upset(
      mat,
      cells,
      name = "Genes"
    )
  )
  
  dev.off()
  
  return(long)
}


############################
# 4. Raw overlap
############################

overlap_plot(
  wilcox_filter,
  "up_in_cfRNA",
  "02_raw_up_in_cfRNA"
)

overlap_plot(
  wilcox_filter,
  "up_in_debris_cell",
  "03_raw_up_in_debris_cell"
)


############################
# 5. Second filter:
# Wilcoxon filtered + >=2 DEG methods
############################

keep <- logical(nrow(wilcox_filter))

for (i in seq_len(nrow(wilcox_filter))) {
  
  row <- wilcox_filter[i]
  
  deg_dir <- ifelse(
    row$direction == "up_in_cfRNA",
    "up",
    "down"
  )
  
  key <- paste(
    row$Cell_type,
    deg_dir,
    sep = "_"
  )
  
  keep[i] <- row$gene %in% support[[key]]
}

filtered <- wilcox_filter[keep]

fwrite(
  filtered,
  file.path(outdir, "04_method_supported_wilcoxon_filtered.tsv"),
  sep = "\t"
)

cat("\nFiltered direction count:\n")
print(table(filtered$direction))


############################
# 6. Filtered overlap
############################

overlap_plot(
  filtered,
  "up_in_cfRNA",
  "05_filtered_up_in_cfRNA"
)

overlap_plot(
  filtered,
  "up_in_debris_cell",
  "06_filtered_up_in_debris_cell"
)

cat("\nFinished.\nOutput dir:\n", outdir, "\n")