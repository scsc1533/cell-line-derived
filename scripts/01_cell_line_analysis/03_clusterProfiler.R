# =========================================================
# Script: GO and KEGG enrichment analysis for up-regulated
#         genes in cfRNA and cell (debris_cell) gene sets
# Purpose:
#   Read up-regulated gene lists (cfRNA and cell),
#   convert gene symbols to ENTREZ IDs, run GO (BP/CC/MF)
#   and KEGG enrichment, save result tables (CSV) and
#   bubble plots (PDF) for each gene set.
# =========================================================

# ===================== Load packages =====================
library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)   # Human (change if using another species)
library(ggplot2)

# ===================== Set paths and parameters =====================
data_dir <- "./07_enrichment_all5"          # Please modify to your actual path
output_dir <- file.path(data_dir, "enrichment_results_all_genes")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

p_cutoff <- 0.05      # p-value cutoff for GO/KEGG
q_cutoff <- 0.2       # q-value cutoff

# ===================== Read data =====================
cfRNA_file <- file.path(data_dir, "up_in_cfRNA_all5_gene_list.txt")
cell_file <- file.path(data_dir, "up_in_debris_cell_all5_gene_list.txt")

df_cfRNA <- fread(cfRNA_file, header = TRUE, data.table = FALSE)
df_cell <- fread(cell_file, header = TRUE, data.table = FALSE)

# Check required columns
if (!"gene" %in% colnames(df_cfRNA)) 
  stop("cfRNA_up_genes.txt file is missing the 'gene' column")
if (!"gene" %in% colnames(df_cell)) 
  stop("cell_up_genes.txt file is missing the 'gene' column")

# ===================== Gene ID conversion function =====================
convert_genes <- function(gene_symbols) {
  gene_symbols <- unique(na.omit(gene_symbols))
  if (length(gene_symbols) == 0) return(NULL)
  suppressMessages({
    mapped <- bitr(gene_symbols, fromType = "SYMBOL", 
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  })
  return(mapped)
}

# ===================== Enrichment analysis + bubble plot function =====================
run_enrichment <- function(gene_list, gene_set_name, sample_name, output_dir) {
  # gene_list: vector of gene symbols (already deduplicated)
  # gene_set_name: "cfRNA" or "cell"
  # sample_name: sample identifier (fixed as "all" in this case)
  
  cat("\n===== Analyzing:", sample_name, "-", gene_set_name, " (gene count:", length(gene_list), ")=====\n")
  
  if (length(gene_list) == 0) {
    cat("  No genes, skipping\n")
    return(NULL)
  }
  
  # Convert IDs
  mapped <- convert_genes(gene_list)
  if (is.null(mapped) || nrow(mapped) == 0) {
    cat("  No valid gene IDs, skipping\n")
    return(NULL)
  }
  entrez_ids <- mapped$ENTREZID
  cat("  Number of valid ENTREZIDs:", length(entrez_ids), "\n")
  
  # Store enrichment objects and result tables
  enrich_objects <- list()
  res_list <- list()
  
  # ---------- GO enrichment (BP, CC, MF) ----------
  go_ontologies <- c("BP", "CC", "MF")
  for (onto in go_ontologies) {
    cat("  Running GO:", onto, "...")
    tryCatch({
      ego <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Hs.eg.db,
                      ont = onto,
                      pAdjustMethod = "BH",
                      pvalueCutoff = p_cutoff,
                      qvalueCutoff = q_cutoff,
                      readable = TRUE)
      if (!is.null(ego) && nrow(ego) > 0) {
        res_df <- as.data.frame(ego)
        res_df$Sample <- sample_name
        res_df$Gene_set <- gene_set_name
        res_df$Ontology <- onto
        res_list[[paste0("GO_", onto)]] <- res_df
        enrich_objects[[paste0("GO_", onto)]] <- ego
        cat(" Significant terms found:", nrow(res_df), "\n")
      } else {
        cat(" No significant terms\n")
      }
    }, error = function(e) {
      cat(" Error:", e$message, "\n")
    })
  }
  
  # ---------- KEGG pathway enrichment ----------
  cat("  Running KEGG...")
  tryCatch({
    ekg <- enrichKEGG(gene = entrez_ids,
                      organism = "hsa",
                      pvalueCutoff = p_cutoff,
                      qvalueCutoff = q_cutoff)
    if (!is.null(ekg) && nrow(ekg) > 0) {
      res_df <- as.data.frame(ekg)
      res_df$Sample <- sample_name
      res_df$Gene_set <- gene_set_name
      res_df$Ontology <- "KEGG"
      res_list[["KEGG"]] <- res_df
      enrich_objects[["KEGG"]] <- ekg
      cat(" Significant terms found:", nrow(res_df), "\n")
    } else {
      cat(" No significant terms\n")
    }
  }, error = function(e) {
    cat(" Error:", e$message, "\n")
  })
  
  # ---------- Save CSV results ----------
  if (length(res_list) > 0) {
    combined <- bind_rows(res_list)
    out_csv <- file.path(output_dir, paste0(gene_set_name, "_enrichment.csv"))
    fwrite(combined, out_csv)
    cat("  Result CSV saved to:", out_csv, "\n")
  } else {
    cat("  No significant enrichment results, skipping plotting\n")
    return(NULL)
  }
  
  # ---------- Draw bubble plots (for each Ontology with significant results) ----------
  for (name in names(enrich_objects)) {
    obj <- enrich_objects[[name]]
    if (is.null(obj) || nrow(obj) == 0) next
    
    # Use clusterProfiler's dotplot
    p <- dotplot(obj, showCategory = 15, orderBy = "p.adjust") +
      ggtitle(paste(sample_name, gene_set_name, name, sep = " | ")) +
      theme(plot.title = element_text(size = 10))
    
    # Save as PDF
    safe_name <- gsub("[^A-Za-z0-9]", "_", paste(sample_name, gene_set_name, name, sep = "_"))
    out_pdf <- file.path(output_dir, paste0(safe_name, "_bubble.pdf"))
    ggsave(out_pdf, p, width = 10, height = 7)
    cat("  Bubble plot saved to:", out_pdf, "\n")
  }
  
  return(invisible(res_list))
}

# ===================== Main program: analyze all genes directly =====================
# Extract genes (deduplicate)
genes_cf <- unique(df_cfRNA$gene)
genes_cell <- unique(df_cell$gene)

cat("\ncfRNA gene count:", length(genes_cf), "; cell gene count:", length(genes_cell), "\n")

# Run enrichment separately
if (length(genes_cf) > 0) {
  run_enrichment(genes_cf, "cfRNA", "all", output_dir)
} else {
  cat("cfRNA file has no valid genes\n")
}

if (length(genes_cell) > 0) {
  run_enrichment(genes_cell, "cell", "all", output_dir)
} else {
  cat("cell file has no valid genes\n")
}

cat("\n========== All analyses completed ==========\n")
cat("All results (CSV + bubble plot PDF) saved in:", output_dir, "\n")