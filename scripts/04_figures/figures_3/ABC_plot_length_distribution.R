###############################################################################
# Housekeeping gene fragment length distribution plotting script
# Goal: plot length distributions of RACK1, EEF1A1, HNRNPA1, RPL23A, RPL31, GAPDH
#       across Components (cfRNA, debris, cell, debris+cell)
# Output: one figure per Component, different Cell in different colors, with error ribbon
# Optimization: filter by sample + gene before melting to greatly reduce memory usage
###############################################################################

# ======================== 0. Load dependencies ========================
library(data.table)
library(ggplot2)
library(scales)

# ======================== 1. Read data ========================
message("[1/6] Reading input files ...")

qc_raw     <- fread("./cellcult_Sample_Information.txt")
length_raw <- fread("./all_samples_length_merged.mlncRNA.txt")

# Clean column names: remove possible UTF-8 BOM and leading/trailing whitespace
setnames(qc_raw, trimws(sub("\ufeff", "", names(qc_raw), fixed = TRUE)))

# Length table columns: first 4 columns are metadata, the rest are length positions 1-100
# If fread did not recognize the header (column names are V1, V2, ...), assign manually
n_len_cols <- ncol(length_raw) - 4L
length_col_names <- c("sample", "Transcript", "Gene", "Type", as.character(seq_len(n_len_cols)))
setnames(length_raw, seq_along(length_raw), length_col_names)

message(sprintf("  Raw length table rows: %d, cols: %d", nrow(length_raw), ncol(length_raw)))
message(sprintf("  Length table columns: %s", paste(head(names(length_raw), 6), collapse = ", ")))

# ======================== 2. Pre-filter by samples + genes ========================
message("[2/6] Pre-filtering by QC samples and housekeeping genes ...")

# Keep only QC == 1 samples
qc_pass <- qc_raw[QC == 1, .(Sample, Cell, Component)]
message(sprintf("  QC-passed samples: %d", nrow(qc_pass)))

# Housekeeping gene list
hk_genes <- c("SPN", "DRAXIN", "NUPR1", "WDR74", "MIR99AHG", "TPT1-AS1")
hk_dt <- data.table(Gene = hk_genes)

# First filter samples via data.table join (avoid %in% type ambiguity on large tables)
length_filt <- length_raw[qc_pass[, .(Sample)], on = .(sample = Sample), nomatch = 0]
# Then filter genes via %chin% (%chin% is data.table's fast %in% optimized for character columns)
length_filt <- length_filt[Gene %chin% hk_genes]

message(sprintf("  Length table rows after pre-filtering: %d (was %d)",
                nrow(length_filt), nrow(length_raw)))

# Check gene coverage
found_genes <- unique(length_filt$Gene)
missing_genes <- setdiff(hk_genes, found_genes)
if (length(missing_genes) > 0) {
  warning("Genes NOT found in data: ", paste(missing_genes, collapse = ", "))
}
message(sprintf("  Genes retained: %s", paste(found_genes, collapse = ", ")))

# ======================== 3. Merge sample metadata ========================
message("[3/6] Merging sample metadata ...")

# Append Cell and Component information via data.table join
length_filt <- qc_pass[length_filt, on = .(Sample = sample), nomatch = 0]
setnames(length_filt, "Sample", "sample")
length_filt[, c("Transcript", "Type") := NULL]  # Transcript information no longer needed

# Cell factor levels (control legend order)
length_filt[, Cell := factor(Cell, levels = c("Hep3B2.1-7", "HepG2", "K562", "HTR-8/SVneo", "HEK293T"))]

# ======================== 4. Reshape and aggregate by gene ========================
message("[4/6] Melting to long format and aggregating ...")

length_cols <- grep("^[0-9]+$", names(length_filt), value = TRUE)
message(sprintf("  Length positions: %d (range %s - %s)",
                length(length_cols), length_cols[1], length_cols[length(length_cols)]))

# Wide -> long (data volume already greatly reduced at this point)
length_long <- melt(length_filt,
                    id.vars       = c("sample", "Gene", "Cell", "Component"),
                    measure.vars  = length_cols,
                    variable.name = "Length",
                    value.name    = "Count")
length_long[, Length := as.integer(as.character(Length))]

# Aggregate all transcripts by sample x Gene, summing counts at each length position
# (since genes were already filtered at the wide stage and each row is one transcript, summing by sample x Gene suffices)
gene_agg <- length_long[, .(Count = sum(Count, na.rm = TRUE)),
                        by = .(sample, Gene, Cell, Component, Length)]

# Compute within-sample x Gene proportion
gene_agg[, Total := sum(Count), by = .(sample, Gene)]
gene_agg[, Proportion := ifelse(Total > 0, Count / Total, 0)]

message(sprintf("  Aggregated rows: %d", nrow(gene_agg)))
rm(length_raw, length_filt, length_long)  # Release intermediate variables
gc()

# ======================== 4.5. Create cell+debris pooled group ========================
message("[4.5/6] Creating cell+debris pooled group ...")

# Merge cell and debris samples' proportions, marking Component uniformly as cell+debris
cd_pooled <- gene_agg[Component %in% c("cell", "debris")]
cd_pooled[, Component := "cell+debris"]
# Do not re-aggregate Count, reuse the already computed Proportion (downstream computes Mean/SE by Cell x Gene x Component)
gene_agg <- rbind(gene_agg, cd_pooled)

message(sprintf("  Aggregated rows (with cell+debris): %d", nrow(gene_agg)))

# ======================== 5. Plotting ========================
message("[5/6] Generating plots ...")

# Theme settings (academic journal style, no background grid lines)
theme_academic <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    strip.background   = element_rect(fill = "grey95", colour = "grey80"),
    strip.text         = element_text(size = 11, face = "bold"),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    axis.text.y        = element_text(size = 9),
    axis.title         = element_text(size = 12),
    legend.position    = "bottom",
    legend.title       = element_text(size = 10),
    legend.text        = element_text(size = 9),
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.margin        = margin(10, 15, 10, 10)
  )

# Pre-compute overall summary statistics (Cell x Gene x Component x Length)
message("  Pre-computing summary statistics ...")
plot_dt_all <- gene_agg[, .(
  Mean = mean(Proportion, na.rm = TRUE),
  SE   = sd(Proportion, na.rm = TRUE) / sqrt(.N),
  N    = .N
), by = .(Cell, Gene, Component, Length)]
plot_dt_all[, `:=`(
  ymin = Mean - SE,
  ymax = Mean + SE
)]
plot_dt_all[N <= 1L, `:=`(ymin = NA_real_, ymax = NA_real_)]

# Plot by gene x component in order (one gene one component per figure)
for (gene in hk_genes) {
  for (comp in sort(unique(plot_dt_all$Component))) {

    sub_data <- plot_dt_all[Gene == gene & Component == comp]

    if (nrow(sub_data) == 0) {
      message(sprintf("  Skipping '%s' x '%s': no data", gene, comp))
      next
    }

    # Dynamically determine x-axis range
    x_min <- sub_data[Mean > 0, min(Length, na.rm = TRUE)]
    x_max <- sub_data[Mean > 0, max(Length, na.rm = TRUE)]
    x_range <- seq(x_min, x_max, by = 1)

    has_se <- sub_data[, any(!is.na(ymin))]

    message(sprintf("  %s | %s: x-axis %d - %d bp, SE: %s",
                    gene, comp, x_min, x_max, has_se))

    p <- ggplot(sub_data,
                aes(x     = Length,
                    y     = Mean,
                    color = Cell,
                    fill  = Cell,
                    group = Cell))

    if (has_se) {
      p <- p + geom_ribbon(data = sub_data[N > 1L],
                           aes(ymin = ymin, ymax = ymax),
                           alpha    = 0.15,
                           colour   = NA)
    }

    p <- p +
      geom_line(linewidth = 0.8, na.rm = TRUE) +

      scale_x_continuous(breaks = x_range,
                         limits = c(min(x_range) - 0.5, max(x_range) + 0.5),
                         expand = c(0, 0)) +

      scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                         labels = label_number(accuracy = 0.001)) +

      labs(x     = "Fragment Length (bp)",
           y     = "Proportion",
           title = sprintf("%s - %s", gene, comp)) +

      theme_academic

    out_name <- sprintf("length_distribution_%s_%s.pdf",
                        gene, gsub("[+ ]", "_", comp))
    ggsave(out_name, plot = p, width = 8, height = 6, device = "pdf")
    message(sprintf("  Saved: %s", out_name))
  }
}

message("[6/6] Done.")