# Plot BP pathway enrichment dotplots for cfRNA and cell groups.
# Input files are expected in the cloud/workspace directory configured below.

required_pkgs <- c("readr", "dplyr", "ggplot2", "stringr", "forcats", "scales")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Please install required R packages first: ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(scales)
})


parse_gene_ratio <- function(x) {
  vapply(strsplit(x, "/", fixed = TRUE), function(z) {
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

# File path settings ---------------------------------------------------------
# Run this script from the project root. If your cloud workspace mounts the
# project at another root, set project_root to that root path.
project_root <- getwd()
base_dir <- "/data/work/01_2603cell_culture/07_figure/01_figure1/07_pathway_enrichment"

input_files <- c(
  cfRNA = file.path(base_dir, "cfRNA_enrichment.csv"),
  cell = file.path(base_dir, "cell_enrichment.csv")
)

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing input file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

output_dir <- file.path(base_dir, "plots")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use "RichFactor" for a pathway-specific enrichment intensity view.
# Change to "GeneRatio_num" if you prefer the conventional clusterProfiler axis.
x_metric <- "RichFactor"

cf_raw <- read_csv(input_files[["cfRNA"]], show_col_types = FALSE)
cell_raw <- read_csv(input_files[["cell"]], show_col_types = FALSE)

cell_modules <- data.frame(
  Description = c(
    "ERAD pathway",
    "response to endoplasmic reticulum stress",
    "protein folding",
    "Golgi vesicle transport",
    "post-Golgi vesicle-mediated transport",
    "glycoprotein metabolic process",
    "protein N-linked glycosylation",
    "mitochondrion organization",
    "oxidative phosphorylation",
    "replication fork processing"
  ),
  Module = c(
    "ER/ERAD",
    "ER/ERAD",
    "ER/ERAD",
    "Golgi transport",
    "Golgi transport",
    "Glycosylation",
    "Glycosylation",
    "Mitochondria",
    "Mitochondria",
    "DNA replication"
  ),
  stringsAsFactors = FALSE
)

cf_plot_data <- cf_raw %>%
  filter(Ontology == "BP") %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(
    Group = "cfRNA",
    pathway_order = row_number(),
    Module = if_else(
      str_detect(Description, regex("complement", ignore_case = TRUE)),
      "Complement signaling",
      "Chemotaxis/migration"
    ),
    GeneRatio_num = parse_gene_ratio(GeneRatio),
    neg_log10_padj = -log10(p.adjust)
  )

cell_plot_data <- cell_raw %>%
  filter(Ontology == "BP") %>%
  inner_join(cell_modules, by = "Description") %>%
  mutate(
    Group = "cell",
    pathway_order = match(Description, cell_modules$Description),
    GeneRatio_num = parse_gene_ratio(GeneRatio),
    neg_log10_padj = -log10(p.adjust)
  ) %>%
  arrange(pathway_order)

missing_cell_terms <- setdiff(cell_modules$Description, cell_plot_data$Description)
if (length(missing_cell_terms) > 0) {
  stop(
    "These selected cell pathways were not found in cell_enrichment.csv:\n",
    paste(missing_cell_terms, collapse = "\n"),
    call. = FALSE
  )
}

selected_pathways <- bind_rows(cf_plot_data, cell_plot_data) %>%
  select(
    Group, Module, ID, Description, GeneRatio, GeneRatio_num, RichFactor,
    FoldEnrichment, zScore, pvalue, p.adjust, qvalue, Count, geneID
  )

write_csv(selected_pathways, file.path(output_dir, "selected_pathways_for_plot.csv"))

module_cols <- c(
  "Chemotaxis/migration" = "#D55E00",
  "Complement signaling" = "#CC79A7",
  "ER/ERAD" = "#0072B2",
  "Golgi transport" = "#009E73",
  "Glycosylation" = "#E69F00",
  "Mitochondria" = "#56B4E9",
  "DNA replication" = "#6A3D9A"
)

x_axis_label <- if (x_metric == "GeneRatio_num") "Gene ratio" else "Rich factor"
get_x_limits <- function(x, lower_mult = 0.18, upper_mult = 0.12) {
  x <- x[is.finite(x)]
  x_range <- range(x)
  x_span <- diff(x_range)
  if (x_span == 0) {
    x_span <- max(abs(x_range[1]), 0.01)
  }
  c(
    max(0, x_range[1] - x_span * lower_mult),
    x_range[2] + x_span * upper_mult
  )
}

prep_single_plot_data <- function(df, order_by = c("p.adjust", "pathway_order")) {
  order_by <- match.arg(order_by)
  if (order_by == "p.adjust") {
    df <- df %>% arrange(desc(p.adjust))
  } else {
    df <- df %>% arrange(desc(pathway_order))
  }
  df %>%
    mutate(Term = fct_inorder(str_wrap(Description, width = 42)))
}

make_dotplot <- function(df, title, order_by = c("p.adjust", "pathway_order")) {
  df <- prep_single_plot_data(df, order_by = match.arg(order_by))
  x_limits <- get_x_limits(df[[x_metric]])

  ggplot(df, aes(x = .data[[x_metric]], y = Term)) +
    geom_point(aes(size = Count, color = Module), alpha = 0.95) +
    scale_color_manual(values = module_cols, drop = FALSE) +
    scale_size_continuous(range = c(3, 8), breaks = pretty_breaks(n = 4)) +
    scale_x_continuous(limits = x_limits, labels = number_format(accuracy = 0.01), expand = expansion(mult = c(0, 0))) +
    labs(
      title = title,
      x = x_axis_label,
      y = NULL,
      color = "Module",
      size = "Gene count"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 10, color = "black"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.height = grid::unit(0.45, "cm"),
      plot.margin = margin(5.5, 12, 5.5, 5.5)
    )
}

cf_plot <- make_dotplot(
  cf_plot_data,
  title = "cfRNA enriched BP pathways",
  order_by = "p.adjust"
)

cell_plot <- make_dotplot(
  cell_plot_data,
  title = "Cell enriched BP pathways",
  order_by = "pathway_order"
)

combined_plot_data <- bind_rows(cf_plot_data, cell_plot_data) %>%
  group_by(Group) %>%
  arrange(desc(pathway_order), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    Term_unique = paste(Group, str_wrap(Description, width = 42), sep = "___"),
    Term_unique = factor(Term_unique, levels = unique(Term_unique))
  )

combined_x_limits <- get_x_limits(combined_plot_data[[x_metric]])

combined_plot <- ggplot(
  combined_plot_data,
  aes(x = .data[[x_metric]], y = Term_unique)
) +
  geom_point(aes(size = Count, color = Module), alpha = 0.95) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_color_manual(values = module_cols, drop = FALSE) +
  scale_size_continuous(range = c(3, 8), breaks = pretty_breaks(n = 4)) +
  scale_x_continuous(limits = combined_x_limits, labels = number_format(accuracy = 0.01), expand = expansion(mult = c(0, 0))) +
  labs(
    x = x_axis_label,
    y = NULL,
    color = "Module",
    size = "Gene count"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92", color = "grey70", linewidth = 0.3),
    strip.text.y = element_text(face = "bold", size = 10, angle = 0),
    axis.text.y = element_text(size = 8.5, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.height = grid::unit(0.45, "cm"),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  )

pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"

ggsave(
  filename = file.path(output_dir, paste0("cfRNA_BP_top10_", x_metric, "_dotplot.pdf")),
  plot = cf_plot,
  width = 6.2,
  height = 4.2,
  device = pdf_device
)
ggsave(
  filename = file.path(output_dir, paste0("cfRNA_BP_top10_", x_metric, "_dotplot.png")),
  plot = cf_plot,
  width = 6.2,
  height = 4.2,
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, paste0("cell_BP_selected10_", x_metric, "_dotplot.pdf")),
  plot = cell_plot,
  width = 6.8,
  height = 4.5,
  device = pdf_device
)
ggsave(
  filename = file.path(output_dir, paste0("cell_BP_selected10_", x_metric, "_dotplot.png")),
  plot = cell_plot,
  width = 6.8,
  height = 4.5,
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, paste0("combined_BP_top10_", x_metric, "_dotplot.pdf")),
  plot = combined_plot,
  width = 7.2,
  height = 7.2,
  device = pdf_device
)
ggsave(
  filename = file.path(output_dir, paste0("combined_BP_top10_", x_metric, "_dotplot.png")),
  plot = combined_plot,
  width = 7.2,
  height = 7.2,
  dpi = 600
)

message("Done. Plots and selected pathway table were saved to: ", output_dir)









