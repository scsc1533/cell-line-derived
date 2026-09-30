#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

# =============================================================================
# length_gdm_lowdim_interpretable_model_6dim_with_independent_validation.R
#
# Goal:
#   Build a GDM vs Healthy classification model based on direction-split
#   low-dimensional interpretable features, and perform external validation
#   in an independent validation set W06.
#
# Training set:
#   SH:
#     /data/work/01_2603cell_culture/07_figure/04_figure4/02_GDM_validation_SH/all
#
# Independent validation set:
#   W06:
#     /data/work/01_2603cell_culture/07_figure/04_figure4/02_GDM_validation_W06
#
# Key principles:
#   1. gene_level_group_tests.tsv is used only in training set SH to select direction genes.
#   2. W06 gene_level_group_tests.tsv is only used for optional direction concordance check,
#      not involved in model feature selection.
#   3. W06 only uses the same batch of genes selected from SH to construct the same low-dim features.
#   4. Model is fit only on SH, then directly predicts W06.
#
# Feature definition:
#   For each target reference, split genes into two classes by direction:
#
#   GDM_closer:
#     gdm_minus_healthy_D_target < 0
#     i.e., D_target is lower in GDM samples, closer to the target cell line reference distribution
#
#   GDM_farther:
#     gdm_minus_healthy_D_target > 0
#     i.e., D_target is higher in GDM samples, farther from the target cell line reference distribution
#
#   At most 6 direction features:
#     1. HTR8_GDM_closer
#     2. HTR8_GDM_farther
#     3. Liver_GDM_closer
#     4. Liver_GDM_farther
#     5. K562_GDM_closer
#     6. K562_GDM_farther
#
# Model:
#   - logistic regression
#   - SH internal evaluation: 5-fold stratified CV, repeated 20 times
#   - W06 external evaluation: independent validation
# =============================================================================


# =============================================================================
# 0. Configuration section
# =============================================================================

training_dataset_label <- "SH"
independent_validation_label <- "W06"

training_validation_dir <- "./02_GDM_validation_SH/all"

independent_validation_dir <- "./02_GDM_validation_W06mz"

output_dir <- "./GDM_SH_with_W06_validation_top3_mz"

training_gene_test_candidates <- c(
  file.path(training_validation_dir, "gene_level_group_tests.tsv"),
  file.path(training_validation_dir, "gene_level_group_tests.tsv.gz")
)

training_gene_distance_candidates <- c(
  file.path(training_validation_dir, "gene_level_distances.tsv"),
  file.path(training_validation_dir, "gene_level_distances.tsv.gz")
)

independent_gene_test_candidates <- c(
  file.path(independent_validation_dir, "gene_level_group_tests.tsv"),
  file.path(independent_validation_dir, "gene_level_group_tests.tsv.gz"),
  file.path(independent_validation_dir, "all", "gene_level_group_tests.tsv"),
  file.path(independent_validation_dir, "all", "gene_level_group_tests.tsv.gz")
)

independent_gene_distance_candidates <- c(
  file.path(independent_validation_dir, "gene_level_distances.tsv"),
  file.path(independent_validation_dir, "gene_level_distances.tsv.gz"),
  file.path(independent_validation_dir, "all", "gene_level_distances.tsv"),
  file.path(independent_validation_dir, "all", "gene_level_distances.tsv.gz")
)

comparison_feature_map <- data.table(
  comparison = c(
    "HTR_8_SVneo_vs_rest",
    "HepG2_Hep3B2.1_7_vs_rest",
    "K562_vs_rest"
  ),
  feature_prefix = c(
    "HTR8",
    "Liver",
    "K562"
  ),
  target_label = c(
    "HTR-8/SVneo",
    "Liver",
    "K562"
  )
)

feature_direction_levels <- c(
  "GDM_closer",
  "GDM_farther"
)

p_cutoff <- 0.05
n_folds <- 5
n_repeats <- 20
random_seed <- 123
classification_threshold <- 0.5


# =============================================================================
# 1. General utility functions
# =============================================================================

resolve_first_existing <- function(paths) {
  
  hit <- paths[file.exists(paths)][1]
  
  if (length(hit) == 0 || is.na(hit)) {
    stop(sprintf(
      "Input file not found: %s",
      paste(paths, collapse = " / ")
    ))
  }
  
  hit
}


make_feature_order <- function(
    comparison_feature_map,
    feature_direction_levels
) {
  
  feature_order_dt <- rbindlist(
    lapply(seq_len(nrow(comparison_feature_map)), function(i) {
      
      data.table(
        comparison = comparison_feature_map$comparison[i],
        feature_prefix = comparison_feature_map$feature_prefix[i],
        target_label = comparison_feature_map$target_label[i],
        feature_direction = feature_direction_levels
      )
    })
  )
  
  feature_order_dt[
    ,
    feature_name_median := paste0(
      "median_D_target_",
      feature_prefix,
      "_",
      feature_direction
    )
  ]
  
  feature_order_dt[
    ,
    feature_name_mean := sub(
      "^median_",
      "mean_",
      feature_name_median
    )
  ]
  
  feature_order_dt
}


make_pretty_feature_label <- function(x) {
  
  y <- x
  
  y <- sub("^median_D_target_", "", y)
  y <- sub("^mean_D_target_", "", y)
  
  y <- sub("^HTR8_", "HTR-8/SVneo\n", y)
  y <- sub("^Liver_", "Liver\n", y)
  y <- sub("^K562_", "K562\n", y)
  
  y <- sub("GDM_closer$", "GDM closer", y)
  y <- sub("GDM_farther$", "GDM farther", y)
  
  y
}


extract_feature_annotation <- function(features) {
  
  x <- data.table(
    feature = features
  )
  
  x[
    ,
    feature_core := sub(
      "^median_D_target_",
      "",
      feature
    )
  ]
  
  x[
    ,
    feature_core := sub(
      "^mean_D_target_",
      "",
      feature_core
    )
  ]
  
  x[
    grepl("^HTR8_", feature_core),
    Target := "HTR-8/SVneo"
  ]
  
  x[
    grepl("^Liver_", feature_core),
    Target := "Liver"
  ]
  
  x[
    grepl("^K562_", feature_core),
    Target := "K562"
  ]
  
  x[
    grepl("GDM_closer$", feature_core),
    Direction := "GDM closer genes"
  ]
  
  x[
    grepl("GDM_farther$", feature_core),
    Direction := "GDM farther genes"
  ]
  
  x[
    ,
    feature_label := make_pretty_feature_label(feature)
  ]
  
  x
}


compute_auc <- function(labels, scores) {
  
  keep <- !is.na(labels) & !is.na(scores)
  labels <- labels[keep]
  scores <- scores[keep]
  
  if (length(unique(labels)) < 2) {
    return(NA_real_)
  }
  
  pos <- scores[labels == 1]
  neg <- scores[labels == 0]
  
  if (length(pos) == 0 || length(neg) == 0) {
    return(NA_real_)
  }
  
  ranks <- rank(
    c(pos, neg),
    ties.method = "average"
  )
  
  n_pos <- length(pos)
  n_neg <- length(neg)
  
  auc <- (
    sum(ranks[seq_len(n_pos)]) -
      n_pos * (n_pos + 1) / 2
  ) / (n_pos * n_neg)
  
  auc
}


compute_metrics <- function(
    labels,
    probs,
    threshold = classification_threshold
) {
  
  keep <- !is.na(labels) & !is.na(probs)
  labels <- labels[keep]
  probs <- probs[keep]
  
  if (length(labels) == 0) {
    return(data.table(
      auc = NA_real_,
      accuracy = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_,
      balanced_accuracy = NA_real_
    ))
  }
  
  pred <- as.integer(probs >= threshold)
  
  tp <- sum(pred == 1 & labels == 1)
  tn <- sum(pred == 0 & labels == 0)
  fp <- sum(pred == 1 & labels == 0)
  fn <- sum(pred == 0 & labels == 1)
  
  sensitivity <- if ((tp + fn) > 0) {
    tp / (tp + fn)
  } else {
    NA_real_
  }
  
  specificity <- if ((tn + fp) > 0) {
    tn / (tn + fp)
  } else {
    NA_real_
  }
  
  accuracy <- mean(pred == labels)
  
  balanced_accuracy <- mean(
    c(sensitivity, specificity),
    na.rm = TRUE
  )
  
  auc <- compute_auc(labels, probs)
  
  data.table(
    auc = auc,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = balanced_accuracy
  )
}


make_stratified_folds <- function(y, k, seed) {
  
  set.seed(seed)
  
  folds <- integer(length(y))
  
  for (cls in sort(unique(y))) {
    
    idx <- which(y == cls)
    idx <- sample(idx)
    
    fold_id <- rep(
      seq_len(k),
      length.out = length(idx)
    )
    
    folds[idx] <- fold_id
  }
  
  folds
}


build_roc_curve <- function(labels, scores) {
  
  keep <- !is.na(labels) & !is.na(scores)
  labels <- labels[keep]
  scores <- scores[keep]
  
  if (length(unique(labels)) < 2) {
    return(data.table(
      threshold = c(Inf, -Inf),
      fpr = c(0, 1),
      tpr = c(0, 1)
    ))
  }
  
  thresholds <- sort(
    unique(scores),
    decreasing = TRUE
  )
  
  roc_list <- lapply(
    thresholds,
    function(thr) {
      
      pred <- as.integer(scores >= thr)
      
      tp <- sum(pred == 1 & labels == 1)
      tn <- sum(pred == 0 & labels == 0)
      fp <- sum(pred == 1 & labels == 0)
      fn <- sum(pred == 0 & labels == 1)
      
      tpr <- if ((tp + fn) > 0) {
        tp / (tp + fn)
      } else {
        0
      }
      
      fpr <- if ((fp + tn) > 0) {
        fp / (fp + tn)
      } else {
        0
      }
      
      data.table(
        threshold = thr,
        fpr = fpr,
        tpr = tpr
      )
    }
  )
  
  roc_dt <- rbindlist(roc_list)
  
  roc_dt <- unique(rbind(
    data.table(
      threshold = Inf,
      fpr = 0,
      tpr = 0
    ),
    roc_dt,
    data.table(
      threshold = -Inf,
      fpr = 1,
      tpr = 1
    )
  ))
  
  setorder(
    roc_dt,
    fpr,
    tpr
  )
  
  roc_dt
}


plot_roc_curve <- function(
    roc_dt,
    auc_value,
    title,
    output_pdf,
    output_png
) {
  
  p <- ggplot(
    roc_dt,
    aes(x = fpr, y = tpr)
  ) +
    geom_line(
      color = "#D55E00",
      linewidth = 1
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    coord_equal() +
    theme_bw(base_size = 12) +
    labs(
      title = sprintf(
        "%s AUC = %.3f",
        title,
        auc_value
      ),
      x = "False Positive Rate",
      y = "True Positive Rate"
    )
  
  ggsave(
    output_pdf,
    p,
    width = 5.5,
    height = 5
  )
  
  ggsave(
    output_png,
    p,
    width = 5.5,
    height = 5,
    dpi = 300
  )
  
  p
}


# =============================================================================
# 2. Gene selection and feature construction functions
# =============================================================================

select_direction_genes_from_training <- function(
    gene_test_dt,
    comparison_feature_map,
    feature_order_dt,
    p_cutoff,
    output_dir
) {
  
  required_gene_test_cols <- c(
    "comparison",
    "Gene",
    "wilcox_p_D_target",
    "gdm_minus_healthy_D_target"
  )
  
  missing_gene_test_cols <- setdiff(
    required_gene_test_cols,
    names(gene_test_dt)
  )
  
  if (length(missing_gene_test_cols) > 0) {
    stop(sprintf(
      "Training set gene_level_group_tests is missing columns: %s",
      paste(missing_gene_test_cols, collapse = ", ")
    ))
  }
  
  selected_gene_dt <- gene_test_dt[
    comparison %in% comparison_feature_map$comparison &
      !is.na(wilcox_p_D_target) &
      wilcox_p_D_target < p_cutoff &
      !is.na(gdm_minus_healthy_D_target) &
      gdm_minus_healthy_D_target != 0,
    .(
      comparison,
      Gene,
      wilcox_p_D_target,
      gdm_minus_healthy_D_target
    )
  ]
  
  if (nrow(selected_gene_dt) == 0) {
    stop(
      "No genes in training set satisfy wilcox_p_D_target < 0.05 with a clear direction."
    )
  }
  
  selected_gene_dt[
    ,
    feature_direction := fifelse(
      gdm_minus_healthy_D_target < 0,
      "GDM_closer",
      "GDM_farther"
    )
  ]
  
  selected_gene_dt <- merge(
    selected_gene_dt,
    comparison_feature_map,
    by = "comparison",
    all.x = TRUE
  )
  
  selected_gene_dt[
    ,
    feature_name := paste0(
      "median_D_target_",
      feature_prefix,
      "_",
      feature_direction
    )
  ]
  
  selected_gene_dt[
    ,
    feature_name_mean := sub(
      "^median_",
      "mean_",
      feature_name
    )
  ]
  
  available_features_median <- unique(
    selected_gene_dt$feature_name
  )
  
  feature_order_dt_available <- feature_order_dt[
    feature_name_median %in% available_features_median
  ]
  
  missing_feature_gene_sets <- setdiff(
    feature_order_dt$feature_name_median,
    available_features_median
  )
  
  if (length(missing_feature_gene_sets) > 0) {
    message(
      ">> The following directions have no significant genes selected in the training set and will be skipped automatically: ",
      paste(missing_feature_gene_sets, collapse = ", ")
    )
  }
  
  if (nrow(feature_order_dt_available) == 0) {
    stop("No available significant genes in any direction, cannot build model.")
  }
  
  needed_features_median <- feature_order_dt_available$feature_name_median
  needed_features_mean <- feature_order_dt_available$feature_name_mean
  
  gene_set_summary <- selected_gene_dt[
    ,
    .(
      n_genes = .N,
      median_effect = median(
        gdm_minus_healthy_D_target,
        na.rm = TRUE
      ),
      min_p = min(
        wilcox_p_D_target,
        na.rm = TRUE
      ),
      max_p = max(
        wilcox_p_D_target,
        na.rm = TRUE
      )
    ),
    by = .(
      comparison,
      target_label,
      feature_direction,
      feature_name
    )
  ]
  
  gene_set_summary[
    ,
    feature_order_tmp := match(
      feature_name,
      needed_features_median
    )
  ]
  
  setorder(
    gene_set_summary,
    feature_order_tmp
  )
  
  gene_set_summary[
    ,
    feature_order_tmp := NULL
  ]
  
  fwrite(
    feature_order_dt_available,
    file.path(
      output_dir,
      "six_feature_definition_available_training_SH.tsv"
    ),
    sep = "\t"
  )
  
  fwrite(
    selected_gene_dt,
    file.path(
      output_dir,
      "selected_genes_for_direction_lowdim_model_training_SH.tsv"
    ),
    sep = "\t"
  )
  
  fwrite(
    gene_set_summary,
    file.path(
      output_dir,
      "gene_set_summary_direction_features_training_SH.tsv"
    ),
    sep = "\t"
  )
  
  list(
    selected_gene_dt = selected_gene_dt,
    gene_set_summary = gene_set_summary,
    feature_order_dt_available = feature_order_dt_available,
    needed_features_median = needed_features_median,
    needed_features_mean = needed_features_mean
  )
}


build_lowdim_feature_matrices <- function(
    gene_distance_dt,
    selected_gene_dt,
    needed_features_median,
    needed_features_mean,
    dataset_label,
    strict_features = TRUE
) {
  
  gene_distance_dt <- copy(gene_distance_dt)
  
  required_distance_cols <- c(
    "comparison",
    "sample",
    "Gene",
    "group",
    "D_target"
  )
  
  missing_distance_cols <- setdiff(
    required_distance_cols,
    names(gene_distance_dt)
  )
  
  if (length(missing_distance_cols) > 0) {
    stop(sprintf(
      "%s gene_level_distances is missing columns: %s",
      dataset_label,
      paste(missing_distance_cols, collapse = ", ")
    ))
  }
  
  gene_distance_dt[
    ,
    group := as.integer(group)
  ]
  
  if (!"group_label" %in% names(gene_distance_dt)) {
    gene_distance_dt[
      ,
      group_label := fifelse(
        group == 0,
        "Healthy",
        "GDM"
      )
    ]
  } else {
    gene_distance_dt[
      ,
      group_label := fifelse(
        group == 0,
        "Healthy",
        "GDM"
      )
    ]
  }
  
  feature_long_dt <- merge(
    gene_distance_dt[
      ,
      .(
        comparison,
        sample,
        Gene,
        group,
        group_label,
        D_target
      )
    ],
    selected_gene_dt[
      ,
      .(
        comparison,
        Gene,
        feature_name
      )
    ],
    by = c(
      "comparison",
      "Gene"
    ),
    all = FALSE
  )
  
  if (nrow(feature_long_dt) == 0) {
    stop(sprintf(
      "%s: the selected training genes were not matched to records in gene_level_distances.",
      dataset_label
    ))
  }
  
  feature_long_dt <- feature_long_dt[
    ,
    .(
      group = first(group),
      group_label = first(group_label),
      n_genes_used = sum(!is.na(D_target)),
      median_D_target = median(
        D_target,
        na.rm = TRUE
      ),
      mean_D_target = mean(
        D_target,
        na.rm = TRUE
      )
    ),
    by = .(
      sample,
      feature_name
    )
  ]
  
  feature_wide_median_dt <- dcast(
    feature_long_dt,
    sample + group + group_label ~ feature_name,
    value.var = "median_D_target"
  )
  
  feature_wide_mean_dt <- dcast(
    feature_long_dt,
    sample + group + group_label ~ feature_name,
    value.var = "mean_D_target"
  )
  
  mean_feature_cols_old <- setdiff(
    names(feature_wide_mean_dt),
    c(
      "sample",
      "group",
      "group_label"
    )
  )
  
  setnames(
    feature_wide_mean_dt,
    old = mean_feature_cols_old,
    new = sub(
      "^median_",
      "mean_",
      mean_feature_cols_old
    )
  )
  
  missing_median_cols <- setdiff(
    needed_features_median,
    names(feature_wide_median_dt)
  )
  
  missing_mean_cols <- setdiff(
    needed_features_mean,
    names(feature_wide_mean_dt)
  )
  
  if (strict_features && length(missing_median_cols) > 0) {
    stop(sprintf(
      "%s is missing median feature columns: %s",
      dataset_label,
      paste(missing_median_cols, collapse = ", ")
    ))
  }
  
  if (strict_features && length(missing_mean_cols) > 0) {
    stop(sprintf(
      "%s is missing mean feature columns: %s",
      dataset_label,
      paste(missing_mean_cols, collapse = ", ")
    ))
  }
  
  if (!strict_features) {
    
    if (length(missing_median_cols) > 0) {
      message(
        ">> ",
        dataset_label,
        " the following median features were not generated and will be skipped: ",
        paste(missing_median_cols, collapse = ", ")
      )
      
      needed_features_median <- setdiff(
        needed_features_median,
        missing_median_cols
      )
    }
    
    if (length(missing_mean_cols) > 0) {
      message(
        ">> ",
        dataset_label,
        " the following mean features were not generated and will be skipped: ",
        paste(missing_mean_cols, collapse = ", ")
      )
      
      needed_features_mean <- setdiff(
        needed_features_mean,
        missing_mean_cols
      )
    }
  }
  
  feature_count_wide <- dcast(
    feature_long_dt,
    sample + group + group_label ~ feature_name,
    value.var = "n_genes_used"
  )
  
  count_feature_cols_old <- setdiff(
    names(feature_count_wide),
    c(
      "sample",
      "group",
      "group_label"
    )
  )
  
  setnames(
    feature_count_wide,
    old = count_feature_cols_old,
    new = paste0(
      "n_genes_",
      count_feature_cols_old
    )
  )
  
  feature_wide_median_dt <- merge(
    feature_wide_median_dt,
    feature_count_wide,
    by = c(
      "sample",
      "group",
      "group_label"
    ),
    all.x = TRUE
  )
  
  feature_wide_mean_dt <- merge(
    feature_wide_mean_dt,
    feature_count_wide,
    by = c(
      "sample",
      "group",
      "group_label"
    ),
    all.x = TRUE
  )
  
  model_dt_median <- feature_wide_median_dt[
    complete.cases(
      feature_wide_median_dt[
        ,
        ..needed_features_median
      ]
    )
  ]
  
  model_dt_mean <- feature_wide_mean_dt[
    complete.cases(
      feature_wide_mean_dt[
        ,
        ..needed_features_mean
      ]
    )
  ]
  
  if (nrow(model_dt_median) == 0) {
    stop(sprintf(
      "%s has no samples with complete median low-dim features.",
      dataset_label
    ))
  }
  
  if (nrow(model_dt_mean) == 0) {
    stop(sprintf(
      "%s has no samples with complete mean low-dim features.",
      dataset_label
    ))
  }
  
  model_dt_median[
    ,
    y := as.integer(group == 1)
  ]
  
  model_dt_mean[
    ,
    y := as.integer(group == 1)
  ]
  
  list(
    feature_long_dt = feature_long_dt,
    model_dt_median = model_dt_median,
    model_dt_mean = model_dt_mean,
    needed_features_median = needed_features_median,
    needed_features_mean = needed_features_mean
  )
}


check_independent_gene_test_concordance <- function(
    independent_gene_test_dt,
    selected_gene_dt,
    comparison_feature_map,
    output_dir,
    validation_label
) {
  
  required_gene_test_cols <- c(
    "comparison",
    "Gene",
    "wilcox_p_D_target",
    "gdm_minus_healthy_D_target"
  )
  
  missing_cols <- setdiff(
    required_gene_test_cols,
    names(independent_gene_test_dt)
  )
  
  if (length(missing_cols) > 0) {
    message(
      ">> Independent validation set gene_level_group_tests is missing required columns, skipping direction concordance check: ",
      paste(missing_cols, collapse = ", ")
    )
    return(NULL)
  }
  
  validation_direction_dt <- independent_gene_test_dt[
    comparison %in% comparison_feature_map$comparison,
    .(
      comparison,
      Gene,
      wilcox_p_D_target,
      gdm_minus_healthy_D_target
    )
  ]
  
  validation_direction_dt[
    ,
    feature_direction_validation := fifelse(
      is.na(gdm_minus_healthy_D_target) |
        gdm_minus_healthy_D_target == 0,
      NA_character_,
      fifelse(
        gdm_minus_healthy_D_target < 0,
        "GDM_closer",
        "GDM_farther"
      )
    )
  ]
  
  training_direction_dt <- selected_gene_dt[
    ,
    .(
      comparison,
      Gene,
      feature_name,
      target_label,
      feature_direction_training = feature_direction,
      training_wilcox_p_D_target = wilcox_p_D_target,
      training_gdm_minus_healthy_D_target = gdm_minus_healthy_D_target
    )
  ]
  
  concordance_dt <- merge(
    training_direction_dt,
    validation_direction_dt,
    by = c(
      "comparison",
      "Gene"
    ),
    all.x = TRUE
  )
  
  setnames(
    concordance_dt,
    old = c(
      "wilcox_p_D_target",
      "gdm_minus_healthy_D_target"
    ),
    new = c(
      "validation_wilcox_p_D_target",
      "validation_gdm_minus_healthy_D_target"
    )
  )
  
  concordance_dt[
    ,
    independent_has_gene_test := !is.na(validation_wilcox_p_D_target)
  ]
  
  concordance_dt[
    ,
    direction_concordant := fifelse(
      independent_has_gene_test &
        !is.na(feature_direction_validation),
      feature_direction_training == feature_direction_validation,
      NA
    )
  ]
  
  concordance_summary_dt <- concordance_dt[
    ,
    .(
      n_training_genes = .N,
      n_genes_with_validation_test = sum(
        independent_has_gene_test,
        na.rm = TRUE
      ),
      n_direction_concordant = sum(
        direction_concordant == TRUE,
        na.rm = TRUE
      ),
      n_direction_discordant = sum(
        direction_concordant == FALSE,
        na.rm = TRUE
      ),
      direction_concordance_rate = ifelse(
        sum(!is.na(direction_concordant)) > 0,
        mean(direction_concordant, na.rm = TRUE),
        NA_real_
      )
    ),
    by = .(
      feature_name,
      target_label,
      feature_direction_training
    )
  ]
  
  fwrite(
    concordance_dt,
    file.path(
      output_dir,
      paste0(
        "independent_validation_",
        validation_label,
        "_gene_test_direction_concordance_detail.tsv"
      )
    ),
    sep = "\t"
  )
  
  fwrite(
    concordance_summary_dt,
    file.path(
      output_dir,
      paste0(
        "independent_validation_",
        validation_label,
        "_gene_test_direction_concordance_summary.tsv"
      )
    ),
    sep = "\t"
  )
  
  list(
    concordance_detail = concordance_dt,
    concordance_summary = concordance_summary_dt
  )
}


# =============================================================================
# 3. Modeling functions
# =============================================================================

run_lowdim_model <- function(
    model_dt,
    needed_features,
    version_label,
    output_dir,
    dataset_label = "SH"
) {
  
  model_dt <- copy(model_dt)
  
  if (length(needed_features) == 0) {
    stop("needed_features is empty, cannot build model.")
  }
  
  formula_str <- paste(
    "y ~",
    paste(needed_features, collapse = " + ")
  )
  
  fit_formula <- as.formula(formula_str)
  
  full_fit <- glm(
    fit_formula,
    data = model_dt,
    family = binomial(),
    control = glm.control(maxit = 100)
  )
  
  coef_mat <- summary(full_fit)$coefficients
  
  coef_dt <- data.table(
    feature = rownames(coef_mat),
    coefficient = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    z_value = coef_mat[, "z value"],
    p_value = coef_mat[, "Pr(>|z|)"]
  )
  
  coef_dt[
    ,
    odds_ratio := exp(coefficient)
  ]
  
  pred_col <- paste0(
    "pred_prob_full_",
    dataset_label,
    "_",
    version_label
  )
  
  model_dt[
    ,
    (pred_col) := as.numeric(
      predict(
        full_fit,
        newdata = model_dt,
        type = "response"
      )
    )
  ]
  
  full_metrics <- compute_metrics(
    model_dt$y,
    model_dt[[pred_col]]
  )
  
  full_metrics[
    ,
    `:=`(
      dataset = dataset_label,
      model = paste0(
        "full_fit_apparent_",
        dataset_label,
        "_",
        version_label
      ),
      n_samples_total = nrow(model_dt),
      n_healthy = sum(model_dt$group == 0),
      n_GDM = sum(model_dt$group == 1)
    )
  ]
  
  cv_pred_list <- vector(
    "list",
    n_repeats
  )
  
  cv_metrics_list <- vector(
    "list",
    n_repeats
  )
  
  for (r in seq_len(n_repeats)) {
    
    folds <- make_stratified_folds(
      model_dt$y,
      k = n_folds,
      seed = random_seed + r
    )
    
    pred_dt <- copy(
      model_dt[
        ,
        .(
          sample,
          group,
          group_label,
          y
        )
      ]
    )
    
    pred_dt[
      ,
      repeat_id := r
    ]
    
    pred_dt[
      ,
      fold_id := folds
    ]
    
    pred_dt[
      ,
      pred_prob := NA_real_
    ]
    
    for (fold in seq_len(n_folds)) {
      
      train_dt <- model_dt[
        folds != fold
      ]
      
      test_dt <- model_dt[
        folds == fold
      ]
      
      if (length(unique(train_dt$y)) < 2) {
        next
      }
      
      fit_cv <- glm(
        fit_formula,
        data = train_dt,
        family = binomial(),
        control = glm.control(maxit = 100)
      )
      
      pred_dt[
        fold_id == fold,
        pred_prob := as.numeric(
          predict(
            fit_cv,
            newdata = test_dt,
            type = "response"
          )
        )
      ]
    }
    
    cv_pred_list[[r]] <- pred_dt
    
    metrics_dt <- compute_metrics(
      pred_dt$y,
      pred_dt$pred_prob
    )
    
    metrics_dt[
      ,
      repeat_id := r
    ]
    
    cv_metrics_list[[r]] <- metrics_dt
  }
  
  cv_pred_dt <- rbindlist(
    cv_pred_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  cv_metrics_dt <- rbindlist(
    cv_metrics_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  cv_summary_dt <- cv_metrics_dt[
    ,
    .(
      auc_mean = mean(auc, na.rm = TRUE),
      auc_sd = sd(auc, na.rm = TRUE),
      accuracy_mean = mean(accuracy, na.rm = TRUE),
      accuracy_sd = sd(accuracy, na.rm = TRUE),
      sensitivity_mean = mean(sensitivity, na.rm = TRUE),
      sensitivity_sd = sd(sensitivity, na.rm = TRUE),
      specificity_mean = mean(specificity, na.rm = TRUE),
      specificity_sd = sd(specificity, na.rm = TRUE),
      balanced_accuracy_mean = mean(
        balanced_accuracy,
        na.rm = TRUE
      ),
      balanced_accuracy_sd = sd(
        balanced_accuracy,
        na.rm = TRUE
      )
    )
  ]
  
  cv_summary_dt[
    ,
    `:=`(
      dataset = dataset_label,
      model = paste0(
        "repeated_cv_",
        dataset_label,
        "_",
        version_label
      ),
      n_repeats = n_repeats,
      n_folds = n_folds
    )
  ]
  
  cv_pred_mean_dt <- cv_pred_dt[
    ,
    .(
      group = first(group),
      group_label = first(group_label),
      y = first(y),
      pred_prob_mean = mean(
        pred_prob,
        na.rm = TRUE
      )
    ),
    by = sample
  ]
  
  roc_dt <- build_roc_curve(
    cv_pred_mean_dt$y,
    cv_pred_mean_dt$pred_prob_mean
  )
  
  roc_auc <- compute_auc(
    cv_pred_mean_dt$y,
    cv_pred_mean_dt$pred_prob_mean
  )
  
  feature_plot_dt <- melt(
    model_dt[
      ,
      c(
        "sample",
        "group_label",
        needed_features
      ),
      with = FALSE
    ],
    id.vars = c(
      "sample",
      "group_label"
    ),
    variable.name = "feature",
    value.name = "value"
  )
  
  feature_plot_dt[
    ,
    feature_label := factor(
      feature,
     