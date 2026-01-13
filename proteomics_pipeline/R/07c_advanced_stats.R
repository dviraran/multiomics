# =============================================================================
# Advanced Statistical Methods for Proteomics
# =============================================================================
# Extended statistical analyses for proteomics data
#
# Key analyses:
# - Mixed effects models (for repeated measures/blocking factors)
# - Robust regression with outlier down-weighting
# - Effect size estimation with bootstrap confidence intervals
# - Power analysis for study design
# - Multiple testing correction comparison
# =============================================================================

# -----------------------------------------------------------------------------
# Main Orchestrator Function
# -----------------------------------------------------------------------------

#' Run Advanced Statistical Analysis
#'
#' Comprehensive advanced statistical methods for proteomics
#'
#' @param normalized_data Normalized protein intensity matrix
#' @param metadata Sample metadata
#' @param da_results Differential abundance results
#' @param config Pipeline configuration list
#' @return List containing advanced statistics results
#' @export
run_advanced_stats <- function(normalized_data, metadata, da_results, config) {
  log_message("=== Running Advanced Statistical Analysis ===")

  # Get config settings
  stats_config <- config$advanced_stats %||% list()

  if (!(stats_config$run_advanced_stats %||% TRUE)) {
    log_message("Advanced statistics disabled in config. Skipping.")
    return(NULL)
  }

  # Set up output directories
  output_dir <- file.path(config$output$output_dir, "advanced_stats")
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  results <- list()

  # 1. Effect size with bootstrap CIs
  if (stats_config$compute_effect_size_ci %||% TRUE) {
    results$effect_sizes <- compute_effect_sizes_with_ci(
      normalized_data = normalized_data,
      metadata = metadata,
      da_results = da_results,
      config = config,
      output_dir = output_dir
    )
  }

  # 2. Mixed effects models (if applicable)
  if (has_random_effects(metadata, config)) {
    results$mixed_effects <- run_mixed_effects_analysis(
      normalized_data = normalized_data,
      metadata = metadata,
      config = config,
      output_dir = output_dir
    )
  }

  # 3. Robust regression
  if (stats_config$run_robust_regression %||% FALSE) {
    results$robust_regression <- run_robust_regression(
      normalized_data = normalized_data,
      metadata = metadata,
      config = config,
      output_dir = output_dir
    )
  }

  # 4. Power analysis
  if (stats_config$run_power_analysis %||% TRUE) {
    results$power_analysis <- run_power_analysis(
      da_results = da_results,
      metadata = metadata,
      config = config,
      output_dir = output_dir
    )
  }

  # 5. Multiple testing comparison
  results$mtc_comparison <- compare_multiple_testing_methods(
    da_results = da_results,
    output_dir = output_dir
  )

  # Generate visualizations
  figures <- generate_advanced_stats_plots(
    results = results,
    da_results = da_results,
    config = config,
    plots_dir = plots_dir
  )
  results$figures <- figures

  # Create summary
  results$summary <- create_advanced_stats_summary(results)

  # Save outputs
  save_advanced_stats_outputs(results, output_dir)

  log_message("Advanced statistical analysis complete!")
  return(results)
}

# -----------------------------------------------------------------------------
# Effect Size with Bootstrap Confidence Intervals
# -----------------------------------------------------------------------------

#' Compute Effect Sizes with Bootstrap CIs
#'
#' @param normalized_data Normalized protein intensities
#' @param metadata Sample metadata
#' @param da_results DA results
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Data frame with effect sizes and CIs
compute_effect_sizes_with_ci <- function(normalized_data, metadata, da_results, config, output_dir) {
  log_message("Computing effect sizes with bootstrap confidence intervals...")

  # Get design info
  condition_col <- config$design$condition_column %||% "condition"
  ref_level <- config$design$reference_level %||% NULL

  if (!condition_col %in% colnames(metadata)) {
    log_message("WARNING: Condition column not found. Skipping effect size calculation.")
    return(NULL)
  }

  # Get conditions
  conditions <- metadata[[condition_col]]
  if (is.null(ref_level)) {
    ref_level <- levels(factor(conditions))[1]
  }

  # Bootstrap parameters
  n_boot <- config$advanced_stats$bootstrap_n %||% 1000
  ci_level <- config$advanced_stats$ci_level %||% 0.95

  log_message("  Bootstrap iterations: ", n_boot)
  log_message("  CI level: ", ci_level)

  # Handle case where da_results might be a list
  if (is.list(da_results) && !is.data.frame(da_results)) {
    if ("table" %in% names(da_results)) {
      da_results <- da_results$table
    } else if (length(da_results) > 0 && is.data.frame(da_results[[1]])) {
      da_results <- da_results[[1]]
    } else {
      log_message("WARNING: Could not extract data frame from da_results. Skipping effect size calculation.")
      return(NULL)
    }
  }

  # Ensure protein_id column exists
  if (!"protein_id" %in% colnames(da_results)) {
    da_results$protein_id <- rownames(da_results)
  }

  # Map column names if needed
  if (!"log2FoldChange" %in% colnames(da_results) && "logFC" %in% colnames(da_results)) {
    da_results$log2FoldChange <- da_results$logFC
  }
  if (!"pvalue" %in% colnames(da_results) && "P.Value" %in% colnames(da_results)) {
    da_results$pvalue <- da_results$P.Value
  }
  if (!"padj" %in% colnames(da_results) && "adj.P.Val" %in% colnames(da_results)) {
    da_results$padj <- da_results$adj.P.Val
  }

  # Get protein list from DA results
  proteins <- da_results$protein_id

  if (is.null(proteins) || length(proteins) == 0) {
    log_message("WARNING: No proteins found in DA results. Skipping effect size calculation.")
    return(NULL)
  }

  # Compute effect sizes in parallel
  if (requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- config$n_cores %||% parallel::detectCores() - 1
    n_cores <- min(n_cores, 8)  # Cap at 8 cores
  } else {
    n_cores <- 1
  }

  log_message("  Processing ", length(proteins), " proteins...")

  effect_results <- lapply(proteins, function(protein) {
    tryCatch({
      compute_protein_effect_size(
        protein = protein,
        expr_data = normalized_data,
        conditions = conditions,
        ref_level = ref_level,
        n_boot = n_boot,
        ci_level = ci_level
      )
    }, error = function(e) {
      data.frame(
        protein_id = protein,
        effect_size = NA,
        effect_size_type = NA,
        ci_lower = NA,
        ci_upper = NA,
        se = NA,
        stringsAsFactors = FALSE
      )
    })
  })

  effect_df <- do.call(rbind, effect_results)

  # Merge with DA results - only include columns that exist
  merge_cols <- intersect(c("protein_id", "log2FoldChange", "pvalue", "padj"), colnames(da_results))
  if (length(merge_cols) > 1) {  # Need at least protein_id + one other
    effect_df <- merge(
      effect_df,
      da_results[, merge_cols, drop = FALSE],
      by = "protein_id",
      all.x = TRUE
    )
  }

  # Determine significance of effect sizes (CI excludes 0)
  effect_df$significant <- effect_df$ci_lower > 0 | effect_df$ci_upper < 0

  # Save results
  write.csv(
    effect_df,
    file.path(output_dir, "effect_sizes_with_ci.csv"),
    row.names = FALSE
  )

  log_message("  Computed effect sizes for ", sum(!is.na(effect_df$effect_size)), " proteins")

  return(effect_df)
}

#' Compute Effect Size for Single Protein
compute_protein_effect_size <- function(protein, expr_data, conditions, ref_level, n_boot, ci_level) {
  # Get expression values
  if (protein %in% rownames(expr_data)) {
    values <- as.numeric(expr_data[protein, ])
  } else {
    return(data.frame(
      protein_id = protein,
      effect_size = NA,
      effect_size_type = NA,
      ci_lower = NA,
      ci_upper = NA,
      se = NA,
      stringsAsFactors = FALSE
    ))
  }

  # Get groups
  unique_conds <- unique(conditions)

  if (length(unique_conds) == 2) {
    # Two groups: compute Cohen's d or Hedge's g
    group1 <- values[conditions == ref_level]
    group2 <- values[conditions != ref_level]

    # Remove NAs
    group1 <- group1[!is.na(group1)]
    group2 <- group2[!is.na(group2)]

    if (length(group1) < 2 || length(group2) < 2) {
      return(data.frame(
        protein_id = protein,
        effect_size = NA,
        effect_size_type = NA,
        ci_lower = NA,
        ci_upper = NA,
        se = NA,
        stringsAsFactors = FALSE
      ))
    }

    # Compute Hedge's g (bias-corrected Cohen's d)
    effect_size <- compute_hedges_g(group1, group2)

    # Bootstrap CI
    boot_effects <- replicate(n_boot, {
      boot_g1 <- sample(group1, replace = TRUE)
      boot_g2 <- sample(group2, replace = TRUE)
      compute_hedges_g(boot_g1, boot_g2)
    })

    boot_effects <- boot_effects[is.finite(boot_effects)]

    if (length(boot_effects) < 10) {
      ci_lower <- NA
      ci_upper <- NA
      se <- NA
    } else {
      alpha <- 1 - ci_level
      ci_lower <- quantile(boot_effects, alpha / 2, na.rm = TRUE)
      ci_upper <- quantile(boot_effects, 1 - alpha / 2, na.rm = TRUE)
      se <- sd(boot_effects, na.rm = TRUE)
    }

    return(data.frame(
      protein_id = protein,
      effect_size = effect_size,
      effect_size_type = "Hedges_g",
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      se = se,
      stringsAsFactors = FALSE
    ))

  } else {
    # Multiple groups: compute eta-squared or partial eta-squared
    df <- data.frame(y = values, group = factor(conditions))
    df <- df[complete.cases(df), ]

    if (nrow(df) < 5) {
      return(data.frame(
        protein_id = protein,
        effect_size = NA,
        effect_size_type = NA,
        ci_lower = NA,
        ci_upper = NA,
        se = NA,
        stringsAsFactors = FALSE
      ))
    }

    # Compute eta-squared
    effect_size <- compute_eta_squared(df$y, df$group)

    # Bootstrap CI
    boot_effects <- replicate(n_boot, {
      boot_idx <- sample(nrow(df), replace = TRUE)
      boot_df <- df[boot_idx, ]
      compute_eta_squared(boot_df$y, boot_df$group)
    })

    boot_effects <- boot_effects[is.finite(boot_effects)]

    if (length(boot_effects) < 10) {
      ci_lower <- NA
      ci_upper <- NA
      se <- NA
    } else {
      alpha <- 1 - ci_level
      ci_lower <- quantile(boot_effects, alpha / 2, na.rm = TRUE)
      ci_upper <- quantile(boot_effects, 1 - alpha / 2, na.rm = TRUE)
      se <- sd(boot_effects, na.rm = TRUE)
    }

    return(data.frame(
      protein_id = protein,
      effect_size = effect_size,
      effect_size_type = "eta_squared",
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      se = se,
      stringsAsFactors = FALSE
    ))
  }
}

#' Compute Hedge's g (Bias-Corrected Cohen's d)
compute_hedges_g <- function(group1, group2) {
  n1 <- length(group1)
  n2 <- length(group2)

  mean1 <- mean(group1, na.rm = TRUE)
  mean2 <- mean(group2, na.rm = TRUE)

  # Pooled SD
  var1 <- var(group1, na.rm = TRUE)
  var2 <- var(group2, na.rm = TRUE)
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

  if (pooled_sd == 0) return(NA)

  # Cohen's d
  d <- (mean2 - mean1) / pooled_sd

  # Hedge's g correction factor
  correction <- 1 - (3 / (4 * (n1 + n2) - 9))

  g <- d * correction

  return(g)
}

#' Compute Eta-Squared
compute_eta_squared <- function(y, group) {
  # ANOVA
  fit <- tryCatch({
    aov(y ~ group)
  }, error = function(e) NULL)

  if (is.null(fit)) return(NA)

  ss <- summary(fit)[[1]]
  ss_between <- ss["group", "Sum Sq"]
  ss_total <- sum(ss[, "Sum Sq"])

  eta_sq <- ss_between / ss_total

  return(eta_sq)
}

# -----------------------------------------------------------------------------
# Mixed Effects Models
# -----------------------------------------------------------------------------

#' Check if Mixed Effects Models are Appropriate
has_random_effects <- function(metadata, config) {
  random_effects_cols <- config$design$random_effects %||% NULL

  if (is.null(random_effects_cols)) {
    # Check for common blocking variables
    potential_random <- c("subject", "patient", "donor", "replicate", "batch", "block")
    random_effects_cols <- intersect(potential_random, colnames(metadata))
  }

  if (length(random_effects_cols) > 0) {
    # Check if any column has multiple observations per level
    for (col in random_effects_cols) {
      if (col %in% colnames(metadata)) {
        counts <- table(metadata[[col]])
        if (any(counts > 1)) {
          return(TRUE)
        }
      }
    }
  }

  return(FALSE)
}

#' Run Mixed Effects Analysis
#'
#' @param normalized_data Normalized protein intensities
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Data frame with mixed effects results
run_mixed_effects_analysis <- function(normalized_data, metadata, config, output_dir) {
  log_message("Running mixed effects models...")

  # Check for lme4
  if (!requireNamespace("lme4", quietly = TRUE)) {
    log_message("WARNING: lme4 package not available. Skipping mixed effects analysis.")
    return(NULL)
  }

  # Get design info
  condition_col <- config$design$condition_column %||% "condition"
  random_effects_cols <- config$design$random_effects %||% c("subject", "patient", "donor")
  random_effects_cols <- intersect(random_effects_cols, colnames(metadata))

  if (length(random_effects_cols) == 0) {
    log_message("  No valid random effects columns found.")
    return(NULL)
  }

  random_col <- random_effects_cols[1]  # Use first available
  log_message("  Fixed effect: ", condition_col)
  log_message("  Random effect: ", random_col)

  # Prepare data
  proteins <- rownames(normalized_data)

  # Fit mixed models
  mixed_results <- lapply(proteins, function(protein) {
    tryCatch({
      fit_mixed_model(
        protein = protein,
        expr_data = normalized_data,
        metadata = metadata,
        condition_col = condition_col,
        random_col = random_col
      )
    }, error = function(e) {
      data.frame(
        protein_id = protein,
        fixed_effect = NA,
        fixed_se = NA,
        fixed_pvalue = NA,
        random_variance = NA,
        residual_variance = NA,
        icc = NA,
        converged = FALSE,
        stringsAsFactors = FALSE
      )
    })
  })

  mixed_df <- do.call(rbind, mixed_results)

  # Adjust p-values
  mixed_df$padj <- p.adjust(mixed_df$fixed_pvalue, method = "BH")

  # Sort by significance
  mixed_df <- mixed_df[order(mixed_df$padj), ]

  # Save results
  write.csv(
    mixed_df,
    file.path(output_dir, "mixed_effects_results.csv"),
    row.names = FALSE
  )

  log_message("  Converged models: ", sum(mixed_df$converged, na.rm = TRUE), "/", nrow(mixed_df))
  log_message("  Significant (FDR < 0.05): ", sum(mixed_df$padj < 0.05, na.rm = TRUE))

  return(mixed_df)
}

#' Fit Mixed Model for Single Protein
fit_mixed_model <- function(protein, expr_data, metadata, condition_col, random_col) {
  # Get expression
  y <- as.numeric(expr_data[protein, ])

  # Build data frame
  df <- data.frame(
    y = y,
    condition = factor(metadata[[condition_col]]),
    random = factor(metadata[[random_col]])
  )
  df <- df[complete.cases(df), ]

  if (nrow(df) < 5 || length(unique(df$random)) < 2) {
    return(data.frame(
      protein_id = protein,
      fixed_effect = NA,
      fixed_se = NA,
      fixed_pvalue = NA,
      random_variance = NA,
      residual_variance = NA,
      icc = NA,
      converged = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  # Fit model
  model <- lme4::lmer(y ~ condition + (1 | random), data = df)

  # Get fixed effects
  fixed_coef <- lme4::fixef(model)
  fixed_se <- sqrt(diag(as.matrix(vcov(model))))

  # Get variance components
  var_comp <- as.data.frame(lme4::VarCorr(model))

  random_var <- var_comp$vcov[var_comp$grp == "random"]
  resid_var <- var_comp$vcov[var_comp$grp == "Residual"]

  # Intraclass correlation
  icc <- random_var / (random_var + resid_var)

  # P-value for fixed effect (using Satterthwaite if available)
  pvalue <- tryCatch({
    if (requireNamespace("lmerTest", quietly = TRUE)) {
      model_test <- lmerTest::lmer(y ~ condition + (1 | random), data = df)
      summary(model_test)$coefficients[2, "Pr(>|t|)"]
    } else {
      # Approximate using z-test
      z <- fixed_coef[2] / fixed_se[2]
      2 * pnorm(-abs(z))
    }
  }, error = function(e) NA)

  return(data.frame(
    protein_id = protein,
    fixed_effect = fixed_coef[2],
    fixed_se = fixed_se[2],
    fixed_pvalue = pvalue,
    random_variance = random_var,
    residual_variance = resid_var,
    icc = icc,
    converged = TRUE,
    stringsAsFactors = FALSE
  ))
}

# -----------------------------------------------------------------------------
# Robust Regression
# -----------------------------------------------------------------------------

#' Run Robust Regression Analysis
#'
#' @param normalized_data Normalized protein intensities
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Data frame with robust regression results
run_robust_regression <- function(normalized_data, metadata, config, output_dir) {
  log_message("Running robust regression analysis...")

  # Check for MASS package
  if (!requireNamespace("MASS", quietly = TRUE)) {
    log_message("WARNING: MASS package not available. Skipping robust regression.")
    return(NULL)
  }

  # Get design info
  condition_col <- config$design$condition_column %||% "condition"

  # Prepare data
  proteins <- rownames(normalized_data)

  # Fit robust models
  robust_results <- lapply(proteins, function(protein) {
    tryCatch({
      fit_robust_regression(
        protein = protein,
        expr_data = normalized_data,
        metadata = metadata,
        condition_col = condition_col
      )
    }, error = function(e) {
      data.frame(
        protein_id = protein,
        robust_coef = NA,
        robust_se = NA,
        robust_pvalue = NA,
        n_outliers = NA,
        stringsAsFactors = FALSE
      )
    })
  })

  robust_df <- do.call(rbind, robust_results)

  # Adjust p-values
  robust_df$padj <- p.adjust(robust_df$robust_pvalue, method = "BH")

  # Save results
  write.csv(
    robust_df,
    file.path(output_dir, "robust_regression_results.csv"),
    row.names = FALSE
  )

  return(robust_df)
}

#' Fit Robust Regression for Single Protein
fit_robust_regression <- function(protein, expr_data, metadata, condition_col) {
  # Get expression
  y <- as.numeric(expr_data[protein, ])

  # Build data frame
  df <- data.frame(
    y = y,
    condition = factor(metadata[[condition_col]])
  )
  df <- df[complete.cases(df), ]

  if (nrow(df) < 5) {
    return(data.frame(
      protein_id = protein,
      robust_coef = NA,
      robust_se = NA,
      robust_pvalue = NA,
      n_outliers = NA,
      stringsAsFactors = FALSE
    ))
  }

  # Fit robust linear model (M-estimation)
  model <- MASS::rlm(y ~ condition, data = df)

  # Get coefficients
  coefs <- summary(model)$coefficients

  if (nrow(coefs) < 2) {
    return(data.frame(
      protein_id = protein,
      robust_coef = NA,
      robust_se = NA,
      robust_pvalue = NA,
      n_outliers = NA,
      stringsAsFactors = FALSE
    ))
  }

  robust_coef <- coefs[2, 1]
  robust_se <- coefs[2, 2]
  t_value <- coefs[2, 3]

  # Approximate p-value
  df_resid <- nrow(df) - 2
  pvalue <- 2 * pt(-abs(t_value), df_resid)

  # Count outliers (weights < 0.9)
  n_outliers <- sum(model$w < 0.9)

  return(data.frame(
    protein_id = protein,
    robust_coef = robust_coef,
    robust_se = robust_se,
    robust_pvalue = pvalue,
    n_outliers = n_outliers,
    stringsAsFactors = FALSE
  ))
}

# -----------------------------------------------------------------------------
# Power Analysis
# -----------------------------------------------------------------------------

#' Run Power Analysis
#'
#' @param da_results DA results
#' @param metadata Sample metadata
#' @param config Pipeline configuration
#' @param output_dir Output directory
#' @return Power analysis results
run_power_analysis <- function(da_results, metadata, config, output_dir) {
  log_message("Running power analysis...")

  # Get current sample size
  condition_col <- config$design$condition_column %||% "condition"
  conditions <- metadata[[condition_col]]
  group_sizes <- table(conditions)

  current_n <- min(group_sizes)

  # Estimate effect sizes from data
  significant <- da_results$padj < 0.05 & !is.na(da_results$padj)
  effect_sizes <- abs(da_results$log2FoldChange[significant])

  if (length(effect_sizes) == 0) {
    log_message("  No significant proteins for power estimation.")
    return(NULL)
  }

  median_effect <- median(effect_sizes, na.rm = TRUE)
  sd_effect <- mad(effect_sizes, na.rm = TRUE)

  log_message("  Current sample size per group: ", current_n)
  log_message("  Median effect size (|log2FC|): ", round(median_effect, 2))

  # Power calculation for different sample sizes
  sample_sizes <- c(3, 5, 8, 10, 15, 20, 30, 50)
  alpha <- 0.05

  power_results <- lapply(sample_sizes, function(n) {
    # Effect sizes to test
    effect_tests <- c(0.5, 1.0, 1.5, 2.0, median_effect)

    powers <- sapply(effect_tests, function(d) {
      calculate_power(n = n, effect_size = d, alpha = alpha)
    })

    data.frame(
      n_per_group = n,
      effect_0.5 = powers[1],
      effect_1.0 = powers[2],
      effect_1.5 = powers[3],
      effect_2.0 = powers[4],
      effect_median = powers[5],
      stringsAsFactors = FALSE
    )
  })

  power_df <- do.call(rbind, power_results)

  # Required sample size for 80% power at different effect sizes
  required_n <- sapply(c(0.5, 1.0, 1.5, 2.0, median_effect), function(d) {
    calculate_required_n(effect_size = d, power = 0.80, alpha = alpha)
  })

  required_df <- data.frame(
    effect_size = c(0.5, 1.0, 1.5, 2.0, median_effect),
    required_n = required_n,
    stringsAsFactors = FALSE
  )

  # Save results
  write.csv(power_df, file.path(output_dir, "power_by_sample_size.csv"), row.names = FALSE)
  write.csv(required_df, file.path(output_dir, "required_sample_sizes.csv"), row.names = FALSE)

  return(list(
    power_table = power_df,
    required_n = required_df,
    current_n = current_n,
    median_effect = median_effect
  ))
}

#' Calculate Power for Two-Sample t-Test
calculate_power <- function(n, effect_size, alpha = 0.05) {
  # Use approximation based on non-central t-distribution
  df <- 2 * n - 2
  ncp <- effect_size * sqrt(n / 2)  # Non-centrality parameter
  t_crit <- qt(1 - alpha / 2, df)

  power <- 1 - pt(t_crit, df, ncp) + pt(-t_crit, df, ncp)

  return(power)
}

#' Calculate Required Sample Size
calculate_required_n <- function(effect_size, power = 0.80, alpha = 0.05) {
  # Binary search for required n
  for (n in seq(3, 500)) {
    current_power <- calculate_power(n, effect_size, alpha)
    if (current_power >= power) {
      return(n)
    }
  }
  return(NA)
}

# -----------------------------------------------------------------------------
# Multiple Testing Correction Comparison
# -----------------------------------------------------------------------------

#' Compare Multiple Testing Correction Methods
#'
#' @param da_results DA results with raw p-values
#' @param output_dir Output directory
#' @return Comparison of different MTC methods
compare_multiple_testing_methods <- function(da_results, output_dir) {
  log_message("Comparing multiple testing correction methods...")

  pvalues <- da_results$pvalue[!is.na(da_results$pvalue)]

  methods <- c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr")

  mtc_comparison <- lapply(methods, function(method) {
    padj <- p.adjust(pvalues, method = method)

    data.frame(
      method = method,
      n_significant_0.01 = sum(padj < 0.01, na.rm = TRUE),
      n_significant_0.05 = sum(padj < 0.05, na.rm = TRUE),
      n_significant_0.10 = sum(padj < 0.10, na.rm = TRUE),
      median_padj = median(padj, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  mtc_df <- do.call(rbind, mtc_comparison)

  # Save results
  write.csv(mtc_df, file.path(output_dir, "mtc_comparison.csv"), row.names = FALSE)

  log_message("  Methods compared: ", paste(methods, collapse = ", "))

  return(mtc_df)
}

# -----------------------------------------------------------------------------
# Visualization Functions
# -----------------------------------------------------------------------------

#' Generate Advanced Statistics Plots
#'
#' @param results All advanced stats results
#' @param da_results DA results
#' @param config Pipeline configuration
#' @param plots_dir Plots output directory
#' @return List of figure paths
generate_advanced_stats_plots <- function(results, da_results, config, plots_dir) {
  log_message("Generating advanced statistics plots...")

  figures <- list()

  # 1. Effect size forest plot
  if (!is.null(results$effect_sizes)) {
    figures$forest <- plot_effect_size_forest(
      effect_df = results$effect_sizes,
      plots_dir = plots_dir
    )
  }

  # 2. Power analysis curves
  if (!is.null(results$power_analysis)) {
    figures$power <- plot_power_curves(
      power_results = results$power_analysis,
      plots_dir = plots_dir
    )
  }

  # 3. MTC comparison
  if (!is.null(results$mtc_comparison)) {
    figures$mtc <- plot_mtc_comparison(
      mtc_df = results$mtc_comparison,
      plots_dir = plots_dir
    )
  }

  # 4. Mixed effects ICC distribution
  if (!is.null(results$mixed_effects)) {
    figures$icc <- plot_icc_distribution(
      mixed_df = results$mixed_effects,
      plots_dir = plots_dir
    )
  }

  return(figures)
}

#' Plot Effect Size Forest Plot
plot_effect_size_forest <- function(effect_df, plots_dir) {
  log_message("  Creating effect size forest plot...")

  # Get significant proteins
  sig_df <- effect_df %>%
    filter(significant == TRUE, !is.na(effect_size)) %>%
    arrange(desc(abs(effect_size))) %>%
    head(30)

  if (nrow(sig_df) == 0) {
    log_message("  No significant effect sizes to plot.")
    return(NULL)
  }

  # Order by effect size
  sig_df$protein_id <- factor(sig_df$protein_id, levels = sig_df$protein_id[order(sig_df$effect_size)])

  p <- ggplot(sig_df, aes(x = effect_size, y = protein_id)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, color = "grey30") +
    geom_point(aes(color = effect_size > 0), size = 2) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                       labels = c("TRUE" = "Positive", "FALSE" = "Negative"),
                       name = "Direction") +
    labs(
      title = "Effect Sizes with 95% Bootstrap CI",
      subtitle = paste0("Top ", nrow(sig_df), " significant proteins (", sig_df$effect_size_type[1], ")"),
      x = "Effect Size",
      y = "Protein"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))

  fig_path <- file.path(plots_dir, "effect_size_forest_plot.png")
  ggsave(fig_path, p, width = 10, height = 12)

  return(fig_path)
}

#' Plot Power Analysis Curves
plot_power_curves <- function(power_results, plots_dir) {
  log_message("  Creating power analysis plots...")

  power_df <- power_results$power_table

  # Pivot longer for plotting
  power_long <- power_df %>%
    tidyr::pivot_longer(
      cols = starts_with("effect_"),
      names_to = "effect_size",
      values_to = "power"
    ) %>%
    mutate(
      effect_size = gsub("effect_", "", effect_size),
      effect_size = ifelse(effect_size == "median",
                           paste0("median (", round(power_results$median_effect, 2), ")"),
                           effect_size)
    )

  p1 <- ggplot(power_long, aes(x = n_per_group, y = power, color = effect_size)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
    geom_vline(xintercept = power_results$current_n, linetype = "dotted", color = "blue") +
    annotate("text", x = power_results$current_n, y = 0.5,
             label = paste0("Current n=", power_results$current_n),
             angle = 90, vjust = -0.5, color = "blue") +
    scale_x_continuous(breaks = unique(power_df$n_per_group)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Statistical Power by Sample Size",
      x = "Sample Size per Group",
      y = "Power",
      color = "Effect Size\n(log2FC)"
    ) +
    theme_minimal()

  # Required sample size plot
  required_df <- power_results$required_n

  p2 <- ggplot(required_df, aes(x = effect_size, y = required_n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = required_n), vjust = -0.5, size = 3) +
    labs(
      title = "Required Sample Size for 80% Power",
      x = "Effect Size (log2FC)",
      y = "Required n per group"
    ) +
    theme_minimal()

  combined <- patchwork::wrap_plots(p1, p2, ncol = 2)

  fig_path <- file.path(plots_dir, "power_analysis.png")
  ggsave(fig_path, combined, width = 14, height = 6)

  return(fig_path)
}

#' Plot MTC Comparison
plot_mtc_comparison <- function(mtc_df, plots_dir) {
  log_message("  Creating MTC comparison plot...")

  # Reshape for plotting
  mtc_long <- mtc_df %>%
    tidyr::pivot_longer(
      cols = starts_with("n_significant"),
      names_to = "threshold",
      values_to = "n_significant"
    ) %>%
    mutate(
      threshold = gsub("n_significant_", "FDR < ", threshold)
    )

  p <- ggplot(mtc_long, aes(x = method, y = n_significant, fill = threshold)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "Multiple Testing Correction Comparison",
      x = "Method",
      y = "Number of Significant Proteins",
      fill = "Threshold"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  fig_path <- file.path(plots_dir, "mtc_comparison.png")
  ggsave(fig_path, p, width = 10, height = 6)

  return(fig_path)
}

#' Plot ICC Distribution
plot_icc_distribution <- function(mixed_df, plots_dir) {
  log_message("  Creating ICC distribution plot...")

  icc_values <- mixed_df$icc[mixed_df$converged & !is.na(mixed_df$icc)]

  if (length(icc_values) == 0) {
    return(NULL)
  }

  p <- ggplot(data.frame(icc = icc_values), aes(x = icc)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = median(icc_values), linetype = "dashed", color = "red") +
    annotate("text", x = median(icc_values), y = Inf,
             label = paste0("Median = ", round(median(icc_values), 2)),
             vjust = 2, hjust = -0.1, color = "red") +
    labs(
      title = "Distribution of Intraclass Correlation Coefficients",
      subtitle = "From mixed effects models",
      x = "ICC",
      y = "Count"
    ) +
    theme_minimal()

  fig_path <- file.path(plots_dir, "icc_distribution.png")
  ggsave(fig_path, p, width = 8, height = 6)

  return(fig_path)
}

# -----------------------------------------------------------------------------
# Output Functions
# -----------------------------------------------------------------------------

#' Create Advanced Stats Summary
create_advanced_stats_summary <- function(results) {
  summary <- list()

  if (!is.null(results$effect_sizes)) {
    summary$n_significant_effects <- sum(results$effect_sizes$significant, na.rm = TRUE)
    summary$median_effect_size <- median(abs(results$effect_sizes$effect_size), na.rm = TRUE)
  }

  if (!is.null(results$mixed_effects)) {
    summary$n_converged_mixed <- sum(results$mixed_effects$converged, na.rm = TRUE)
    summary$median_icc <- median(results$mixed_effects$icc, na.rm = TRUE)
    summary$n_significant_mixed <- sum(results$mixed_effects$padj < 0.05, na.rm = TRUE)
  }

  if (!is.null(results$power_analysis)) {
    summary$current_power <- results$power_analysis$power_table$effect_median[
      results$power_analysis$power_table$n_per_group == results$power_analysis$current_n
    ]
  }

  return(summary)
}

#' Save Advanced Stats Outputs
save_advanced_stats_outputs <- function(results, output_dir) {
  log_message("Saving advanced statistics outputs...")

  # Save full results as RDS
  saveRDS(results, file.path(output_dir, "advanced_stats_results.rds"))

  log_message("  Outputs saved to: ", output_dir)
}
