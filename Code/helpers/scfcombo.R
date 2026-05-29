# scfcombo.R — R translation of scfcombo.ado
# Handles SCF multiple imputation + bootstrap variance estimation
#
# The SCF has 5 implicates (imp 1-5) and provides bootstrap replicate weights
# (wt1b1-wt1b999) and replicate draw counts (mm1-mm999) for variance estimation.
#
# Total variance = sampling_variance + (1 + 1/n_implicates) * imputation_variance

scfcombo <- function(formula, data, command = "lm", imps = 5, reps = 200,
                     weight_var = NULL, imp_var = "imp") {
  # Run the estimation command on each implicate and collect coefficients
  # Then run bootstrap replicates on implicate 1 to get sampling variance
  # Combine per Rubin's rules

  formula <- as.formula(formula)

  # Step 1: Imputation variance
  # Run model on each implicate, collect coefficients
  coef_list <- list()

  for (j in 1:imps) {
    sub <- data[data[[imp_var]] == j, ]

    if (!is.null(weight_var)) {
      fit <- do.call(command, list(formula = formula, data = sub, weights = sub[[weight_var]]))
    } else {
      fit <- do.call(command, list(formula = formula, data = sub))
    }

    coef_list[[j]] <- coef(fit)
  }

  # Average coefficients across implicates
  coef_mat <- do.call(rbind, coef_list)
  coef_avg <- colMeans(coef_mat, na.rm = TRUE)

  # Imputation variance (between-imputation)
  imp_var_mat <- var(coef_mat)

  # Step 2: Sampling variance via bootstrap (on first implicate)
  sub1 <- data[data[[imp_var]] == 1, ]

  boot_coefs <- list()
  mm_cols <- paste0("mm", 1:reps)
  wt_cols <- paste0("wt1b", 1:reps)

  # Check which bootstrap columns exist
  available_mm <- intersect(mm_cols, names(sub1))
  available_wt <- intersect(wt_cols, names(sub1))
  n_boot <- min(length(available_mm), length(available_wt), reps)

  if (n_boot > 0) {
    for (b in 1:n_boot) {
      mm_col <- available_mm[b]
      wt_col <- available_wt[b]

      # Create bootstrap replicate
      boot_data <- sub1[!is.na(sub1[[mm_col]]) & sub1[[mm_col]] < 1000, ]
      if (nrow(boot_data) == 0) next

      # Expand by mm counts
      boot_data <- boot_data[rep(seq_len(nrow(boot_data)), boot_data[[mm_col]]), ]

      tryCatch({
        if (!is.null(weight_var)) {
          fit_b <- do.call(command, list(formula = formula, data = boot_data,
                                         weights = boot_data[[wt_col]]))
        } else {
          fit_b <- do.call(command, list(formula = formula, data = boot_data))
        }
        boot_coefs[[length(boot_coefs) + 1]] <- coef(fit_b)
      }, error = function(e) NULL)
    }

    if (length(boot_coefs) > 1) {
      boot_mat <- do.call(rbind, boot_coefs)
      sampling_var_mat <- var(boot_mat)
    } else {
      sampling_var_mat <- matrix(0, length(coef_avg), length(coef_avg))
    }
  } else {
    # No bootstrap weights available; use model-based SEs
    if (!is.null(weight_var)) {
      fit_base <- do.call(command, list(formula = formula, data = sub1, weights = sub1[[weight_var]]))
    } else {
      fit_base <- do.call(command, list(formula = formula, data = sub1))
    }
    sampling_var_mat <- vcov(fit_base)
  }

  # Step 3: Combine variances (Rubin's rules)
  total_var <- sampling_var_mat + (1 + 1/imps) * imp_var_mat

  se <- sqrt(diag(total_var))
  t_stat <- coef_avg / se
  p_val <- 2 * pt(-abs(t_stat), df = Inf)

  result <- list(
    coefficients = coef_avg,
    se = se,
    vcov = total_var,
    t = t_stat,
    p = p_val,
    imp_variance = imp_var_mat,
    sampling_variance = sampling_var_mat,
    n_implicates = imps,
    n_boot_reps = n_boot
  )

  class(result) <- "scfcombo"
  return(result)
}

print.scfcombo <- function(x, ...) {
  cat("SCF Combined Estimates (", x$n_implicates, " implicates, ",
      x$n_boot_reps, " bootstrap reps)\n\n", sep = "")
  tab <- data.frame(
    Estimate = x$coefficients,
    Std.Error = x$se,
    t.value = x$t,
    p.value = x$p
  )
  print(round(tab, 4))
}
