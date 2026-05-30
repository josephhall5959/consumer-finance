# 21_nested_logit.R — Nested logit demand model with counterfactuals

cat("Running nested logit model...\n")

source(file.path(ROOT, "Code", "helpers", "nlfunctions.R"))

alpha <- -4

# Market structure calibration
oo <- 0.92  # Outside option share
nf <- 7     # Total firms
ts <- 0.47  # Top firm conditional PQ share
na_firms <- 2  # Already-accepting firms
as_share <- 0.17 / na_firms  # Accepting firms conditional share
ns <- 1     # Switching firms
ss <- 0.20 / ns  # Switching firm share
no <- nf - 1 - na_firms - ns  # Other firms

pq_0 <- c(
  ts * (1 - oo),
  rep((1 - ts - ss * ns - as_share * na_firms) * (1 - oo) / no, no),
  rep(as_share * (1 - oo), na_firms),
  rep(ss * (1 - oo), ns)
)

# Consumer usage
u_0 <- 0.37
u_1 <- u_0 + 0.055

# Pre-shock structure
a_0 <- c(0, 0, 0, 0, 1, 1, 0) * u_0
g_0 <- cbind(c(1, 2, 3, 4, 5, 5, 6), c(1, 1, 1, 1, 1, 1, 1))

# Post-shock structure
a_1 <- c(0, 0, 0, 0, 1, 1, 1) * u_1
g_1 <- cbind(c(1, 2, 3, 4, 5, 5, 5), c(1, 1, 1, 1, 1, 1, 1))

# Optimization objective
to_optimize <- function(free_params, vector = FALSE) {
  sigma <- c(free_params[1], 0)
  zeta <- free_params[2]

  c_0 <- rep(1, 7) + a_0 * zeta
  c_1_local <- rep(1, 7) + a_1 * zeta

  p_delta_0 <- c(rep(2, 7), c(1, rep(0, 6)))

  pq_moments <- function(params) {
    p <- params[1:7]
    delta_0 <- params[8:14]
    delta <- delta_0 + alpha * p
    s <- f.share(delta, sigma, g_0)$sj
    s_cond <- f.share(delta, sigma, g_0)$sj_hg
    pq <- s * p
    foc <- s + (p - c_0) * (alpha / (1 - sigma[1])) * s * (1 - sigma[1] * s_cond - (1 - sigma[1]) * s)
    return(c(pq - pq_0, foc))
  }

  p_delta_1 <- nleqslv(p_delta_0, pq_moments)
  delta_0 <- p_delta_1$x[8:14]
  p_0 <- p_delta_1$x[1:7]
  s_0 <- f.share(delta_0 + alpha * p_0, sigma, g_0)$sj

  p_moments <- function(params) {
    p <- params[1:7]
    delta <- delta_0 + alpha * p
    s <- f.share(delta, sigma, g_1)$sj
    s_cond <- f.share(delta, sigma, g_1)$sj_hg
    foc <- s + (p - c_1_local) * (alpha / (1 - sigma[1])) * s * (1 - sigma[1] * s_cond - (1 - sigma[1]) * s)
    return(foc)
  }

  p_1 <- nleqslv(p_0, p_moments)$x
  s_1 <- f.share(delta_0 + alpha * p_1, sigma, g_1)$sj

  delta_t <- mean(s_1 * p_1 / (s_0 * p_0) - 1)
  delta_a <- mean(s_1[5:6] * p_1[5:6] / (s_0[5:6] * p_0[5:6])) -
    mean(s_1[-c(5:6)] * p_1[-c(5:6)] / (s_0[-c(5:6)] * p_0[-c(5:6)]))

  if (vector) return(c(delta_t, delta_a))

  errors <- c(delta_t, delta_a) - target_moments
  return(errors %*% errors)
}

# Target moments and their standard errors flow from the reduced-form sales
# DiD regressions estimated in 13_did_marquette.R (written to smm_moments.csv):
#   Moment 1 = average sales change from Marquette (treated x post)
#   Moment 2 = accept-vs-not sales differential    (treated x post x card)
# This makes structural-parameter uncertainty propagate directly from the
# regression SEs through the SMM to SE(sigma).
smm_file <- file.path(DATA_GEN, "smm_moments.csv")
if (file.exists(smm_file)) {
  smm_in <- fread(smm_file)
  target_moments <- as.numeric(smm_in$estimate)
  target_se <- as.numeric(smm_in$se)
  # Full 2x2 covariance of the moments (off-diagonal from cov12 if present).
  cov12 <- if ("cov12" %in% names(smm_in)) as.numeric(smm_in$cov12[1]) else 0
  Sigma_m <- matrix(c(target_se[1]^2, cov12, cov12, target_se[2]^2), 2, 2)
  cat(sprintf("  Loaded target moments from smm_moments.csv: (%.4f, %.4f), SE (%.4f, %.4f), corr %.3f\n",
              target_moments[1], target_moments[2], target_se[1], target_se[2],
              cov12 / (target_se[1] * target_se[2])))
} else {
  warning("smm_moments.csv not found; using fallback hardcoded moments.")
  target_moments <- c(0.0027, -0.048)
  target_se <- c(0.005, 0.02)
  Sigma_m <- diag(target_se^2)
}

# --- Helper: estimate for a given alpha and target_moments ---
estimate_model <- function(alpha_val, tm, verbose = FALSE) {
  # Override alpha in the to_optimize closure
  env <- environment(to_optimize)
  old_alpha <- env$alpha
  env$alpha <- alpha_val

  P <- tryCatch(
    optim(c(0, 0), to_optimize),
    error = function(e) list(par = c(NA, NA), convergence = 1)
  )

  env$alpha <- old_alpha  # restore
  if (verbose && P$convergence == 0) {
    cat("  alpha =", alpha_val, ": sigma =", P$par[1], ", zeta =", P$par[2], "\n")
  }
  P$par
}

# Optimize at baseline alpha = -4
sigma_0 <- 0
zeta_0 <- 0

P <- optim(c(sigma_0, zeta_0), to_optimize)
params <- P$par
cat("  Estimated parameters: sigma =", params[1], ", zeta =", params[2], "\n")
cat("  Cost inclusive of usage:", params[2] * 0.37, "\n")
cat("  Model moments:", to_optimize(params, vector = TRUE), "\n")
cat("  Target moments:", target_moments, "\n")

# --- A3: Delta-method standard errors ---
# Because the SMM is exactly identified (2 moments -> 2 parameters), the model
# matches the targets exactly and theta(m) is the inverse of the moment map
# M(theta).  We propagate the regression moments' covariance Sigma_m through
# this map: Var(theta) = J Sigma_m J', with J = (dM/dtheta)^{-1} = dtheta/dm.
# The Jacobian dM/dtheta is computed by central finite differences around the
# estimate -- a few model solves rather than thousands of bootstrap re-solves.
cat("  Computing delta-method standard errors...\n")

model_moments <- function(theta) to_optimize(theta, vector = TRUE)

h_fd <- 1e-5
K <- matrix(0, 2, 2)                       # K[, j] = dM / d theta_j
for (j in 1:2) {
  tp <- params; tm <- params
  tp[j] <- tp[j] + h_fd; tm[j] <- tm[j] - h_fd
  K[, j] <- (model_moments(tp) - model_moments(tm)) / (2 * h_fd)
}
J <- solve(K)                              # d theta / d m
V_theta <- J %*% Sigma_m %*% t(J)          # full param covariance (off-diag kept)
delta_se <- sqrt(diag(V_theta))

cat("  Delta-method SE: sigma =", delta_se[1], ", zeta =", delta_se[2], "\n")
cat("  Delta-method SE of cost:", delta_se[2] * 0.37, "\n")
cat(sprintf("  corr(sigma, zeta) = %.3f\n",
            V_theta[1, 2] / sqrt(V_theta[1, 1] * V_theta[2, 2])))
cat(sprintf("  sigma = %.4f (SE %.4f, t = %.2f);  zeta = %.4f (SE %.4f, t = %.2f)\n",
            params[1], delta_se[1], params[1] / delta_se[1],
            params[2], delta_se[2], params[2] / delta_se[2]))

# Calibration-results table with delta-method SEs (read by the paper, Table 5).
# Reproduces the full Parameters + Target Moments table entirely from the model
# objects so nothing is hardcoded in main.tex.
tryCatch({
  mm <- model_moments(params)                      # model-implied moments at solution
  cost_reduction <- params[2] * u_0                # implied cost reduction (zeta * u_0)

  # value/percent formatters
  fmt_se   <- function(v, s) sprintf("$%.3f$ (%.3f)", v, s)     # signed value (s.e.)
  fmt_pct  <- function(x) sprintf("$%+g\\%%$", round(x * 100, 1))
  fmt_pcts <- function(x, s) sprintf("$%+g\\%%$ (%g\\%%)",
                                     round(x * 100, 1), round(s * 100, 1))

  firm_store_share <- 1 - as_share * na_firms      # firms' store-card conditional share
  d_store <- -ss * ns                              # store-card share change (calibration)
  d_bank  <- u_1 - u_0                             # bank-card share change (calibration)

  est_lines <- c(
    "\\begin{tabular}{lrc}",
    "\\multicolumn{3}{c}{\\textbf{Parameters}} \\\\",
    "Parameter & Value & Method \\\\",
    "\\hline",
    sprintf("$\\alpha$ & $%g$ & Literature \\\\", alpha),
    sprintf("$\\sigma$ & %s & Matched \\\\", fmt_se(params[1], delta_se[1])),
    sprintf("$\\zeta$ & %s & Matched \\\\", fmt_se(params[2], delta_se[2])),
    sprintf("N Firms & %d & Observed\\\\", as.integer(nf)),
    "N Bank Networks & 1 & Observed\\\\",
    sprintf("Top Firm Share & %g\\%% & Observed\\\\", round(ts * 100)),
    sprintf("Firm Store Card Share & %g\\%% & Observed\\\\", round(firm_store_share * 100)),
    sprintf("Consumer Bank Card Share & %g\\%% & Observed", round(u_0 * 100)),
    "\\\\ \\multicolumn{3}{c}{\\textbf{Target Moments}} \\\\",
    "Moment & Model & Data \\\\",
    "\\hline",
    sprintf("$\\Delta$ Store Card Share & %s & %s \\\\", fmt_pct(d_store), fmt_pct(d_store)),
    sprintf("$\\Delta$ Bank Card Share & %s & %s \\\\", fmt_pct(d_bank), fmt_pct(d_bank)),
    sprintf("Avg $\\Delta$ Sales ($Treated\\times Post$) & %s & %s \\\\",
            fmt_pct(mm[1]), fmt_pcts(target_moments[1], target_se[1])),
    sprintf("Accept Differential ($\\times Accept$) & %s & %s",
            fmt_pct(mm[2]), fmt_pcts(target_moments[2], target_se[2])),
    "\\end{tabular}"
  )
  writeLines(est_lines, file.path(OUTPUT_TAB, "estimation_results.tex"))
  cat("  Created estimation_results.tex\n")

  # Compact slide version: calibrated parameters (slides.tex)
  slide_lines <- c(
    "\\begin{tabular}{lc}",
    "\\toprule",
    "Parameter & Value \\\\",
    "\\midrule",
    sprintf("$\\alpha$ (price sensitivity) & $%g$ \\\\", alpha),
    sprintf("$\\sigma$ (card differentiation) & %.3f (%.3f) \\\\", params[1], delta_se[1]),
    sprintf("$\\zeta$ (cost change) & $%.3f$ (%.3f) \\\\", params[2], delta_se[2]),
    "\\midrule",
    sprintf("Implied cost reduction & %g\\%% \\\\", round(abs(cost_reduction) * 100, 1)),
    "\\bottomrule",
    "\\end{tabular}"
  )
  writeLines(slide_lines, file.path(OUTPUT_TAB, "estimation_results_slide.tex"))
  cat("  Created estimation_results_slide.tex\n")
}, error = function(e) cat("  Estimation table output failed:", e$message, "\n"))

# --- A4: Sensitivity to alpha ---
cat("  Computing sensitivity to alpha...\n")
alpha_grid <- c(-2, -3, -4, -5, -6, -0.258)
sensitivity_results <- data.frame(
  alpha = numeric(),
  sigma = numeric(),
  zeta = numeric(),
  cost_reduction = numeric(),
  price_change = numeric()
)

for (a_val in alpha_grid) {
  # Create a modified objective with different alpha
  to_opt_alpha <- function(free_params, vector = FALSE) {
    sigma_l <- c(free_params[1], 0)
    zeta_l <- free_params[2]
    c_0_l <- rep(1, 7) + a_0 * zeta_l
    c_1_l <- rep(1, 7) + a_1 * zeta_l

    p_delta_0_l <- c(rep(2, 7), c(1, rep(0, 6)))

    pq_moments_l <- function(params_l) {
      p <- params_l[1:7]
      d0 <- params_l[8:14]
      delta <- d0 + a_val * p
      s <- f.share(delta, sigma_l, g_0)$sj
      s_cond <- f.share(delta, sigma_l, g_0)$sj_hg
      pq <- s * p
      foc <- s + (p - c_0_l) * (a_val / (1 - sigma_l[1])) * s *
        (1 - sigma_l[1] * s_cond - (1 - sigma_l[1]) * s)
      return(c(pq - pq_0, foc))
    }

    sol <- nleqslv(p_delta_0_l, pq_moments_l)
    d0 <- sol$x[8:14]
    p0 <- sol$x[1:7]
    s0 <- f.share(d0 + a_val * p0, sigma_l, g_0)$sj

    p_moments_l <- function(p) {
      delta <- d0 + a_val * p
      s <- f.share(delta, sigma_l, g_1)$sj
      s_cond <- f.share(delta, sigma_l, g_1)$sj_hg
      foc <- s + (p - c_1_l) * (a_val / (1 - sigma_l[1])) * s *
        (1 - sigma_l[1] * s_cond - (1 - sigma_l[1]) * s)
      return(foc)
    }

    p1 <- nleqslv(p0, p_moments_l)$x
    s1 <- f.share(d0 + a_val * p1, sigma_l, g_1)$sj

    delta_t <- mean(s1 * p1 / (s0 * p0) - 1)
    delta_a <- mean(s1[5:6] * p1[5:6] / (s0[5:6] * p0[5:6])) -
      mean(s1[-c(5:6)] * p1[-c(5:6)] / (s0[-c(5:6)] * p0[-c(5:6)]))

    if (vector) return(c(delta_t, delta_a))
    errors <- c(delta_t, delta_a) - target_moments
    return(errors %*% errors)
  }

  P_a <- tryCatch(optim(c(0, 0), to_opt_alpha),
                  error = function(e) list(par = c(NA, NA), convergence = 1))

  if (!any(is.na(P_a$par))) {
    sensitivity_results <- rbind(sensitivity_results, data.frame(
      alpha = a_val,
      sigma = round(P_a$par[1], 4),
      zeta = round(P_a$par[2], 4),
      cost_reduction = round(P_a$par[2] * u_0, 4),
      price_change = NA
    ))
  }
}

cat("\nSensitivity to alpha:\n")
print(sensitivity_results)

# Output sensitivity table
tryCatch({
  tex_out_sa <- sensitivity_results %>%
    dplyr::select(-price_change) %>%
    kbl(format = "latex", booktabs = TRUE, digits = 4,
        col.names = c("$\\alpha$", "$\\sigma$", "$\\zeta$",
                       "Cost reduction"),
        escape = FALSE) %>%
    kable_classic() %>%
    as.character()
  tex_lines_sa <- strsplit(tex_out_sa, "\n")[[1]]
  tex_lines_sa <- tex_lines_sa[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption", tex_lines_sa)]
  writeLines(tex_lines_sa, file.path(OUTPUT_TAB, "sensitivity_alpha.tex"))
  cat("  Created sensitivity_alpha.tex\n")
}, error = function(e) cat("  Sensitivity table output failed:", e$message, "\n"))

# Counterfactual analysis
sigma <- c(params[1], 0)
zeta <- params[2]
c_0 <- rep(1, 7) + a_0 * zeta

p_delta_0 <- c(rep(2, 7), c(1, rep(0, 6)))

pq_moments <- function(params_inner) {
  p <- params_inner[1:7]
  delta_0 <- params_inner[8:14]
  delta <- delta_0 + alpha * p
  s <- f.share(delta, sigma, g_0)$sj
  s_cond <- f.share(delta, sigma, g_0)$sj_hg
  pq <- s * p
  foc <- s + (p - c_0) * (alpha / (1 - sigma[1])) * s * (1 - sigma[1] * s_cond - (1 - sigma[1]) * s)
  return(c(pq - pq_0, foc))
}

p_delta_1 <- nleqslv(p_delta_0, pq_moments)
delta_0 <- p_delta_1$x[8:14]
p_0 <- p_delta_1$x[1:7]
s_0 <- f.share(delta_0 + alpha * p_0, sigma, g_0)$sj

# Baseline values
b0 <- mean(p_0)
b1 <- mean(c_0)
b2 <- as.numeric(s_0 %*% (p_0 - c_0))
b2b <- mean(s_0 * (p_0 - c_0))
b3 <- 1 / exp(1 - sum(s_0))
b4 <- as.numeric(s_0 %*% (p_0 - c_0)) + 1 / exp(1 - sum(s_0))

# Run counterfactuals and collect results
counterfactual_results <- list()

# Helper to solve and report
solve_counterfactual <- function(a_cf, g_cf, label) {
  c_cf <- rep(1, 7) + a_cf * zeta

  p_moments_cf <- function(p) {
    delta <- delta_0 + alpha * p
    s <- f.share(delta, sigma, g_cf)$sj
    s_cond <- f.share(delta, sigma, g_cf)$sj_hg
    foc <- s + (p - c_cf) * (alpha / (1 - sigma[1])) * s * (1 - sigma[1] * s_cond - (1 - sigma[1]) * s)
    return(foc)
  }

  p_cf <- nleqslv(p_0, p_moments_cf)$x
  s_cf <- f.share(delta_0 + alpha * p_cf, sigma, g_cf)$sj

  data.frame(
    Scenario = label,
    Price_level = mean(p_cf),
    Cost_level = mean(c_cf),
    Price_change = mean(p_cf) / b0 - 1,
    Cost_change = mean(c_cf) / b1 - 1,
    Profit_change = as.numeric(s_cf %*% (p_cf - c_cf)) / b2 - 1,
    CS_change = (1 / exp(1 - sum(s_cf))) / b3 - 1,
    Welfare_change = (as.numeric(s_cf %*% (p_cf - c_cf)) + 1 / exp(1 - sum(s_cf))) / b4 - 1
  )
}

cf_results <- rbind(
  solve_counterfactual(a_1, g_1, "Post-Marquette"),
  solve_counterfactual(c(0, 0, 0, 0, 1, 1, 0) * u_1,
                       cbind(c(1, 2, 3, 4, 5, 5, 5), c(1, 1, 1, 1, 1, 1, 1)),
                       "Change nests only"),
  solve_counterfactual(c(0, 0, 0, 0, 1, 1, 1) * u_1,
                       cbind(c(1, 2, 3, 4, 5, 5, 6), c(1, 1, 1, 1, 1, 1, 1)),
                       "Change costs only"),
  solve_counterfactual(c(0, 0, 0, 0, 0, 0, 0),
                       cbind(c(1, 2, 3, 4, 5, 6, 7), c(1, 1, 1, 1, 1, 1, 1)),
                       "No acceptance"),
  solve_counterfactual(c(1, 1, 1, 1, 1, 1, 1),
                       cbind(c(1, 1, 1, 1, 1, 1, 1), c(1, 1, 1, 1, 1, 1, 1)),
                       "Homogeneous network")
)

cat("\nCounterfactual Results:\n")
print(cf_results)

# --- Decomposition table (Price / Cost / Markup) read by the paper ---
# Express changes as basis-point movements of the index level relative to the
# baseline (treated = Post-Marquette; "No cost changes" = nests only; "Only cost
# changes" = costs only).  Markup level = price level - cost level.
tryCatch({
  mk_base <- b0 - b1
  bp <- function(x) sprintf("%+dbp", round(x * 1e4))
  get_row <- function(lbl) cf_results[cf_results$Scenario == lbl, ]
  tr <- get_row("Post-Marquette")
  nc <- get_row("Change nests only")
  oc <- get_row("Change costs only")
  d_price <- function(r) r$Price_level - b0
  d_cost  <- function(r) r$Cost_level  - b1
  d_mk    <- function(r) (r$Price_level - r$Cost_level) - mk_base
  cf_lines <- c(
    "\\begin{tabular}{c|cccc}",
    "Variable & Baseline & Treated & No Cost & Only Cost \\\\",
    " & & (Marquette) & Changes & Changes \\\\",
    "\\hline",
    " & (1) & (2) & (3) & (4) \\\\",
    sprintf("Price Index & %.2f & %s & %s & %s \\\\", b0,
            bp(d_price(tr)), bp(d_price(nc)), bp(d_price(oc))),
    sprintf("Cost Index & %.2f & %s & %s & %s \\\\", b1,
            bp(d_cost(tr)), bp(d_cost(nc)), bp(d_cost(oc))),
    sprintf("Markup & %.2f & %s & %s & %s \\\\", mk_base,
            bp(d_mk(tr)), bp(d_mk(nc)), bp(d_mk(oc))),
    "\\end{tabular}"
  )
  writeLines(cf_lines, file.path(OUTPUT_TAB, "counterfactuals.tex"))
  cat("  Created counterfactuals.tex\n")
  cat(sprintf("  Decomposition: total price %s = cost %s + markup %s\n",
              bp(d_price(tr)), bp(d_cost(tr)), bp(d_mk(tr))))

  # Compact slide version: cost/markup decomposition (slides.tex)
  cf_slide <- c(
    "\\begin{tabular}{lcc}",
    "\\toprule",
    " & Cost & Markup \\\\",
    "\\midrule",
    sprintf("Baseline $\\to$ Marquette & %s & %s \\\\", bp(d_cost(tr)), bp(d_mk(tr))),
    sprintf("Cost only & %s & %s \\\\", bp(d_cost(oc)), bp(d_mk(oc))),
    sprintf("Competition only & %s & %s \\\\", bp(d_cost(nc)), bp(d_mk(nc))),
    "\\bottomrule",
    "\\end{tabular}"
  )
  writeLines(cf_slide, file.path(OUTPUT_TAB, "counterfactuals_slide.tex"))
  cat("  Created counterfactuals_slide.tex\n")
}, error = function(e) cat("  Counterfactual table output failed:", e$message, "\n"))

cat("Nested logit model complete.\n")
