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

target_moments <- c(0.0027, -0.048)
# Standard errors of target moments (from reduced-form regressions)
target_se <- c(0.005, 0.02)

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

# --- A3: Bootstrap standard errors ---
cat("  Computing bootstrap standard errors...\n")
n_boot <- 50  # Increase to 1000 for final version
boot_params <- matrix(NA, nrow = n_boot, ncol = 2)

# Create a version of the objective that takes target_moments as argument
to_optimize_tm <- function(free_params, tm) {
  result <- to_optimize(free_params, vector = TRUE)
  errors <- result - tm
  return(errors %*% errors)
}

set.seed(42)
for (b in 1:n_boot) {
  # Perturb target moments by their standard errors
  tm_b <- target_moments + rnorm(2) * target_se

  P_b <- tryCatch(
    optim(params, to_optimize_tm, tm = tm_b),
    error = function(e) list(par = c(NA, NA))
  )
  boot_params[b, ] <- P_b$par
}

boot_se <- apply(boot_params, 2, sd, na.rm = TRUE)
cat("  Bootstrap SE: sigma =", boot_se[1], ", zeta =", boot_se[2], "\n")
cat("  Bootstrap SE of cost:", boot_se[2] * 0.37, "\n")

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

cat("Nested logit model complete.\n")
