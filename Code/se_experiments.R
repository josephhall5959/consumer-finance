# se_experiments.R — explore avenues to tighten SE(sigma) / SE(cost change)
# Each avenue re-estimates the reduced-form sales DiD, recomputes the 2 SMM
# target moments + their 2x2 covariance, then propagates through the structural
# delta method to report SE(sigma), SE(zeta), SE(cost).
# Run standalone:  Rscript Code/se_experiments.R

suppressMessages(source("Code/00_setup.R"))
suppressMessages({
  library(data.table); library(fixest)
  library(marginaleffects)
})
source(file.path(ROOT, "Code", "helpers", "nlfunctions.R"))

# ----------------------------------------------------------------------------
# 1. Structural delta-method machinery (alpha = -4, same calibration as script 21)
# ----------------------------------------------------------------------------
alpha <- -4
oo <- 0.92; nf <- 7; ts <- 0.47; na_firms <- 2
as_share <- 0.17 / na_firms; ns <- 1; ss <- 0.20 / ns
no <- nf - 1 - na_firms - ns
pq_0 <- c(ts * (1 - oo),
          rep((1 - ts - ss * ns - as_share * na_firms) * (1 - oo) / no, no),
          rep(as_share * (1 - oo), na_firms),
          rep(ss * (1 - oo), ns))
u_0 <- 0.37; u_1 <- u_0 + 0.055
a_0 <- c(0,0,0,0,1,1,0) * u_0
g_0 <- cbind(c(1,2,3,4,5,5,6), rep(1,7))
a_1 <- c(0,0,0,0,1,1,1) * u_1
g_1 <- cbind(c(1,2,3,4,5,5,5), rep(1,7))

model_moments <- function(theta) {
  sigma <- c(theta[1], 0); zeta <- theta[2]
  c_0 <- rep(1,7) + a_0 * zeta
  c_1 <- rep(1,7) + a_1 * zeta
  p_delta_0 <- c(rep(2,7), c(1, rep(0,6)))
  pq_moments <- function(pp) {
    p <- pp[1:7]; d0 <- pp[8:14]; delta <- d0 + alpha * p
    s <- f.share(delta, sigma, g_0)$sj; sc <- f.share(delta, sigma, g_0)$sj_hg
    foc <- s + (p - c_0) * (alpha/(1-sigma[1])) * s * (1 - sigma[1]*sc - (1-sigma[1])*s)
    c(s*p - pq_0, foc)
  }
  sol <- nleqslv(p_delta_0, pq_moments)
  d0 <- sol$x[8:14]; p_0 <- sol$x[1:7]
  s_0 <- f.share(d0 + alpha*p_0, sigma, g_0)$sj
  p_moments <- function(p) {
    delta <- d0 + alpha*p
    s <- f.share(delta, sigma, g_1)$sj; sc <- f.share(delta, sigma, g_1)$sj_hg
    s + (p - c_1) * (alpha/(1-sigma[1])) * s * (1 - sigma[1]*sc - (1-sigma[1])*s)
  }
  p_1 <- nleqslv(p_0, p_moments)$x
  s_1 <- f.share(d0 + alpha*p_1, sigma, g_1)$sj
  dt <- mean(s_1*p_1/(s_0*p_0) - 1)
  da <- mean(s_1[5:6]*p_1[5:6]/(s_0[5:6]*p_0[5:6])) -
        mean(s_1[-c(5:6)]*p_1[-c(5:6)]/(s_0[-c(5:6)]*p_0[-c(5:6)]))
  c(dt, da)
}

structural_se <- function(tm, Sigma_m) {
  obj <- function(fp) { e <- model_moments(fp) - tm; as.numeric(e %*% e) }
  P <- optim(c(0,0), obj)
  th <- P$par
  h <- 1e-5; K <- matrix(0,2,2)
  for (j in 1:2) {
    tp <- th; tmn <- th; tp[j] <- tp[j]+h; tmn[j] <- tmn[j]-h
    K[,j] <- (model_moments(tp) - model_moments(tmn)) / (2*h)
  }
  J <- solve(K); V <- J %*% Sigma_m %*% t(J)
  se <- sqrt(diag(V))
  list(sigma=th[1], zeta=th[2], cost=th[2]*u_0,
       se_sigma=se[1], se_zeta=se[2], se_cost=se[2]*u_0,
       t_sigma=th[1]/se[1], t_cost=(th[2]*u_0)/(se[2]*u_0))
}

# Map regression coefs -> (moment1 = avg effect, moment2 = differential) for the
# binary treated_post*card spec, returning moments + full 2x2 covariance.
moments_from_binary <- function(m, p_card) {
  b <- coef(m); V <- vcov(m); nm <- names(b)
  i_tp  <- which(nm == "treated_post")
  i_tpc <- which(nm == "treated_post:card")
  L <- matrix(0, 2, length(b))
  L[1, i_tp] <- 1; L[1, i_tpc] <- p_card
  L[2, i_tpc] <- 1
  list(m = as.numeric(L %*% b), S = L %*% V %*% t(L),
       b_tp = b[i_tp], b_tpc = b[i_tpc])
}

# ----------------------------------------------------------------------------
# 2. Load data once
# ----------------------------------------------------------------------------
cat("Loading retail.csv ...\n")
rs <- fread(file.path(DATA_GEN, "retail.csv"),
            select = c("STATE","YEAR","SIC1","DUNSNO","log_sls","card",
                       "treated","r_77","log_emp","single","SLS"))
rs <- rs[!is.na(log_sls) & is.finite(log_sls) & !is.na(treated)]
rs[, post := as.integer(YEAR >= 1978)]
rs[, treated_post := treated * post]
rs[, gap := pmax(0, 18 - r_77)]          # usury intensity (binding tightness)
rs[, gap_post := gap * post]
cat(sprintf("  N = %d firm-years; card share (treated-post) = %.3f\n",
            nrow(rs), rs[treated_post==1, mean(card)]))

res_file <- file.path(DATA_GEN, "se_experiments.csv")
if (file.exists(res_file)) file.remove(res_file)
results <- list()
add_res <- function(label, mo) {
  r <- structural_se(mo$m, mo$S)
  row <- data.frame(
    avenue = label,
    moment1 = mo$m[1], se_m1 = sqrt(mo$S[1,1]),
    moment2 = mo$m[2], se_m2 = sqrt(mo$S[2,2]),
    sigma = r$sigma, se_sigma = r$se_sigma, t_sigma = r$t_sigma,
    cost = r$cost, se_cost = r$se_cost, t_cost = r$t_cost
  )
  results[[label]] <<- row
  fwrite(row, res_file, append = file.exists(res_file))   # incremental save
  cat(sprintf("[%s] m1=%.4f(%.4f) m2=%.4f(%.4f) | sigma=%.4f(se %.4f,t %.2f) cost=%.4f(se %.4f,t %.2f)\n",
              label, mo$m[1], sqrt(mo$S[1,1]), mo$m[2], sqrt(mo$S[2,2]),
              r$sigma, r$se_sigma, r$t_sigma, r$cost, r$se_cost, r$t_cost))
}

# ----------------------------------------------------------------------------
# 3. VERIFY the averaging concern
# ----------------------------------------------------------------------------
cat("\n=== Verify moment 1 = AVERAGE partial effect (not non-accepter effect) ===\n")
m0 <- feols(log_sls ~ treated_post * card | DUNSNO + YEAR, data = rs, cluster = ~STATE)
p_card0 <- rs[treated_post==1, mean(card)]
mo0 <- moments_from_binary(m0, p_card0)
cat(sprintf("  b(treated_post) [card==0, NON-accepters] = %.4f\n", mo0$b_tp))
cat(sprintf("  b(treated_post:card) [differential]      = %.4f\n", mo0$b_tpc))
cat(sprintf("  hand-built avg = b_tp + p_card*b_tpc      = %.4f\n", mo0$m[1]))
# Independent cross-check WITHOUT per-row marginaleffects (too heavy on 13M rows):
# refit with card centered so the treated_post main effect *is* the avg effect.
rs[, card_c := card - p_card0]
mchk <- feols(log_sls ~ treated_post * card_c | DUNSNO + YEAR, data = rs, cluster = ~STATE)
cat(sprintf("  cross-check (card centered) b(treated_post) = %.4f (se %.4f)\n",
            coef(mchk)["treated_post"], se(mchk)["treated_post"]))
rs[, card_c := NULL]; rm(mchk); gc()

# ----------------------------------------------------------------------------
# 4. AVENUES  (one model at a time; drop + gc after each to stay within RAM)
# ----------------------------------------------------------------------------
cat("\n=== Avenues ===\n")
setFixest_nthreads(2)

# V0 baseline
add_res("V0_baseline", mo0); rm(m0); gc()

# V1 add SIC1 x YEAR FE (industry-year shocks)
m <- feols(log_sls ~ treated_post * card | DUNSNO + YEAR + SIC1^YEAR,
           data = rs, cluster = ~STATE)
add_res("V1_sicXyear_FE", moments_from_binary(m, p_card0)); rm(m); gc()

# V2 two-way cluster STATE + SIC1
m <- feols(log_sls ~ treated_post * card | DUNSNO + YEAR, data = rs,
           cluster = ~STATE + SIC1)
add_res("V2_cluster_state_sic", moments_from_binary(m, p_card0)); rm(m); gc()

# V3 two-way cluster STATE + YEAR
m <- feols(log_sls ~ treated_post * card | DUNSNO + YEAR, data = rs,
           cluster = ~STATE + YEAR)
add_res("V3_cluster_state_year", moments_from_binary(m, p_card0)); rm(m); gc()

# V4 trim log_sls to 1-99 pct (heavy tails)
qlo <- rs[, quantile(log_sls, 0.01)]; qhi <- rs[, quantile(log_sls, 0.99)]
keep <- rs$log_sls >= qlo & rs$log_sls <= qhi
pc_trim <- rs[keep][treated_post==1, mean(card)]
m <- feols(log_sls ~ treated_post * card | DUNSNO + YEAR, data = rs[keep],
           cluster = ~STATE)
add_res("V4_trim_1_99", moments_from_binary(m, pc_trim)); rm(m); gc()

# V5 weight by baseline firm size (downweight tiny noisy firms)
rs[, w := as.numeric(SLS)]; rs[!is.finite(w) | w <= 0, w := NA]
pc_w <- rs[treated_post==1 & !is.na(w), weighted.mean(card, w)]
m <- feols(log_sls ~ treated_post * card | DUNSNO + YEAR,
           data = rs, weights = ~w, cluster = ~STATE)
add_res("V5_weight_sales", moments_from_binary(m, pc_w)); rm(m); gc()

# V6 continuous treatment intensity (gap), rescaled to mean treated gap
m <- feols(log_sls ~ gap_post * card | DUNSNO + YEAR, data = rs, cluster = ~STATE)
{
  b <- coef(m); V <- vcov(m); nm <- names(b)
  i_g  <- which(nm == "gap_post"); i_gc <- which(nm == "gap_post:card")
  mean_gap <- rs[treated==1, mean(gap)]            # avg binding gap among treated
  p_card_g <- rs[treated_post==1, mean(card)]
  L <- matrix(0, 2, length(b))
  L[1, i_g] <- mean_gap; L[1, i_gc] <- mean_gap * p_card_g
  L[2, i_gc] <- mean_gap
  add_res("V6_continuous_gap", list(m = as.numeric(L %*% b), S = L %*% V %*% t(L)))
}
rm(m); gc()

# V7 add log_emp + single controls (residual variance)
m <- feols(log_sls ~ treated_post * card + log_emp + single | DUNSNO + YEAR,
           data = rs, cluster = ~STATE)
add_res("V7_firm_controls", moments_from_binary(m, p_card0)); rm(m); gc()

# V8 combine structural FE + controls + trim
pc8 <- rs[keep][treated_post==1, mean(card)]
m <- feols(log_sls ~ treated_post * card + log_emp + single | DUNSNO + YEAR + SIC1^YEAR,
           data = rs[keep], cluster = ~STATE)
add_res("V8_combo", moments_from_binary(m, pc8)); rm(m); gc()

# ----------------------------------------------------------------------------
# 5. Report
# ----------------------------------------------------------------------------
tab <- rbindlist(results)
cat("\n================ SUMMARY (sorted by SE(cost)) ================\n")
print(tab[order(se_cost),
          .(avenue, sigma=round(sigma,4), se_sigma=round(se_sigma,4), t_sigma=round(t_sigma,2),
            cost=round(cost,4), se_cost=round(se_cost,4), t_cost=round(t_cost,2),
            se_m1=round(se_m1,4), se_m2=round(se_m2,4))])
fwrite(tab, file.path(DATA_GEN, "se_experiments.csv"))
cat("\nSaved Data/Generated/se_experiments.csv\n")
