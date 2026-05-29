# 22_macro_model.R — General equilibrium macro model

cat("Solving macro model...\n")

sigma_m <- 2
phi <- 0.5
lbar <- 3
theta_j <- c(1, 1, 1.2)
theta_0 <- c(1, 1, 0.8)
alpha_m <- c(1, 1, 1)
z <- c(1, 1, 1)
r <- 0.05

s_0 <- rep(1, 10)

equilibrium <- function(s, theta = theta_0) {
  theta_f <- theta[1]
  theta_h <- theta[2]
  theta_b <- theta[3]
  theta_c <- theta_h * (1 - theta_b)
  theta_j_hat <- (theta_j - c(0, 0, theta_c)) * theta_f

  p <- c(1, s[1:2])
  l <- s[3:5]
  cc <- s[6:8]
  w <- s[9]
  y <- s[10]

  eq1_3 <- l - (p * z * (1 - phi) / w)^(1 / phi) * theta_j_hat
  eq4_6 <- cc - y * alpha_m * p^(-sigma_m) / sum(alpha_m * p^sigma_m)
  eq7_8 <- cc[2:3] - z[2:3] * theta_j_hat[2:3]^phi * l[2:3]^(1 - phi)
  eq9 <- lbar - sum(l)
  eq10 <- y - theta_h / (1 + r) - w * lbar - sum(p * cc)

  return(c(eq1_3, eq4_6, eq7_8, eq9, eq10))
}

Value <- nleqslv(s_0, equilibrium)$x

Variable <- c("Non-tradable price", "Card-adopting retail price",
              "Tradable labor", "Non-tradable labor", "Card-adopting retail labor",
              "Tradable consumption", "Non-tradable consumption",
              "Card-adopting retail consumption", "Nominal wage", "Nominal income")

sum_table <- data.frame(Variable, Value)

sum_table %>%
  kbl(format = "latex", booktabs = TRUE) %>%
  kable_classic() %>%
  save_kable(file = file.path(OUTPUT_TAB, "model_equilibrium.tex"),
             self_contained = FALSE)

eq_theta <- function(theta) {
  tosolve <- function(s) equilibrium(s, theta)
  return(nleqslv(s_0, tosolve)$x)
}

# Comparative statics: vary theta_f
vary_param <- function(param_name, param_seq, theta_idx) {
  results <- data.frame(
    param = param_seq,
    p_n = 0, p_c = 0,
    c_t = 0, c_n = 0, c_c = 0,
    y_t = 0, y_n = 0, y_c = 0,
    share_t = 0, share_n = 0, share_c = 0, share_cn = 0
  )

  for (i in 1:nrow(results)) {
    theta_1 <- theta_0
    theta_1[theta_idx] <- results$param[i]

    s_1 <- eq_theta(theta_1)

    results$p_n[i] <- s_1[1]
    results$p_c[i] <- s_1[2]
    results$c_t[i] <- s_1[6]
    results$c_n[i] <- s_1[7]
    results$c_c[i] <- s_1[8]

    l <- s_1[3:5]
    theta_f <- theta_1[1]
    theta_h <- theta_1[2]
    theta_b <- theta_1[3]
    theta_c_local <- theta_h * (1 - theta_b)
    theta_j_hat <- (theta_j - c(0, 0, theta_c_local)) * theta_f

    y <- z * theta_j_hat^phi * l^(1 - phi)
    results$y_t[i] <- y[1]
    results$y_n[i] <- y[2]
    results$y_c[i] <- y[3]

    p <- c(1, s_1[1:2])
    cc <- s_1[6:8]
    results$share_t[i] <- cc[1] / sum(cc * p)
    results$share_n[i] <- cc[2] / sum(cc * p)
    results$share_c[i] <- cc[3] / sum(cc * p)
    results$share_cn[i] <- results$share_c[i] / results$share_n[i]
  }

  results$price_level <- (results$c_t + results$c_n * results$p_n + results$c_c * results$p_c) /
    (results$c_t + results$c_n + results$c_c)

  return(results)
}

# Plot for each parameter
for (param_info in list(
  list(name = "theta_f", seq = seq(0.8, 1.0, 0.04), idx = 1, xlab = "Credit to firms"),
  list(name = "theta_h", seq = seq(0.8, 1.0, 0.04), idx = 2, xlab = "Credit to households"),
  list(name = "theta_b", seq = seq(0.6, 0.8, 0.04), idx = 3, xlab = "Credit intermediated by banks")
)) {
  results <- vary_param(param_info$name, param_info$seq, param_info$idx)

  pdf(file.path(OUTPUT_FIG, paste0(param_info$name, ".pdf")), height = 6, width = 9)
  print(
    ggplot(data = results, aes(x = param, y = y_t, color = "Tradable output")) +
      geom_line() +
      geom_line(aes(y = price_level, color = "Price level")) +
      geom_line(aes(y = share_cn, color = "Card share of non-tradables")) +
      scale_y_continuous(name = "Value", limits = c(0.5, 1.5)) +
      xlab(param_info$xlab) +
      theme_classic() +
      theme(legend.position = "bottom")
  )
  dev.off()
}

cat("Macro model complete.\n")
