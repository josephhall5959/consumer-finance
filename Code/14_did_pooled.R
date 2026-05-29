# 14_did_pooled.R — DiD using pooled CPS/PSID microdata

cat("Running pooled microdata DiD...\n")

pretty_did <- function(did, ref_yr, ylab = "Outcome", title = "", breaks = 1) {
  name <- names(summary(did)$coefficients)
  beta <- summary(did)$coefficients
  se <- summary(did)$se
  toplot <- data.frame(name, beta, se) %>%
    mutate(year = readr::parse_number(name))
  zero <- data.frame(name = "_", year = ref_yr, beta = 0, se = 0)
  toplot <- rbind(toplot, zero)
  start <- min(toplot$year)
  end <- max(toplot$year)

  ggplot(data = toplot, aes(x = year, y = beta)) +
    geom_line() +
    geom_ribbon(aes(ymin = beta - 2 * se, ymax = beta + 2 * se), alpha = 0.2) +
    theme_classic() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = ref_yr, linetype = "dashed") +
    scale_x_continuous(name = "Year", breaks = seq(start, end, breaks)) +
    scale_y_continuous(name = ylab) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
}

pooled_file <- file.path(DATA_GEN, "pooled_microdata.csv")
if (!file.exists(pooled_file)) pooled_file <- file.path(DATA_RAW, "pooled_microdata.csv")

if (!file.exists(pooled_file) || file.info(pooled_file)$size < 100) {
  cat("  Pooled microdata not available. Skipping.\n")
} else {
  pooled <- fread(pooled_file)

  # Add treatment variable if not present
  if (!"treated" %in% names(pooled) && "r_77" %in% names(pooled)) {
    pooled[, treated := ifelse(as.numeric(r_77) < 18, 1, 0)]
  }

  if ("treated" %in% names(pooled)) {
    # Bankruptcy DiD (PSID only)
    psid <- pooled[source == "PSID"]
    if (nrow(psid) > 0 && "bankrupt" %in% names(psid)) {
      tryCatch({
        est <- feols(bankrupt ~ i(year, treated, 1977) | year + state,
                     data = psid, cluster = ~state)
        p <- pretty_did(est, 1977, ylab = "Bankruptcy rate", title = "PSID: Bankruptcy")
        ggsave(file.path(OUTPUT_FIG, "bankruptcy.pdf"), p, width = 9, height = 6)
      }, error = function(e) cat("  Bankruptcy DiD failed:", e$message, "\n"))
    }

    # Self-employment DiD (CPS)
    cps <- pooled[source == "CPS"]
    if (nrow(cps) > 0 && "self.employed" %in% names(cps)) {
      tryCatch({
        est <- feols(self.employed ~ i(year, treated, 1977) | year + fipstate,
                     data = cps, cluster = ~fipstate)
        p <- pretty_did(est, 1977, ylab = "Self-employment rate", title = "CPS: Self-employment")
        ggsave(file.path(OUTPUT_FIG, "self_employment.pdf"), p, width = 9, height = 6)
      }, error = function(e) cat("  Self-employment DiD failed:", e$message, "\n"))
    }
  }
}

cat("Pooled microdata DiD complete.\n")
