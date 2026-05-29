# 15_analyze_scranton.R — Scranton, PA credit card adoption analysis

cat("Analyzing Scranton credit card adoption...\n")

scranton_file <- file.path(DATA_RAW, "scranton_71.csv")
if (!file.exists(scranton_file)) {
  stop("scranton_71.csv not found at: ", scranton_file)
}

scranton_71 <- fread(scranton_file)
cat(sprintf("  Loaded %d firms\n", nrow(scranton_71)))

# The 'card' column indicates credit card acceptance (1=accepts, 0=does not)
if (!"card" %in% names(scranton_71)) {
  stop("scranton_71.csv missing 'card' column")
}

scranton_71[, SIC1_2dig := substr(as.character(SIC1), 1, 2)]
scranton_71[, industry := fcase(
  SIC1_2dig == "52", "Hardware",
  SIC1_2dig == "53", "Department Stores",
  SIC1_2dig == "54", "Food",
  SIC1_2dig == "55", "Auto",
  SIC1_2dig == "56", "Clothing",
  SIC1_2dig == "57", "Furnishing",
  SIC1_2dig == "58", "Restaurants",
  SIC1_2dig == "59", "Misc."
)]

format_pct <- function(x) paste0(round(100 * x, 0), "%")

# Table: Adoption by industry
adoption_ind <- scranton_71[!is.na(industry), .(
  acceptance = format_pct(mean(card, na.rm = TRUE)),
  n = .N
), by = industry][order(-acceptance)]

adoption_ind %>%
  kbl(format = "latex", booktabs = TRUE,
      caption = "Credit Card Acceptance by Industry, Scranton 1971") %>%
  kable_classic() %>%
  save_kable(file = file.path(OUTPUT_TAB, "adoption_industries.tex"),
             self_contained = FALSE)
cat("  Created adoption_industries.tex\n")

# Regression: acceptance on firm characteristics (logit with marginal effects)
tryCatch({
  logit1 <- feglm(card ~ log_sls + single | SIC1_2dig,
                  data = scranton_71[SLSCODE == 0],
                  family = binomial(link = "logit"))
  # Average marginal effects (vcov=FALSE for FE models)
  ame1 <- tryCatch(
    marginaleffects::avg_slopes(logit1, vcov = FALSE),
    error = function(e2) { cat("  AME computation note:", e2$message, "\n"); NULL }
  )
  if (!is.null(ame1)) {
    cat("  Scranton logit AME:\n")
    print(ame1)
  }
}, error = function(e) {
  cat("  Logit failed, falling back to OLS:", e$message, "\n")
  logit1 <- feols(card ~ log_sls + single | SIC1_2dig,
                  data = scranton_71[SLSCODE == 0])
})

adoption_dict <- c(log_sls = "log(Sales)",
                    single = "Single-establishment",
                    age = "Firm age",
                    acceptance = "Card acceptance",
                    fe = "Industry FE")
tex_out <- capture.output(etable(logit1, tex = TRUE, dict = adoption_dict))
writeLines(tex_out, file.path(OUTPUT_TAB, "adoption_firms.tex"))
cat("  Created adoption_firms.tex\n")

rm(scranton_71); gc()
cat("Scranton analysis complete.\n")
