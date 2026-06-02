# yp_visualizations.R — Six Yellow Pages data visualizations for review
# Generates figures and a compiled LaTeX review document

cat("Generating Yellow Pages visualizations...\n")

library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(scales)
library(sf)

if (!exists("DATA_RAW")) DATA_RAW <- file.path(getwd(), "Data", "Raw")
OUTPUT <- if (exists("OUTPUT_FIG")) file.path(OUTPUT_FIG, "yp_review") else "Output/Figures/yp_review"
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)

# Figures 1, 4 and 6 use the proprietary Dun & Bradstreet match
# (db_phonebooks_merged.csv); figures 2, 3 and 5 are built entirely from the
# public Library of Congress phonebook files. Skip the D&B panels gracefully
# when the proprietary file is absent so the public figures still build.
DB_PB <- file.path(DATA_RAW, "DB", "db_phonebooks_merged.csv")
HAVE_DB <- file.exists(DB_PB)
if (!HAVE_DB) {
  cat("  NOTE: db_phonebooks_merged.csv (proprietary D&B) not found —\n")
  cat("        skipping D&B-dependent figures 1, 4, 6.\n")
}

# ============================================================
# 1. Card Plan Market Shares
# ============================================================
cat("  1. Card plan market shares...\n")

if (HAVE_DB) {
db <- fread(file.path(DATA_RAW, "DB", "db_phonebooks_merged.csv"),
            select = c("offers_credit", "phonebook_credit_plan", "year"))
acceptors <- db[offers_credit != 0 & offers_credit != "0"]

plan_counts <- acceptors[, .N, by = phonebook_credit_plan]
# Group small plans
plan_counts[, plan := fcase(
  grepl("Master Charge", phonebook_credit_plan), "Master Charge",
  grepl("Bankamericard", phonebook_credit_plan) & !grepl("Pittsburgh", phonebook_credit_plan), "Bankamericard",
  grepl("Pittsburgh", phonebook_credit_plan), "Pittsburgh Nat'l\nBankamericard",
  grepl("NAC", phonebook_credit_plan), "NAC Charge Plan",
  default = "Other"
)]
plan_agg <- plan_counts[, .(N = sum(N)), by = plan][order(-N)]
plan_agg[, pct := N / sum(N)]
plan_agg[, label := paste0(N, " (", round(100 * pct), "%)")]

p1 <- ggplot(plan_agg, aes(x = reorder(plan, N), y = N, fill = plan)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Merchant-year listings")
ggsave(file.path(OUTPUT, "1_plan_shares.pdf"), p1, width = 8, height = 4)
cat("    Created 1_plan_shares.pdf\n")
}  # end HAVE_DB (figure 1)

# ============================================================
# 2. Card Introduction Map (state-level)
# ============================================================
cat("  2. Card introduction map...\n")

intro <- tryCatch(read_excel(file.path(DATA_RAW, "Yellow_Pages", "Introduction Dates.xlsx")),
                  error = function(e) NULL)

if (!is.null(intro) && "Card_intro" %in% names(intro)) {
  intro_dt <- as.data.table(intro)
  intro_dt[, card_year := as.numeric(substr(Card_intro, 1, 4))]

  # State-level: earliest introduction
  state_intro <- intro_dt[!is.na(card_year) & !is.na(State),
                           .(intro_year = min(card_year, na.rm = TRUE)),
                           by = State]

  # Use built-in US state map from maps package
  us_states <- map_data("state")
  state_intro[, state_lower := tolower(state.name[match(State, state.abb)])]

  map_data_merged <- merge(us_states, state_intro,
                            by.x = "region", by.y = "state_lower", all.x = TRUE)

  p2 <- ggplot(map_data_merged, aes(x = long, y = lat, group = group, fill = intro_year)) +
    geom_polygon(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(name = "Earliest card\nintroduction year",
                          na.value = "grey90",
                          breaks = seq(1958, 1972, 2)) +
    coord_map("albers", lat0 = 30, lat1 = 40) +
    theme_void(base_size = 13) +
    theme(legend.position = "bottom",
          legend.key.width = unit(2, "cm"))
  ggsave(file.path(OUTPUT, "2_intro_map.pdf"), p2, width = 9, height = 5.5)
  cat("    Created 2_intro_map.pdf\n")
} else {
  cat("    Skipping map: Introduction Dates.xlsx not available\n")
}

# ============================================================
# 3. Phonebook Coverage Over Time by State
# ============================================================
cat("  3. Phonebook coverage by state over time...\n")

pages <- fread(file.path(DATA_RAW, "Yellow_Pages", "pages_include_pagenumbers.csv"))
pages_inc <- pages[include == 1 & yellow == 1]

# Extract state from title
pages_inc[, state_name := str_extract(title, "(?<=Image 1 of )\\w+")]

books_by_year_state <- pages_inc[!is.na(year) & year >= 1950 & year <= 1980,
                                  .N, by = .(year, state_name)]

p3 <- ggplot(books_by_year_state, aes(x = year, y = N, fill = state_name)) +
  geom_col() +
  scale_fill_brewer(palette = "Set3", name = "State") +
  theme_classic(base_size = 13) +
  labs(x = "Year", y = "Number of Yellow Pages phonebooks") +
  scale_x_continuous(breaks = seq(1950, 1980, 5)) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))
ggsave(file.path(OUTPUT, "3_coverage_by_state.pdf"), p3, width = 9, height = 5)
cat("    Created 3_coverage_by_state.pdf\n")

# ============================================================
# 4. Acceptance by Industry × Firm Size
# ============================================================
cat("  4. Acceptance by industry and firm size...\n")

if (HAVE_DB) {
db2 <- fread(file.path(DATA_RAW, "DB", "db_phonebooks_merged.csv"),
             select = c("covered_by_phonebooks", "offers_credit", "sic1", "sls", "year"))
in_pb <- db2[covered_by_phonebooks != "0" & covered_by_phonebooks != ""]
in_pb[, acceptance := fifelse(offers_credit != 0 & offers_credit != "0", 1L, 0L)]
in_pb[, sls_num := as.numeric(sls)]
in_pb[, SIC1_2dig := substr(as.character(sic1), 1, 2)]
in_pb[, industry := fcase(
  SIC1_2dig == "52", "Hardware",
  SIC1_2dig == "53", "Dept Stores",
  SIC1_2dig == "54", "Food",
  SIC1_2dig == "55", "Auto",
  SIC1_2dig == "56", "Clothing",
  SIC1_2dig == "57", "Furnishing",
  SIC1_2dig == "58", "Restaurants",
  SIC1_2dig == "59", "Misc. Retail",
  default = "Non-retail"
)]

# Sales quartiles (among firms with sales data)
in_pb_sales <- in_pb[!is.na(sls_num) & sls_num > 0 & industry != "Non-retail"]
in_pb_sales[, size_q := cut(sls_num,
                              breaks = quantile(sls_num, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                              labels = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"),
                              include.lowest = TRUE)]

heat_data <- in_pb_sales[!is.na(size_q), .(
  acceptance_rate = mean(acceptance, na.rm = TRUE),
  n = .N
), by = .(industry, size_q)]

p4 <- ggplot(heat_data[n >= 10], aes(x = size_q, y = industry, fill = acceptance_rate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(round(acceptance_rate * 100, 1), "%\n(n=", n, ")")),
            size = 2.8) +
  scale_fill_gradient(low = "white", high = "steelblue",
                       name = "Acceptance\nrate",
                       labels = percent_format()) +
  theme_minimal(base_size = 13) +
  labs(x = "Sales quartile", y = NULL) +
  theme(panel.grid = element_blank())
ggsave(file.path(OUTPUT, "4_acceptance_heatmap.pdf"), p4, width = 8, height = 5)
cat("    Created 4_acceptance_heatmap.pdf\n")

rm(db2, in_pb, in_pb_sales); gc()
}  # end HAVE_DB (figure 4)

# ============================================================
# 5. Phonebook Size Trends
# ============================================================
cat("  5. Phonebook size trends...\n")

pages_inc[, pages_num := as.numeric(pages)]
size_ts <- pages_inc[!is.na(pages_num) & year >= 1950 & year <= 1980,
                      .(mean_pages = mean(pages_num, na.rm = TRUE),
                        median_pages = median(pages_num, na.rm = TRUE),
                        n_books = .N),
                      by = year]

p5 <- ggplot(size_ts, aes(x = year)) +
  geom_line(aes(y = median_pages), linewidth = 0.8) +
  geom_point(aes(y = median_pages, size = n_books), alpha = 0.6) +
  scale_size_continuous(name = "N phonebooks", range = c(1, 5)) +
  theme_classic(base_size = 13) +
  labs(x = "Year", y = "Median pages per phonebook") +
  scale_x_continuous(breaks = seq(1950, 1980, 5))
ggsave(file.path(OUTPUT, "5_phonebook_size.pdf"), p5, width = 8, height = 5)
cat("    Created 5_phonebook_size.pdf\n")

# ============================================================
# 6. Geographic Concentration of Early Adopters
# ============================================================
cat("  6. Geographic concentration of early adopters...\n")

if (HAVE_DB) {
db3 <- fread(file.path(DATA_RAW, "DB", "db_phonebooks_merged.csv"),
             select = c("offers_credit", "city", "state", "year", "sic1",
                         "covered_by_phonebooks"))
acceptors2 <- db3[offers_credit != 0 & offers_credit != "0"]

# By city: count accepting merchants and total firms in those cities
accepting_cities <- acceptors2[, .(n_accepting = .N,
                                    first_year = min(year)),
                                 by = .(city, state)]

# Get total firms in those cities
city_totals <- db3[covered_by_phonebooks != "0" & covered_by_phonebooks != "",
                    .N, by = .(city, state)]
setnames(city_totals, "N", "n_total")

city_data <- merge(accepting_cities, city_totals, by = c("city", "state"), all.x = TRUE)
city_data[, acceptance_rate := n_accepting / n_total]
city_data <- city_data[order(-n_accepting)]

# Top 15 cities bar chart
top15 <- head(city_data, 15)
top15[, city_label := paste0(city, ", ", state)]

p6 <- ggplot(top15, aes(x = reorder(city_label, n_accepting), y = n_accepting)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = paste0("(", first_year, ")")), hjust = -0.1, size = 3) +
  coord_flip() +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Card-accepting merchant-years") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave(file.path(OUTPUT, "6_geographic_concentration.pdf"), p6, width = 8, height = 5)
cat("    Created 6_geographic_concentration.pdf\n")

rm(db3, acceptors2); gc()
}  # end HAVE_DB (figure 6)

# Clean up any D&B objects that were created
rm(list = intersect(c("db", "acceptors"), ls())); gc()

if (HAVE_DB) cat("All 6 visualizations generated.\n") else cat("Public visualizations (2, 3, 5) generated; D&B figures skipped.\n")
