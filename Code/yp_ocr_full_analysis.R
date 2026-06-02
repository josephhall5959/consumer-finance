# yp_ocr_full_analysis.R — Analyze the FULL OCR corpus (merged + page-level)
# 110 merged book-level CSVs + 58 additional page-level books = 165 total

cat("Analyzing full Yellow Pages OCR corpus...\n")

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(mapproj)

# Public OCR corpus location. The corpus is distributed via Zenodo (see the
# repo README / fetch_data.sh) and extracts into Data/Raw/Yellow_Pages/OCR/,
# preserving the "Merged Files/" and "CSV_MANUAL/CSV_MANUAL/" subfolders.
# Override the root with the YP_OCR_DIR environment variable if needed.
.data_raw <- if (exists("DATA_RAW")) DATA_RAW else file.path(getwd(), "Data", "Raw")
.ocr_root <- Sys.getenv("YP_OCR_DIR",
                        unset = file.path(.data_raw, "Yellow_Pages", "OCR"))
OCR_MERGED <- file.path(.ocr_root, "Merged Files")
OCR_PAGES  <- file.path(.ocr_root, "CSV_MANUAL", "CSV_MANUAL")
OUTPUT     <- if (exists("OUTPUT_FIG")) file.path(OUTPUT_FIG, "yp_review") else "Output/Figures/yp_review"
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(OCR_MERGED) && !dir.exists(OCR_PAGES)) {
  cat(sprintf(paste0(
    "  NOTE: OCR corpus not found (looked in %s).\n",
    "  Run `bash fetch_data.sh` to download it from Zenodo, then re-run.\n",
    "  Skipping full-corpus analysis (committed ocr_full_summary.tex is retained).\n"),
    .ocr_root))
  quit(save = "no", status = 0)
}

# ============================================================
# Step 1: Scan merged books (110)
# ============================================================
cat("  Scanning 110 merged books...\n")
merged_files <- list.files(OCR_MERGED, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)

scan_merged <- function(f) {
  tryCatch({
    d <- fread(f, select = c("state", "city_county", "year", "page",
                              "is_title", "title_corrected", "has.number"),
               showProgress = FALSE)
    book_id <- gsub(".*book_([0-9]+).*", "\\1", basename(f))
    titles <- d[is_title == 1 & !is.na(title_corrected) & nchar(title_corrected) > 1]
    cats <- unique(gsub("\\s*\\(Cont'd\\)$|\\s*-\\s*$|\\s*\\(Cont.d\\)$", "",
                        titles$title_corrected, ignore.case = TRUE))
    has_credit <- any(grepl("credit card|charge plan", cats, ignore.case = TRUE))
    data.table(book_id = book_id, source = "merged", state = d$state[1],
               city = d$city_county[1], year = as.integer(d$year[1]),
               n_pages = length(unique(d$page)),
               n_listings = sum(d$has.number == 1, na.rm = TRUE),
               n_categories = length(cats), has_credit_section = has_credit)
  }, error = function(e) NULL)
}
merged_summ <- rbindlist(lapply(merged_files, scan_merged))
cat(sprintf("    Done: %d books\n", nrow(merged_summ)))

# ============================================================
# Step 2: Scan page-level books (165, includes overlap)
# ============================================================
cat("  Scanning page-level books...\n")
page_dirs <- list.dirs(OCR_PAGES, recursive = FALSE)

scan_pages <- function(dir_path) {
  tryCatch({
    bid <- basename(dir_path)
    csvs <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)
    csvs <- csvs[file.size(csvs) > 0]
    if (length(csvs) == 0) return(NULL)

    # Page-level CSVs lack is_title/title_corrected — only raw OCR
    # Read a sample (first + every 10th) for efficiency
    sample_idx <- unique(c(1, seq(1, length(csvs), by = 10), length(csvs)))
    d <- rbindlist(lapply(csvs[sample_idx], function(f) {
      tryCatch(fread(f, select = c("state", "city_county", "year", "page", "has.number"),
                     showProgress = FALSE),
               error = function(e) NULL)
    }), fill = TRUE)

    if (is.null(d) || nrow(d) == 0) return(NULL)

    scale_factor <- length(csvs) / length(sample_idx)
    yr <- suppressWarnings(as.integer(d$year[1]))
    data.table(book_id = bid, source = "page_level", state = d$state[1],
               city = d$city_county[1], year = yr,
               n_pages = length(csvs),
               n_listings = as.integer(sum(d$has.number == 1, na.rm = TRUE) * scale_factor),
               n_categories = NA_integer_, has_credit_section = NA)
  }, error = function(e) { cat("    Error:", bid, e$message, "\n"); NULL })
}
page_summ <- rbindlist(lapply(page_dirs, scan_pages))
cat(sprintf("    Done: %d books\n", nrow(page_summ)))

# ============================================================
# Step 3: Combine, dedup by book_id (prefer merged if available)
# ============================================================
# Normalize book_ids: merged files have numeric IDs, page-level have "book_NNNNN"
merged_summ[, book_id_norm := gsub("^0+", "", book_id)]
if (nrow(page_summ) > 0) {
  page_summ[, book_id_norm := gsub("^0+", "", gsub("book_", "", book_id))]
  # Take page-level books NOT in merged
  new_books <- page_summ[!book_id_norm %in% merged_summ$book_id_norm]
} else {
  new_books <- data.table()
}
cat(sprintf("  New books from page-level: %d\n", nrow(new_books)))

all_books <- rbind(merged_summ, new_books, fill = TRUE)
all_books <- all_books[!is.na(year) & year > 1900]
setorder(all_books, state, city, year)

cat(sprintf("\n  === FULL CORPUS ===\n"))
cat(sprintf("  Total phonebooks: %d\n", nrow(all_books)))
cat(sprintf("  States: %d (%s)\n", length(unique(all_books$state)),
            paste(unique(all_books$state[!is.na(all_books$state)]), collapse = ", ")))
cat(sprintf("  Cities: %d\n", length(unique(all_books$city))))
cat(sprintf("  Year range: %d-%d\n", min(all_books$year), max(all_books$year)))
cat(sprintf("  Total pages: %s\n", format(sum(all_books$n_pages), big.mark = ",")))
cat(sprintf("  Total listings: %s\n", format(sum(all_books$n_listings), big.mark = ",")))
cat(sprintf("  Books with credit card sections: %d\n", sum(all_books$has_credit_section, na.rm = TRUE)))

# ============================================================
# Figure A: Full corpus — books by state and year
# ============================================================
cat("\n  Creating figures...\n")

p_a <- ggplot(all_books, aes(x = factor(year), fill = state)) +
  geom_bar() +
  theme_classic(base_size = 12) +
  labs(x = "Year", y = "Phonebooks in OCR corpus", fill = "State") +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(OUTPUT, "A_full_corpus_by_year.pdf"), p_a, width = 10, height = 5)
cat("    A_full_corpus_by_year.pdf\n")

# ============================================================
# Figure B: Geographic coverage — city × year panel
# ============================================================
# Create a city-year indicator matrix
city_year <- all_books[, .(n_books = .N, has_credit = any(has_credit_section, na.rm = TRUE)),
                        by = .(state, city, year)]
city_year[, city_label := paste0(city, ", ", state)]

# Order cities by state then first year
city_order <- city_year[, .(first_year = min(year)), by = city_label][order(first_year)]
city_year[, city_label := factor(city_label, levels = rev(city_order$city_label))]

p_b <- ggplot(city_year, aes(x = year, y = city_label)) +
  geom_tile(aes(fill = has_credit), color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                     labels = c("No credit section", "Has credit card section"),
                     name = NULL) +
  theme_minimal(base_size = 10) +
  labs(x = "Year", y = NULL) +
  scale_x_continuous(breaks = seq(1950, 1980, 5)) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 6))
ggsave(file.path(OUTPUT, "B_city_year_coverage.pdf"), p_b, width = 10, height = 12)
cat("    B_city_year_coverage.pdf\n")

# ============================================================
# Figure C: California deep dive — BankAmericard's home state
# ============================================================
ca_books <- all_books[state == "California"]
if (nrow(ca_books) > 0) {
  ca_city_year <- ca_books[, .(n_books = .N, has_credit = any(has_credit_section, na.rm = TRUE),
                                n_listings = sum(n_listings)),
                            by = .(city, year)]

  p_c <- ggplot(ca_city_year, aes(x = year, y = reorder(city, year, FUN = min))) +
    geom_point(aes(size = n_listings, color = has_credit), alpha = 0.7) +
    scale_size_continuous(name = "Business\nlistings", range = c(1, 8),
                           labels = comma) +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                        labels = c("No credit section", "Credit card section"),
                        name = NULL) +
    theme_classic(base_size = 12) +
    labs(x = "Year", y = NULL,
         title = "California: BankAmericard's Home State") +
    scale_x_continuous(breaks = seq(1950, 1980, 5)) +
    theme(legend.position = "right")
  ggsave(file.path(OUTPUT, "C_california_deep_dive.pdf"), p_c, width = 10, height = 7)
  cat("    C_california_deep_dive.pdf\n")
}

# ============================================================
# Figure D: When did credit card sections first appear?
# ============================================================
credit_intro <- all_books[has_credit_section == TRUE,
                           .(first_credit_year = min(year)),
                           by = .(state, city)]
credit_intro <- credit_intro[order(first_credit_year)]

if (nrow(credit_intro) > 0) {
  credit_intro[, city_label := paste0(city, ", ", state)]
  p_d <- ggplot(credit_intro, aes(x = first_credit_year,
                                    y = reorder(city_label, -first_credit_year))) +
    geom_point(size = 3, color = "firebrick") +
    geom_segment(aes(xend = first_credit_year - 0.3, yend = city_label),
                 color = "firebrick", linewidth = 0.8) +
    theme_classic(base_size = 12) +
    labs(x = "Year of first observed credit card section",
         y = NULL) +
    scale_x_continuous(breaks = seq(1955, 1975, 2))
  ggsave(file.path(OUTPUT, "D_credit_section_intro.pdf"), p_d, width = 8, height = 6)
  cat("    D_credit_section_intro.pdf\n")
}

# ============================================================
# Figure E: Listings density — comparing cities
# ============================================================
city_stats <- all_books[, .(
  n_books = .N,
  total_listings = sum(n_listings),
  avg_listings = mean(n_listings),
  total_pages = sum(n_pages),
  year_range = paste0(min(year), "-", max(year))
), by = .(state, city)]
city_stats[, city_label := paste0(city, ", ", state)]

p_e <- ggplot(city_stats[n_books >= 2],
              aes(x = reorder(city_label, avg_listings), y = avg_listings)) +
  geom_col(aes(fill = state), alpha = 0.8) +
  geom_text(aes(label = n_books), hjust = -0.3, size = 3) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", name = "State") +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Avg. business listings per phonebook",
       caption = "Numbers show count of phonebooks per city") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
ggsave(file.path(OUTPUT, "E_listings_by_city.pdf"), p_e, width = 9, height = 8)
cat("    E_listings_by_city.pdf\n")

# ============================================================
# Figure F: Category richness over time
# ============================================================
p_f <- ggplot(all_books[!is.na(n_categories)], aes(x = year, y = n_categories)) +
  geom_point(aes(color = state), alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", color = "black", linewidth = 0.8, se = TRUE) +
  theme_classic(base_size = 13) +
  labs(x = "Year", y = "Unique business categories per phonebook",
       color = "State") +
  scale_color_brewer(palette = "Set2")
ggsave(file.path(OUTPUT, "F_categories_over_time.pdf"), p_f, width = 8, height = 5)
cat("    F_categories_over_time.pdf\n")

# ============================================================
# Save full corpus summary
# ============================================================
fwrite(all_books, file.path(OUTPUT, "full_corpus_summary.csv"))

# Updated summary table
ocr_full_summary <- data.frame(
  Metric = c(
    "Phonebooks in OCR corpus",
    "--- from merged book-level files",
    "--- from page-level files (new)",
    "Total pages scanned",
    "Business listings (with phone numbers)",
    "Unique business categories",
    "Books with credit card sections",
    "States covered",
    "Cities covered",
    "Year range"
  ),
  Value = c(
    nrow(all_books),
    nrow(merged_summ[!is.na(year) & year > 1900]),
    nrow(new_books[!is.na(year) & year > 1900]),
    format(sum(all_books$n_pages), big.mark = ","),
    format(sum(all_books$n_listings), big.mark = ","),
    format(sum(all_books$n_categories, na.rm = TRUE), big.mark = ","),
    sum(all_books$has_credit_section, na.rm = TRUE),
    length(unique(all_books$state[!is.na(all_books$state)])),
    length(unique(all_books$city)),
    paste0(min(all_books$year), "-", max(all_books$year))
  )
)

library(kableExtra)
tex_out <- ocr_full_summary %>%
  kbl(format = "latex", booktabs = TRUE, col.names = c("", "")) %>%
  kable_classic() %>%
  as.character()
tex_lines <- strsplit(tex_out, "\n")[[1]]
tex_lines <- tex_lines[!grepl("^\\\\begin\\{table|^\\\\end\\{table|^\\\\centering$|^\\\\caption",
                               tex_lines)]
writeLines(tex_lines, file.path(OUTPUT, "ocr_full_summary.tex"))
cat("  Created ocr_full_summary.tex\n")

cat("\nFull corpus analysis complete.\n")
