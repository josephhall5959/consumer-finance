# US Telephone Directory OCR Corpus — Data Dictionary

This documents the machine-readable Yellow Pages corpus released with
**"Credit Cards and Retail Firms: Historical Evidence from the United States"**
(Joseph P. Hall, Georgia Institute of Technology).

## Provenance & license

- **Source scans:** U.S. Telephone Directories collection, Library of Congress
  (<https://www.loc.gov/collections/united-states-telephone-directory-collection/>).
  The scans are in the **public domain**.
- **OCR processing:** performed by Hi-Tech Inc. on the Library of Congress scans
  (statements of work in `01-MSA - Hitech INC.pdf`, `02-SOW - Hitech INC*.pdf`).
- **Compilation & code:** Joseph P. Hall. Released under **CC-BY-4.0** — reuse
  freely with attribution.
- **Not included:** the proprietary Dun & Bradstreet establishment microdata used
  to attach sales to listings. That data must be licensed separately and is never
  redistributed in this package.

## How to obtain

Run `bash fetch_data.sh` from the repository root. It downloads the corpus from
Zenodo and unpacks it into `Data/Raw/Yellow_Pages/OCR/`. Override the location
with the `YP_OCR_DIR` environment variable if you store it elsewhere.

## Layout

```
Data/Raw/Yellow_Pages/OCR/
├── Merged Files/            # 110 merged, book-level CSVs (one file per phonebook)
└── CSV_MANUAL/CSV_MANUAL/   # page-level CSVs, one subdirectory per phonebook
```

There are **two layers** with different schemas.

### Layer 1 — Merged book-level files (`Merged Files/*.csv`)

The analytically richest layer: one row per OCR'd text token, with bounding
boxes and automated category-header detection.

| Column | Type | Description |
|--------|------|-------------|
| `index` | int | Row index within the book |
| `set` | str | OCR batch / set identifier |
| `state` | str | State of the phonebook (from LoC metadata) |
| `city_county` | str | City or county the directory covers |
| `year` | int | Directory year |
| `page` | int | Page number the token appears on |
| `is_title` | 0/1 | 1 if the token was detected as a Yellow Pages category header |
| `title_corrected` | str | Cleaned category-header text (when `is_title == 1`) |
| `has.number` | 0/1 | 1 if the entry contains a telephone number (proxies a business listing) |
| `text` | str | Raw OCR text of the token/entry |
| `body.height` | num | OCR'd glyph height (proxy for font size / display prominence) |

Counts: ~3.7M tokens across 110 books, 13,548 pages, 607,344 listings
(`has.number == 1`), 17,995 distinct `title_corrected` category headers.

### Layer 2 — Page-level files (`CSV_MANUAL/CSV_MANUAL/<book>/*.csv`)

Raw OCR output only (no category-header detection). One CSV per scanned page.

| Column | Type | Description |
|--------|------|-------------|
| `state` | str | State of the phonebook |
| `city_county` | str | City or county covered |
| `year` | int | Directory year |
| `page` | int | Page number |
| `has.number` | 0/1 | 1 if the entry contains a telephone number |

Combining both layers (de-duplicating books that appear in both, preferring the
merged file) gives the full corpus: **168 phonebooks, 5 states, 44 cities,
1950–1978, ~2.1M business listings, 49,447 unique category headers, 13 books
with a detected "Credit Card & Other Credit Plans" section.**

## Companion manual files (already in `Data/Raw/Yellow_Pages/`)

These small hand-curated files are tracked in git directly:

| File | Contents |
|------|----------|
| `Merchant List.xlsx` | The 510 hand-transcribed card-accepting merchant-year records |
| `Introduction Dates.xlsx` | Earliest observed card listing per geography |
| `pages_include_pagenumbers.csv` | Phonebook inclusion flags + page counts |
| `index_manual*.csv` | Manually transcribed directory index (covered municipalities) |
| `titles_list.txt` | Flat list of detected category headers |

## Which scripts consume what

| Script | Reads | Public-only? |
|--------|-------|--------------|
| `Code/yp_ocr_analysis.R` | `Merged Files/` | ✅ yes |
| `Code/yp_ocr_full_analysis.R` | `Merged Files/` + `CSV_MANUAL/` | ✅ yes |
| `Code/yp_visualizations.R` | manual files (figs 2,3,5) + D&B match (figs 1,4,6) | partial |
| `Code/16_analyze_phonebooks.R` | `db_phonebooks_merged.csv` (D&B) | ❌ needs D&B |

The OCR-corpus scripts skip gracefully (exit 0 with a NOTE) if the corpus is not
present, so the rest of the pipeline still runs.
