# CLAUDE.md — Guide for Claude Code (and other agents)

This file orients an AI coding assistant working in this replication package.
If you are a researcher who just cloned this repo and asked Claude to "reproduce
this paper" or "build on this data," start here.

## What this repo is

Replication package for **"Credit Cards and Retail Firms: Historical Evidence
from the United States"** by Joseph P. Hall (Georgia Institute of Technology).
It reproduces every figure, table, the paper PDF, and the slide deck from raw
data via an R pipeline.

## One-prompt quickstart

```bash
# 1. (Optional but recommended) download the public Yellow Pages OCR corpus
bash fetch_data.sh

# 2. Run the whole thing: installs R packages, runs the pipeline, builds PDFs
bash run.sh
```

Outputs land in `Output/Figures/`, `Output/Tables/`, `Paper/main.pdf`, and
`Paper/slides.pdf`.

## Data availability — read this before reproducing

Datasets fall into three buckets. **The pipeline degrades gracefully**: scripts
whose inputs are missing print a NOTE and skip rather than halting the run.

| Bucket | What | How to get it |
|--------|------|---------------|
| ✅ **Public, shipped here** | Author compilations (usury rates), small generated CSVs, all `Output/` artifacts | already in the repo |
| ✅ **Public, fetched** | Yellow Pages / US Telephone Directory OCR corpus (Library of Congress, public domain) | `bash fetch_data.sh` (Zenodo) |
| ⬇️ **Public, download yourself** | SCF, CBP, CPS, PSID, BEA, FRED, QFR | from each agency; see `README.md` "Data Sources" |
| 🔒 **Proprietary, NOT included** | Dun & Bradstreet establishment microdata; Compustat | license via D&B / WRDS |

**The headline reduced-form and structural results in the paper depend on the
proprietary D&B match, which is not redistributable.** What you *can* fully
reproduce from public data alone: the entire Yellow Pages / US Telephone
Directory corpus analysis (`Code/yp_ocr_analysis.R`, `Code/yp_ocr_full_analysis.R`,
and the public panels of `Code/yp_visualizations.R`). This is the data-collection
public good released with the paper — see `Data/Raw/Yellow_Pages/OCR_CORPUS_DATA_DICTIONARY.md`.

## Repository map

- `Code/00_setup.R` — package install + path globals (`DATA_RAW`, `DATA_GEN`,
  `OUTPUT_FIG`, `OUTPUT_TAB`). Every script is sourced after this.
- `Code/run_all.R` — master driver. Phase 1 = data cleaning (01–09),
  Phase 2 = analysis (10–22), Phase 3 = Yellow Pages corpus (public).
- `Code/yp_*.R` — Yellow Pages corpus characterization (Phase 3).
- `Data/Raw/` — source data (mostly git-ignored; see buckets above).
- `Data/Generated/` — intermediate + final datasets (pipeline writes here).
- `Output/Figures`, `Output/Tables` — generated artifacts the paper `\input`s.
- `Paper/main.tex`, `Paper/slides.tex` — paper and BoC slide deck.

## House rules (important — these reflect the author's standing preferences)

1. **No hardcoded numbers/tables in the paper.** Everything in `main.tex` /
   `slides.tex` must be an `\input{../Output/Tables/*.tex}` or `\includegraphics`
   of a script-generated artifact. The lone exception is tables lifted verbatim
   from published literature (clearly attributed).
2. **Never patch an output as a one-off.** If you tweak a figure or table, edit
   the generating `Code/*.R` script and regenerate — do not hand-edit the PDF/TeX
   output, or the next `run_all.R` will silently revert it.
3. **Never commit or redistribute D&B or other proprietary data.** It is listed
   under "NEVER share" in `.gitignore`. Only public, LoC-derived data is released.
4. **Building LaTeX:** always pass `</dev/null` to `latexmk` so it can't hang at
   the `?` error prompt:
   `cd Paper && latexmk -pdf -interaction=nonstopmode main.tex </dev/null`
5. After an interrupted `latexmk` run, biblatex aux can corrupt; fix with
   `latexmk -C` then rebuild.

## Extending the corpus (for the next researcher)

The OCR corpus is two layers (merged book-level with bounding boxes +
category-header flags; page-level raw OCR). Schema and provenance are in
`Data/Raw/Yellow_Pages/OCR_CORPUS_DATA_DICTIONARY.md`. To analyze it, model your
script on `Code/yp_ocr_analysis.R`: read from the setup-derived `DATA_RAW` path
(or the `YP_OCR_DIR` env var), write figures to `Output/Figures/yp_review/`, and
add the script to Phase 3 of `Code/run_all.R`.
