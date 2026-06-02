#!/usr/bin/env bash
# ============================================================================
# fetch_data.sh — Download the public Yellow Pages / US Telephone Directory
# OCR corpus that accompanies "Credit Cards and Retail Firms."
#
# The corpus is built from the Library of Congress US Telephone Directories
# collection (public domain) and is archived on Zenodo. This script downloads
# it and unpacks it into Data/Raw/Yellow_Pages/OCR/ in the layout the pipeline
# expects:
#
#   Data/Raw/Yellow_Pages/OCR/
#   ├── Merged Files/                 # 110 merged book-level OCR CSVs
#   └── CSV_MANUAL/CSV_MANUAL/        # page-level OCR CSVs (one dir per book)
#
# Usage:
#   bash fetch_data.sh                # download from Zenodo and extract
#   YP_OCR_ARCHIVE=/path/to/corpus.zip bash fetch_data.sh   # use a local zip
#
# NOTE: This corpus is the *public* data only. The proprietary Dun & Bradstreet
# microdata is NOT distributed and must be licensed separately (see README).
# ============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="$ROOT/Data/Raw/Yellow_Pages/OCR"
MERGED="$DEST/Merged Files"

# ---------------------------------------------------------------------------
# Zenodo source. Replace ZENODO_URL with the direct archive URL once the
# deposit is published (see ZENODO_DEPOSIT.md). It can also be overridden via
# the environment, e.g.  ZENODO_URL=https://zenodo.org/records/XXXX/files/yp_ocr_corpus.zip
# ---------------------------------------------------------------------------
# Zenodo concept DOI 10.5281/zenodo.20513427; this pins the published version
# record (20513428) so replication fetches the exact archive used in the paper.
ZENODO_URL="${ZENODO_URL:-https://zenodo.org/records/20513428/files/yp_ocr_corpus.zip?download=1}"

# ---------------------------------------------------------------------------
# Already present? Nothing to do.
# ---------------------------------------------------------------------------
if [ -d "$MERGED" ] && [ -n "$(find "$MERGED" -maxdepth 1 -iname '*.csv' -print -quit 2>/dev/null)" ]; then
  echo "[fetch_data] OCR corpus already present at: $DEST"
  echo "[fetch_data] Nothing to do. (Delete the folder to re-download.)"
  exit 0
fi

mkdir -p "$DEST"
# Keep the temp dir on the same filesystem as the destination so the final
# move is a fast, warning-free rename rather than a cross-device copy.
TMP="$(mktemp -d "$ROOT/.fetch_tmp.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
ARCHIVE="$TMP/yp_ocr_corpus.zip"

# ---------------------------------------------------------------------------
# Obtain the archive: local override > Zenodo download.
# ---------------------------------------------------------------------------
if [ -n "${YP_OCR_ARCHIVE:-}" ]; then
  echo "[fetch_data] Using local archive: $YP_OCR_ARCHIVE"
  ARCHIVE="$YP_OCR_ARCHIVE"
elif [ "$ZENODO_URL" = "REPLACE_WITH_ZENODO_ARCHIVE_URL" ]; then
  cat >&2 <<'EOF'
[fetch_data] ERROR: No data source configured.

The Zenodo URL has not been set yet. Do one of:
  1. Publish the corpus on Zenodo (see ZENODO_DEPOSIT.md), then edit the
     ZENODO_URL line at the top of this script, or run:
       ZENODO_URL=<direct-archive-url> bash fetch_data.sh
  2. Point at a local copy of the corpus zip:
       YP_OCR_ARCHIVE=/path/to/yp_ocr_corpus.zip bash fetch_data.sh
EOF
  exit 1
else
  echo "[fetch_data] Downloading corpus from: $ZENODO_URL"
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail -o "$ARCHIVE" "$ZENODO_URL"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$ARCHIVE" "$ZENODO_URL"
  else
    echo "[fetch_data] ERROR: need curl or wget to download." >&2
    exit 1
  fi
fi

# ---------------------------------------------------------------------------
# Extract, then normalize the layout regardless of how the zip nests things:
# locate the "Merged Files" and "CSV_MANUAL" folders and move them into place.
# ---------------------------------------------------------------------------
echo "[fetch_data] Extracting..."
EXDIR="$TMP/extracted"
mkdir -p "$EXDIR"
unzip -q "$ARCHIVE" -d "$EXDIR"

found_merged="$(find "$EXDIR" -type d -iname 'Merged Files' -print -quit || true)"
# The page-level data nests as CSV_MANUAL/CSV_MANUAL/<book>/...; target the
# INNER directory (the one holding the per-book folders), not the outer wrapper.
found_pages="$(find "$EXDIR" -type d -path '*/CSV_MANUAL/CSV_MANUAL' -print -quit || true)"
if [ -z "$found_pages" ]; then
  # Fallback for an archive that stores a single, un-doubled CSV_MANUAL.
  found_pages="$(find "$EXDIR" -type d -iname 'CSV_MANUAL' -print -quit || true)"
fi

if [ -n "$found_merged" ]; then
  rm -rf "$MERGED"
  mv "$found_merged" "$MERGED"
  echo "[fetch_data] Installed Merged Files/ ($(find "$MERGED" -maxdepth 1 -iname '*.csv' | wc -l | tr -d ' ') books)"
fi
if [ -n "$found_pages" ]; then
  # Reproduce the CSV_MANUAL/CSV_MANUAL/ layout the R scripts expect.
  rm -rf "$DEST/CSV_MANUAL"
  mkdir -p "$DEST/CSV_MANUAL"
  mv "$found_pages" "$DEST/CSV_MANUAL/CSV_MANUAL"
  echo "[fetch_data] Installed CSV_MANUAL/CSV_MANUAL/ ($(find "$DEST/CSV_MANUAL/CSV_MANUAL" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ') books)"
fi

if [ ! -d "$MERGED" ] && [ ! -d "$DEST/CSV_MANUAL/CSV_MANUAL" ]; then
  echo "[fetch_data] ERROR: could not find 'Merged Files' or 'CSV_MANUAL' in the archive." >&2
  echo "[fetch_data] Extracted tree was:" >&2
  find "$EXDIR" -maxdepth 2 -type d >&2
  exit 1
fi

echo "[fetch_data] Done. Corpus is at: $DEST"
echo "[fetch_data] You can now run: Rscript Code/run_all.R  (or: bash run.sh)"
