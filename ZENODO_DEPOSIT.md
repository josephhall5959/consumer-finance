# Zenodo deposit checklist (Yellow Pages OCR corpus)

This is the step-by-step for publishing the public Yellow Pages corpus so the
paper can cite a DOI and `fetch_data.sh` can download it. Metadata is
pre-filled in `.zenodo.json`.

## What to upload — ALREADY BUILT

The upload bundle has been assembled for you at **`/workspace/zenodo_upload/`**:

| File | Size | Notes |
|------|------|-------|
| `yp_ocr_corpus.zip` | ~609 MB | the corpus (verified: 110/110 merged books, 165/165 page-level book dirs, integrity OK) |
| `yp_ocr_corpus.zip.sha256` | — | checksum for your records / Zenodo integrity |
| `Merchant List.xlsx` | ~720 KB | the 510 hand-transcribed acceptance records |
| `OCR_CORPUS_DATA_DICTIONARY.md` | ~4 KB | schema + provenance |
| `.zenodo.json` | — | deposit metadata (reference copy) |

The zip's internal layout is exactly what `fetch_data.sh` expects:

```
yp_ocr_corpus.zip
├── Merged Files/                 # 110 merged book-level OCR CSVs
├── CSV_MANUAL/CSV_MANUAL/        # page-level OCR CSVs, one dir per book (165)
└── OCR_CORPUS_DATA_DICTIONARY.md
```

To rebuild it from scratch (only if needed — e.g. the corpus changes), the
source lives unzipped at `/workspace/Input/_extracted/`
(`Yellow_Pages_OCR/Merged Files` and `CSV_MANUAL/CSV_MANUAL`):

```bash
cd "/workspace/Input/_extracted/Yellow_Pages_OCR" && \
  zip -r -q /workspace/yp_ocr_corpus.zip "Merged Files" -x "*.DS_Store" "*__MACOSX*"
cd "/workspace/Input/_extracted" && \
  zip -r -q /workspace/yp_ocr_corpus.zip "CSV_MANUAL/CSV_MANUAL" -x "*.DS_Store" "*__MACOSX*"
cd /workspace/Replication/Data/Raw/Yellow_Pages && \
  zip -q /workspace/yp_ocr_corpus.zip "OCR_CORPUS_DATA_DICTIONARY.md"
```

**Upload to Zenodo:** `yp_ocr_corpus.zip` plus, as separate files for
discoverability, `OCR_CORPUS_DATA_DICTIONARY.md` and `Merchant List.xlsx`.

**Do NOT upload** anything derived from Dun & Bradstreet
(`db_phonebooks_merged.csv`, `retail.csv`, `db_70-85.csv`).

## Steps

1. Sign in at <https://zenodo.org> (ORCID or GitHub login).
2. **New upload** → drag in `yp_ocr_corpus.zip` + the two companion files.
3. Zenodo can ingest `.zenodo.json` automatically if you connect the GitHub repo
   and cut a release; otherwise copy the fields from `.zenodo.json` by hand
   (title, description, dataset type, CC-BY-4.0, creator, keywords).
4. **Reserve a DOI** before publishing (button in the form) so you can cite it in
   the paper.
5. **Publish**. Note the record number and the direct file URL, which looks like:
   `https://zenodo.org/records/<RECORD_ID>/files/yp_ocr_corpus.zip`

## After publishing — wire up the DOI

1. In `fetch_data.sh`, set `ZENODO_URL` to the direct archive URL above.
2. In `Paper/main.tex`, replace `[Zenodo DOI to be inserted upon deposit]`
   (Appendix \ref{sec:construction}) with the reserved DOI, e.g.
   `\href{https://doi.org/10.5281/zenodo.XXXXXXX}{10.5281/zenodo.XXXXXXX}`.
3. Add a data-availability footnote/citation if the journal requires one.
4. Rebuild the paper (`cd Paper && latexmk -pdf -interaction=nonstopmode main.tex </dev/null`).

## License note

Source scans are public domain (Library of Congress). The compiled corpus and
code are released CC-BY-4.0. Confirm this is the license you want before
publishing — CC-BY-4.0 requires attribution; CC0 would waive even that.
