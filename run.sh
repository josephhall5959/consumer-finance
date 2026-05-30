#!/bin/bash
set -e
cd "$(dirname "$0")"

echo "=== Credit Cards and Retail Firms — Replication Package ==="
echo ""

echo "=== Step 1: Installing R packages ==="
Rscript Code/00_setup.R

echo ""
echo "=== Step 2: Running analysis pipeline ==="
Rscript Code/run_all.R

echo ""
echo "=== Step 3: Compiling paper ==="
cd Paper
pdflatex -interaction=nonstopmode main.tex || true
bibtex main || true
pdflatex -interaction=nonstopmode main.tex || true
pdflatex -interaction=nonstopmode main.tex || true

echo ""
echo "=== Step 4: Compiling slides ==="
pdflatex -interaction=nonstopmode slides.tex || true
pdflatex -interaction=nonstopmode slides.tex || true
cd ..

echo ""
echo "=== Done! ==="
echo "  Paper:  Paper/main.pdf"
echo "  Slides: Paper/slides.pdf"
echo ""
echo "  Figures: Output/Figures/"
echo "  Tables:  Output/Tables/"
