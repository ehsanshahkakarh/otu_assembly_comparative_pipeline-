#!/bin/bash
# 18S Census Parse Cleanup Script
# Date: 2026-03-03
# Purpose: Reorganize 18S_censusparse to match clean ncbi_parse structure

set -e  # Exit on error

echo "================================================================================"
echo "18S CENSUS PARSE CLEANUP"
echo "================================================================================"
echo ""

# Get the script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Working directory: $SCRIPT_DIR"
echo ""

# Step 1: Create archive directory
echo "[1/9] Creating archive directory..."
mkdir -p 18S_censusparse_old
echo "✓ Created 18S_censusparse_old/"

# Step 2: Move archives from py_18S
echo ""
echo "[2/9] Moving archive directories..."
cd py_18S
if [ -d "archive" ]; then
    mv archive ../18S_censusparse_old/
    echo "✓ Moved archive/"
fi
if [ -d "archive_current_src" ]; then
    mv archive_current_src ../18S_censusparse_old/
    echo "✓ Moved archive_current_src/"
fi
if [ -d "sanity_checks" ]; then
    mv sanity_checks ../18S_censusparse_old/
    echo "✓ Moved sanity_checks/"
fi

# Step 3: Consolidate documentation
echo ""
echo "[3/9] Consolidating documentation..."
if [ -f "AI_CACHE_IMPLEMENTATION_SUMMARY.md" ]; then
    mv AI_CACHE_IMPLEMENTATION_SUMMARY.md docs/
    echo "✓ Moved AI_CACHE_IMPLEMENTATION_SUMMARY.md"
fi
if [ -f "AI_RESOLUTION_WORKFLOW.md" ]; then
    mv AI_RESOLUTION_WORKFLOW.md docs/
    echo "✓ Moved AI_RESOLUTION_WORKFLOW.md"
fi
if [ -f "PIPELINE_STRUCTURE.md" ]; then
    mv PIPELINE_STRUCTURE.md docs/
    echo "✓ Moved PIPELINE_STRUCTURE.md"
fi
if [ -f "QUICKSTART.md" ]; then
    mv QUICKSTART.md docs/
    echo "✓ Moved QUICKSTART.md"
fi
if [ -f "CLEANUP_PLAN.md" ]; then
    mv CLEANUP_PLAN.md docs/
    echo "✓ Moved CLEANUP_PLAN.md"
fi

# Step 4: Remove redundant scripts
echo ""
echo "[4/9] Removing redundant runner scripts..."
rm -f run_complete_pipeline.py
rm -f run_division_context_adder.py
rm -f run_systematic_resolver.py
rm -f run_taxonkit_parser.py
rm -f fix_corrupted_csv.py
echo "✓ Removed old run_* scripts"

# Step 5: Clean cache directory
echo ""
echo "[5/9] Cleaning cache directory..."
cd cache
rm -f ai_cache_run*.log
rm -f ai_resolution_review_*.txt
echo "✓ Removed old cache logs"
cd ..

# Step 6: Move __pycache__ to archive
echo ""
echo "[6/9] Archiving __pycache__..."
if [ -d "__pycache__" ]; then
    mv __pycache__ ../18S_censusparse_old/
    echo "✓ Moved __pycache__/"
fi

# Step 7: Reorganize parent directory
echo ""
echo "[7/9] Reorganizing parent directory..."
cd ..

# Rename csv_outputs to output
if [ -d "csv_outputs" ]; then
    mv csv_outputs output
    echo "✓ Renamed csv_outputs/ → output/"
fi

# Step 8: Flatten py_18S structure
echo ""
echo "[8/9] Flattening directory structure..."

# Move key files and directories up
mv py_18S/src .
mv py_18S/cache .
mv py_18S/logs .
mv py_18S/docs .
mv py_18S/run_18S_pipeline.py .
mv py_18S/example_ai_resolutions.json .

# Move README to docs
if [ -f "py_18S/README.md" ]; then
    mv py_18S/README.md docs/PIPELINE_README.md
    echo "✓ Moved README.md → docs/PIPELINE_README.md"
fi

# Remove empty py_18S directory
if [ -d "py_18S" ] && [ -z "$(ls -A py_18S)" ]; then
    rmdir py_18S
    echo "✓ Removed empty py_18S/"
fi

# Step 9: Create archive README
echo ""
echo "[9/9] Creating archive README..."
cat > 18S_censusparse_old/README.md << 'EOF'
# 18S Census Parse - Archived Files

This directory contains archived files from the 18S census parse reorganization on 2026-03-03.

## Contents

- `archive/` - Old experimental code and legacy parsers
- `archive_current_src/` - Previous version of source code
- `sanity_checks/` - Old sanity check scripts
- `old_docs/` - Archived documentation

## Current Code

The current, clean 18S census parse code is in:
`00_gaps_taxonomic/parse_repaa_table/18S_censusparse/`

## Reason for Archiving

These files were archived during a cleanup to:
1. Flatten the directory structure (remove nested py_18S/)
2. Consolidate documentation
3. Remove redundant scripts
4. Create a clean, maintainable structure

## Date Archived
2026-03-03
EOF
echo "✓ Created archive README"

echo ""
echo "================================================================================"
echo "CLEANUP COMPLETE!"
echo "================================================================================"
echo ""
echo "New structure:"
echo "  18S_censusparse/"
echo "    ├── README.md"
echo "    ├── run_18S_pipeline.py"
echo "    ├── example_ai_resolutions.json"
echo "    ├── src/          (17 modules)"
echo "    ├── cache/        (ai_resolutions.json)"
echo "    ├── logs/         (5 log files)"
echo "    ├── output/       (3 CSV files)"
echo "    ├── metadata/     (input data)"
echo "    └── docs/         (consolidated docs)"
echo ""
echo "Archive created:"
echo "  18S_censusparse_old/"
echo ""
echo "Next step: Move 18S_censusparse_old/ to old_pipeline/"
echo "  mv 18S_censusparse_old ../old_pipeline/"
echo ""

