#!/bin/bash
# Run All Validation Scripts for New Merger Pipeline
# ===================================================
# 
# This script runs all validation checks in sequence:
# 1. Sanity check (compare with old merger)
# 2. Comprehensive validation (data quality)
# 3. Dashboard (summary view)
#
# Usage: ./run_all_validations.sh

set -e  # Exit on error

echo "================================================================================"
echo "NEW MERGER PIPELINE - COMPREHENSIVE VALIDATION"
echo "================================================================================"
echo ""
echo "This will run all validation scripts in sequence:"
echo "  1. Sanity check (compare with old merger)"
echo "  2. Data quality validation"
echo "  3. Summary dashboard"
echo ""
echo "================================================================================"
echo ""

# Check if we're in the right directory
if [ ! -f "run_18s_ncbi_merger.py" ]; then
    echo "❌ Error: Must run from new_merger/ directory"
    echo "   cd to new_merger/ and try again"
    exit 1
fi

# Check if output files exist
if [ ! -f "outputs/18s_ncbi_merged_division.csv" ]; then
    echo "⚠️  Warning: 18S output files not found"
    echo "   Run: python run_18s_ncbi_merger.py"
    echo ""
fi

if [ ! -f "outputs/16s_ncbi_merged_division.csv" ]; then
    echo "⚠️  Warning: 16S output files not found"
    echo "   Run: python run_16s_ncbi_merger.py"
    echo ""
fi

echo "================================================================================"
echo "STEP 1/3: SANITY CHECK (Compare with old merger)"
echo "================================================================================"
echo ""

python utils/sanity_check_merger.py

echo ""
echo "================================================================================"
echo "STEP 2/3: DATA QUALITY VALIDATION"
echo "================================================================================"
echo ""

python utils/validate_outputs.py

echo ""
echo "================================================================================"
echo "STEP 3/3: SUMMARY DASHBOARD"
echo "================================================================================"
echo ""

python utils/validation_dashboard.py

echo ""
echo "================================================================================"
echo "VALIDATION COMPLETE"
echo "================================================================================"
echo ""
echo "✅ All validation scripts completed successfully!"
echo ""
echo "📁 Generated reports:"
echo "   - sanity_check/18s_sanity_check_*.txt"
echo "   - sanity_check/16s_sanity_check_*.txt"
echo "   - validation_reports/validation_report_*.txt"
echo ""
echo "📊 Review the reports to verify merger outputs are correct."
echo ""
echo "================================================================================"

