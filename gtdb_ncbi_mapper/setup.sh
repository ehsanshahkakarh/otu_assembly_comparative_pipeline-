#!/bin/bash
# GTDB-NCBI Mapper Setup Script
# ==============================
# This script sets up the required conda environments for the pipeline

set -e  # Exit on error

echo "========================================================================"
echo "GTDB-NCBI Mapper - Environment Setup"
echo "========================================================================"
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Please install Miniconda or Anaconda first."
    echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✓ Conda found: $(conda --version)"
echo ""

# Create main pipeline environment
echo "Creating main pipeline environment (gtdb_ncbi_mapper)..."
if conda env list | grep -q "^gtdb_ncbi_mapper "; then
    echo "  Environment 'gtdb_ncbi_mapper' already exists."
    read -p "  Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n gtdb_ncbi_mapper -y
        conda env create -f environment.yml
    fi
else
    conda env create -f environment.yml
fi
echo "✓ Main environment created"
echo ""

# Create NCBI datasets environment
echo "Creating NCBI datasets environment (ncbi_datasets)..."
if conda env list | grep -q "^ncbi_datasets "; then
    echo "  Environment 'ncbi_datasets' already exists."
    read -p "  Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n ncbi_datasets -y
        conda env create -f ncbi_datasets_env.yml
    fi
else
    conda env create -f ncbi_datasets_env.yml
fi
echo "✓ NCBI datasets environment created"
echo ""

# Verify installations
echo "Verifying installations..."
echo ""

echo "1. Main environment (gtdb_ncbi_mapper):"
conda run -n gtdb_ncbi_mapper python -c "import pandas; import yaml; import tqdm; print('  ✓ All packages installed')"

echo ""
echo "2. NCBI datasets environment (ncbi_datasets):"
if conda run -n ncbi_datasets datasets --version &> /dev/null; then
    echo "  ✓ NCBI datasets tool installed"
else
    echo "  ⚠ WARNING: NCBI datasets tool not found"
    echo "  You may need to install it manually:"
    echo "  https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
fi

echo ""
echo "========================================================================"
echo "✓ Setup Complete!"
echo "========================================================================"
echo ""
echo "To use the pipeline:"
echo "  1. Activate the environment:"
echo "     conda activate gtdb_ncbi_mapper"
echo ""
echo "  2. Run the merger:"
echo "     python merge_gtdb_ncbi_accessions.py"
echo ""
echo "  3. Generate mapping tables:"
echo "     python create_gtdb_ncbi_map.py"
echo ""
echo "  4. (Optional) Check unmapped status:"
echo "     python check_unmapped_status.py archaea"
echo ""
echo "For more information, see README.md"
echo ""

