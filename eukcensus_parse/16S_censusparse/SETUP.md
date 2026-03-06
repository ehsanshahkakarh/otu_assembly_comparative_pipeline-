# 16S Census Parser - Setup Guide

## Quick Start (Recommended)

### Option 1: Conda Environment (Easiest)

```bash
# Create the environment from the provided file
conda env create -f environment.yml

# Activate the environment
conda activate 16s_censusparse

# Run the parser
cd py_16S
python run_16S_parser.py
```

### Option 2: Manual Installation

```bash
# Install Python dependencies
pip install -r requirements.txt

# Install taxonkit (choose one method):

# Method A: Using conda
conda install -c bioconda taxonkit

# Method B: Download binary from GitHub
# Visit: https://github.com/shenwei356/taxonkit/releases
# Download the appropriate binary for your system
# Add to PATH

# Method C: Build from source (requires Go)
go install github.com/shenwei356/taxonkit@latest
```

## System Requirements

- **Python**: 3.8 or higher
- **Memory**: 4GB RAM minimum (8GB recommended for large datasets)
- **Disk Space**: ~50GB for NCBI taxonomy database
- **OS**: Linux, macOS, or Windows (with WSL)

## External Dependencies

### taxonkit (Required)

The parser requires `taxonkit` for NCBI taxonomy lookups.

**Installation verification:**
```bash
taxonkit version
```

**NCBI Taxonomy Database:**

Taxonkit needs the NCBI taxonomy database. Download it:

```bash
# Create directory for taxonomy data
mkdir -p ~/.taxonkit

# Download and extract NCBI taxonomy
cd ~/.taxonkit
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

# Set environment variable (add to ~/.bashrc or ~/.zshrc)
export TAXONKIT_DB=~/.taxonkit
```

## Input Data

The parser expects the input file at:
```
metadata/eukcensus_16S.clusters.97.tsv
```

**Required columns:**
- `centroid`: OTU cluster centroid ID
- `members`: Number of sequences in cluster
- `size`: Total sequence abundance
- `phylum`: Phylum/division name
- `familiy`: Family name (note: typo in original data)
- `genus`: Genus name

## Running the Parser

### Basic Usage

```bash
cd py_16S
python run_16S_parser.py
```

### Custom Input File

```bash
python run_16S_parser.py path/to/custom_input.tsv
```

### Custom Output Prefix

```bash
python run_16S_parser.py custom_input.tsv custom_output
```

## Output Files

The parser generates the following files:

### CSV Outputs (in `csv_16S/`)
- `eukcensus16S_by_division.csv` - Phylum-level aggregation
- `eukcensus16S_by_family.csv` - Family-level aggregation
- `eukcensus16S_by_genus.csv` - Genus-level aggregation

### Log Files (in `py_16S/logs/`)
- `eukcensus16S_processing.log` - Processing log
- `eukcensus16S_comprehensive_unmapped.log` - Unmapped taxa analysis

## Troubleshooting

### "taxonkit: command not found"

**Solution:** Install taxonkit using one of the methods above and ensure it's in your PATH.

### "NCBI taxonomy database not found"

**Solution:** Download the taxonomy database and set the `TAXONKIT_DB` environment variable.

### "Memory error" or slow performance

**Solution:** The parser uses chunked reading for large files. If you still encounter issues:
- Increase available RAM
- Reduce chunk size in `run_census_parser.py` (line 83)

### Import errors

**Solution:** Ensure you're running from the correct directory and all dependencies are installed:
```bash
cd py_16S
python -c "import pandas, tqdm; print('Dependencies OK')"
```

## Performance

**Expected processing time:**
- Small datasets (<10K OTUs): ~30 seconds
- Medium datasets (10K-100K OTUs): 1-5 minutes
- Large datasets (>100K OTUs): 5-15 minutes

**Optimization tips:**
- Use SSD storage for faster I/O
- Ensure taxonkit database is on local disk (not network drive)
- Use the conda environment for optimized dependencies

## Next Steps

After successful setup and execution:

1. Check output files in `csv_16S/`
2. Review logs in `py_16S/logs/`
3. Use outputs for downstream analysis (database merging, visualization)

For detailed documentation, see:
- `README.md` - Pipeline overview
- `REPRODUCIBILITY_GUIDE.md` - Reproducibility analysis
- `RESOLUTION_SYSTEM_EXPLAINED.md` - Resolution system details

