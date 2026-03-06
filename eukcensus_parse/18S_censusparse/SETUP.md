# Setup Instructions - 18S Census Parser

## Prerequisites

- **Python:** 3.7 or higher
- **Package Manager:** conda (recommended) or pip
- **Disk Space:** ~500 MB (for data and outputs)
- **RAM:** 4 GB minimum, 8 GB recommended

---

## Installation

### Option 1: Conda (Recommended)

This will install Python, pandas, and taxonkit automatically:

```bash
# Create conda environment from file
conda env create -f environment.yml

# Activate the environment
conda activate 18s_censusparse

# Verify taxonkit is installed
taxonkit version
```

### Option 2: pip + Manual taxonkit Installation

If you prefer pip or don't have conda:

```bash
# Install Python dependencies
pip install -r requirements.txt

# Install taxonkit manually
# For Linux:
wget https://github.com/shenwei356/taxonkit/releases/download/v0.15.1/taxonkit_linux_amd64.tar.gz
tar -xzf taxonkit_linux_amd64.tar.gz
sudo mv taxonkit /usr/local/bin/

# For macOS:
wget https://github.com/shenwei356/taxonkit/releases/download/v0.15.1/taxonkit_darwin_amd64.tar.gz
tar -xzf taxonkit_darwin_amd64.tar.gz
sudo mv taxonkit /usr/local/bin/

# Verify installation
taxonkit version
```

---

## Input Data

### Download EukCensus 18S Data

The pipeline requires the EukCensus 18S cluster file:

```bash
# If you have the file locally, copy it to metadata/
cp /path/to/eukcensus_18S.clusters.97.tsv metadata/

# Or download from EukCensus repository (if publicly available)
# wget <URL> -O metadata/eukcensus_18S.clusters.97.tsv
```

**File details:**
- **Filename:** `eukcensus_18S.clusters.97.tsv`
- **Size:** ~56 MB
- **Location:** Must be in `metadata/` directory
- **Format:** Tab-separated values with columns: cluster_id, size, taxonomy

---

## Directory Structure

The pipeline expects this structure (will auto-create output directories):

```
18S_censusparse/
├── metadata/
│   └── eukcensus_18S.clusters.97.tsv  ← INPUT FILE (you provide)
├── src/                                ← Python modules
├── output/                             ← Auto-created
├── logs/                               ← Auto-created
├── cache/                              ← Auto-created
├── run_18S_pipeline.py                 ← Main script
├── requirements.txt
├── environment.yml
└── SETUP.md (this file)
```

---

## Running the Pipeline

### Quick Start (All Steps)

```bash
# Run complete pipeline
python run_18S_pipeline.py --all
```

### Individual Steps (Advanced)

```bash
# Step 1: Taxonkit parser (NCBI lookups)
python run_18S_pipeline.py --step taxonkit

# Step 2: Systematic resolver (custom resolutions)
python run_18S_pipeline.py --step resolve

# Step 3: AI cache (optional, for research)
python run_18S_pipeline.py --step ai-cache
```

---

## Expected Outputs

After successful run, you should see:

```
output/
├── eukcensus_18S_by_division.csv  ← Division-level taxonomy
├── eukcensus_18S_by_family.csv    ← Family-level taxonomy
└── eukcensus_18S_by_genus.csv     ← Genus-level taxonomy

logs/
├── eukcensus_18S_unmapped_final.log        ← Taxa not in NCBI
├── eukcensus_taxonkit_only_unmapped_from_taxonkit.log
└── systematic_resolver.log
```

---

## Verification

Check that outputs were created successfully:

```bash
# Check output files exist
ls -lh output/eukcensus_18S_by_*.csv

# Check number of lines (should be ~70,000+ each)
wc -l output/eukcensus_18S_by_*.csv

# View first few lines
head output/eukcensus_18S_by_division.csv
```

---

## Troubleshooting

### Error: "taxonkit: command not found"

**Solution:** Install taxonkit (see Option 2 above) or use conda environment

### Error: "No module named 'pandas'"

**Solution:** Install pandas: `pip install pandas`

### Error: "FileNotFoundError: metadata/eukcensus_18S.clusters.97.tsv"

**Solution:** Download/copy the input file to `metadata/` directory

### Error: "Permission denied" when writing outputs

**Solution:** Ensure you have write permissions in the directory

---

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **Python** | 3.7 | 3.9+ |
| **RAM** | 4 GB | 8 GB |
| **Disk Space** | 500 MB | 1 GB |
| **OS** | Linux, macOS, Windows | Linux |
| **taxonkit** | 0.15.0 | Latest |

---

## Testing Your Installation

Run this quick test:

```bash
# Test taxonkit
echo "Homo sapiens" | taxonkit name2taxid

# Expected output:
# Homo sapiens    9606

# Test Python imports
python -c "import pandas; print('pandas:', pandas.__version__)"

# Expected output:
# pandas: 1.x.x
```

If both tests pass, you're ready to run the pipeline! ✅

---

## Getting Help

- **Documentation:** See `README.md` for pipeline details
- **Reproducibility:** See `REPRODUCIBILITY_GUIDE.md` for portability info
- **Architecture:** See `docs/ARCHITECTURE_README.md` for technical details

---

## Quick Reference

```bash
# Setup (one time)
conda env create -f environment.yml
conda activate 18s_censusparse

# Run (every time)
python run_18S_pipeline.py --all

# Check outputs
ls -lh output/
```

That's it! 🎉

