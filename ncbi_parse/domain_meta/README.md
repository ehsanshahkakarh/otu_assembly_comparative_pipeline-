# Domain Splitter

Splits species-grouped NCBI assembly data into domain-specific CSV files.

## Overview

This tool takes the output from the main `ncbi_parse` pipeline and splits it by taxonomic domain:
- **Bacteria**
- **Archaea**
- **Eukaryota**
- **Viruses**
- **Unknown** (metagenomes, environmental samples, unclassified)

## Usage

```bash
cd 00_gaps_taxonomic/parse_repaa_table/ncbi_parse/domain_meta
python3 run_domain_splitter.py
```

## Input

The script automatically finds the most recent `species_grouped_*.csv` file from:
```
../output/species_grouped_YYYYMMDD_HHMMSS.csv
```

## Output

### Domain-Specific CSV Files

Located in `output/`:
- `bacteria_species.csv` - All bacterial species
- `archaea_species.csv` - All archaeal species
- `eukaryota_species.csv` - All eukaryotic species
- `viruses_species.csv` - All viral species
- `unknown_species.csv` - Metagenomes and unclassified entries
- `domain_summary_YYYYMMDD_HHMMSS.csv` - Summary statistics

### Columns in Domain CSV Files

Same as main pipeline output, plus:
- `domain` - Taxonomic domain (Bacteria, Archaea, Eukaryota, Viruses, Unknown)

All other columns from the species_grouped output are preserved:
- `species_taxid`
- `species_name`
- `total_genome_count`
- `isolate_genome_count`
- `uncultured_genome_count`
- `isolate_genome_percentage`
- `species_genome_percentage`
- `lineage`
- `lineage_ranks`
- `lineage_taxids`

## Domain Extraction Logic

The domain is extracted from the `lineage` field:

1. **Viruses**: First element is "Viruses"
2. **Bacteria/Archaea/Eukaryota**: Second element (after "cellular organisms")
3. **Unknown**: Missing lineage or unclassified entries

## Latest Run (Mar 2, 2026)

**Input:** `species_grouped_20260301_214126.csv`  
**Total species:** 162,731

**Results:**
- Bacteria: 99,085 species (2,930,149 genomes)
- Viruses: 35,847 species (256,154 genomes)
- Eukaryota: 24,211 species (59,130 genomes)
- Archaea: 3,433 species (34,539 genomes)
- Unknown: 155 species (35,849 genomes - metagenomes/environmental)

**Success rate:** 99.9% (only 155 legitimately unclassified entries)

## Directory Structure

```
domain_meta/
├── README.md                    # This file
├── ANALYSIS_SUMMARY.md          # Detailed analysis of old vs new results
├── run_domain_splitter.py       # Main script
├── src/
│   ├── __init__.py
│   └── domain_splitter.py       # Domain splitting logic
├── logs/
│   └── domain_splitter_*.log    # Execution logs
└── output/
    ├── bacteria_species.csv
    ├── archaea_species.csv
    ├── eukaryota_species.csv
    ├── viruses_species.csv
    ├── unknown_species.csv
    └── domain_summary_*.csv
```

## Integration

Use the domain-specific CSV files for downstream analysis:

```python
import pandas as pd

# Load bacteria-only data
bacteria = pd.read_csv('output/bacteria_species.csv')

# Filter for high-quality isolate data
high_quality = bacteria[bacteria['isolate_genome_percentage'] > 50]

# Get top species by genome count
top_species = bacteria.nlargest(100, 'total_genome_count')
```

## Notes

- The script automatically finds the latest species_grouped file
- Sorts each domain CSV by total_genome_count (descending)
- Generates summary statistics for each domain
- Logs all operations to `logs/domain_splitter_*.log`

## See Also

- `ANALYSIS_SUMMARY.md` - Detailed comparison of old vs new results
- `../README.md` - Main ncbi_parse pipeline documentation

