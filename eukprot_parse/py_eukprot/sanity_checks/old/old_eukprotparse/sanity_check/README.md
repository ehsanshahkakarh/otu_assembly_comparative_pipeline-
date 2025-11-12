# EukProt Taxonomy Sanity Check

This directory contains scripts and data for verifying the taxonomic classifications in the EukProt database.

## Directory Structure

```
sanity_check/
├── README.md                      # This file
├── Eukaryota_groups.py            # Script to analyze eukaryotic taxonomic groups
├── taxon_names/                   # Directory containing taxonomic name files
│   ├── create_eukprot_taxon_files.py  # Script to create taxonomic name files with taxids
│   ├── verify_eukprot_taxonomy.py     # Script to verify taxonomic classifications
│   ├── eukaryotic_phyla.csv           # Phylum names with taxids
│   ├── eukaryotic_families.csv        # Family names with taxids
│   ├── eukaryotic_genera.csv          # Genus names with taxids
│   ├── logs/                          # Directory for log files
│   └── temp_taxa_batch_*.txt          # Temporary files for batch processing
```

## Workflow

1. **Create Taxonomic Name Files**:
   - Extract taxonomic names from EukProt data
   - Look up NCBI taxids for each name
   - Create CSV files with taxon_name and taxid columns

2. **Verify Taxonomic Classifications**:
   - Check if each taxon belongs to the expected domain (Eukaryota)
   - Verify that each taxon is classified at the expected rank (phylum, family, or genus)
   - Identify and flag Candidatus taxa (which should be excluded)
   - Generate verification logs and summary statistics

## Scripts

### `create_eukprot_taxon_files.py`

This script creates CSV files with taxonomic names and their corresponding NCBI taxids.

```bash
# Run the script
cd taxon_names
python create_eukprot_taxon_files.py
```

The script:
- Reads taxonomic names from the EukProt data files
- Uses taxonkit to look up NCBI taxids for each name
- Creates CSV files with taxon_name and taxid columns
- Skips Candidatus taxa
- Generates log files with statistics

### `verify_eukprot_taxonomy.py`

This script verifies the taxonomic classifications in the EukProt database.

```bash
# Run the script
cd taxon_names
python verify_eukprot_taxonomy.py
```

The script:
- Reads taxonomic names from the CSV files
- Uses taxonkit to get lineage information for each taxon
- Checks if each taxon belongs to the expected domain (Eukaryota)
- Verifies that each taxon is classified at the expected rank
- Identifies and flags Candidatus taxa
- Generates verification logs and summary statistics

Optional arguments:
- `--taxdump PATH`: Path to NCBI taxdump directory
- `--output-dir PATH`: Directory for output files

### `Eukaryota_groups.py`

This script analyzes the taxonomic groups in the EukProt database.

```bash
# Run the script
python Eukaryota_groups.py
```

The script:
- Analyzes the distribution of taxonomic groups in the EukProt database
- Identifies the most common phyla, families, and genera
- Generates statistics on taxonomic coverage
- Creates visualizations of taxonomic distributions

## Output Files

### Taxonomic Name Files

- `eukaryotic_phyla.csv`: Phylum names with taxids
- `eukaryotic_families.csv`: Family names with taxids
- `eukaryotic_genera.csv`: Genus names with taxids

Format:
```
taxon_name,taxid
Ascomycota,4890
Basidiomycota,5204
...
```

### Verification Logs

- `eukaryotic_phyla_verification.log`: Verification log for phyla
- `eukaryotic_families_verification.log`: Verification log for families
- `eukaryotic_genera_verification.log`: Verification log for genera

Format:
```
Taxon    Issue
Taxon1   Domain mismatch: Not in Eukaryota
Taxon2   Rank missing: No phylum classification
Taxon3   Candidatus taxon should be excluded
...
```

### Summary File

- `verification_summary.txt`: Summary of verification results

Format:
```
File                    Total Taxa    Mismatches    Percentage Correct
eukaryotic_phyla.csv    100           5             95.00%
eukaryotic_families.csv 500           25            95.00%
eukaryotic_genera.csv   1000          50            95.00%
```

## Integration with Other Databases

This sanity check is part of a larger system that also processes:

1. **NCBI Database**: Contains both prokaryotic and eukaryotic data
2. **GTDB Database**: Contains prokaryotic data

The verified taxonomic data will be used for:
- Cross-database taxonomic comparisons
- Identification of taxonomic inconsistencies
- Analysis of genome representation across databases
- Creation of unified taxonomic reference datasets
