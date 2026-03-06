
# Species-Level NCBI Assembly Statistics

## Overview

This directory contains a reproducible pipeline for generating species-level statistics from NCBI assembly data by grouping genomes by their `species_taxid` values. This method provides comprehensive species-level analysis with accurate isolate vs uncultured genome counts and full taxonomic information.

## Key Features

- **True Species-Level Grouping**: Uses NCBI's official `species_taxid` for biologically accurate species identification
- **Unbiased Genome Classification**: Classifies ALL genomes as isolate or uncultured with no filtering or selection bias
- **Accurate Statistical Counts**: Provides true counts and percentages of isolate vs uncultured genomes per species
- **Comprehensive Taxonomic Information**: Integrates taxonkit to provide full lineage information for each species
- **Reproducible Results**: Consistent methodology that can be easily replicated
- **Complete Data Retention**: Processes all genomes from the assembly file for accurate statistics

## Files

### Main Entry Point

**`run_species_grouper.py`** - Main pipeline script
   - Orchestrates the complete species-level analysis workflow
   - Loads NCBI assembly data
   - Classifies genomes as isolate or uncultured
   - Groups by species_taxid and calculates statistics
   - Enriches with taxonomic lineage information via taxonkit
   - Outputs comprehensive species-level statistics

### Source Modules (`src/`)

1. **`data_loader.py`** - Loads NCBI assembly summary data
2. **`genome_classifier.py`** - Classifies genomes as isolate or uncultured based on organism_name and excluded_from_refseq
3. **`species_grouper.py`** - Groups genomes by species_taxid and calculates counts
4. **`lineage_enricher.py`** - Adds taxonomic lineage information using taxonkit
5. **`percentage_calculator.py`** - Calculates isolate percentages and species percentages
6. **`output_writer.py`** - Saves species statistics to CSV
7. **`taxid_fallback_resolver.py`** - Resolves missing lineages using alternative taxids (optional)

### Configuration

**`config.py`** - Centralized configuration for paths, filters, and parameters

## Usage

### Basic Usage

```bash
# Run the main pipeline (processes full dataset)
python run_species_grouper.py

# Run with sample data for testing
python run_species_grouper.py --sample 10000

# Specify custom output directory
python run_species_grouper.py --output-dir custom_output/
```

## Output Files

### Main Output

**`output/species_grouped_YYYYMMDD_HHMMSS.csv`** - Species-level statistics with full lineage information

**Columns:**
- `species_taxid` - NCBI species-level taxonomic ID
- `species_name` - Species name extracted from lineage
- `total_genome_count` - Total number of genomes for this species in the assembly file
- `isolate_genome_count` - Number of genomes classified as isolates
- `uncultured_genome_count` - Number of genomes classified as uncultured
- `isolate_genome_percentage` - Percentage of genomes that are isolates (isolate_count / total_count * 100)
- `species_genome_percentage` - Percentage of total dataset this species represents
- `lineage` - Full taxonomic lineage (semicolon-separated)
- `lineage_ranks` - Taxonomic ranks for each lineage level
- `lineage_taxids` - Taxonomic IDs for each lineage level

### Log Files

**`logs/species_parser_YYYYMMDD_HHMMSS.log`** - Detailed processing log with statistics and progress information

## Method Advantages

### Compared to Previous Methods

1. **Biological Accuracy**: Uses NCBI's official species_taxid instead of parsing organism names
2. **Consistency**: Eliminates variability from strain naming conventions
3. **Completeness**: 100% coverage - processes ALL genomes in the assembly file
4. **Reproducibility**: Same results every time, regardless of data order
5. **Efficiency**: Direct grouping without complex string parsing
6. **Unbiased Statistics**: No filtering or selection - provides true counts of isolate vs uncultured genomes

### Data Quality Benefits

- **No Strain Confusion**: Different strains of same species are properly grouped together
- **Taxonomic Precision**: Official NCBI taxonomic assignments via species_taxid
- **Accurate Metrics**: True counts and percentages of isolate vs uncultured genomes per species
- **Comprehensive Lineage**: Full taxonomic hierarchy via taxonkit integration
- **Complete Data Retention**: All genomes counted for accurate statistical representation

## Requirements

### Software Dependencies

- Python 3.7+
- pandas
- taxonkit (must be installed and in PATH)

### Python Packages

```bash
pip install pandas
```

Or use the provided requirements.txt:

```bash
pip install -r requirements.txt
```

### System Requirements

- taxonkit must be installed and accessible via command line
- Sufficient memory for processing large datasets (recommend 8GB+ RAM for full NCBI dataset)

## Data Requirements

The pipeline expects to find the NCBI assembly summary file:
- `00assembly_summary_genbank.txt`

The script will automatically search for this file in the parent metadata directory:
- `../00assembly_summary_genbank.txt` (relative to ncbi_parse directory)

You can download the latest version from NCBI:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
```

## Performance Notes

- Processing full NCBI dataset (~3M+ records): ~10-20 minutes depending on system
- Memory usage: ~2-4GB for full dataset
- Taxonkit lineage retrieval: Processes all unique species_taxid values
- Output file size: ~5-10MB for species-level statistics

## Troubleshooting

### Common Issues

1. **taxonkit not found**: Ensure taxonkit is installed and in your PATH
   ```bash
   # Test taxonkit installation
   taxonkit version
   ```

2. **Memory errors**: Use `--sample` for testing with smaller datasets
   ```bash
   python run_species_grouper.py --sample 10000
   ```

3. **File not found**: Check that `00assembly_summary_genbank.txt` is in the parent metadata directory
   ```bash
   ls ../00assembly_summary_genbank.txt
   ```

4. **Missing lineages**: Ensure taxonkit has access to NCBI taxonomy database
   ```bash
   # Download/update taxonkit database
   taxonkit create-taxdump
   ```

## Pipeline Workflow

The pipeline executes the following steps:

1. **Load Assembly Data** - Reads NCBI assembly_summary_genbank.txt file
2. **Classify Genome Sources** - Classifies each genome as 'isolate' or 'uncultured' based on:
   - organism_name patterns (uncultured, metagenome, environmental, etc.)
   - excluded_from_refseq notes
   - Special handling for enrichment cultures with strain names
3. **Group by Species** - Groups all genomes by species_taxid and calculates:
   - Total genome count per species
   - Isolate genome count per species
   - Uncultured genome count per species
   - Isolate percentage per species
4. **Enrich with Lineage** - Adds full taxonomic lineage using taxonkit
5. **Calculate Percentages** - Adds species_genome_percentage (% of total dataset)
6. **Save Output** - Writes species-level statistics to timestamped CSV file

## Integration with Downstream Analysis

This output can be used for:

1. **Database Merging**: Merge with 16S/18S census data using species_taxid
2. **Gap Analysis**: Identify underrepresented taxonomic groups
3. **Isolate Availability**: Determine which species have cultured isolates available
4. **Coverage Assessment**: Evaluate genome coverage across taxonomic hierarchy
5. **Priority Setting**: Identify species/lineages for targeted sequencing

## Example Output

Sample rows from `species_grouped_YYYYMMDD_HHMMSS.csv`:

```csv
species_taxid,species_name,total_genome_count,isolate_genome_count,uncultured_genome_count,isolate_genome_percentage,species_genome_percentage,lineage,lineage_ranks,lineage_taxids
562,Escherichia coli,45230,42150,3080,93.19,1.52,Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli,superkingdom;phylum;class;order;family;genus;species,2;1224;1236;91347;543;561;562
1280,Staphylococcus aureus,23456,22890,566,97.59,0.79,Bacteria;Bacillota;Bacilli;Bacillales;Staphylococcaceae;Staphylococcus;Staphylococcus aureus,superkingdom;phylum;class;order;family;genus;species,2;1239;91061;1385;90964;1279;1280
```

## Genome Classification Logic

The pipeline classifies genomes as **isolate** or **uncultured** using the following logic:

1. **Default**: All genomes start as 'isolate'
2. **Check organism_name** for uncultured patterns:
   - uncultured, environmental, metagenome, unclassified
   - unknown, unidentified, mixed culture, enrichment culture
   - derived from metagenome, metagenome-assembled, MAG
   - single amplified genome, SAG, environmental sample
3. **Check excluded_from_refseq** column for same patterns
4. **Special case**: Enrichment cultures with strain names remain as 'isolate'

This ensures accurate classification while handling edge cases appropriately.

---

*Created: 2025-01-27*
*Updated: 2025-03-03 - Clarified statistics-focused approach, removed isolate-preferential selection references*
*Part of the NCBI assembly parsing pipeline*
