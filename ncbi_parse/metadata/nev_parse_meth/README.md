# Species Taxid Grouper - New Reproducible Method

## Overview

This directory contains a new, more reproducible approach to parsing NCBI assembly data by grouping identical `species_taxid` values. This method provides true species-level analysis with isolate-preferential selection and comprehensive taxonomic information.

## Key Features

- **True Species-Level Grouping**: Uses NCBI's official `species_taxid` for biologically accurate species identification
- **Isolate-Preferential Selection**: Prioritizes isolate genomes over uncultured when both are available for a species
- **Comprehensive Taxonomic Information**: Integrates taxonkit to provide full lineage information
- **Reproducible Results**: Consistent methodology that can be easily replicated
- **Quality Assurance**: Built-in verification and analysis tools
- **Progress Tracking**: Uses tqdm for progress visualization

## Files

### Core Scripts

1. **`species_taxid_grouper.py`** - Main parser script
   - Groups NCBI assembly data by species_taxid
   - Creates species-level statistics and isolate-preferred subsets
   - Integrates taxonomic lineage information via taxonkit

2. **`species_analysis_tools.py`** - Analysis and verification companion
   - Provides taxonomic distribution analysis
   - Analyzes isolate vs uncultured patterns
   - Identifies potential data quality issues
   - Compares with existing datasets

3. **`README.md`** - This documentation file

## Usage

### Basic Usage

```bash
# Run the main grouper script
python species_taxid_grouper.py

# Run with sample data for testing
python species_taxid_grouper.py --sample-size 10000

# Specify output directory
python species_taxid_grouper.py --output-dir results/
```

### Analysis Tools

```bash
# Run all analyses
python species_analysis_tools.py --input-dir results/

# Run specific analysis
python species_analysis_tools.py --analysis taxonomy
python species_analysis_tools.py --analysis isolates
python species_analysis_tools.py --analysis diversity
python species_analysis_tools.py --analysis issues

# Compare with existing data
python species_analysis_tools.py --comparison-file ../taxid_species_grouped.csv
```

## Output Files

### From species_taxid_grouper.py

1. **`species_grouped_statistics.csv`** - Complete species-level statistics
   - Columns: species_taxid, representative_organism_name, representative_taxid, total_genomes, isolate_count, uncultured_count, domain, kingdom, phylum, class, order, family, genus, species

2. **`species_isolate_preferred_subset.csv`** - One genome per species (isolate-preferred)
   - Contains all original columns from assembly file
   - One representative genome per species_taxid
   - Prioritizes isolate genomes when available

### From species_analysis_tools.py

1. **`species_analysis_report_YYYYMMDD_HHMMSS.txt`** - Comprehensive analysis report
   - Dataset overview and statistics
   - Top species by genome count
   - Quality assessment results

## Method Advantages

### Compared to Previous Methods

1. **Biological Accuracy**: Uses NCBI's official species_taxid instead of parsing organism names
2. **Consistency**: Eliminates variability from strain naming conventions
3. **Completeness**: 100% coverage of species_taxid in NCBI data
4. **Reproducibility**: Same results every time, regardless of data order
5. **Efficiency**: Direct grouping without complex string parsing

### Data Quality Benefits

- **No Strain Confusion**: Different strains of same species are properly grouped
- **Taxonomic Precision**: Official NCBI taxonomic assignments
- **Isolate Preference**: Prioritizes cultured isolates for downstream analysis
- **Comprehensive Lineage**: Full taxonomic hierarchy via taxonkit integration

## Requirements

### Software Dependencies

- Python 3.7+
- pandas
- numpy
- tqdm
- taxonkit (must be installed and in PATH)

### Python Packages

```bash
pip install pandas numpy tqdm matplotlib
```

### System Requirements

- taxonkit must be installed and accessible via command line
- Sufficient memory for processing large datasets (recommend 8GB+ RAM)

## Data Requirements

The scripts expect to find the NCBI assembly summary file:
- `00assembly_summary_genbank.txt`

The script will automatically search for this file in common locations:
- Current directory
- Parent directory
- `../metadata/` directory

## Performance Notes

- Processing full NCBI dataset (~3M records): ~10-15 minutes
- Memory usage: ~2-4GB for full dataset
- Taxonkit lineage retrieval: Processed in batches of 1000 for efficiency

## Troubleshooting

### Common Issues

1. **taxonkit not found**: Ensure taxonkit is installed and in your PATH
2. **Memory errors**: Use `--sample-size` for testing with smaller datasets
3. **File not found**: Check that `00assembly_summary_genbank.txt` is accessible

### Performance Optimization

- Use SSD storage for faster I/O
- Increase batch size for taxonkit if you have more memory
- Consider using `--sample-size` for development and testing

## Integration with Existing Pipeline

This method can replace existing parsing approaches while maintaining compatibility:

1. Output formats are similar to existing parsers
2. Column names follow established conventions
3. Can be integrated into existing merger scripts
4. Provides both detailed statistics and subset data

## Future Enhancements

Potential improvements for future versions:

1. **Parallel Processing**: Multi-threading for taxonkit calls
2. **Caching**: Cache taxonkit results for repeated runs
3. **Visualization**: Built-in plotting capabilities
4. **Export Formats**: Additional output formats (JSON, HDF5)
5. **Database Integration**: Direct database connectivity options

---

*Created: 2025-01-27*
*Part of the metagenomics_toolkit project*
