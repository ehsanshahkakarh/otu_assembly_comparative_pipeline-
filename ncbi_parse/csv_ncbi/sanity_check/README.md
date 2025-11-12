# NCBI Taxonomy Sanity Check Scripts

## Purpose
This directory contains scripts for validating and analyzing NCBI taxonomy data, focusing on extracting and organizing taxonomic groups across different domains (Eukaryota, Bacteria, Archaea, and Viruses). Each script processes NCBI taxonomy CSV files to create clean, deduplicated lists of taxa at different taxonomic levels (phylum, family, genus) with their corresponding NCBI taxonomy IDs.

## File Structure

### Input Files
Located in the parent directory (`csv_ncbi/`):
- `ncbi_genus_counts.csv`: Contains genus-level taxonomy data
- `ncbi_family_counts.csv`: Contains family-level taxonomy data
- `ncbi_phylum_counts.csv`: Contains phylum-level taxonomy data

### Output Files
Generated in the `taxon_names/` subdirectory:
- Domain-specific CSV files for each taxonomic level:
  - **Eukaryota**: `eukaryotic_genera.csv`, `eukaryotic_families.csv`, `eukaryotic_phyla.csv`
  - **Bacteria**: `bacterial_genera.csv`, `bacterial_families.csv`, `bacterial_phyla.csv`
  - **Archaea**: `archaeal_genera.csv`, `archaeal_families.csv`, `archaeal_phyla.csv`
  - **Viruses**: `viral_genera.csv`, `viral_families.csv`, `viral_phyla.csv`

## Scripts Overview

The directory contains four main scripts, each focusing on a specific domain:

1. **`Eukaryota_groups.py`**: Extracts eukaryotic taxonomic groups
2. **`Bacteria_groups.py`**: Extracts bacterial taxonomic groups
3. **`Archaea_groups.py`**: Extracts archaeal taxonomic groups
4. **`Virus_groups.py`**: Extracts viral taxonomic groups

All scripts follow the same structure and workflow but filter for different domains.

## Code Documentation

### 1. Imports and Setup
```python
import pandas as pd
from pathlib import Path
```
- Uses pandas for data manipulation
- Uses pathlib for cross-platform path handling

### 2. Main Functions
Each script has a main function that orchestrates the process of extracting taxonomic data:
- `process_taxonomic_groups()` in Eukaryota_groups.py
- `process_bacterial_groups()` in Bacteria_groups.py
- `process_archaeal_groups()` in Archaea_groups.py
- `process_viral_groups()` in Virus_groups.py

#### Path Setup
```python
script_dir = Path(__file__).resolve().parent
parent_dir = script_dir.parent
output_dir = script_dir / "taxon_names"
```
- Determines script location
- Sets up paths for input and output files
- Creates output directory if it doesn't exist

#### File Definitions
```python
input_files = {
    'genus': parent_dir / "ncbi_genus_counts.csv",
    'family': parent_dir / "ncbi_family_counts.csv",
    'phylum': parent_dir / "ncbi_phylum_counts.csv"
}

output_files = {
    'genus': output_dir / "domain_genera.csv",
    'family': output_dir / "domain_families.csv",
    'phylum': output_dir / "domain_phyla.csv"
}
```
- Maps taxonomic levels to their corresponding input/output files
- Uses consistent naming convention for easy maintenance
- Output filenames are prefixed with the domain name (e.g., "eukaryotic_", "bacterial_")

#### Processing Loop
For each taxonomic level (genus, family, phylum):

1. **Data Loading**:
   - Reads CSV file using pandas
   - Validates required columns ('domain', taxonomic level, 'taxid')

2. **Data Filtering**:
   - Filters for specific domain (case insensitive, handles whitespace)
   - For viruses, accepts both 'virus' and 'viruses' as domain values

3. **Data Processing**:
   - Extracts unique taxa names and their corresponding taxids
   - Handles taxids stored as lists by extracting the first taxid
   - Creates a clean data structure with 'taxon_name' and 'taxid' columns
   - Removes duplicates and sorts alphabetically by taxon name

4. **Output Generation**:
   - Saves to CSV file with headers
   - Includes both taxon name and taxid in a tabular format
   - Preserves the relationship between taxa and their taxonomy IDs

### 3. Error Handling
The scripts include comprehensive error handling for:
- File not found errors
- Missing required columns
- Data processing errors
- Provides clear error messages with emoji indicators

## Usage
Run any script from the command line:
```bash
python Eukaryota_groups.py
python Bacteria_groups.py
python Archaea_groups.py
python Virus_groups.py
```

## Output Format
Each output CSV file contains:
- Header row with column names: 'taxon_name', 'taxid'
- One row per unique taxon
- Alphabetically sorted by taxon name
- No duplicates
- Properly formatted CSV with comma separators

## Dependencies
- Python 3.x
- pandas
- pathlib (standard library)

## Terminal Output
The scripts use emoji indicators for different types of messages:
- 📥 Processing start
- 🔍 Processing specific level
- ✅ Successful completion
- ❌ Error messages

Terminal output shows only filenames (not full paths) for cleaner display:
```
📥 Processing archaeal taxonomic groups...

🔍 Processing genus level...
  Input: ncbi_genus_counts.csv
  Output: archaeal_genera.csv
✅ Saved 123 unique Archaeal genus with taxids to: archaeal_genera.csv
```

## Configuration and Customization

### Adding a New Domain
To add support for a new domain:
1. Create a new script following the pattern of existing scripts
2. Modify the domain filtering condition to match the new domain
3. Update output filenames to reflect the new domain

### Changing Output Format
To modify the output format:
1. Edit the output file extension in the `output_files` dictionary
2. Adjust the data processing section to format data as needed
3. Modify the `to_csv()` call with appropriate parameters

### Adding New Taxonomic Levels
To support additional taxonomic levels:
1. Add the new level to the `input_files` and `output_files` dictionaries
2. Ensure the input CSV files contain the required columns
3. No other code changes are needed as the processing loop handles all levels

## Notes
- All scripts are case-insensitive when filtering for domains
- Whitespace is handled automatically in domain names
- Output files are overwritten if they already exist
- The scripts create the output directory if it doesn't exist
- Taxids are extracted from the input files, which may contain them as lists
- The first taxid is used when multiple are available for a single taxon