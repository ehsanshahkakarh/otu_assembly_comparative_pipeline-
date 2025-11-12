# GTDB CSV Data Directory

This directory contains processed GTDB (Genome Taxonomy Database) data in CSV format, organized by taxonomic ranks and including both count summaries and detailed accession information.

## Main Data Files

### Count Files
These files contain taxonomic counts at different ranks:

1. `gtdb_phylum_counts.csv` (11MB)
   - Columns:
     - `phylum`: Taxonomic phylum name
     - `gtdb_genome_count`: Number of genomes in that phylum
   - Contains comprehensive list of bacterial and archaeal phyla
   - Example: Pseudomonadota (214,930 genomes), Bacillota (82,709 genomes)

2. `gtdb_family_counts.csv` (11MB)
   - Columns:
     - `family`: Taxonomic family name
     - `gtdb_genome_count`: Number of genomes in that family
   - Contains family-level taxonomic classifications

3. `gtdb_genus_counts.csv` (11MB)
   - Columns:
     - `genus`: Taxonomic genus name
     - `gtdb_genome_count`: Number of genomes in that genus
   - Contains genus-level taxonomic classifications

### Accession Files
These files contain the same taxonomic information plus accession numbers:

1. `gtdb_phylum_with_accessions.csv` (22MB)
   - Columns:
     - `phylum`: Taxonomic phylum name
     - `accession_clean`: Cleaned accession number
     - `domain`: Domain classification (Bacteria/Archaea)
     - Additional metadata columns

2. `gtdb_family_with_accessions.csv` (23MB)
   - Columns:
     - `family`: Taxonomic family name
     - `accession_clean`: Cleaned accession number
     - `domain`: Domain classification
     - Additional metadata columns

3. `gtdb_genus_with_accessions.csv` (21MB)
   - Columns:
     - `genus`: Taxonomic genus name
     - `accession_clean`: Cleaned accession number
     - `domain`: Domain classification
     - Additional metadata columns

## Metadata Libraries

Located in the `mtdata_libraries/` subdirectory:

1. `gtdb_accession_phylum_table.csv` (17MB)
   - Maps accessions to phyla
   - Used for quick lookup of phylum assignments

2. `gtdb_phylum_counts_polars_fixed.csv`
   - Processed version of phylum counts
   - Optimized for Polars framework
   - Contains same data as `gtdb_phylum_counts.csv` but in a different format

3. `gtdb_phylum_counts_dask.csv`
   - Processed version of phylum counts
   - Optimized for Dask framework
   - Contains same data as `gtdb_phylum_counts.csv` but in a different format

## Logs

Located in the `logs/` subdirectory:
- Contains timestamped log files from parsing operations
- Format: `gtdb_phylum_parse_YYYYMMDD_HHMMSS.log`
- Records processing steps, errors, and statistics

## Usage Notes

1. **Count Files**
   - Use for quick taxonomic distribution analysis
   - Suitable for summary statistics and visualization
   - Smaller file size, faster to process

2. **Accession Files**
   - Use when detailed genome information is needed
   - Required for mapping genomes to taxonomic classifications
   - Contains additional metadata for each genome

3. **Metadata Libraries**
   - Use for specific processing frameworks (Polars/Dask)
   - Optimized for performance in respective environments
   - Contains processed versions of the main data files

## Data Processing

The files in this directory are generated from raw GTDB data and processed to:
- Clean and standardize taxonomic names
- Remove version numbers from accessions
- Add domain classifications
- Generate count summaries
- Create optimized versions for different processing frameworks

## File Size Considerations

- Count files are approximately 11MB each
- Accession files are approximately 21-23MB each
- Total directory size is around 100MB
- Consider using appropriate tools for handling large CSV files (e.g., Polars, Dask, or pandas with chunking) 