# Test Species-Based Merger

## Overview

This is a **test implementation** of the merger workflow that uses the **correct approach**:

1. **Census files are the BACKBONE** - Use census taxa as search terms
2. **Build synonym dictionary** from NCBI `names.dmp` for flexible name matching
3. **Search species_grouped_*.csv** for species belonging to each census taxon
4. **Aggregate** matched species genome counts

## Why This Approach?

The census files contain taxa that may not align with true NCBI taxonomic ranks (e.g., environmental clades, paraphyletic groups). Therefore, we cannot simply extract phylum/family/genus from species lineages and aggregate. Instead, we must:

- Use census taxa as the search backbone
- Search the comprehensive NCBI species data for each census taxon
- Handle taxonomic name changes and synonyms using `names.dmp`
- Aggregate genome counts from all matching species

## Modules

### `species_data_loader.py`
Loads the `species_grouped_*.csv` file (162K species with full lineage information).

### `census_synonym_builder.py`
Builds a comprehensive synonym dictionary for each census taxon using NCBI `names.dmp`:
- Maps each census name to all possible matching names
- Handles taxonomic name changes (e.g., Euryarchaeota → Methanobacteriota)
- Uses the pre-built synonym dictionary (1.7M+ mappings)

### `species_searcher.py`
Searches `species_grouped_*.csv` for species belonging to a census taxon using hierarchical matching:

1. **Direct taxid match**: Census taxid in `lineage_taxids`
2. **Name match**: Any possible name in `lineage` string
3. **Rank extraction**: (Disabled for performance) Extract taxon at specified rank

### `species_aggregator.py`
Aggregates matched species to calculate:
- Total genome count
- Isolate genome count
- Uncultured genome count
- Isolate percentage
- Species count

## Test Runner

`test_species_merger.py` - Main test script

**Usage:**
```bash
# Test division level (18S)
python3 test_species_merger.py --level division --census-type 18S

# Test family level (18S)
python3 test_species_merger.py --level family --census-type 18S

# Test genus level (18S)
python3 test_species_merger.py --level genus --census-type 18S

# Test with debug logging
python3 test_species_merger.py --level division --census-type 18S --debug
```

## Test Results (Division Level, 18S)

**Performance:**
- Completed in ~42 seconds
- Processed 162,731 species
- ~4,000 species/second

**Matching:**
- Matched: 14/22 census divisions (63.6%)
- Census-only: 8 divisions (no NCBI genomes)

**Top Matches:**
- Eukaryota.U.division: 58,939 genomes (24,044 species)
- Opisthokonta: 47,927 genomes (19,748 species)
- Streptophyta: 7,831 genomes (3,039 species)

**Census-only taxa:**
- TSAR.U.division, Archaeplastida.U.division
- Centroplasthelida, Rigifilida, Picozoa
- Hemimastigophora, Nibbleridia, Nebulidia

## Output

**Location:** `test_outputs/test_18s_species_merged_division_TIMESTAMP.csv`

**Columns:**
- `division` - Census taxon name
- `census_taxid` - Census taxon ID (clean integer string)
- `census_otu_count` - Census OTU count
- `census_size_count` - Census size count
- `ncbi_genome_count` - Total NCBI genomes
- `ncbi_isolate_count` - Isolate genomes
- `ncbi_uncultured_count` - Uncultured genomes
- `ncbi_isolate_percentage` - Percentage of isolate genomes
- `ncbi_species_count` - Number of species
- **`novelty_factor`** - census_otu_count / ncbi_species_count (environmental diversity vs genomic coverage)
- **`overrepresentation_factor`** - ncbi_species_count / census_otu_count (database bias)
- `match_status` - `matched` or `census_only`

### Key Metrics Interpretation

**Novelty Factor** = census_otu_count / ncbi_species_count
- **High (>10)**: Many environmental OTUs, few genomes → **Priority for sequencing**
- **Low (<2)**: Well-represented in databases → Good genomic coverage

**Overrepresentation Factor** = ncbi_species_count / census_otu_count
- **High (>1)**: More genomes than environmental OTUs → Database bias toward cultured taxa
- **Low (<1)**: Fewer genomes than environmental OTUs → Underrepresented

## Next Steps

1. ✅ Test division level (DONE)
2. ⏳ Test family level
3. ⏳ Test genus level
4. ⏳ Compare with current merger output
5. ⏳ Validate correctness
6. ⏳ Integrate into main merger workflow

## Performance Notes

- **Optimized**: Disabled slow row-by-row rank extraction (redundant with name matching)
- **Vectorized**: Uses pandas string operations for taxid/name matching
- **Progress tracking**: Logs every 10 taxa processed
- **Scalable**: Can process 162K species in ~42 seconds

## Dependencies

- pandas
- NCBI `names.dmp` file (via `synonym_dict.py`)
- `species_grouped_*.csv` file from `nev_parse_meth`
- Census files from `18S_censusparse` or `16S_censusparse`

