# Comprehensive Visual: All Filtrations and Settings Summary

## 📋 Overview
This document summarizes all filtrations, thresholds, and settings used in the mega comprehensive stacked visual script based on our conversation history and script analysis.

## 🔍 Data Filtrations Applied

### 1. **Basic Data Quality Filters**
- **Census OTU Count**: Must be > 0
- **NCBI Species Count**: Must be > 0
- **Domain Filtering**: 
  - 16S data: Filter by domain ("Bacteria" or "Archaea")
  - 18S data: All Eukaryota data included

### 2. **Taxonomic Factor Thresholds**
- **Novelty Factor Threshold**: > 1.0
  - Only taxa with novelty factor above 1.0 are considered for highlighting
  - Represents taxa that are more novel/uncultured in census vs NCBI
- **Overrepresentation Factor Threshold**: > 1.0
  - Only taxa with coverage factor above 1.0 are considered for highlighting
  - Represents taxa that are overrepresented in NCBI vs census

### 3. **Top Taxa Selection**
- **Top Novelty Taxa**: Maximum 10 per dataset
  - Selected from taxa above novelty threshold
  - Sorted by novelty factor (descending)
  - If fewer than 10 meet criteria, only those are selected
- **Top Coverage Taxa**: Maximum 10 per dataset
  - Selected from taxa above coverage threshold  
  - Sorted by overrepresentation factor (descending)
  - If fewer than 10 meet criteria, only those are selected

### 4. **Column Cleanup Filters**
**Duplicate columns removed during processing:**
- `Census_OTU_Count` (duplicate of `census_otu_count`)
- `NCBI_Species_Count` (duplicate of `ncbi_species_count`)
- `NCBI_Genome_Count` (duplicate of `ncbi_genome_count`)
- `Isolate_Count` (duplicate of `isolate_count`)
- `Novelty_Ratio` (duplicate of `novelty_factor`)
- `Coverage_Factor` (duplicate of `overrepresentation_factor`)
- `Domain` (redundant uppercase version)
- `Level` (obvious from filename)
- `coverage_percentage` (not needed)

## 📊 Data Organization Rules

### 1. **Source Data Sorting Order**
1. **Novelty Taxa First**: `Is_Top_Novelty = TRUE`
   - Sorted by novelty factor (descending)
2. **Coverage Taxa Second**: `Is_Top_Coverage = TRUE`  
   - Sorted by overrepresentation factor (descending)
3. **Other Taxa Last**: All remaining taxa
   - Sorted alphabetically by taxon name

### 2. **Circle Size Calculation**
Based on isolate percentage (isolate_count / ncbi_genome_count * 100):
- **0% isolates**: Circle size = 25
- **<10% isolates**: Circle size = 20  
- **<50% isolates**: Circle size = 15
- **≥50% isolates**: Circle size = 10

*Inverted scale: Larger circles = poorly cultured organisms*

## 🎨 Visual Settings

### 1. **Plot Dimensions**
- **Width**: 54 inches
- **Height**: 20 inches
- **DPI**: 300
- **Layout**: 2 rows × 3 columns (phylum/family × Bacteria/Archaea/Eukaryota)

### 2. **Axis Settings**
- **X-axis**: Log scale, limits [1, 10000] (census OTU count)
- **Y-axis**: Log scale, limits [1, 10000] (NCBI species count)

### 3. **Circle Aesthetics**
- **Shape**: 21 (filled circle with border)
- **Size Range**: [10, 22] in plots
- **Stroke**: 0.6 thickness
- **Alpha**: 0.9 transparency
- **Background Alpha**: 0.3

### 4. **Text and Annotations**
- **Text Size**: 11
- **Novelty Annotations**: Nudged upward (nudge_y = 0.2)
- **Coverage Annotations**: Nudged downward (nudge_y = -0.2)
- **Random Seeds**: 123 (novelty), 456 (coverage) for reproducible positioning

## 🎯 Current Data Status

### **Final Source Data Structure (22 columns):**
1. `Taxon` - taxonomic name
2. `census_otu_count` - OTU count from census
3. `census_size_count` - size count from census
4. `otu_percentage` - OTU percentage
5. `size_percentage` - size percentage
6. `ncbi_genome_count` - genome count from NCBI
7. `ncbi_species_count` - species count from NCBI
8. `isolate_count` - isolate count
9. `genome_pct_db` - genome percentage in database
10. `species_pct` - species percentage
11. `isolate_percentage` - isolate percentage
12. `novelty_factor` - novelty ratio
13. `overrepresentation_factor` - coverage factor
14. `direct_matches` - direct matches
15. `lineage_matches` - lineage matches
16. `total_matches` - total matches
17. `match_status` - match status
18. `domain` - domain (Eukaryota/Bacteria/Archaea)
19. `Circle_Size` - calculated circle size for visualization
20. `Division` - taxonomic division/phylum
21. `Is_Top_Novelty` - boolean flag for top novelty taxa
22. `Is_Top_Coverage` - boolean flag for top coverage taxa

### **Current Flag Counts:**
| File | Novelty TRUE | Coverage TRUE |
|------|-------------|---------------|
| Archaea_family_source_data.csv | 10 | 5 |
| Archaea_phylum_source_data.csv | 4 | 0 |
| Bacteria_family_source_data.csv | 10 | 10 |
| Bacteria_phylum_source_data.csv | 10 | 7 |
| Eukaryota_family_source_data.csv | 10 | 10 |
| Eukaryota_phylum_source_data.csv | 10 | 2 |

## 🔧 Configuration Management

All settings are now centralized in:
- **`config/comprehensive_visual_config.yaml`** - Main configuration file
- **`config/FILTRATIONS_AND_SETTINGS_SUMMARY.md`** - This documentation

## 🚨 Key Issues Identified

1. **Library Compatibility**: R library conflicts with GLIBCXX versions
2. **Inconsistent Flag Counts**: Some datasets have fewer than 10 taxa meeting thresholds
3. **Column Structure Variations**: Different files had different column counts before cleanup

## ✅ Solutions Implemented

1. **Removed dplyr dependencies** where possible
2. **Implemented proper threshold-based filtering** (only taxa meeting criteria are flagged)
3. **Standardized column structure** across all source data files
4. **Created centralized configuration** for easy management
