# AI-Assisted Taxonomic Resolution Workflow

This document describes the workflow for using AI to resolve unmapped taxonomic names in the 18S census data.

## Overview

The workflow uses a **cache-first approach** to resolve taxonomic names that are not found in NCBI taxonomy:

1. **Check AI Cache** - Look for previously researched resolutions
2. **Check Manual Database** - Look in `known_parents.py` 
3. **AI Research** - Generate prompts for AI to research unmapped taxa
4. **Cache Results** - Store findings for future runs

## Files

- `src/ai_resolution_cache.py` - Cache management class
- `src/cached_ai_resolver.py` - Main resolver with AI integration
- `run_ai_cache_simple.py` - Main script to run resolution
- `add_resolution_to_cache.py` - Interactive tool to add resolutions
- `batch_import_resolutions.py` - Batch import from JSON
- `cache/ai_resolutions.json` - Persistent cache file
- `example_ai_resolutions.json` - Example AI research results

## Workflow Steps

### Step 1: Run Initial Resolution

```bash
cd 00_gaps_taxonomic/parse_repaa_table/18S_censusparse/py_18S
python run_ai_cache_simple.py
```

This will:
- Extract unmapped families from the taxonkit log
- Check the AI cache for existing resolutions
- Check the manual database (`known_parents.py`)
- Generate AI research prompts for remaining unmapped taxa
- Create a review log in `cache/ai_resolution_review_*.txt`

### Step 2: Research Unmapped Taxa

The script will output prompts like this for each unmapped taxon:

```
🤖 AI SEARCH NEEDED: Dino-Group-II.U.family

TAXONOMIC RESEARCH REQUEST
==========================
Taxon Name: Dino-Group-II.U.family
Rank: family

TASK:
Please search for information about this taxonomic name and identify:
1. What is the parent taxon?
2. What is the NCBI TaxID of the parent taxon?
3. What is the full taxonomic lineage?
```

Use an AI assistant (like Claude, ChatGPT, or web search) to research each taxon.

### Step 3: Add Resolutions to Cache

#### Option A: Interactive Entry (One at a time)

```bash
python add_resolution_to_cache.py
```

Follow the prompts to enter:
- Taxon name
- Parent NCBI TaxID
- Parent taxon name
- Rank (family/genus)
- Confidence (0.0-1.0)
- Research notes
- Validation status

#### Option B: Batch Import (Multiple at once)

1. Create a JSON file with your research results (see `example_ai_resolutions.json`)
2. Import the batch:

```bash
python batch_import_resolutions.py your_resolutions.json
```

### Step 4: Re-run Resolution

After adding resolutions to the cache, re-run the script:

```bash
python run_ai_cache_simple.py
```

Now the cached resolutions will be used automatically!

### Step 5: Review and Validate

Check the review log:

```bash
cat cache/ai_resolution_review_*.txt
```

The log shows:
- ✅ **Validated resolutions** - Ready to use
- ⚠️ **Unvalidated resolutions** - Need human review
- 🔍 **Not resolved** - Still need research

### Step 6: Validate Resolutions

To mark a resolution as validated, you can:

1. Edit the cache file directly (`cache/ai_resolutions.json`)
2. Set `"validated": true` for correct resolutions
3. Or use the Python API:

```python
from pathlib import Path
from src.ai_resolution_cache import AIResolutionCache

cache = AIResolutionCache(Path("cache/ai_resolutions.json"))
cache.validate_resolution("Dino-Group-II.U.family", validator="your_name")
```

## JSON Format for Batch Import

```json
{
  "resolutions": [
    {
      "taxon_name": "MAST-12",
      "parent_taxid": "33634",
      "parent_name": "Stramenopiles",
      "confidence": 0.9,
      "research_notes": "Marine Stramenopile clade",
      "rank": "family",
      "validated": false
    }
  ]
}
```

## Cache File Structure

The cache file (`cache/ai_resolutions.json`) contains:

```json
{
  "metadata": {
    "created": "2026-03-03T10:00:00",
    "last_updated": "2026-03-03T12:00:00",
    "total_resolutions": 25,
    "version": "1.0"
  },
  "resolutions": {
    "Dino-Group-II.U.family": {
      "parent_taxid": "2864",
      "parent_name": "Dinophyceae",
      "lineage": "Eukaryota;Alveolata;Dinoflagellata;Dinophyceae;Dino-Group-II.U.family",
      "lineage_ranks": "superkingdom;kingdom;phylum;class;family",
      "lineage_taxids": "2759;33630;2864;2864;NA",
      "rank": "family",
      "source": "AI-web-search",
      "confidence": 0.95,
      "research_date": "2026-03-03",
      "research_notes": "Well-documented environmental clade",
      "validated": true,
      "validator": "human",
      "validation_date": "2026-03-03"
    }
  }
}
```

## Tips for AI Research

When researching taxa:

1. **Search terms**: Use taxon name + "taxonomy" + "18S rRNA"
2. **Look for**:
   - Scientific publications
   - NCBI Taxonomy Browser
   - SILVA, PR2, or other rRNA databases
   - Environmental clade documentation
3. **For environmental clades** (e.g., "MAST-12", "Dino-Group-II"):
   - Find the broader taxonomic group they belong to
   - Look for phylogenetic placement studies
4. **Confidence levels**:
   - 0.9-1.0: Well-documented, clear parent
   - 0.7-0.9: Good evidence, some uncertainty
   - 0.5-0.7: Tentative placement
   - <0.5: Highly uncertain

## Integration with Pipeline

Once you have validated resolutions, they can be:

1. **Used directly** - The cache is checked on every run
2. **Exported to known_parents.py** - For permanent inclusion
3. **Applied to CSV files** - Via the systematic resolver

## Troubleshooting

**Problem**: Taxonkit fails to get lineage for parent TaxID

**Solution**: Verify the TaxID exists in NCBI:
```bash
echo "2864" | taxonkit lineage -R -t
```

**Problem**: Cache not being used

**Solution**: Check cache file exists and is valid JSON:
```bash
cat cache/ai_resolutions.json | python -m json.tool
```

**Problem**: Duplicate entries

**Solution**: The script will skip taxa already in cache. To overwrite, delete the entry from the cache file first.

