# AI-Assisted Taxonomic Resolution - Quick Start Guide

## 🚀 Get Started in 5 Minutes

### Step 1: Test the System

```bash
cd 00_gaps_taxonomic/parse_repaa_table/18S_censusparse/py_18S
python test_ai_cache.py
```

Expected output:
```
✅ Resolver initialized
Total resolved: 3/3
Cache statistics: 13 resolutions
```

### Step 2: Import Example Resolutions

```bash
python batch_import_resolutions.py example_ai_resolutions.json
```

This imports 13 pre-researched taxonomic resolutions into the cache.

### Step 3: Run the Main Workflow

```bash
python run_ai_cache_simple.py
```

This will:
- Extract unmapped families from the taxonkit log
- Check AI cache and manual database
- Generate research prompts for remaining taxa
- Create a review log

### Step 4: Review Results

```bash
ls -lh cache/
cat cache/ai_resolution_review_*.txt
```

The review log shows:
- ✅ Validated resolutions (ready to use)
- ⚠️ Unvalidated resolutions (need review)
- 🔍 Not resolved (need research)

## 📝 Adding Your Own Resolutions

### Option A: Interactive Entry (Recommended for 1-5 taxa)

```bash
python add_resolution_to_cache.py
```

Follow the prompts:
1. Enter taxon name
2. Enter parent NCBI TaxID
3. Enter parent taxon name
4. Enter confidence (0.0-1.0)
5. Enter research notes

### Option B: Batch Import (Recommended for 5+ taxa)

1. Create a JSON file (e.g., `my_resolutions.json`):

```json
{
  "resolutions": [
    {
      "taxon_name": "Your-Taxon-Name",
      "parent_taxid": "12345",
      "parent_name": "Parent Taxon",
      "confidence": 0.9,
      "research_notes": "Found in literature XYZ",
      "rank": "family",
      "validated": false
    }
  ]
}
```

2. Import:

```bash
python batch_import_resolutions.py my_resolutions.json
```

## ✅ Validating Resolutions

After reviewing resolutions, mark them as validated:

```bash
python validate_resolutions.py
```

This interactive tool lets you:
- Review each unvalidated resolution
- Validate correct ones
- Skip uncertain ones
- Track who validated what

## 🔍 How to Research a Taxon

When you see a research prompt like:

```
🤖 AI SEARCH NEEDED: Dino-Group-II.U.family
```

1. **Search Google Scholar** or **PubMed**:
   - "Dino-Group-II taxonomy"
   - "Dino-Group-II 18S rRNA"
   - "Dino-Group-II phylogeny"

2. **Check NCBI Taxonomy Browser**:
   - https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

3. **Look for**:
   - What broader group does this belong to?
   - What is the NCBI TaxID of the parent?
   - Is this an environmental clade or formal taxon?

4. **Add to cache** using one of the methods above

## 📊 Checking Progress

```bash
# View cache statistics
python -c "
from pathlib import Path
from src.ai_resolution_cache import AIResolutionCache
cache = AIResolutionCache(Path('cache/ai_resolutions.json'))
stats = cache.get_statistics()
print(f'Total: {stats[\"total_resolutions\"]}')
print(f'Validated: {stats[\"validated\"]}')
print(f'Pending: {stats[\"unvalidated\"]}')
"
```

## 🎯 Common Workflows

### Workflow 1: Research One Taxon

```bash
# 1. Note the taxon name from run_ai_cache_simple.py output
# 2. Research it (Google Scholar, NCBI, etc.)
# 3. Add to cache
python add_resolution_to_cache.py
# 4. Re-run to verify
python run_ai_cache_simple.py
```

### Workflow 2: Batch Research with AI

```bash
# 1. Run main script to get research prompts
python run_ai_cache_simple.py > research_prompts.txt

# 2. Use AI (ChatGPT, Claude, etc.) to research all taxa
#    Copy prompts from research_prompts.txt

# 3. AI generates JSON with results
#    Save as ai_research_results.json

# 4. Import batch
python batch_import_resolutions.py ai_research_results.json

# 5. Validate
python validate_resolutions.py
```

### Workflow 3: Incremental Progress

```bash
# Day 1: Research 5 taxa
python add_resolution_to_cache.py
# (add 5 resolutions)

# Day 2: Research 5 more
python add_resolution_to_cache.py
# (add 5 more)

# Day 3: Validate all
python validate_resolutions.py
```

## 🐛 Troubleshooting

**Problem**: "Cache file not found"
```bash
# Solution: Create cache directory
mkdir -p cache
python test_ai_cache.py
```

**Problem**: "Failed to get lineage for parent TaxID"
```bash
# Solution: Verify TaxID exists in NCBI
echo "YOUR_TAXID" | taxonkit lineage -R -t
```

**Problem**: "Import failed"
```bash
# Solution: Validate JSON syntax
cat your_file.json | python -m json.tool
```

## 📚 Files Reference

| File | Purpose |
|------|---------|
| `test_ai_cache.py` | Quick system test |
| `run_ai_cache_simple.py` | Main workflow |
| `add_resolution_to_cache.py` | Interactive entry |
| `batch_import_resolutions.py` | Batch import |
| `validate_resolutions.py` | Validation tool |
| `cache/ai_resolutions.json` | Persistent cache |
| `example_ai_resolutions.json` | Example data |

## 🎓 Next Steps

1. ✅ Test the system (`test_ai_cache.py`)
2. ✅ Import examples (`batch_import_resolutions.py`)
3. 📝 Research unmapped taxa
4. ✅ Validate resolutions (`validate_resolutions.py`)
5. 🔄 Integrate with pipeline

## 💡 Tips

- Start with high-confidence, well-documented taxa
- Use confidence scores to prioritize validation
- Keep research notes detailed for future reference
- Validate in batches after researching multiple taxa
- The cache prevents duplicate work - use it!

## 📖 Full Documentation

- `AI_RESOLUTION_WORKFLOW.md` - Complete workflow guide
- `AI_CACHE_IMPLEMENTATION_SUMMARY.md` - Technical details

