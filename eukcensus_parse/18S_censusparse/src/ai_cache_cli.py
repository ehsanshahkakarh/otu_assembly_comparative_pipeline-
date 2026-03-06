#!/usr/bin/env python3
"""
AI Cache Command-Line Interface
================================

Unified CLI for managing AI-researched taxonomic resolutions.

Commands:
  - add: Add a single resolution interactively
  - import: Batch import from JSON file
  - validate: Validate resolutions interactively
  - stats: Show cache statistics
  - export: Export validated resolutions
"""

import sys
import json
import argparse
from pathlib import Path
import logging

# Handle imports for both module and script execution
try:
    from .ai_resolution_cache import AIResolutionCache, get_default_cache_path
    from .cached_ai_resolver import CachedAIResolver
except ImportError:
    # Running as script, add src to path
    sys.path.insert(0, str(Path(__file__).parent))
    from ai_resolution_cache import AIResolutionCache, get_default_cache_path
    from cached_ai_resolver import CachedAIResolver

logger = logging.getLogger(__name__)


def cmd_add(cache_file: Path):
    """Interactive command to add a single resolution"""
    resolver = CachedAIResolver(cache_file)
    
    print("\n" + "="*80)
    print("ADD RESOLUTION TO CACHE")
    print("="*80)
    
    taxon_name = input("\nTaxon name: ").strip()
    if not taxon_name:
        print("❌ Taxon name required")
        return
    
    if resolver.cache.has_resolution(taxon_name):
        print(f"⚠️  '{taxon_name}' already in cache!")
        if input("Overwrite? (yes/no): ").strip().lower() not in ['yes', 'y']:
            return
    
    parent_taxid = input("Parent NCBI TaxID: ").strip()
    parent_name = input("Parent taxon name: ").strip()
    rank = input("Rank (family/genus) [family]: ").strip() or "family"
    
    try:
        confidence = float(input("Confidence (0.0-1.0) [0.8]: ").strip() or "0.8")
    except ValueError:
        print("❌ Invalid confidence value")
        return
    
    print("\nResearch notes (press Enter twice when done):")
    notes_lines = []
    while True:
        line = input()
        if not line:
            break
        notes_lines.append(line)
    research_notes = " ".join(notes_lines) or f"Manually added resolution for {taxon_name}"
    
    validated = input("\nMark as validated? (yes/no) [no]: ").strip().lower() in ['yes', 'y']
    
    success = resolver.add_manual_resolution(
        taxon_name=taxon_name,
        parent_taxid=parent_taxid,
        parent_name=parent_name,
        confidence=confidence,
        research_notes=research_notes,
        rank=rank,
        validated=validated
    )
    
    if success:
        print(f"\n✅ Successfully added: {taxon_name} → {parent_name}")
    else:
        print(f"\n❌ Failed to add resolution")


def cmd_import(cache_file: Path, input_file: Path):
    """Batch import resolutions from JSON file"""
    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return
    
    resolver = CachedAIResolver(cache_file)
    
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    resolutions = data.get('resolutions', [])
    logger.info(f"Importing {len(resolutions)} resolutions from {input_file}")
    
    success_count = 0
    skip_count = 0
    fail_count = 0
    
    for i, res in enumerate(resolutions, 1):
        taxon_name = res.get('taxon_name')
        if not taxon_name:
            skip_count += 1
            continue
        
        if resolver.cache.has_resolution(taxon_name):
            logger.info(f"[{i}/{len(resolutions)}] Skipping {taxon_name} (already in cache)")
            skip_count += 1
            continue
        
        logger.info(f"[{i}/{len(resolutions)}] Importing: {taxon_name}")
        
        success = resolver.add_manual_resolution(
            taxon_name=taxon_name,
            parent_taxid=res.get('parent_taxid', ''),
            parent_name=res.get('parent_name', ''),
            confidence=float(res.get('confidence', 0.5)),
            research_notes=res.get('research_notes', ''),
            rank=res.get('rank', 'family'),
            validated=res.get('validated', False)
        )
        
        if success:
            success_count += 1
        else:
            fail_count += 1
    
    logger.info(f"\n✅ Imported: {success_count}, ⏭️  Skipped: {skip_count}, ❌ Failed: {fail_count}")


def cmd_stats(cache_file: Path):
    """Show cache statistics"""
    cache = AIResolutionCache(cache_file)
    stats = cache.get_statistics()

    print("\n" + "="*80)
    print("AI CACHE STATISTICS")
    print("="*80)
    print(f"Cache file: {cache_file}")
    print(f"Total resolutions: {stats['total_resolutions']}")
    print(f"  ✅ Validated: {stats['validated']}")
    print(f"  ⚠️  Unvalidated: {stats['unvalidated']}")
    print(f"Average confidence: {stats['average_confidence']:.2f}")
    print(f"Last updated: {stats['last_updated']}")
    print("="*80)


def cmd_validate(cache_file: Path, validator: str):
    """Interactively validate unvalidated resolutions"""
    cache = AIResolutionCache(cache_file)
    unvalidated = cache.get_unvalidated()

    if not unvalidated:
        print("✅ All resolutions are validated!")
        return

    print(f"\n📋 Found {len(unvalidated)} unvalidated resolutions\n")

    for i, taxon_name in enumerate(unvalidated, 1):
        resolution = cache.get_resolution(taxon_name)

        print("="*80)
        print(f"[{i}/{len(unvalidated)}] {taxon_name}")
        print("="*80)
        print(f"Parent TaxID:  {resolution['parent_taxid']}")
        print(f"Parent Name:   {resolution['parent_name']}")
        print(f"Rank:          {resolution['rank']}")
        print(f"Confidence:    {resolution['confidence']:.2f}")
        print(f"Research Date: {resolution['research_date']}")
        print(f"\nLineage:\n{resolution['lineage']}")
        if resolution.get('research_notes'):
            print(f"\nNotes: {resolution['research_notes']}")
        print("="*80)

        response = input("\nValidate this resolution? (y/n/q to quit): ").strip().lower()

        if response == 'q':
            print("Validation session ended.")
            break
        elif response == 'y':
            cache.validate_resolution(taxon_name, validator)
            print(f"✅ Validated: {taxon_name}\n")
        else:
            print(f"⏭️  Skipped: {taxon_name}\n")

    # Show final stats
    stats = cache.get_statistics()
    print("\n" + "="*80)
    print(f"Validation complete: {stats['validated']}/{stats['total_resolutions']} validated")
    print("="*80)


def cmd_list(cache_file: Path, unvalidated_only: bool, validated_only: bool):
    """List resolutions"""
    cache = AIResolutionCache(cache_file)

    if unvalidated_only:
        taxa = cache.get_unvalidated()
        title = "UNVALIDATED RESOLUTIONS"
    elif validated_only:
        taxa = cache.get_validated()
        title = "VALIDATED RESOLUTIONS"
    else:
        taxa = list(cache.cache["resolutions"].keys())
        title = "ALL RESOLUTIONS"

    print("="*80)
    print(title)
    print("="*80)

    for taxon_name in sorted(taxa):
        resolution = cache.get_resolution(taxon_name)
        status = "✅" if resolution.get('validated') else "⏳"
        print(f"{status} {taxon_name:40s} → {resolution['parent_name']:30s} (TaxID: {resolution['parent_taxid']})")

    print("="*80)
    print(f"Total: {len(taxa)}")
    print("="*80)


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description='AI Cache Management CLI',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Add command
    subparsers.add_parser('add', help='Add a single resolution interactively')
    
    # Import command
    import_parser = subparsers.add_parser('import', help='Batch import from JSON')
    import_parser.add_argument('file', type=Path, help='JSON file to import')
    
    # Stats command
    subparsers.add_parser('stats', help='Show cache statistics')

    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Interactively validate unvalidated resolutions')
    validate_parser.add_argument('--validator', default='human', help='Name of validator (default: human)')

    # List command
    list_parser = subparsers.add_parser('list', help='List resolutions')
    list_parser.add_argument('--unvalidated', action='store_true', help='Show only unvalidated')
    list_parser.add_argument('--validated', action='store_true', help='Show only validated')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return

    cache_file = get_default_cache_path()
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    if args.command == 'add':
        cmd_add(cache_file)
    elif args.command == 'import':
        cmd_import(cache_file, args.file)
    elif args.command == 'stats':
        cmd_stats(cache_file)
    elif args.command == 'validate':
        cmd_validate(cache_file, args.validator)
    elif args.command == 'list':
        cmd_list(cache_file, args.unvalidated, args.validated)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
    main()

