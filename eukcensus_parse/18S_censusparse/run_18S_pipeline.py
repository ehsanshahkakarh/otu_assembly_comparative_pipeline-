#!/usr/bin/env python3
"""
18S Census Parser Pipeline Orchestrator

Orchestrates the complete 18S parsing workflow:
1. Clean taxonkit parser (NCBI-only lookups)
2. Systematic resolver (custom resolutions for unmapped families)
3. AI cache management (research and cache unmapped taxa)
4. Merge results (apply resolutions to CSV files)
"""

import sys
import argparse
import logging
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))


def setup_logging():
    """Set up basic logging for the pipeline."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )


def run_taxonkit_parser():
    """Run the clean taxonkit parser."""
    logging.info("\n" + "=" * 80)
    logging.info("STEP 1: Running Clean Taxonkit Parser")
    logging.info("=" * 80)

    from src.pipeline_taxonkit import run_taxonkit_pipeline
    run_taxonkit_pipeline()


def run_systematic_resolver():
    """Run the systematic family resolver."""
    logging.info("\n" + "=" * 80)
    logging.info("STEP 2: Running Systematic Family Resolver")
    logging.info("=" * 80)

    from src.pipeline_resolver import run_resolver_pipeline
    run_resolver_pipeline()


def run_ai_cache():
    """Run AI cache resolution workflow."""
    logging.info("\n" + "=" * 80)
    logging.info("STEP 3: AI Cache Resolution")
    logging.info("=" * 80)

    from src.cached_ai_resolver import CachedAIResolver
    from src.ai_resolution_cache import get_default_cache_path

    # Setup paths
    base_dir = Path(__file__).parent
    cache_file = get_default_cache_path()
    log_dir = base_dir / "logs"
    unmapped_log = log_dir / "eukcensus_taxonkit_only_unmapped_from_taxonkit.log"

    if not unmapped_log.exists():
        logging.warning(f"Unmapped log not found: {unmapped_log}")
        logging.info("Run 'taxonkit' step first to generate unmapped log")
        return

    # Extract unmapped families
    unmapped_families = []
    with open(unmapped_log, 'r') as f:
        in_family_section = False
        for line in f:
            if '=== FAMILY LEVEL UNMAPPED NAMES ===' in line:
                in_family_section = True
                continue
            elif '===' in line and in_family_section:
                break
            elif in_family_section and line.startswith('FAMILY |'):
                parts = line.split('|')
                if len(parts) >= 2:
                    family_name = parts[1].strip()
                    if family_name and family_name != 'original_name':
                        unmapped_families.append(family_name)

    if not unmapped_families:
        logging.info("No unmapped families found")
        return

    logging.info(f"Found {len(unmapped_families)} unmapped families")

    # Run AI resolver
    resolver = CachedAIResolver(cache_file)
    resolutions = resolver.resolve_unmapped_taxa(unmapped_families, rank="family")

    # Print summary
    stats = resolver.cache.get_statistics()
    logging.info(f"\nCache Statistics:")
    logging.info(f"  Total resolutions: {stats['total_resolutions']}")
    logging.info(f"  Validated: {stats['validated']}")
    logging.info(f"  Unvalidated: {stats['unvalidated']}")

    logging.info(f"\nResolved {len(resolutions)}/{len(unmapped_families)} families")
    logging.info(f"Cache file: {cache_file}")


def main():
    """Main pipeline orchestrator."""
    parser = argparse.ArgumentParser(
        description='18S Census Parser Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline (all steps)
  python run_18S_pipeline.py --all

  # Run individual steps
  python run_18S_pipeline.py --step taxonkit
  python run_18S_pipeline.py --step resolve
  python run_18S_pipeline.py --step ai-cache

  # Run taxonkit parser and AI cache
  python run_18S_pipeline.py --step taxonkit --step ai-cache

  # Run just AI cache (after taxonkit has been run)
  python run_18S_pipeline.py --step ai-cache
        """
    )
    
    parser.add_argument(
        '--step',
        action='append',
        choices=['taxonkit', 'resolve', 'ai-cache'],
        help='Run specific pipeline step(s). Can be specified multiple times.'
    )

    parser.add_argument(
        '--all',
        action='store_true',
        help='Run all pipeline steps in sequence'
    )
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging()
    
    logging.info("=" * 80)
    logging.info("18S CENSUS PARSER PIPELINE")
    logging.info("=" * 80)
    
    # Determine which steps to run
    steps_to_run = []

    if args.all:
        steps_to_run = ['taxonkit', 'resolve', 'ai-cache']
    elif args.step:
        steps_to_run = args.step
    else:
        # Default: run all steps
        logging.info("No steps specified, running all steps by default")
        steps_to_run = ['taxonkit', 'resolve']

    logging.info(f"Pipeline steps to run: {', '.join(steps_to_run)}")

    # Run steps in order
    try:
        if 'taxonkit' in steps_to_run:
            run_taxonkit_parser()

        if 'resolve' in steps_to_run:
            run_systematic_resolver()

        if 'ai-cache' in steps_to_run:
            run_ai_cache()
        
        logging.info("\n" + "=" * 80)
        logging.info("PIPELINE COMPLETE")
        logging.info("=" * 80)
        logging.info("Final output files:")
        logging.info("  - csv_outputs/eukcensus_18S_by_division.csv")
        logging.info("  - csv_outputs/eukcensus_18S_by_family.csv")
        logging.info("  - csv_outputs/eukcensus_18S_by_genus.csv")
        logging.info("  - logs/eukcensus_18S_unmapped_final.log")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

