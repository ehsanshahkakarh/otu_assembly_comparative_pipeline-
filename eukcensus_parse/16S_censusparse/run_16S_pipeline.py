#!/usr/bin/env python3
"""
16S Census Parser - Main Pipeline Entry Point
==============================================

Simple wrapper script to run the 16S census parser with the new modular structure.

Usage:
    # Run with default input file
    python run_16S_pipeline.py
    
    # Run with custom input file
    python run_16S_pipeline.py custom_input.tsv
    
    # Run with custom input and output prefix
    python run_16S_pipeline.py custom_input.tsv custom_output

This script processes EukCensus 16S cluster data and generates taxonomic summaries
at division (phylum), family, and genus levels with NCBI taxonomy integration.

Output files:
    - output/eukcensus16S_by_division.csv
    - output/eukcensus16S_by_family.csv
    - output/eukcensus16S_by_genus.csv
    - logs/eukcensus16S_processing.log
    - logs/eukcensus16S_comprehensive_unmapped.log
"""

import sys
from pathlib import Path

# Add the parent directory to the path so we can import src as a package
sys.path.insert(0, str(Path(__file__).parent))

# Import and run the main parser
from src.run_census_parser import main

if __name__ == "__main__":
    print("=" * 80)
    print("16S CENSUS PARSER")
    print("=" * 80)
    print()
    
    # Run the parser
    main()

