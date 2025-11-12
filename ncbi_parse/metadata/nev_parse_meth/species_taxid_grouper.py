#!/usr/bin/env python3
"""
Species Taxid Grouper - New Reproducible Method
Created: 2025-01-27

This script groups NCBI assembly data by identical species_taxid values to create
a reproducible, species-level dataset. It implements isolate-preferential selection
and generates comprehensive taxonomic information using taxonkit.

Key Features:
- Groups by species_taxid for true biological species-level analysis
- Isolate-preferential selection (chooses isolate genomes over uncultured when available)
- Comprehensive taxonomic lineage information via taxonkit
- Multiple output formats for different analysis needs
- Progress tracking with tqdm
- Detailed statistics and verification

Usage:
    python species_taxid_grouper.py [--sample-size N] [--output-dir DIR]
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys
import subprocess
import tempfile
import os
from datetime import datetime
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

def find_assembly_file():
    """Find the NCBI assembly summary file."""
    possible_paths = [
        Path("../00assembly_summary_genbank.txt"),
        Path("00assembly_summary_genbank.txt"),
        Path("metadata/00assembly_summary_genbank.txt"),
        Path("../metadata/00assembly_summary_genbank.txt")
    ]
    
    for path in possible_paths:
        if path.exists():
            return path
    
    print("❌ Could not find 00assembly_summary_genbank.txt")
    print("   Searched in:", [str(p) for p in possible_paths])
    sys.exit(1)

def classify_genome_source(organism_name):
    """Classify genome as isolate or uncultured based on organism name."""
    if pd.isna(organism_name):
        return 'unknown'
    
    name_lower = str(organism_name).lower()
    uncultured_indicators = [
        'uncultured', 'unculture', 'environmental', 'metagenome',
        'unidentified', 'unknown', 'candidate', 'endosymbiont'
    ]
    
    for indicator in uncultured_indicators:
        if indicator in name_lower:
            return 'uncultured'
    
    return 'isolate'

def get_lineages_from_taxids(taxids, batch_size=1000):
    """Get taxonomic lineages for species_taxids using taxonkit."""
    if not taxids:
        return {}
    
    print(f"🔍 Getting lineages for {len(taxids)} species taxids...")
    lineage_data = {}
    
    # Process in batches to avoid command line length limits
    for i in tqdm(range(0, len(taxids), batch_size), desc="Processing taxid batches"):
        batch_taxids = taxids[i:i+batch_size]
        
        # Create temporary file for this batch
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_file:
            temp_filename = temp_file.name
            for taxid in batch_taxids:
                temp_file.write(f"{taxid}\n")
        
        try:
            # Use taxonkit lineage with ranks and taxids
            result = subprocess.run([
                "taxonkit", "lineage", "-R", "-t", temp_filename
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    if not line.strip():
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        taxid = parts[0].strip()
                        lineage = parts[1].strip()
                        lineage_taxids = parts[2].strip()
                        lineage_ranks = parts[3].strip()
                        
                        if lineage and lineage != taxid:
                            lineage_data[taxid] = {
                                'lineage': lineage,
                                'lineage_ranks': lineage_ranks,
                                'lineage_taxids': lineage_taxids
                            }
            else:
                print(f"⚠️  Taxonkit warning for batch {i//batch_size + 1}: {result.stderr}")
                
        except Exception as e:
            print(f"❌ Error processing batch {i//batch_size + 1}: {e}")
        finally:
            try:
                os.unlink(temp_filename)
            except:
                pass
    
    print(f"✅ Retrieved lineages for {len(lineage_data)} species")
    return lineage_data

def parse_lineage_to_columns(lineage_info):
    """Parse lineage information into taxonomic rank columns."""
    if not lineage_info:
        return {}
    
    lineage = lineage_info.get('lineage', '')
    ranks = lineage_info.get('lineage_ranks', '')
    taxids = lineage_info.get('lineage_taxids', '')
    
    # Initialize all ranks
    rank_data = {
        'domain': '', 'kingdom': '', 'phylum': '', 'class': '', 
        'order': '', 'family': '', 'genus': '', 'species': ''
    }
    
    if not lineage or not ranks:
        return rank_data
    
    lineage_parts = lineage.split(';')
    rank_parts = ranks.split(';')
    
    # Map ranks to values
    for lineage_part, rank_part in zip(lineage_parts, rank_parts):
        lineage_part = lineage_part.strip()
        rank_part = rank_part.strip().lower()
        
        if rank_part in rank_data and lineage_part:
            rank_data[rank_part] = lineage_part
    
    return rank_data

def create_species_grouped_dataset(df, output_dir):
    """Create the main species-grouped dataset with isolate preference."""
    print("\n🧬 Creating species-grouped dataset...")
    
    # Add genome source classification
    print("📋 Classifying genome sources...")
    df['genome_source'] = df['organism_name'].apply(classify_genome_source)
    
    # Group by species_taxid and get counts
    print("📊 Grouping by species_taxid...")
    species_stats = df.groupby('species_taxid').agg({
        'assembly_accession': 'count',
        'organism_name': 'first',  # Get representative organism name
        'taxid': 'first',  # Get representative taxid
        'genome_source': lambda x: (x == 'isolate').sum()  # Count isolates
    }).rename(columns={
        'assembly_accession': 'total_genomes',
        'organism_name': 'representative_organism_name',
        'taxid': 'representative_taxid',
        'genome_source': 'isolate_count'
    })
    
    species_stats['uncultured_count'] = species_stats['total_genomes'] - species_stats['isolate_count']
    species_stats.reset_index(inplace=True)
    
    print(f"✅ Found {len(species_stats)} unique species")
    print(f"   Total genomes: {species_stats['total_genomes'].sum():,}")
    print(f"   Species with isolates: {(species_stats['isolate_count'] > 0).sum():,}")
    print(f"   Species with only uncultured: {(species_stats['isolate_count'] == 0).sum():,}")
    
    # Get lineage information
    unique_species_taxids = species_stats['species_taxid'].astype(str).tolist()
    lineage_data = get_lineages_from_taxids(unique_species_taxids)
    
    # Add taxonomic columns
    print("🏷️  Adding taxonomic information...")
    for rank in ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        species_stats[rank] = ''
    
    for idx, row in tqdm(species_stats.iterrows(), total=len(species_stats), desc="Processing lineages"):
        species_taxid_str = str(row['species_taxid'])
        if species_taxid_str in lineage_data:
            rank_data = parse_lineage_to_columns(lineage_data[species_taxid_str])
            for rank, value in rank_data.items():
                species_stats.at[idx, rank] = value
    
    # Save species statistics
    output_file = output_dir / "species_grouped_statistics.csv"
    species_stats.to_csv(output_file, index=False)
    print(f"💾 Saved species statistics: {output_file}")
    
    return species_stats

def create_isolate_preferred_subset(df, species_stats, output_dir):
    """Create isolate-preferred subset (one genome per species)."""
    print("\n🎯 Creating isolate-preferred subset...")
    
    # Separate isolate and uncultured genomes
    isolate_df = df[df['genome_source'] == 'isolate'].copy()
    uncultured_df = df[df['genome_source'] == 'uncultured'].copy()
    
    # Get one isolate per species (if available)
    isolate_subset = isolate_df.drop_duplicates(subset=['species_taxid'], keep='first')
    
    # Get species that don't have isolates
    species_with_isolates = set(isolate_subset['species_taxid'])
    uncultured_only = uncultured_df[~uncultured_df['species_taxid'].isin(species_with_isolates)]
    uncultured_subset = uncultured_only.drop_duplicates(subset=['species_taxid'], keep='first')
    
    # Combine isolate-preferred subset
    species_subset = pd.concat([isolate_subset, uncultured_subset], ignore_index=True)
    
    print(f"✅ Isolate-preferred subset: {len(species_subset):,} genomes")
    print(f"   - Isolate genomes: {len(isolate_subset):,}")
    print(f"   - Uncultured genomes: {len(uncultured_subset):,}")
    print(f"   - Species coverage: {species_subset['species_taxid'].nunique():,}")
    
    # Save subset
    output_file = output_dir / "species_isolate_preferred_subset.csv"
    species_subset.to_csv(output_file, index=False)
    print(f"💾 Saved isolate-preferred subset: {output_file}")
    
    return species_subset

def main():
    parser = argparse.ArgumentParser(description="Group NCBI data by species_taxid")
    parser.add_argument("--sample-size", type=int, help="Process only N records for testing")
    parser.add_argument("--output-dir", type=str, default=".", help="Output directory")
    args = parser.parse_args()
    
    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("🧬 Species Taxid Grouper - New Reproducible Method")
    print("=" * 60)
    print(f"📅 Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Find and load assembly file
    assembly_file = find_assembly_file()
    print(f"📂 Loading: {assembly_file}")
    
    # Load data
    try:
        df = pd.read_csv(assembly_file, sep='\t', low_memory=False)
        if args.sample_size:
            df = df.head(args.sample_size)
            print(f"🔬 Using sample size: {len(df):,} records")
        else:
            print(f"📊 Total records: {len(df):,}")
    except Exception as e:
        print(f"❌ Error loading file: {e}")
        sys.exit(1)
    
    # Verify required columns
    required_cols = ['assembly_accession', 'taxid', 'species_taxid', 'organism_name']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"❌ Missing required columns: {missing_cols}")
        sys.exit(1)
    
    # Add genome source classification
    df['genome_source'] = df['organism_name'].apply(classify_genome_source)
    
    # Create outputs
    species_stats = create_species_grouped_dataset(df, output_dir)
    species_subset = create_isolate_preferred_subset(df, species_stats, output_dir)
    
    # Generate summary report
    print("\n📋 Summary Report")
    print("=" * 40)
    print(f"Total input records: {len(df):,}")
    print(f"Unique species (species_taxid): {df['species_taxid'].nunique():,}")
    print(f"Isolate genomes: {(df['genome_source'] == 'isolate').sum():,}")
    print(f"Uncultured genomes: {(df['genome_source'] == 'uncultured').sum():,}")
    print(f"Species with isolates: {(species_stats['isolate_count'] > 0).sum():,}")
    print(f"Species with only uncultured: {(species_stats['isolate_count'] == 0).sum():,}")
    
    print(f"\n✅ Processing completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"📁 Output files saved in: {output_dir}")

if __name__ == "__main__":
    main()
