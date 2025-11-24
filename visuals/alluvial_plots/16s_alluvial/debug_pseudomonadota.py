#!/usr/bin/env python3

import pandas as pd
import sys

print("=== Debugging Pseudomonadota in 16S Bacteria Data ===\n")

# Load the merged data
try:
    data_path = "../../../Eukcensus_merge/16s_merged/csv_results/16s_ncbi_merged_clean_phylum.csv"
    data = pd.read_csv(data_path)
    print(f"Loaded data: {len(data)} rows")
except Exception as e:
    print(f"Error loading data: {e}")
    sys.exit(1)

# Filter for bacteria
bacteria = data[data['domain'] == 'Bacteria']
print(f"Bacteria entries: {len(bacteria)} rows")

# Check if Pseudomonadota exists
pseudo = bacteria[bacteria['phylum'] == 'Pseudomonadota']
print(f"Pseudomonadota entries: {len(pseudo)} rows")

if len(pseudo) == 0:
    print("❌ ERROR: No Pseudomonadota found in bacteria data!")
    print("\nAvailable bacterial phyla:")
    phyla = bacteria['phylum'].value_counts().head(15)
    for phylum, count in phyla.items():
        print(f"  {phylum}: {count} entries")
    sys.exit(1)

# Check Pseudomonadota data
print(f"\n=== Pseudomonadota Data ===")
print(f"Entries: {len(pseudo)}")
print(f"Match status: {pseudo['match_status'].value_counts().to_dict()}")

# Calculate totals
bacteria_totals = {
    'genomes': bacteria['ncbi_genome_count'].sum(),
    'species': bacteria['ncbi_species_count'].sum(),
    'otus': bacteria['census_otu_count'].sum(),
    'sequences': bacteria['census_size_count'].sum()
}

pseudo_totals = {
    'genomes': pseudo['ncbi_genome_count'].sum(),
    'species': pseudo['ncbi_species_count'].sum(),
    'otus': pseudo['census_otu_count'].sum(),
    'sequences': pseudo['census_size_count'].sum()
}

print(f"\n=== Bacteria Totals ===")
for key, value in bacteria_totals.items():
    print(f"{key.capitalize()}: {value:,}")

print(f"\n=== Pseudomonadota Totals ===")
for key, value in pseudo_totals.items():
    print(f"{key.capitalize()}: {value:,}")

print(f"\n=== Pseudomonadota Percentages ===")
for key in bacteria_totals:
    if bacteria_totals[key] > 0:
        pct = (pseudo_totals[key] / bacteria_totals[key]) * 100
        print(f"{key.capitalize()}: {pct:.2f}%")

# Calculate total representation
pseudo_total_rep = sum([
    (pseudo_totals['genomes'] / bacteria_totals['genomes']) * 100,
    (pseudo_totals['species'] / bacteria_totals['species']) * 100,
    (pseudo_totals['otus'] / bacteria_totals['otus']) * 100,
    (pseudo_totals['sequences'] / bacteria_totals['sequences']) * 100
])

print(f"\nPseudomonadota Total Representation: {pseudo_total_rep:.2f}%")

# Check top bacteria by total representation
print(f"\n=== Top 15 Bacteria by Total Representation ===")
bacteria_summary = bacteria.groupby('phylum').agg({
    'ncbi_genome_count': 'sum',
    'ncbi_species_count': 'sum', 
    'census_otu_count': 'sum',
    'census_size_count': 'sum'
}).reset_index()

bacteria_summary['genome_pct'] = (bacteria_summary['ncbi_genome_count'] / bacteria_totals['genomes']) * 100
bacteria_summary['species_pct'] = (bacteria_summary['ncbi_species_count'] / bacteria_totals['species']) * 100
bacteria_summary['otu_pct'] = (bacteria_summary['census_otu_count'] / bacteria_totals['otus']) * 100
bacteria_summary['seq_pct'] = (bacteria_summary['census_size_count'] / bacteria_totals['sequences']) * 100
bacteria_summary['total_rep'] = bacteria_summary['genome_pct'] + bacteria_summary['species_pct'] + bacteria_summary['otu_pct'] + bacteria_summary['seq_pct']

top_bacteria = bacteria_summary.sort_values('total_rep', ascending=False).head(15)

for i, row in top_bacteria.iterrows():
    print(f"{row.name+1:2d}. {row['phylum']:<20} | G: {row['genome_pct']:5.1f}% | S: {row['species_pct']:5.1f}% | O: {row['otu_pct']:5.1f}% | Seq: {row['seq_pct']:5.1f}% | Total: {row['total_rep']:6.1f}%")

print("\n=== Analysis Complete ===")
