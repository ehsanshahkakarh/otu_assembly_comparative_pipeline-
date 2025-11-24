#!/usr/bin/env python3

import pandas as pd
import re

print("=== Diagnosing Pseudomonadota Visibility Issue ===\n")

# Check if the flow annotations file has the right structure
try:
    flow_file = "alluvial_16s_bacteria_pct_flow_annotations.tsv"
    flows = pd.read_csv(flow_file, sep='\t')
    print(f"Flow annotations loaded: {len(flows)} rows")
    
    # Check Pseudomonadota entries
    pseudo_flows = flows[flows['Taxon'] == 'Pseudomonadota']
    print(f"Pseudomonadota flow entries: {len(pseudo_flows)}")
    
    if len(pseudo_flows) > 0:
        print("\nPseudomonadota flows:")
        for _, row in pseudo_flows.iterrows():
            print(f"  {row['Node']}: {row['Percentage']:.2f}% (width: {row['Flow_Width']:.2f})")
        
        # Check if nodes 1 and 4 have substantial values
        node1 = pseudo_flows[pseudo_flows['Node'] == 'Genbank_Genome_%']
        node4 = pseudo_flows[pseudo_flows['Node'] == 'Genbank_Species_%']
        
        if len(node1) > 0 and len(node4) > 0:
            node1_pct = node1.iloc[0]['Percentage']
            node4_pct = node4.iloc[0]['Percentage']
            print(f"\nNode 1 (Genbank Genomes): {node1_pct:.2f}%")
            print(f"Node 4 (Genbank Species): {node4_pct:.2f}%")
            
            if node1_pct > 10 and node4_pct > 10:
                print("✅ Both nodes have substantial values - should be visible!")
            else:
                print("❌ One or both nodes have low values")
        else:
            print("❌ Missing node 1 or node 4 data")
    else:
        print("❌ No Pseudomonadota found in flow annotations!")
        
except Exception as e:
    print(f"Error reading flow annotations: {e}")

# Check the summary file
try:
    summary_file = "alluvial_16s_bacteria_pct_summary.tsv"
    summary = pd.read_csv(summary_file, sep='\t')
    print(f"\nSummary file loaded: {len(summary)} rows")
    
    pseudo_summary = summary[summary['Phylum'] == 'Pseudomonadota']
    if len(pseudo_summary) > 0:
        print("Pseudomonadota in summary:")
        for col in summary.columns:
            if col != 'Phylum':
                print(f"  {col}: {pseudo_summary.iloc[0][col]}")
    else:
        print("❌ Pseudomonadota not in summary file!")
        print("Available phyla in summary:")
        for phylum in summary['Phylum'].head(10):
            print(f"  {phylum}")
            
except Exception as e:
    print(f"Error reading summary: {e}")

# Check if there are any numbered prefixes in the data
try:
    # Look for patterns like "1. Pseudomonadota" in flow annotations
    numbered_pseudo = flows[flows['Taxon'].str.contains(r'\d+\.\s*Pseudomonadota', na=False)]
    if len(numbered_pseudo) > 0:
        print(f"\n⚠️  Found numbered Pseudomonadota entries: {len(numbered_pseudo)}")
        for taxon in numbered_pseudo['Taxon'].unique():
            print(f"  {taxon}")
    else:
        print("\n✅ No numbered prefixes found in flow annotations")
        
except Exception as e:
    print(f"Error checking numbered prefixes: {e}")

print("\n=== Diagnosis Complete ===")
