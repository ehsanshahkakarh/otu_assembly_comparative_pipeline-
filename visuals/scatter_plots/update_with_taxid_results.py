#!/usr/bin/env python3
"""
Update Overrepresented Files with TaxID-Based Results
====================================================
Replace the existing overrepresented files with TaxID-based extraction results.
"""

import subprocess
import os
from datetime import datetime

def run_taxid_extraction(gene, category):
    """Run TaxID-based extraction and save to family_diversity_results."""
    
    output_dir = "family_diversity_results"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create filename matching the existing format
    filename = f"{gene}_{category}_family_diversity.txt"
    filepath = os.path.join(output_dir, filename)
    
    # Run the TaxID-based extraction command
    cmd = ["python", "extract_ncbi_taxid_diversity.py", "--gene", gene, "--category", category, "--output", filepath]
    
    print(f"🔍 Running TaxID-based extraction: {' '.join(cmd)}")
    print(f"📁 Overwriting: {filepath}")
    
    try:
        # Run command and capture output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Write output to file
        with open(filepath, 'w') as f:
            f.write(f"# NCBI TaxID-Based Overrepresented Taxonomic Diversity Analysis\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Command: {' '.join(cmd)}\n")
            f.write("# " + "="*70 + "\n\n")
            f.write(result.stdout)
            
            if result.stderr:
                f.write("\n\n# STDERR:\n")
                f.write(result.stderr)
        
        print(f"✅ Success: Updated {filepath}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running command: {e}")
        print(f"   stdout: {e.stdout}")
        print(f"   stderr: {e.stderr}")
        return False

def main():
    """Update overrepresented family diversity files with TaxID-based results."""
    print("🚀 UPDATING OVERREPRESENTED FILES WITH TAXID-BASED RESULTS")
    print("="*70)
    print("📁 Target directory: family_diversity_results/")
    
    # Define overrepresented extractions to update
    extractions = [
        ("16S", "overrepresented"),
        ("18S", "overrepresented"),
    ]
    
    successful = 0
    failed = 0
    
    for gene, category in extractions:
        print(f"\n{'='*60}")
        success = run_taxid_extraction(gene, category)
        if success:
            successful += 1
        else:
            failed += 1
    
    print(f"\n🎯 FINAL SUMMARY:")
    print(f"✅ Successfully updated: {successful}")
    print(f"❌ Failed updates: {failed}")
    
    # List all files in family_diversity_results
    output_dir = "family_diversity_results"
    if os.path.exists(output_dir):
        files = [f for f in os.listdir(output_dir) if f.endswith('.txt')]
        files.sort()
        
        print(f"\n📋 Files in {output_dir}/:")
        for f in files:
            filepath = os.path.join(output_dir, f)
            size = os.path.getsize(filepath)
            mtime = os.path.getmtime(filepath)
            mtime_str = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
            
            # Mark updated files
            if f.endswith('overrepresented_family_diversity.txt'):
                status = "🔄 TAXID-BASED"
            else:
                status = "📄 Census-based"
            
            print(f"   {status} {f} ({size:,} bytes, {mtime_str})")

if __name__ == "__main__":
    main()
