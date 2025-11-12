#!/usr/bin/env python3
"""
Family NCBI Workflow Runner

This script runs both family NCBI scripts in sequence:
1. Family accessions builder (generates classified accession files)
2. Family parser (generates count files)

Both scripts use species_taxid for true species-level subsetting.

Usage:
    python run_family_workflow.py [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
"""

import subprocess
import sys
import argparse
from pathlib import Path
from datetime import datetime

def check_prerequisites():
    """Check if required scripts exist."""
    required_scripts = [
        "family_ncbi_accessions_new.py",
        "family_ncbi_parser_new.py"
    ]
    
    missing_scripts = []
    for script in required_scripts:
        if not Path(script).exists():
            missing_scripts.append(script)
    
    if missing_scripts:
        print(f"❌ Missing required scripts: {', '.join(missing_scripts)}")
        return False
    
    print(f"✅ All required scripts found")
    return True

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Run family parser workflow in correct order")
    parser.add_argument("--output-dir", help="Output directory for both scripts")
    parser.add_argument("--input-dir", help="Input directory for both scripts")
    
    args = parser.parse_args()
    
    print("🧬 FAMILY NCBI WORKFLOW")
    print("=" * 60)
    print("This workflow runs both family scripts:")
    print("1. Accessions script (generates classified accession files)")
    print("2. Parser script (generates count files)")
    print("\n🧬 Species Subsetting Approach:")
    print("Both scripts now use species_taxid for true species-level subsetting")
    print("This is the most biologically accurate approach that:")
    print("- Uses NCBI taxonomy to identify true species")
    print("- Groups all strains of the same species together")
    print("- Avoids issues with strain naming conventions")
    
    # Check prerequisites
    if not check_prerequisites():
        return 1
    
    # Prepare common arguments
    common_args = []
    if args.input_dir:
        common_args.extend(["--input-dir", args.input_dir])
    if args.output_dir:
        common_args.extend(["--output-dir", args.output_dir])
    
    start_time = datetime.now()
    
    # Step 1: Run family accessions script
    print(f"\n🚀 Step 1: Running family accessions builder...")
    print(f"⏰ Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    accessions_cmd = ["python", "family_ncbi_accessions_new.py"] + common_args
    print(f"📝 Command: {' '.join(accessions_cmd)}")
    
    try:
        result = subprocess.run(accessions_cmd, check=True)
        print(f"✅ Family accessions builder completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"❌ Family accessions builder failed with return code {e.returncode}")
        return 1
    except Exception as e:
        print(f"❌ Error running family accessions builder: {e}")
        return 1
    
    # Step 2: Run family parser script
    print(f"\n🚀 Step 2: Running family parser...")
    print(f"⏰ Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    parser_cmd = ["python", "family_ncbi_parser_new.py"] + common_args
    print(f"📝 Command: {' '.join(parser_cmd)}")
    
    try:
        result = subprocess.run(parser_cmd, check=True)
        print(f"✅ Family parser completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"❌ Family parser failed with return code {e.returncode}")
        return 1
    except Exception as e:
        print(f"❌ Error running family parser: {e}")
        return 1
    
    # Summary
    end_time = datetime.now()
    duration = end_time - start_time
    
    print(f"\n🎉 FAMILY WORKFLOW COMPLETED SUCCESSFULLY!")
    print("=" * 60)
    print(f"⏰ Total duration: {duration}")
    print(f"📂 Output files should be in the specified output directory")
    print(f"\n📊 Generated files:")
    print(f"   - ncbi_family_with_accessions_classified.csv (total genome accessions)")
    print(f"   - ncbi_family_species_subset_with_accessions_classified.csv (species subset accessions)")
    print(f"   - ncbi_family_counts.csv (total genome counts)")
    print(f"   - ncbi_family_species_counts.csv (species subset counts)")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
