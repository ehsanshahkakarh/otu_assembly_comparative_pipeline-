#!/usr/bin/env python3
"""
Genus Workflow Runner

This script runs the genus parser and accessions scripts in the correct order.
The accessions script must run first to generate the classified file that the parser needs.

Usage:
    python run_genus_workflow.py [--output-dir OUTPUT_DIR]
"""

import subprocess
import argparse
import sys
from pathlib import Path

def run_script(script_name: str, description: str, args: list = None) -> bool:
    """Run a script and return success status."""
    print(f"\n🚀 {description}")
    print("=" * 60)
    
    cmd = ["python", script_name]
    if args:
        cmd.extend(args)
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print(f"✅ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with return code {e.returncode}")
        return False
    except FileNotFoundError:
        print(f"❌ Script not found: {script_name}")
        return False

def check_prerequisites():
    """Check if required scripts exist."""
    required_scripts = [
        "genus_ncbi_accessions_new.py",
        "genus_ncbi_parser_new.py"
    ]
    
    missing_scripts = []
    for script in required_scripts:
        if not Path(script).exists():
            missing_scripts.append(script)
    
    if missing_scripts:
        print(f"❌ Missing required scripts:")
        for script in missing_scripts:
            print(f"   - {script}")
        return False
    
    print(f"✅ All required scripts found")
    return True

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Run genus parser workflow in correct order")
    parser.add_argument("--output-dir", help="Output directory for both scripts")
    parser.add_argument("--input-dir", help="Input directory for both scripts")
    
    args = parser.parse_args()
    
    print("🧬 GENUS NCBI WORKFLOW")
    print("=" * 60)
    print("This workflow runs both genus scripts:")
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
        sys.exit(1)
    
    # Prepare arguments for scripts
    script_args = []
    if args.output_dir:
        script_args.extend(["--output-dir", args.output_dir])
    if args.input_dir:
        script_args.extend(["--input-dir", args.input_dir])
    
    # Step 1: Run accessions script first
    accessions_success = run_script(
        "genus_ncbi_accessions_new.py",
        "Step 1: Building genus accessions with classification",
        script_args
    )
    
    if not accessions_success:
        print(f"\n❌ Accessions script failed - cannot proceed with parser")
        sys.exit(1)
    
    # Step 2: Run parser script (depends on classified file from step 1)
    parser_success = run_script(
        "genus_ncbi_parser_new.py", 
        "Step 2: Generating genus counts using classified data",
        script_args
    )
    
    # Summary
    print(f"\n🎉 WORKFLOW SUMMARY")
    print("=" * 60)
    print(f"   Accessions: {'✅ Success' if accessions_success else '❌ Failed'}")
    print(f"   Parser: {'✅ Success' if parser_success else '❌ Failed'}")
    
    if accessions_success and parser_success:
        print(f"\n✅ Genus workflow completed successfully!")
        
        # Show output files
        output_dir = Path(args.output_dir) if args.output_dir else Path(__file__).resolve().parent.parent
        expected_files = [
            "ncbi_genus_counts.csv",
            "ncbi_genus_species_counts.csv", 
            "ncbi_genus_with_accessions_classified.csv",
            "ncbi_genus_species_subset_with_accessions_classified.csv"
        ]
        
        print(f"\n📁 Expected output files in {output_dir}:")
        for filename in expected_files:
            file_path = output_dir / filename
            if file_path.exists():
                size_info = f"({file_path.stat().st_size / (1024*1024):.1f} MB)" if file_path.stat().st_size > 1024*1024 else f"({file_path.stat().st_size / 1024:.1f} KB)"
                print(f"   ✅ {filename} {size_info}")
            else:
                print(f"   ❌ {filename} - not found")
        
        sys.exit(0)
    else:
        print(f"\n❌ Workflow failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()
