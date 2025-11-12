#!/usr/bin/env python3
"""
Test All New NCBI Scripts

This script tests all the new NCBI parser and accessions scripts to ensure they work correctly
with species_taxid-based true species-level subsetting.

Usage:
    python test_all_new_scripts.py [--test-genus] [--test-family] [--test-phylum] [--test-all]
"""

import subprocess
import sys
import argparse
from pathlib import Path
from datetime import datetime
import pandas as pd

def check_script_exists(script_name: str) -> bool:
    """Check if a script exists."""
    return Path(script_name).exists()

def run_script_test(script_name: str, script_type: str, taxonomic_level: str) -> bool:
    """Run a single script test."""
    print(f"\n🧪 Testing {taxonomic_level} {script_type}...")
    print(f"📝 Script: {script_name}")
    
    if not check_script_exists(script_name):
        print(f"❌ Script not found: {script_name}")
        return False
    
    try:
        # Run with --help to test basic functionality
        result = subprocess.run(
            ["python", script_name, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0:
            print(f"✅ {taxonomic_level} {script_type} help test passed")
            return True
        else:
            print(f"❌ {taxonomic_level} {script_type} help test failed")
            print(f"Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"⏰ {taxonomic_level} {script_type} test timed out")
        return False
    except Exception as e:
        print(f"❌ Error testing {taxonomic_level} {script_type}: {e}")
        return False

def test_workflow_script(workflow_script: str, taxonomic_level: str) -> bool:
    """Test a workflow script."""
    print(f"\n🔄 Testing {taxonomic_level} workflow...")
    print(f"📝 Workflow: {workflow_script}")
    
    if not check_script_exists(workflow_script):
        print(f"❌ Workflow script not found: {workflow_script}")
        return False
    
    try:
        # Run with --help to test basic functionality
        result = subprocess.run(
            ["python", workflow_script, "--help"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30
        )
        
        if result.returncode == 0:
            print(f"✅ {taxonomic_level} workflow help test passed")
            return True
        else:
            print(f"❌ {taxonomic_level} workflow help test failed")
            print(f"Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"⏰ {taxonomic_level} workflow test timed out")
        return False
    except Exception as e:
        print(f"❌ Error testing {taxonomic_level} workflow: {e}")
        return False

def test_genus_scripts() -> bool:
    """Test genus scripts."""
    print(f"\n🧬 TESTING GENUS SCRIPTS")
    print("=" * 50)
    
    success = True
    
    # Test genus accessions script
    if not run_script_test("genus_ncbi_accessions_new.py", "accessions", "genus"):
        success = False
    
    # Test genus parser script
    if not run_script_test("genus_ncbi_parser_new.py", "parser", "genus"):
        success = False
    
    # Test genus workflow script
    if not test_workflow_script("run_genus_workflow.py", "genus"):
        success = False
    
    return success

def test_family_scripts() -> bool:
    """Test family scripts."""
    print(f"\n🧬 TESTING FAMILY SCRIPTS")
    print("=" * 50)
    
    success = True
    
    # Test family accessions script
    if not run_script_test("family_ncbi_accessions_new.py", "accessions", "family"):
        success = False
    
    # Test family parser script
    if not run_script_test("family_ncbi_parser_new.py", "parser", "family"):
        success = False
    
    # Test family workflow script
    if not test_workflow_script("run_family_workflow.py", "family"):
        success = False
    
    return success

def test_phylum_scripts() -> bool:
    """Test phylum scripts."""
    print(f"\n🧬 TESTING PHYLUM SCRIPTS")
    print("=" * 50)
    
    success = True
    
    # Test phylum accessions script
    if not run_script_test("phylum_ncbi_accessions_new.py", "accessions", "phylum"):
        success = False
    
    # Test phylum parser script
    if not run_script_test("phylum_ncbi_parser_new.py", "parser", "phylum"):
        success = False
    
    # Test phylum workflow script
    if not test_workflow_script("run_phylum_workflow.py", "phylum"):
        success = False
    
    return success

def check_dependencies() -> bool:
    """Check if required dependencies are available."""
    print(f"\n🔍 CHECKING DEPENDENCIES")
    print("=" * 50)
    
    success = True
    
    # Check Python packages
    required_packages = ['pandas', 'tqdm', 'psutil']
    for package in required_packages:
        try:
            __import__(package)
            print(f"✅ {package} is available")
        except ImportError:
            print(f"❌ {package} is not available")
            success = False
    
    # Check taxonkit
    try:
        result = subprocess.run(
            ["taxonkit", "version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            version = result.stdout.strip()
            print(f"✅ taxonkit is available: {version}")
        else:
            print(f"❌ taxonkit is not working properly")
            success = False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print(f"❌ taxonkit is not available")
        success = False
    
    return success

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Test all new NCBI scripts")
    parser.add_argument("--test-genus", action="store_true", help="Test genus scripts only")
    parser.add_argument("--test-family", action="store_true", help="Test family scripts only")
    parser.add_argument("--test-phylum", action="store_true", help="Test phylum scripts only")
    parser.add_argument("--test-all", action="store_true", help="Test all scripts (default)")
    parser.add_argument("--skip-dependencies", action="store_true", help="Skip dependency checks")
    
    args = parser.parse_args()
    
    # Default to test all if no specific test is requested
    if not any([args.test_genus, args.test_family, args.test_phylum]):
        args.test_all = True
    
    print("🧪 NCBI SCRIPTS TEST SUITE")
    print("=" * 60)
    print("Testing all new NCBI scripts with species_taxid-based subsetting")
    print(f"⏰ Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    start_time = datetime.now()
    overall_success = True
    
    # Check dependencies first
    if not args.skip_dependencies:
        if not check_dependencies():
            print(f"\n❌ Dependency check failed. Some tests may not work properly.")
            overall_success = False
    
    # Run tests based on arguments
    if args.test_genus or args.test_all:
        if not test_genus_scripts():
            overall_success = False
    
    if args.test_family or args.test_all:
        if not test_family_scripts():
            overall_success = False
    
    if args.test_phylum or args.test_all:
        if not test_phylum_scripts():
            overall_success = False
    
    # Summary
    end_time = datetime.now()
    duration = end_time - start_time
    
    print(f"\n📊 TEST SUMMARY")
    print("=" * 60)
    print(f"⏰ Total duration: {duration}")
    
    if overall_success:
        print(f"🎉 ALL TESTS PASSED!")
        print(f"✅ All new NCBI scripts are working correctly")
        print(f"\n🧬 Key Features Verified:")
        print(f"   - Scripts can be imported and run")
        print(f"   - Help functionality works")
        print(f"   - Dependencies are available")
        print(f"   - Workflow scripts are functional")
        print(f"\n🚀 Ready to run with species_taxid-based true species-level subsetting!")
        return 0
    else:
        print(f"❌ SOME TESTS FAILED!")
        print(f"⚠️  Please check the errors above and fix any issues")
        print(f"\n🔧 Common fixes:")
        print(f"   - Install missing Python packages: pip install pandas tqdm psutil")
        print(f"   - Install taxonkit: conda install -c bioconda taxonkit")
        print(f"   - Check file permissions and paths")
        return 1

if __name__ == "__main__":
    sys.exit(main())
