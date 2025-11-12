import subprocess
import pandas as pd
from pathlib import Path
import sys
from datetime import datetime
import logging

def setup_logging(log_dir):
    """Set up logging configuration"""
    log_file = log_dir / f"verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return log_file

def verify_mapping(mapping_file, rank, taxdump_dir=None):
    """
    Verify taxonomic mappings using taxonkit
    
    Args:
        mapping_file (Path): Path to the mapping CSV file
        rank (str): Taxonomic rank to verify (genus/family)
        taxdump_dir (Path): Optional path to taxdump directory
    """
    logging.info(f"Verifying {rank} mappings from {mapping_file}")
    
    # Read mapping file
    try:
        df = pd.read_csv(mapping_file)
        if 'taxid' not in df.columns or rank not in df.columns:
            raise ValueError(f"Input file must contain 'taxid' and '{rank}' columns")
    except Exception as e:
        logging.error(f"Error reading mapping file: {e}")
        return None

    # Prepare taxonkit command
    cmd_lineage = ["taxonkit", "lineage"]
    if taxdump_dir:
        cmd_lineage += ["--data-dir", str(taxdump_dir)]

    # Get taxids
    taxid_list = df['taxid'].astype(str).tolist()
    logging.info(f"Found {len(taxid_list)} taxids to verify")

    # Call taxonkit lineage
    try:
        result = subprocess.run(cmd_lineage, input="\n".join(taxid_list),
                              text=True, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error calling TaxonKit: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Unexpected error running TaxonKit: {str(e)}")
        return None

    # Process results
    discrepancies = []
    lineage_lines = result.stdout.strip().split("\n")
    
    for line in lineage_lines:
        try:
            taxid, lineage = line.split("\t")
            taxid = int(taxid)
            
            # Get expected mapping
            expected = df[df['taxid'] == taxid][rank].iloc[0]
            
            # Parse lineage
            ranks = lineage.split(";")
            found = None
            
            # Look for rank in lineage
            rank_patterns = {
                'genus': ['genus', 'g__'],
                'family': ['family', 'f__']
            }
            
            for r in ranks:
                for pattern in rank_patterns[rank]:
                    if r.lower().startswith(f"{pattern}:"):
                        found = r.split(":", 1)[1].strip()
                        break
                if found:
                    break
            
            # Check for discrepancy
            if found and found != expected:
                discrepancies.append({
                    'taxid': taxid,
                    'expected': expected,
                    'found': found,
                    'lineage': lineage
                })
                logging.warning(f"Discrepancy found for taxid {taxid}:")
                logging.warning(f"  Expected {rank}: {expected}")
                logging.warning(f"  Found {rank}: {found}")
                logging.warning(f"  Full lineage: {lineage}")
                
        except Exception as e:
            logging.error(f"Error processing taxid {taxid}: {str(e)}")
            continue

    # Save discrepancies to file
    if discrepancies:
        output_file = log_dir / f"{rank}_discrepancies.csv"
        pd.DataFrame(discrepancies).to_csv(output_file, index=False)
        logging.info(f"Found {len(discrepancies)} discrepancies. Saved to {output_file}")
    else:
        logging.info(f"No discrepancies found for {rank} mappings")

    return discrepancies

def main():
    # Set up paths
    script_dir = Path(__file__).resolve().parent
    log_dir = script_dir / "error_log"
    taxdump_dir = script_dir / "taxdump_ncbi"
    
    # Create log directory if it doesn't exist
    log_dir.mkdir(exist_ok=True)
    
    # Set up logging
    log_file = setup_logging(log_dir)
    logging.info("Starting taxonomic mapping verification")
    
    # Verify family mappings
    family_file = script_dir / "taxid_to_family.csv"
    if family_file.exists():
        family_discrepancies = verify_mapping(family_file, "family", taxdump_dir)
    else:
        logging.error(f"Family mapping file not found: {family_file}")
    
    # Verify genus mappings
    genus_file = script_dir / "taxid_to_genus.csv"
    if genus_file.exists():
        genus_discrepancies = verify_mapping(genus_file, "genus", taxdump_dir)
    else:
        logging.error(f"Genus mapping file not found: {genus_file}")
    
    logging.info("Verification complete")
    logging.info(f"Log file: {log_file}")

if __name__ == "__main__":
    main() 