import subprocess
import pandas as pd
from pathlib import Path
import sys
import logging
from datetime import datetime

def setup_logging(log_dir):
    """Set up logging configuration"""
    log_file = log_dir / f"taxonomic_verification_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return log_file

def verify_taxonomic_names(input_file, expected_domain, expected_rank, taxdump_dir=None):
    """Verify taxonomic names using taxonkit"""
    try:
        # Read the input file
        with open(input_file, 'r') as f:
            names = [line.strip() for line in f if line.strip()]
        
        if not names:
            logging.warning(f"No names found in {input_file}")
            return pd.DataFrame(), pd.DataFrame()
        
        # Prepare taxonkit command
        cmd = ["taxonkit", "name2taxid"]
        if taxdump_dir:
            cmd.extend(["--data-dir", str(taxdump_dir)])
        
        # Run taxonkit to get taxids
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Send names to taxonkit
        stdout, stderr = process.communicate(input='\n'.join(names))
        
        if stderr:
            logging.error(f"Error running taxonkit: {stderr}")
            return pd.DataFrame(), pd.DataFrame()
        
        # Parse taxonkit output
        taxid_data = []
        for line in stdout.strip().split('\n'):
            if line:
                name, taxid = line.split('\t')
                taxid_data.append({'name': name, 'taxid': taxid})
        
        if not taxid_data:
            logging.warning(f"No taxids found for names in {input_file}")
            return pd.DataFrame(), pd.DataFrame()
        
        # Get lineage information
        taxids = [d['taxid'] for d in taxid_data]
        cmd = ["taxonkit", "lineage", "-i", "2"]
        if taxdump_dir:
            cmd.extend(["--data-dir", str(taxdump_dir)])
        
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        stdout, stderr = process.communicate(input='\n'.join(taxids))
        
        if stderr:
            logging.error(f"Error running taxonkit lineage: {stderr}")
            return pd.DataFrame(), pd.DataFrame()
        
        # Parse lineage data
        results = []
        errors = []
        
        for line in stdout.strip().split('\n'):
            if not line:
                continue
                
            taxid, lineage = line.split('\t')
            name = next((d['name'] for d in taxid_data if d['taxid'] == taxid), None)
            
            # Parse lineage string
            ranks = {}
            for item in lineage.split(';'):
                if ':' in item:
                    rank, value = item.split(':')
                    ranks[rank.strip()] = value.strip()
            
            # Check domain and rank
            domain = ranks.get('domain', '').lower()
            rank = ranks.get(expected_rank, '').lower()
            
            is_valid = True
            error_reasons = []
            
            if domain != expected_domain.lower():
                is_valid = False
                error_reasons.append(f"Expected domain {expected_domain}, got {domain}")
            
            if not rank:
                is_valid = False
                error_reasons.append(f"No {expected_rank} found in lineage")
            
            result = {
                'name': name,
                'taxid': taxid,
                'domain': domain,
                'rank': rank,
                'is_valid': is_valid,
                'error_reasons': '; '.join(error_reasons) if error_reasons else ''
            }
            
            if is_valid:
                results.append(result)
            else:
                errors.append(result)
        
        return pd.DataFrame(results), pd.DataFrame(errors)
    
    except Exception as e:
        logging.error(f"Error verifying taxonomic names: {e}")
        return pd.DataFrame(), pd.DataFrame()

def main():
    # Set up paths
    script_dir = Path(__file__).resolve().parent
    csv_dir = script_dir / "csv"
    output_dir = script_dir / "error_log"
    taxdump_dir = script_dir / "taxdump_ncbi"
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(exist_ok=True)
    
    # Set up logging
    log_file = setup_logging(output_dir)
    logging.info("Starting taxonomic name verification")
    
    # Define files to verify
    files_to_verify = [
        ("archaeal_families.txt", "archaea", "family"),
        ("bacterial_families.txt", "bacteria", "family"),
        ("archaeal_genera.txt", "archaea", "genus"),
        ("bacterial_genera.txt", "bacteria", "genus")
    ]
    
    # Process each file
    for filename, domain, rank in files_to_verify:
        input_file = csv_dir / filename
        if not input_file.exists():
            logging.warning(f"File not found: {input_file}")
            continue
        
        logging.info(f"Verifying {filename}...")
        valid_df, error_df = verify_taxonomic_names(input_file, domain, rank, taxdump_dir)
        
        # Save results
        if not valid_df.empty:
            valid_file = output_dir / f"valid_{filename.replace('.txt', '.csv')}"
            valid_df.to_csv(valid_file, index=False)
            logging.info(f"Valid names saved to: {valid_file}")
        
        if not error_df.empty:
            error_file = output_dir / f"invalid_{filename.replace('.txt', '.csv')}"
            error_df.to_csv(error_file, index=False)
            logging.info(f"Invalid names saved to: {error_file}")
        
        # Log statistics
        total = len(valid_df) + len(error_df)
        if total > 0:
            valid_count = len(valid_df)
            error_count = len(error_df)
            logging.info(f"Total names processed: {total}")
            logging.info(f"Valid names: {valid_count}")
            logging.info(f"Invalid names: {error_count}")
            logging.info(f"Validation success rate: {(valid_count/total)*100:.2f}%")
    
    logging.info(f"Verification complete. Log file: {log_file}")

if __name__ == "__main__":
    main() 