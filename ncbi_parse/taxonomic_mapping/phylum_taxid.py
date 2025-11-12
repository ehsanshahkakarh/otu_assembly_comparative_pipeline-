import pandas as pd
import numpy as np
import sys
from datetime import datetime
import logging
from pathlib import Path
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import lru_cache
import time

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler()
        ]
    )

def load_taxonomy_files():
    """Load taxonomy files from taxdump directory"""
    try:
        # Get script directory and set up paths
        script_dir = Path(__file__).resolve().parent
        taxdump_dir = script_dir / "taxdump_ncbi"

        # Check if taxdump directory exists
        if not taxdump_dir.exists():
            raise FileNotFoundError(f"Taxdump directory not found: {taxdump_dir}")

        nodes_file = taxdump_dir / "nodes.dmp"
        names_file = taxdump_dir / "names.dmp"

        # Check if required files exist
        if not nodes_file.exists():
            raise FileNotFoundError(f"Nodes file not found: {nodes_file}")
        if not names_file.exists():
            raise FileNotFoundError(f"Names file not found: {names_file}")

        logging.info(f"Loading taxonomy files from: {taxdump_dir}")

        # Load nodes file with progress indication
        with tqdm(desc="Loading nodes.dmp", unit="MB") as pbar:
            nodes_df = pd.read_csv(nodes_file, sep="|", header=None)
            pbar.update(1)

        # Load names file with progress indication
        with tqdm(desc="Loading names.dmp", unit="MB") as pbar:
            names_df = pd.read_csv(names_file, sep="|", header=None)
            pbar.update(1)

        # Clean up the data
        nodes_df = nodes_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
        names_df = names_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

        return nodes_df, names_df
    except Exception as e:
        logging.error(f"Error loading taxonomy files: {e}")
        sys.exit(1)

def build_taxonomy_maps(nodes_df, names_df):
    """Build taxonomy maps for efficient lookup"""
    try:
        logging.info(f"📊 Building taxonomy maps from {len(nodes_df):,} nodes and {len(names_df):,} names")

        # Create taxid to rank map
        with tqdm(desc="Building taxid→rank map", unit="entries") as pbar:
            taxid_to_rank = dict(zip(nodes_df[0], nodes_df[2]))
            pbar.update(len(nodes_df))

        # Create taxid to parent map
        with tqdm(desc="Building taxid→parent map", unit="entries") as pbar:
            taxid_to_parent = dict(zip(nodes_df[0], nodes_df[1]))
            pbar.update(len(nodes_df))

        # Create taxid to name map (scientific names only)
        with tqdm(desc="Building taxid→name map", unit="entries") as pbar:
            scientific_names = names_df[names_df[3] == "scientific name"]
            taxid_to_name = dict(zip(scientific_names[0], scientific_names[1]))
            pbar.update(len(scientific_names))

        logging.info(f"✅ Built maps: {len(taxid_to_rank):,} ranks, {len(taxid_to_name):,} names")
        return taxid_to_rank, taxid_to_parent, taxid_to_name
    except Exception as e:
        logging.error(f"Error building taxonomy maps: {e}")
        sys.exit(1)

# Global variables for multiprocessing
taxid_to_rank_global = None
taxid_to_parent_global = None
taxid_to_name_global = None

def init_worker(taxid_to_rank, taxid_to_parent, taxid_to_name):
    """Initialize worker process with shared data"""
    global taxid_to_rank_global, taxid_to_parent_global, taxid_to_name_global
    taxid_to_rank_global = taxid_to_rank
    taxid_to_parent_global = taxid_to_parent
    taxid_to_name_global = taxid_to_name

@lru_cache(maxsize=100000)
def get_full_lineage_cached(taxid):
    """Get full lineage for a taxid with caching - returns dict of rank:name"""
    try:
        lineage = {}
        current_taxid = taxid
        visited = set()  # Prevent infinite loops

        while current_taxid != 1 and current_taxid not in visited:  # 1 is the root node
            visited.add(current_taxid)
            rank = taxid_to_rank_global.get(current_taxid)
            name = taxid_to_name_global.get(current_taxid)

            if rank and name:
                lineage[rank] = name

            current_taxid = taxid_to_parent_global.get(current_taxid)
            if current_taxid is None:
                break

        return lineage
    except Exception as e:
        logging.error(f"Error getting lineage for taxid {taxid}: {e}")
        return {}

def get_rank_target(taxid, target_rank, taxid_to_rank, taxid_to_parent, taxid_to_name):
    """Get the target rank for a given taxid (legacy function for compatibility)"""
    try:
        current_taxid = taxid
        while current_taxid != 1:  # 1 is the root node
            if taxid_to_rank.get(current_taxid) == target_rank:
                return taxid_to_name.get(current_taxid)
            current_taxid = taxid_to_parent.get(current_taxid)
        return None
    except Exception as e:
        logging.error(f"Error getting {target_rank} for taxid {taxid}: {e}")
        return None

def detect_domain_from_lineage(lineage):
    """Detect domain from lineage, with special handling for viruses"""
    # First try standard domain/superkingdom
    domain = lineage.get("domain", lineage.get("superkingdom"))
    if domain:
        return domain

    # Check for viral indicators in the lineage
    viral_indicators = [
        'viruses', 'virus', 'viral', 'viridae', 'viricota', 'viricetes',
        'orthornavirae', 'bamfordvirae', 'varidnaviria', 'riboviria',
        'duplodnaviria', 'monodnaviria', 'adnaviria', 'ribozyviria'
    ]

    # Check all ranks in lineage for viral indicators
    for rank, name in lineage.items():
        if name and any(indicator in name.lower() for indicator in viral_indicators):
            return "Viruses"

    # Check for eukaryotic indicators
    eukaryotic_indicators = [
        'eukaryota', 'eukarya', 'fungi', 'metazoa', 'viridiplantae',
        'stramenopiles', 'alveolata', 'rhizaria', 'excavata', 'amoebozoa'
    ]

    for rank, name in lineage.items():
        if name and any(indicator in name.lower() for indicator in eukaryotic_indicators):
            return "Eukaryota"

    # If we have cellular organisms, it's likely prokaryotic
    if 'cellular organisms' in [name.lower() for name in lineage.values() if name]:
        # Check for archaea indicators
        archaea_indicators = ['archaea', 'archaeal', 'methanobrevibacter', 'methanococcus']
        for rank, name in lineage.items():
            if name and any(indicator in name.lower() for indicator in archaea_indicators):
                return "Archaea"
        # Default to Bacteria for cellular organisms
        return "Bacteria"

    return "Unknown"

def process_taxid_chunk(taxid_chunk, target_rank):
    """Process a chunk of taxids using cached lineage lookup"""
    results = []
    unmapped = []

    for taxid in taxid_chunk:
        lineage = get_full_lineage_cached(taxid)
        target_name = lineage.get(target_rank)
        # Use improved domain detection
        domain = detect_domain_from_lineage(lineage)

        if target_name:
            results.append({
                'taxid': taxid,
                target_rank: target_name,
                'domain': domain
            })
        else:
            unmapped.append({
                'taxid': taxid,
                'reason': f'No {target_rank} found in lineage'
            })

    return results, unmapped

def process_assembly_file_parallel(assembly_file, taxid_to_rank, taxid_to_parent, taxid_to_name, target_rank, chunk_size=50000, n_processes=None):
    """Process assembly file with parallel processing and chunking"""
    try:
        # Convert to Path object for better handling
        assembly_path = Path(assembly_file)

        # If file doesn't exist, try metadata folder
        if not assembly_path.exists():
            script_dir = Path(__file__).resolve().parent
            metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
            if metadata_assembly.exists():
                assembly_path = metadata_assembly
                logging.info(f"📁 Using assembly file from metadata folder: {assembly_path}")
            else:
                raise FileNotFoundError(f"Assembly file not found: {assembly_file}")

        # Read assembly file with header handling
        logging.info(f"📖 Reading assembly file: {assembly_path}")
        df = pd.read_csv(assembly_path, sep="\t", low_memory=False)
        if 'taxid' not in df.columns:
            # Try reading with skiprows=1 if taxid column not found
            df = pd.read_csv(assembly_path, sep="\t", low_memory=False, skiprows=1)
            if 'taxid' not in df.columns:
                raise ValueError("Could not find 'taxid' column in assembly file")

        logging.info(f"📊 Processing {len(df):,} assembly entries")

        # Get unique taxids to avoid redundant processing
        unique_taxids = df['taxid'].unique()
        logging.info(f"🔍 Found {len(unique_taxids):,} unique taxids (reduced from {len(df):,})")

        # Set up multiprocessing
        if n_processes is None:
            n_processes = min(cpu_count(), 8)  # Cap at 8 to avoid overwhelming system

        logging.info(f"🚀 Using {n_processes} processes with chunk size {chunk_size:,}")

        # Split unique taxids into chunks
        taxid_chunks = [unique_taxids[i:i + chunk_size] for i in range(0, len(unique_taxids), chunk_size)]
        logging.info(f"📦 Split into {len(taxid_chunks)} chunks")

        # Process chunks in parallel
        start_time = time.time()
        all_results = []
        all_unmapped = []

        with Pool(processes=n_processes, initializer=init_worker,
                 initargs=(taxid_to_rank, taxid_to_parent, taxid_to_name)) as pool:

            # Create tasks for each chunk
            tasks = [(chunk, target_rank) for chunk in taxid_chunks]

            # Process with progress bar
            with tqdm(total=len(tasks), desc="Processing chunks", unit="chunks") as pbar:
                for results, unmapped in pool.starmap(process_taxid_chunk, tasks):
                    all_results.extend(results)
                    all_unmapped.extend(unmapped)
                    pbar.update(1)

        processing_time = time.time() - start_time
        logging.info(f"⏱️  Parallel processing completed in {processing_time:.2f} seconds")

        # Create result DataFrames
        result_df = pd.DataFrame(all_results)
        unmapped_df = pd.DataFrame(all_unmapped)

        # Expand results back to original dataframe size by merging
        logging.info("🔄 Expanding results to match original dataframe...")
        final_df = df[['taxid']].merge(result_df, on='taxid', how='left')

        # Separate mapped and unmapped
        mapped_mask = final_df[target_rank].notna()
        final_mapped = final_df[mapped_mask][['taxid', target_rank, 'domain']].copy()
        final_unmapped_taxids = final_df[~mapped_mask]['taxid'].unique()

        final_unmapped_df = pd.DataFrame({
            'taxid': final_unmapped_taxids,
            'reason': f'No {target_rank} found in lineage'
        })

        return final_mapped, final_unmapped_df

    except Exception as e:
        logging.error(f"Error processing assembly file: {e}")
        sys.exit(1)

def process_assembly_file(assembly_file, taxid_to_rank, taxid_to_parent, taxid_to_name, target_rank):
    """Process assembly file and map taxids to target rank (legacy sequential version)"""
    try:
        # Convert to Path object for better handling
        assembly_path = Path(assembly_file)

        # If file doesn't exist, try metadata folder
        if not assembly_path.exists():
            script_dir = Path(__file__).resolve().parent
            metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
            if metadata_assembly.exists():
                assembly_path = metadata_assembly
                logging.info(f"📁 Using assembly file from metadata folder: {assembly_path}")
            else:
                raise FileNotFoundError(f"Assembly file not found: {assembly_file}")

        # Read assembly file with header handling
        logging.info(f"📖 Reading assembly file: {assembly_path}")
        df = pd.read_csv(assembly_path, sep="\t", low_memory=False)
        if 'taxid' not in df.columns:
            # Try reading with skiprows=1 if taxid column not found
            df = pd.read_csv(assembly_path, sep="\t", low_memory=False, skiprows=1)
            if 'taxid' not in df.columns:
                raise ValueError("Could not find 'taxid' column in assembly file")

        logging.info(f"📊 Processing {len(df):,} assembly entries")

        # Initialize lists for results
        taxids = []
        ranks = []
        domains = []
        unmapped_taxids = []

        # Process each taxid with progress bar
        for taxid in tqdm(df['taxid'], desc=f"Mapping taxids to {target_rank}", unit="taxids"):
            rank = get_rank_target(taxid, target_rank, taxid_to_rank, taxid_to_parent, taxid_to_name)
            if rank:
                taxids.append(taxid)
                ranks.append(rank)
                # Use improved domain detection
                lineage = get_full_lineage_cached(taxid)
                domain = detect_domain_from_lineage(lineage)
                domains.append(domain)
            else:
                unmapped_taxids.append(taxid)

        # Create result DataFrame
        result_df = pd.DataFrame({
            'taxid': taxids,
            target_rank: ranks,
            'domain': domains
        })

        # Create unmapped DataFrame
        unmapped_df = pd.DataFrame({
            'taxid': unmapped_taxids,
            'reason': f'No {target_rank} found in lineage'
        })

        return result_df, unmapped_df
    except Exception as e:
        logging.error(f"Error processing assembly file: {e}")
        sys.exit(1)

def main():
    # Set up logging
    setup_logging()
    logging.info("🚀 Starting optimized phylum mapping process")

    # Set up paths
    script_dir = Path(__file__).resolve().parent

    # Try metadata folder first, then script directory
    metadata_assembly = script_dir.parent / "metadata" / "00assembly_summary_genbank.txt"
    if metadata_assembly.exists():
        assembly_file = metadata_assembly
        logging.info(f"📁 Using assembly file from metadata folder: {assembly_file}")
    else:
        assembly_file = script_dir / "00assembly_summary_genbank.txt"

    # Load taxonomy files
    logging.info("Loading taxonomy files...")
    start_time = time.time()
    nodes_df, names_df = load_taxonomy_files()
    load_time = time.time() - start_time
    logging.info(f"⏱️  Taxonomy files loaded in {load_time:.2f} seconds")

    # Build taxonomy maps
    logging.info("Building taxonomy maps...")
    start_time = time.time()
    taxid_to_rank, taxid_to_parent, taxid_to_name = build_taxonomy_maps(nodes_df, names_df)
    build_time = time.time() - start_time
    logging.info(f"⏱️  Taxonomy maps built in {build_time:.2f} seconds")

    # Ask user for processing method
    use_parallel = True  # Default to parallel processing

    if use_parallel:
        logging.info("🚀 Using parallel processing (recommended for large datasets)")
        # Process assembly file with parallel processing
        start_time = time.time()
        result_df, unmapped_df = process_assembly_file_parallel(
            assembly_file, taxid_to_rank, taxid_to_parent, taxid_to_name, "phylum",
            chunk_size=50000, n_processes=None
        )
        process_time = time.time() - start_time
        logging.info(f"⏱️  Parallel processing completed in {process_time:.2f} seconds")
    else:
        logging.info("🐌 Using sequential processing")
        # Process assembly file sequentially (legacy method)
        start_time = time.time()
        result_df, unmapped_df = process_assembly_file(
            assembly_file, taxid_to_rank, taxid_to_parent, taxid_to_name, "phylum"
        )
        process_time = time.time() - start_time
        logging.info(f"⏱️  Sequential processing completed in {process_time:.2f} seconds")

    # Count Candidatus taxa (but preserve them in output)
    candidatus_count = 0
    if len(result_df) > 0:
        candidatus_count = result_df["phylum"].str.contains("Candidatus", na=False).sum()

    # Set up output paths - output file in script directory, unmapped in error_log
    script_dir = Path(__file__).resolve().parent
    output_dir = script_dir / "error_log"
    output_dir.mkdir(exist_ok=True)

    output_file = script_dir / "taxid_to_phylum.csv"
    unmapped_file = output_dir / "unmapped_taxids_phylum.csv"

    # Save results (preserving all taxa including Candidatus)
    logging.info("💾 Saving results...")
    result_df.to_csv(output_file, index=False)
    unmapped_df.to_csv(unmapped_file, index=False)

    # Log statistics
    total_taxids = len(result_df) + len(unmapped_df)
    mapped_taxids = len(result_df)
    unmapped_taxids = len(unmapped_df)

    logging.info("=" * 60)
    logging.info("📊 FINAL STATISTICS")
    logging.info("=" * 60)
    logging.info(f"Total taxids processed: {total_taxids:,}")
    logging.info(f"Successfully mapped taxids: {mapped_taxids:,}")
    logging.info(f"Unmapped taxids: {unmapped_taxids:,}")
    logging.info(f"Candidatus taxa found and preserved: {candidatus_count:,}")
    logging.info(f"Mapping success rate: {(mapped_taxids/total_taxids)*100:.2f}%")
    logging.info(f"Processing speed: {total_taxids/process_time:.0f} taxids/second")
    logging.info(f"Results saved to: {output_file}")
    logging.info(f"Unmapped taxids saved to: {unmapped_file}")
    logging.info("=" * 60)

if __name__ == "__main__":
    main()
