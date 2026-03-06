"""
Organelle Detection and Handling for 16S Census Parser
=======================================================

Handles detection and processing of organellar sequences (chloroplast,
mitochondria, plastid, apicoplast) with vectorized operations for performance.
"""

import pandas as pd
import subprocess
import os


def detect_organelle_type(taxon_name):
    """
    Detect the type of organelle from taxon name.
    
    Args:
        taxon_name: The taxon name to analyze
        
    Returns:
        Tuple of (is_organelle, organelle_type, cleaned_name)
    """
    organelle_patterns = {
        'chloroplast': ['.Chloroplast', ':plas.Chloroplast', '.plas.Chloroplast', 'chloroplast'],
        'mitochondria': ['.Mitochondria', ':mito.Mitochondria', '.mito.Mitochondria', 'mitochondria'],
        'plastid': ['.Plastid', ':plas.Plastid', '.plas.Plastid', 'plastid'],
        'apicoplast': ['.Apicoplast', ':api.Apicoplast', '.api.Apicoplast', 'apicoplast']
    }
    
    taxon_lower = taxon_name.lower()
    
    for organelle_type, patterns in organelle_patterns.items():
        for pattern in patterns:
            if pattern.lower() in taxon_lower:
                # Extract the base name before the organelle indicator
                if '.' in taxon_name:
                    base_name = taxon_name.split('.')[0]
                elif ':' in taxon_name:
                    # Handle cases like "Vitis_vinifera:plas.Chloroplast"
                    base_name = taxon_name.split(':')[0]
                else:
                    base_name = taxon_name
                
                # Further clean any remaining organelle indicators
                if ':plas' in base_name:
                    base_name = base_name.replace(':plas', '')
                if ':mito' in base_name:
                    base_name = base_name.replace(':mito', '')
                if ':api' in base_name:
                    base_name = base_name.replace(':api', '')
                
                return True, organelle_type, base_name
    
    return False, None, taxon_name


def vectorized_organelle_detection(taxon_names):
    """
    Vectorized organelle detection and host extraction using pandas operations.
    Much faster than individual name processing.

    Args:
        taxon_names: List of taxon names

    Returns:
        tuple: (organellar_names_dict, all_lookup_names)
            - organellar_names_dict: {original_name: host_organism_name}
            - all_lookup_names: list of names to use for taxid lookup
    """
    # Convert to pandas Series for vectorized operations
    names_series = pd.Series(taxon_names)

    # Define organelle patterns
    organelle_patterns = ['.Chloroplast', '.Mitochondria', '.Apicoplast', '.Plastid',
                         ':Chloroplast', ':Mitochondria', ':Apicoplast', ':Plastid']

    # Vectorized organelle detection
    is_organellar = names_series.str.contains('|'.join([p.replace('.', r'\.').replace(':', r':') for p in organelle_patterns]), na=False)

    organellar_names = {}
    lookup_names = []

    if is_organellar.any():
        # Process organellar sequences
        organellar_series = names_series[is_organellar]

        # Vectorized host extraction
        for pattern in organelle_patterns:
            mask = organellar_series.str.contains(pattern.replace('.', r'\.').replace(':', r':'), na=False)
            if mask.any():
                # Extract host organisms for this pattern
                hosts = organellar_series[mask].str.split(pattern, expand=True)[0].str.replace('_', ' ').str.strip()

                # Map original names to host names
                for orig_name, host_name in zip(organellar_series[mask].index, hosts):
                    original_name = names_series.iloc[orig_name]
                    if host_name:
                        organellar_names[original_name] = host_name
                        lookup_names.append(host_name)
                    else:
                        lookup_names.append(original_name)

    # Add non-organellar names
    non_organellar_names = names_series[~is_organellar].tolist()
    lookup_names.extend(non_organellar_names)

    print(f"🧬 Vectorized organelle detection: {len(organellar_names)} organellar sequences found")
    if len(organellar_names) > 0:
        print(f"⚡ PERFORMANCE BOOST: Avoided {len(organellar_names) * 2} individual subprocess calls!")

    return organellar_names, lookup_names


def recover_organelle_taxonomy(taxon_name, target_rank):
    """
    Recover taxonomic information from organelle sequences by extracting host organism info.

    NOTE: This function is currently DISABLED in the main processing workflow to prevent
    performance bottlenecks. It's kept for potential future use or debugging.
    Organellar sequences are now handled by vectorized_organelle_detection() instead.

    Args:
        taxon_name: Original taxonomic name (e.g., "Prasiola_crispa.Chloroplast")
        target_rank: Target rank (phylum, family, genus)

    Returns:
        Tuple of (host_species_name, appropriate_rank_name) or (None, None) if recovery fails
    """
    # Quick check if this is an organelle sequence - return early if not
    organelle_patterns = ['.Chloroplast', '.Mitochondria', '.Apicoplast', '.Plastid', ':Chloroplast', ':Mitochondria', ':Apicoplast', ':Plastid']

    host_species = None
    for pattern in organelle_patterns:
        if pattern in taxon_name:
            # Extract host species name
            host_species = taxon_name.split(pattern)[0]
            break

    # Early return if no organelle pattern found - this avoids expensive subprocess calls
    if not host_species:
        return None, None

    # Clean the host species name
    host_species = host_species.replace('_', ' ').strip()

    # Get taxid and lineage for host species
    try:
        # Use taxonkit to get taxid for host species
        env = os.environ.copy()
        result = subprocess.run(
            ["taxonkit", "name2taxid"],
            input=host_species,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
            cwd="."
        )

        if result.returncode == 0 and result.stdout.strip():
            parts = result.stdout.strip().split('\t')
            if len(parts) >= 2 and parts[1] != "0":
                host_taxid = parts[1]

                # Get lineage for host taxid
                lineage_result = subprocess.run(
                    ["taxonkit", "lineage", "-R"],
                    input=host_taxid,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=env,
                    cwd="."
                )

                if lineage_result.returncode == 0 and lineage_result.stdout.strip():
                    lineage_parts = lineage_result.stdout.strip().split('\t')
                    if len(lineage_parts) >= 4:  # taxid, lineage, lineage_taxids, lineage_ranks
                        lineage = lineage_parts[1]
                        lineage_ranks = lineage_parts[3]

                        # Extract appropriate rank from lineage
                        lineage_list = lineage.split(';')
                        ranks_list = lineage_ranks.split(';')

                        # Find the target rank in the lineage
                        for i, rank in enumerate(ranks_list):
                            if rank.lower() == target_rank.lower():
                                if i < len(lineage_list):
                                    return host_species, lineage_list[i]

                        # If exact rank not found, use fallback logic
                        if target_rank == "genus" and len(lineage_list) > 0:
                            # For genus, try to extract genus from species name
                            genus = host_species.split()[0] if ' ' in host_species else host_species
                            return host_species, genus
                        elif target_rank == "family" and "family" in [r.lower() for r in ranks_list]:
                            family_idx = [r.lower() for r in ranks_list].index("family")
                            if family_idx < len(lineage_list):
                                return host_species, lineage_list[family_idx]
                        elif target_rank == "phylum" and "phylum" in [r.lower() for r in ranks_list]:
                            phylum_idx = [r.lower() for r in ranks_list].index("phylum")
                            if phylum_idx < len(lineage_list):
                                return host_species, lineage_list[phylum_idx]

    except Exception as e:
        # If organelle recovery fails, return None
        pass

    return None, None

