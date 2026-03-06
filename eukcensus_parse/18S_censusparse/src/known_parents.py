#!/usr/bin/env python3
"""
Known Parent Taxids Database for Systematic Resolution

This module contains the curated database of taxonomic names (families and genera)
that are not in NCBI taxonomy, mapped to their known parent taxids. The systematic
resolver uses these mappings to generate complete lineages.
"""

# Database of taxonomic names mapped to their parent taxids
# Format: {taxon_name: (parent_taxid, parent_name, notes, rank)}
KNOWN_PARENTS = {
    # ========== FAMILY LEVEL ==========

    # Dinoflagellate families - all map to Dinophyceae (taxid: 2864)
    "Dino-Group-II.U.family": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II-Clade-7": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II-Clade-10-and-11": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II-Clade-6": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-I.U.family": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-I-Clade-5": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-I-Clade-4": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-I-Clade-1": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II-Clade-3": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II-Clade-1": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II-Clade-14": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II-Clade-21": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),
    "Dino-Group-II_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "family"),

    # Stramenopile families - MAST clades map to Stramenopiles (taxid: 33634)
    "MAST-12": ("33634", "Stramenopiles", "Marine Stramenopile clade", "family"),
    "MAST-3": ("33634", "Stramenopiles", "Marine Stramenopile clade", "family"),

    # Amoebozoa families
    "Vermamoebidae": ("554915", "Echinamoebida", "Amoebozoa family", "family"),

    # Excavata families
    "Neobodonidae": ("5752", "Kinetoplastida", "Excavata family", "family"),

    # Ciliate families
    "Tholoniidae": ("5932", "Litostomatea", "Ciliate family", "family"),
    "Ophryoglenida": ("5932", "Litostomatea", "Ciliate order treated as family", "family"),

    # Arthropod families
    "Maxillopoda": ("6657", "Crustacea", "Crustacean subclass", "family"),

    # Oomycete families
    "Haliphthorales": ("4762", "Oomycetes", "Oomycete order", "family"),

    # Land plant families - map to Embryophyta (taxid: 3193)
    "Embryophyceae_XX": ("3193", "Embryophyta", "Archaic name for land plants (Embryophyta)", "family"),

    # High-level taxonomic groups (families)
    "Gyrista.U.family": ("2698737", "Sar", "High-level SAR clade (Gyrista not in NCBI)", "family"),
    "TSAR.U.family": ("2698737", "Sar", "High-level SAR clade", "family"),
    "Archaeplastida.U.family": ("33090", "Viridiplantae", "Archaeplastida → Viridiplantae", "family"),

    # Amoebozoa families
    "Archamoebea.U.family": ("2605435", "Evosea", "Archamoebae group in Evosea", "family"),

    # Apicomplexa families
    "Coccidiomorphea.U.family": ("5794", "Apicomplexa", "Coccidiomorphea group in Apicomplexa", "family"),

    # Oomycete families (additional)
    "Peronosporomycetes_X.U.family": ("4762", "Oomycetes", "Peronosporomycetes class", "family"),

    # Foraminifera families
    "Allogromiina": ("29178", "Foraminifera", "Foraminifera order", "family"),

    # Dinoflagellate families (additional)
    "Nolandellidae": ("2864", "Dinophyceae", "Dinoflagellate family", "family"),

    # Excavata families (additional)
    "Eupetalomonads_X": ("33682", "Euglenozoa", "Euglenozoa group", "family"),

    # Misspelled families (have correct NCBI taxids)
    "Massiteriidae": ("2802177", "Massisteriidae", "Misspelling of Massisteriidae (Cercozoa)", "family"),
    "Skeletonemaceae": ("33848", "Skeletonemataceae", "Misspelling of Skeletonemataceae (diatom)", "family"),

    # Green algae families
    "Ulvales-relatives_X": ("3113", "Ulvales", "Ulvales order (green algae)", "family"),

    # Cercozoa families
    "Paracercomonadidae": ("45108", "Cercomonadidae", "Variant/synonym of Cercomonadidae", "family"),
    "Filoretidae_X": ("136419", "Cercozoa", "Filosa group in Cercozoa", "family"),
    "Filosa-Sarcomonadea.U.family": ("136419", "Cercozoa", "Sarcomonadea group in Cercozoa", "family"),

    # Fungi families
    "Chytridiomycetaceae": ("451435", "Chytridiomycetes", "Chytrid fungi class", "family"),

    # Haptista families
    "Choanocystidae": ("299888", "Choanocystis", "Centroplasthelida genus", "family"),

    # Opisthokonta families
    "Nucleohelea": ("154967", "Nuclearia", "Rotosphaerida genus", "family"),

    # Rhizaria families (Endomyxa)
    "Paradinidae": ("492102", "Paradinium", "Ascetosporea genus", "family"),
    "Paradinida_X": ("2686029", "Paradiniida", "Ascetosporea order", "family"),

    # Stramenopile families (Labyrinthulomycetes)
    "Diplophrydaceae": ("1287658", "Amphitremida", "Labyrinthulomycetes order", "family"),

    # Euglenozoa families
    "Allobodonidae": ("2571230", "Allobodo", "Kinetoplastida genus", "family"),

    # Apicomplexa families (additional)
    "Gregarinomorphea.U.family": ("5794", "Apicomplexa", "Gregarinomorphea group in Apicomplexa", "family"),

    # ========== GENUS LEVEL ==========

    # Dinoflagellate genera - all map to Dinophyceae (taxid: 2864)
    "Dino-Group-II.U.genus": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II-Clade-7_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II-Clade-10-and-11_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II-Clade-6_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-I.U.genus": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-I-Clade-5_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-I-Clade-4_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-I-Clade-1_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II-Clade-3_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II-Clade-1_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II-Clade-14_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II-Clade-21_X": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),
    "Dino-Group-II_XX": ("2864", "Dinophyceae", "Environmental dinoflagellate clade", "genus"),

    # Stramenopile genera - MAST clades map to Stramenopiles (taxid: 33634)
    "MAST-12C": ("33634", "Stramenopiles", "Marine Stramenopile clade", "genus"),
    "MAST-3E": ("33634", "Stramenopiles", "Marine Stramenopile clade", "genus"),

    # Arthropod genera
    "Maxillopoda.U.genus": ("6657", "Crustacea", "Crustacean subclass", "genus"),
    "Maxillopoda_X": ("6657", "Crustacea", "Crustacean subclass", "genus"),

    # Rhizaria genera - map to Thecofilosea (taxid: 1004930)
    "Filosa-Thecofilosea.U.genus": ("1004930", "Thecofilosea", "Cercozoan class", "genus"),

    # Land plant genera - map to Embryophyta (taxid: 3193)
    "Embryophyceae_XX.U.genus": ("3193", "Embryophyta", "Archaic name for land plants (Embryophyta)", "genus"),

    # High-level taxonomic groups (genera)
    "Gyrista.U.genus": ("2698737", "Sar", "High-level SAR clade (Gyrista not in NCBI)", "genus"),
    "TSAR.U.genus": ("2698737", "Sar", "High-level SAR clade", "genus"),
    "Archaeplastida.U.genus": ("33090", "Viridiplantae", "Archaeplastida → Viridiplantae", "genus"),

    # Amoebozoa genera
    "Archamoebea.U.genus": ("2605435", "Evosea", "Archamoebae group in Evosea", "genus"),

    # Apicomplexa genera
    "Coccidiomorphea.U.genus": ("5794", "Apicomplexa", "Coccidiomorphea group in Apicomplexa", "genus"),

    # Oomycete genera (additional)
    "Peronosporomycetes_X.U.genus": ("4762", "Oomycetes", "Peronosporomycetes class", "genus"),

    # Foraminifera genera
    "Allogromiina": ("29178", "Foraminifera", "Foraminifera order", "genus"),

    # Dinoflagellate genera (additional)
    "Nolandellidae_X": ("2864", "Dinophyceae", "Dinoflagellate family", "genus"),

    # Excavata genera (additional)
    "Eupetalomonads_XX": ("33682", "Euglenozoa", "Euglenozoa group", "genus"),

    # Cercozoa genera (additional)
    "Filosa-Sarcomonadea.U.genus": ("136419", "Cercozoa", "Sarcomonadea group in Cercozoa", "genus"),

    # Rhizaria genera (Endomyxa)
    "Paradinida_XX": ("2686029", "Paradiniida", "Ascetosporea order", "genus"),

    # Apicomplexa genera (additional)
    "Gregarinomorphea.U.genus": ("5794", "Apicomplexa", "Gregarinomorphea group in Apicomplexa", "genus"),
    "Eimeriida.U.genus": ("5794", "Apicomplexa", "Eimeriida order in Apicomplexa", "genus"),

    # Radiolaria genera
    "Acanthopodida.U.genus": ("65574", "Acantharea", "Acantharea/Radiolaria group", "genus"),

    # Ciliate genera
    "Hyplorynchus": ("5989", "Haptorida", "Litostomatea ciliate order", "genus"),

    # Dinoflagellate genera (parasitic)
    "Corallicolla": ("88547", "Syndiniales", "Parasitic dinoflagellate order", "genus"),

    # Alveolata genera
    "Alphamonadea_ALPH1": ("2605705", "Alphamonas", "Colpodellida genus", "genus"),
}


def get_parent_info(taxon_name: str):
    """
    Get parent taxid information for a taxonomic name.

    Args:
        taxon_name: Name of the taxon (family or genus) to look up

    Returns:
        Tuple of (parent_taxid, parent_name, notes, rank) or None if not found
    """
    return KNOWN_PARENTS.get(taxon_name)


def get_all_taxa():
    """
    Get list of all taxa in the known parents database.

    Returns:
        List of taxon names
    """
    return list(KNOWN_PARENTS.keys())


def get_taxa_by_rank(rank: str):
    """
    Get all taxa of a specific rank.

    Args:
        rank: Taxonomic rank ('family' or 'genus')

    Returns:
        List of taxon names
    """
    return [
        taxon for taxon, (_, _, _, r) in KNOWN_PARENTS.items()
        if r == rank
    ]


def get_taxa_by_parent(parent_taxid: str):
    """
    Get all taxa that map to a specific parent taxid.

    Args:
        parent_taxid: Parent taxid to search for

    Returns:
        List of taxon names
    """
    return [
        taxon for taxon, (pid, _, _, _) in KNOWN_PARENTS.items()
        if pid == parent_taxid
    ]


def get_statistics():
    """
    Get statistics about the known parents database.

    Returns:
        Dictionary with statistics
    """
    unique_parents = len(set(pid for pid, _, _, _ in KNOWN_PARENTS.values()))

    # Count by rank
    rank_counts = {}
    for taxon, (_, _, _, rank) in KNOWN_PARENTS.items():
        if rank not in rank_counts:
            rank_counts[rank] = 0
        rank_counts[rank] += 1

    # Count by parent
    parent_counts = {}
    for taxon, (pid, pname, _, _) in KNOWN_PARENTS.items():
        if pname not in parent_counts:
            parent_counts[pname] = 0
        parent_counts[pname] += 1
    
    return {
        'total_taxa': len(KNOWN_PARENTS),
        'unique_parents': unique_parents,
        'taxa_by_rank': rank_counts,
        'taxa_by_parent': parent_counts
    }


if __name__ == "__main__":
    # Print database statistics
    stats = get_statistics()
    print("Known Parents Database Statistics")
    print("=" * 60)
    print(f"Total taxa: {stats['total_taxa']}")
    print(f"Unique parent taxa: {stats['unique_parents']}")
    print("\nTaxa by rank:")
    for rank, count in sorted(stats['taxa_by_rank'].items()):
        print(f"  {rank}: {count} taxa")
    print("\nTaxa by parent taxon:")
    for parent, count in sorted(stats['taxa_by_parent'].items(), key=lambda x: -x[1]):
        print(f"  {parent}: {count} taxa")

    print("\n" + "=" * 60)
    print("All taxa in database:")
    for taxon in sorted(KNOWN_PARENTS.keys()):
        pid, pname, notes, rank = KNOWN_PARENTS[taxon]
        print(f"  [{rank:6s}] {taxon} → {pname} (taxid: {pid})")

