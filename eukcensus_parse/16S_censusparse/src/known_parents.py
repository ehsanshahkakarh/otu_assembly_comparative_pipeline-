#!/usr/bin/env python3
"""
Known Parent Taxids Database for 16S Systematic Resolution
===========================================================

This module contains the curated database of prokaryotic taxonomic names
(phyla, families, and genera) that are not in NCBI taxonomy, mapped to their
known parent taxids. The systematic resolver uses these mappings to generate
complete lineages.

Focus: Bacteria and Archaea (16S rRNA)
- Asgard archaea (Lokiarchaeia, etc.)
- Candidate Phyla Radiation (CPR)
- GTDB-specific names not in NCBI
- Misspellings and alternative names
"""

# Database of taxonomic names mapped to their parent taxids
# Format: {taxon_name: (parent_taxid, parent_name, notes, rank)}
KNOWN_PARENTS = {
    # ========== PHYLUM LEVEL ==========
    
    # Asgard Archaea - Lokiarchaeia is now part of Asgardarchaeota
    # Asgardarchaeota (taxid: 1935183) is the parent phylum
    "Lokiarchaeia": ("1935183", "Asgardarchaeota", "Asgard archaea phylum, includes Loki/Thor/Odin/Heimdall", "phylum"),
    
    # Cyanobacteria/Melainabacteria group
    "Cyanobacteriota/Melainabacteria": ("1798711", "Cyanobacteriota", "Combined Cyanobacteria and Melainabacteria", "phylum"),
    
    # ========== FAMILY LEVEL ==========
    
    # CPR (Candidate Phyla Radiation) families
    # Microgenomatia is part of CPR, maps to Patescibacteria (taxid: 1783273)
    "Microgenomatia": ("1783273", "Patescibacteria", "CPR group, ultra-small bacteria", "family"),
    
    # ABY1 - Environmental bacterial clade, part of Patescibacteria
    "ABY": ("1783273", "Patescibacteria", "CPR environmental clade", "family"),
    
    # WWE3 - Environmental bacterial clade, part of Cloacimonadota
    "WWE": ("2138240", "Cloacimonadota", "WWE3 environmental clade", "family"),
    
    # Lokiarchaeia at family level
    "Lokiarchaeia": ("1935183", "Asgardarchaeota", "Asgard archaea", "family"),
    
    # Spirosomaceae - Actually exists in NCBI (taxid: 89373) but may have spelling variants
    "Spirosomaceae": ("89373", "Spirosomaceae", "Bacteroidota family", "family"),
    
    # Abitibacteriaceae - Part of Fibrobacterota
    "Abitibacteriaceae": ("2058777", "Chitinivibrionia", "Fibrobacterota family", "family"),
    
    # CPR2 - Another CPR clade
    "CPR": ("1783273", "Patescibacteria", "Candidate Phyla Radiation", "family"),
    
    # Kazania - Environmental clade, part of Verrucomicrobiota
    "Kazania": ("74201", "Verrucomicrobiota", "Environmental verrucomicrobial clade", "family"),
    
    # Hydrogenimonaceae - Part of Nitrospirota
    "Hydrogenimonaceae": ("40117", "Nitrospirota", "Nitrospira family", "family"),

    # MD2896-B216 - Environmental clade
    "MD2896-B": ("1783272", "Bacillati", "Environmental bacterial clade", "family"),

    # Procabacteriaceae - Candidatus Procabacter family, Chlamydiae-related
    "Procabacteriaceae": ("809", "Chlamydiaceae", "Family for Candidatus Procabacter, amoeba endosymbionts", "family"),
    
    # ========== GENUS LEVEL ==========
    
    # Microgenomatia genus
    "Microgenomatia": ("1783273", "Patescibacteria", "CPR ultra-small bacteria", "genus"),
    
    # Candidatus genera - Many are valid in NCBI, but some need mapping
    # Gracilibacteria (Candidatus) - Part of Patescibacteria
    "Candidatus Gracilibacteria": ("1783273", "Patescibacteria", "CPR superphylum", "genus"),
    
    # ABY1 genus
    "ABY": ("1783273", "Patescibacteria", "CPR environmental clade", "genus"),
    
    # Marinimicrobia (Candidatus) - Part of Marinimicrobia phylum (taxid: 1117647)
    "Candidatus Marinimicrobia": ("1117647", "Marinimicrobia", "Marine bacterial phylum", "genus"),
    
    # Latescibacteria (Candidatus) - Part of Latescibacterota (taxid: 1850252)
    "Candidatus Latescibacteria": ("1850252", "Latescibacterota", "Bacterial phylum", "genus"),
    
    # Aenigmarchaeota (Candidatus) - Archaeal phylum (taxid: 1935204)
    "Candidatus Aenigmarchaeota": ("1935204", "Aenigmarchaeota", "DPANN archaea", "genus"),
    
    # NC10 - Methylomirabilis-related bacteria
    "candidate division NC": ("1843080", "Methylomirabilota", "Methane-oxidizing bacteria", "genus"),
    
    # WWE3 genus
    "WWE": ("2138240", "Cloacimonadota", "WWE3 environmental clade", "genus"),
    
    # Lokiarchaeia genus
    "Lokiarchaeia": ("1935183", "Asgardarchaeota", "Asgard archaea", "genus"),
    
    # Hydrogenedentes (Candidatus) - Part of Hydrogenedentes phylum (taxid: 1784270)
    "Candidatus Hydrogenedentes": ("1784270", "Hydrogenedentota", "Hydrogen-utilizing bacteria", "genus"),
    
    # Diapherotrites (Candidatus) - DPANN archaea (taxid: 1802293)
    "Candidatus Diapherotrites": ("1802293", "Diapherotrites", "DPANN archaea", "genus"),
    
    # Margulisbacteria (Candidatus) - CPR group
    "Candidatus Margulisbacteria": ("1783273", "Patescibacteria", "CPR group", "genus"),
    
    # Tectomicrobia (Candidatus) - Part of Tectomicrobia phylum (taxid: 1850253)
    "Candidatus Tectomicrobia": ("1850253", "Tectomicrobiota", "Bacterial phylum", "genus"),
    
    # Altiarchaeales (Candidatus) - Archaeal order
    "Candidatus Altiarchaeales": ("28890", "Halobacteria", "Halophilic archaea", "genus"),
    
    # Cloacimonetes (Candidatus) - Part of Cloacimonadota
    "Candidatus Cloacimonetes": ("2138240", "Cloacimonadota", "Anaerobic bacteria", "genus"),
    
    # Spirosomaceae genus
    "Spirosomaceae": ("89373", "Spirosomaceae", "Bacteroidota family", "genus"),
    
    # CPR2 genus
    "CPR": ("1783273", "Patescibacteria", "Candidate Phyla Radiation", "genus"),
    
    # Schekmanbacteria (Candidatus) - CPR group
    "Candidatus Schekmanbacteria": ("1783273", "Patescibacteria", "CPR group", "genus"),
    
    # Nanohaloarchaeota (Candidatus) - DPANN archaea (taxid: 1935205)
    "Candidatus Nanohaloarchaeota": ("1935205", "Nanohaloarchaeota", "DPANN archaea", "genus"),
    
    # Caenarcaniphilales (Candidatus) - Archaeal order
    "Candidatus Caenarcaniphilales": ("28890", "Halobacteria", "Halophilic archaea", "genus"),
    
    # Kazania genus
    "Kazania": ("74201", "Verrucomicrobiota", "Environmental verrucomicrobial clade", "genus"),
    
    # Microgenomates (Candidatus) - CPR group
    "Candidatus Microgenomates": ("1783273", "Patescibacteria", "CPR group", "genus"),
    
    # MD2896-B216 genus
    "MD2896-B": ("1783272", "Bacillati", "Environmental bacterial clade", "genus"),
    
    # Tetrasphaera - Actually exists in NCBI (taxid: 92827)
    "Tetrasphaera": ("92827", "Tetrasphaera", "Actinobacteria genus", "genus"),
    
    # Catenococcus - Actually exists in NCBI (taxid: 1264)
    "Catenococcus": ("1264", "Catenococcus", "Firmicutes genus", "genus"),
    
    # Additional Candidatus genera
    "Candidatus Edwardsbacteria": ("1783273", "Patescibacteria", "CPR group", "genus"),
    "Candidatus Firestonebacteria": ("1783273", "Patescibacteria", "CPR group", "genus"),
    "Candidatus Fervidibacteria": ("1783273", "Patescibacteria", "CPR group", "genus"),
    "Candidatus Desantisbacteria": ("1783273", "Patescibacteria", "CPR group", "genus"),
    "Candidatus Omnitrophica": ("1850254", "Omnitrophota", "Bacterial phylum", "genus"),
    "Candidatus Aerophobetes": ("1783273", "Patescibacteria", "CPR group", "genus"),
    "Candidatus Calescamantes": ("1783273", "Patescibacteria", "CPR group", "genus"),

    # Reclassified genera (2024 taxonomic updates)
    # Coenonia was reclassified to Allocoenonia (Oren & Molinari Novoa 2024)
    "Coenonia": ("78328", "Allocoenonia", "Reclassified to Allocoenonia in 2024, Flavobacteriaceae family", "genus"),

    # Marine Alphaproteobacteria
    # Yangia - Marine bacteria, Rhodobacteraceae (now Paracoccaceae)
    "Yangia": ("31989", "Paracoccaceae", "Marine Alphaproteobacteria, formerly Rhodobacteraceae", "genus"),

    # Chlamydiae-related genera
    # Procabacter - Candidatus genus, amoeba endosymbiont
    "Procabacter": ("809", "Chlamydiaceae", "Candidatus Procabacter, amoeba endosymbiont, Chlamydiae-related", "genus"),
}


def get_parent_info(taxon_name: str):
    """
    Get parent taxid information for a taxonomic name.

    Args:
        taxon_name: Name of the taxon to look up

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
        rank: Taxonomic rank ('phylum', 'family', or 'genus')

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
    print("Known Parents Database Statistics (16S - Prokaryotes)")
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

