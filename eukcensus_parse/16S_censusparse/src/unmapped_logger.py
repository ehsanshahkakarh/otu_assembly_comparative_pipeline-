"""
Unmapped Taxa Logging for 16S Census Parser
===========================================

Creates comprehensive logs of unmapped taxonomic names with
pattern analysis and recommendations.
"""

from datetime import datetime
from pathlib import Path


def create_comprehensive_unmapped_log(phylum_data, family_data, genus_data,
                                    phylum_to_taxid, family_to_taxid, genus_to_taxid,
                                    taxid_to_lineage, log_dir, output_prefix):
    """
    Create a comprehensive log of all unmapped taxonomic names with enhanced analysis.

    Args:
        phylum_data, family_data, genus_data: Data dictionaries for each rank
        phylum_to_taxid, family_to_taxid, genus_to_taxid: Taxid mapping dictionaries
        taxid_to_lineage: Lineage information dictionary
        log_dir: Directory to store log files
        output_prefix: Prefix for output files
    """
    log_file = log_dir / f"{output_prefix}_comprehensive_unmapped.log"
    print(f"📝 Creating enhanced comprehensive unmapped log...")

    with open(log_file, 'w') as f:
        f.write("# Enhanced Comprehensive Unmapped Names Log - 16S Census Parser with Fallback Analysis\n")
        f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# This log contains all taxonomic names that failed to get NCBI taxids or lineages\n")
        f.write("# Enhanced with fallback strategy analysis and improved pattern recognition\n")
        f.write("# Format: Rank | Original_Name | Cleaned_Name | Appropriate_Name | OTU_Count | Size_Count | Taxid | Reason\n\n")

        # Summary statistics
        f.write("=== SUMMARY STATISTICS ===\n")

        # Calculate statistics for each rank
        rank_stats = {}
        for rank_name, data_dict, taxid_dict in [
            ('phylum', phylum_data, phylum_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            total = len(data_dict)
            mapped = len([t for t in taxid_dict.values() if t != 'NA'])
            unmapped = total - mapped
            unmapped_pct = (unmapped / total * 100) if total > 0 else 0

            rank_stats[rank_name] = {
                'total': total,
                'mapped': mapped,
                'unmapped': unmapped,
                'unmapped_pct': unmapped_pct
            }

            f.write(f"{rank_name.capitalize()} Statistics:\n")
            f.write(f"  Total: {total}\n")
            f.write(f"  Mapped: {mapped} ({(mapped/total*100):.1f}%)\n")
            f.write(f"  Unmapped: {unmapped} ({unmapped_pct:.1f}%)\n\n")

        # Overall statistics
        total_all = sum(stats['total'] for stats in rank_stats.values())
        mapped_all = sum(stats['mapped'] for stats in rank_stats.values())
        unmapped_all = total_all - mapped_all
        overall_unmapped_pct = (unmapped_all / total_all * 100) if total_all > 0 else 0

        f.write(f"Overall Statistics:\n")
        f.write(f"  Total entries: {total_all}\n")
        f.write(f"  Successfully mapped: {mapped_all} ({(mapped_all/total_all*100):.1f}%)\n")
        f.write(f"  Failed to map: {unmapped_all} ({overall_unmapped_pct:.1f}%)\n\n")

        # Detailed unmapped entries by rank
        for rank_name, data_dict, taxid_dict in [
            ('phylum', phylum_data, phylum_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            f.write(f"=== {rank_name.upper()} LEVEL UNMAPPED NAMES ===\n")

            unmapped_entries = []
            for orig_name, data in data_dict.items():
                taxid = taxid_dict.get(orig_name, "NA")
                if taxid == "NA" or (taxid != "NA" and taxid not in taxid_to_lineage):
                    # Determine reason for failure
                    if taxid == "NA":
                        reason = "NO_TAXID_FOUND"
                    else:
                        reason = "TAXID_NO_LINEAGE"

                    unmapped_entries.append({
                        'original_name': orig_name,
                        'cleaned_name': data.get('cleaned_name', ''),
                        'appropriate_name': data.get('appropriate_name', ''),
                        'otu_count': data['otu_count'],
                        'size_count': data.get('size_count', 0),
                        'taxid': taxid,
                        'reason': reason
                    })

            f.write(f"Total unmapped {rank_name} entries: {len(unmapped_entries)}\n\n")

            # Sort by OTU count (descending)
            unmapped_entries.sort(key=lambda x: x['otu_count'], reverse=True)

            for entry in unmapped_entries:
                f.write(f"{rank_name.upper()} | {entry['original_name']} | {entry['cleaned_name']} | {entry['appropriate_name']} | {entry['otu_count']} | {entry['size_count']} | {entry['taxid']} | {entry['reason']}\n")

            f.write(f"\n")

        # Pattern analysis
        f.write("=== PATTERN ANALYSIS ===\n")

        # Collect all unmapped names for pattern analysis
        all_unmapped = []
        for rank_name, data_dict, taxid_dict in [
            ('phylum', phylum_data, phylum_to_taxid),
            ('family', family_data, family_to_taxid),
            ('genus', genus_data, genus_to_taxid)
        ]:
            for orig_name, data in data_dict.items():
                taxid = taxid_dict.get(orig_name, "NA")
                if taxid == "NA" or (taxid != "NA" and taxid not in taxid_to_lineage):
                    all_unmapped.append({
                        'rank': rank_name,
                        'name': orig_name,
                        'otu_count': data['otu_count'],
                        'size_count': data.get('size_count', 0)
                    })

        # Analyze patterns
        patterns = {
            'Candidatus': [n for n in all_unmapped if 'Candidatus' in n['name']],
            'Organelles': [n for n in all_unmapped if any(org in n['name'] for org in ['Mitochondria', 'Chloroplast', 'Apicoplast', 'Plastid', ':plas', ':mito', ':api'])],
            'Uncultured': [n for n in all_unmapped if 'uncultured' in n['name'].lower()],
            'Environmental': [n for n in all_unmapped if any(env in n['name'].lower() for env in ['metagenome', 'environmental', 'marine'])],
            'Species_level': [n for n in all_unmapped if '_' in n['name'] and '.' not in n['name']],
            'Modern_phyla': [n for n in all_unmapped if n['name'].endswith('ota') or n['name'].endswith('eia')],
            'Candidate_division': [n for n in all_unmapped if 'candidate division' in n['name'].lower()]
        }

        for pattern_name, pattern_names in patterns.items():
            if pattern_names:
                count = len(pattern_names)
                total_otu_occurrences = sum(n['otu_count'] for n in pattern_names)
                total_size_occurrences = sum(n['size_count'] for n in pattern_names)
                f.write(f"{pattern_name.replace('_', ' ').title()}: {count} names ({total_otu_occurrences} total OTU occurrences, {total_size_occurrences} total sequence occurrences)\n")

                # Show top 5 most frequent for each pattern
                top_names = sorted(pattern_names, key=lambda x: x['otu_count'], reverse=True)[:5]
                for name_info in top_names:
                    f.write(f"  {name_info['name']} ({name_info['rank']}) - {name_info['otu_count']} OTU occurrences, {name_info['size_count']} sequence occurrences\n")
                f.write("\n")

        f.write("=== RECOMMENDATIONS ===\n")
        f.write("1. Organelle sequences (mitochondria, chloroplast, plastid, apicoplast) are now mapped to host organisms\n")
        f.write("2. Original organellar names are preserved to prevent merging with non-organellar OTUs\n")
        f.write("3. Candidatus taxa are now properly preserved in NCBI taxonomy (no filtering)\n")
        f.write("4. Modern bacterial phylum names may need mapping to older NCBI names\n")
        f.write("5. Environmental/uncultured samples may not have valid NCBI taxids\n")
        f.write("6. Species-level entries need genus extraction for genus-level parsing\n")
        f.write("7. Consider using GTDB taxonomy database for better bacterial coverage\n")

    print(f"✅ Comprehensive unmapped log written to {log_file}")
    return log_file

