#!/usr/bin/env bash
# Parse lineages from 18S NCBI taxonkit files and organize by phyla
# Usage: ./parse_lineages_to_phyla.sh

shopt -s nocasematch  # case-insensitive matching

# Read phyla names from the phyla file
declare -A phyla_names
phyla_file="18s_ncbi_taxonkit_phyla.txt"

if [[ ! -r "$phyla_file" ]]; then
    echo "Error: Cannot read $phyla_file" >&2
    exit 1
fi

echo "Reading phyla names from $phyla_file..."

# Extract phyla names (skip header lines)
while IFS=$'\t' read -r taxon_name taxid lineage; do
    # Skip header and separator lines
    if [[ "$taxon_name" =~ ^(Taxon_Name|---|=|18S|Generated|Total|Taxa) ]]; then
        continue
    fi
    
    # Skip empty lines
    if [[ -z "$taxon_name" ]]; then
        continue
    fi
    
    # Store phylum name
    phyla_names["$taxon_name"]=1
    echo "  Found phylum: $taxon_name"
done < "$phyla_file"

echo "Found ${#phyla_names[@]} phyla"
echo

# Function to find phylum in lineage
find_phylum_in_lineage() {
    local lineage="$1"
    for phylum in "${!phyla_names[@]}"; do
        if [[ "$lineage" =~ $phylum ]]; then
            echo "$phylum"
            return 0
        fi
    done
    echo "Unknown"
    return 1
}

# Process families and genera files
files_to_process=("18s_ncbi_taxonkit_families.txt" "18s_ncbi_taxonkit_genera.txt")
output_file="18s_ncbi_hierarchical_by_phyla.txt"

# Create associative arrays to store taxa by phylum
declare -A phyla_families
declare -A phyla_genera

echo "Processing taxonkit files..."

for file in "${files_to_process[@]}"; do
    if [[ ! -r "$file" ]]; then
        echo "Warning: Cannot read $file, skipping..." >&2
        continue
    fi
    
    echo "Processing $file..."
    
    # Determine if this is families or genera file
    if [[ "$file" =~ families ]]; then
        taxa_type="FAMILY"
        declare -n current_array=phyla_families
    else
        taxa_type="GENUS"
        declare -n current_array=phyla_genera
    fi
    
    while IFS=$'\t' read -r taxon_name taxid lineage; do
        # Skip header and separator lines
        if [[ "$taxon_name" =~ ^(Taxon_Name|---|=|18S|Generated|Total|Taxa) ]]; then
            continue
        fi
        
        # Skip empty lines
        if [[ -z "$taxon_name" || -z "$lineage" ]]; then
            continue
        fi
        
        # Find which phylum this taxon belongs to
        phylum=$(find_phylum_in_lineage "$lineage")
        
        # Add to the appropriate phylum group
        if [[ -n "${current_array[$phylum]}" ]]; then
            current_array[$phylum]+=$'\n'"$taxon_name"
        else
            current_array[$phylum]="$taxon_name"
        fi
        
    done < "$file"
    
    # Unset the nameref
    unset -n current_array
done

echo "Writing organized results to $output_file..."

# Write organized output
{
    echo "18S NCBI Taxa Organized by Phyla"
    echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Organized from taxonkit lineage data"
    echo "========================================================================"
    echo
    
    # Sort phyla by name
    for phylum in $(printf '%s\n' "${!phyla_names[@]}" | sort); do
        echo "PHYLUM: $phylum"
        echo "------------------------------------------------------------"
        
        # Count families and genera for this phylum
        family_count=0
        genus_count=0
        
        if [[ -n "${phyla_families[$phylum]}" ]]; then
            family_count=$(echo "${phyla_families[$phylum]}" | wc -l)
        fi
        
        if [[ -n "${phyla_genera[$phylum]}" ]]; then
            genus_count=$(echo "${phyla_genera[$phylum]}" | wc -l)
        fi
        
        echo "Families: $family_count | Genera: $genus_count"
        echo
        
        # List families
        if [[ -n "${phyla_families[$phylum]}" ]]; then
            echo "  FAMILIES ($family_count):"
            echo "  --------------------------------------------------"
            echo "${phyla_families[$phylum]}" | sort | sed 's/^/    /'
            echo
        else
            echo "  No families found for $phylum"
            echo
        fi
        
        # List genera
        if [[ -n "${phyla_genera[$phylum]}" ]]; then
            echo "  GENERA ($genus_count):"
            echo "  --------------------------------------------------"
            echo "${phyla_genera[$phylum]}" | sort | sed 's/^/    /'
            echo
        else
            echo "  No genera found for $phylum"
            echo
        fi
        
        echo "========================================================================"
        echo
    done
    
    # Handle unknown phylum taxa
    if [[ -n "${phyla_families[Unknown]}" || -n "${phyla_genera[Unknown]}" ]]; then
        echo "PHYLUM: Unknown (no matching phylum found in lineage)"
        echo "------------------------------------------------------------"
        
        if [[ -n "${phyla_families[Unknown]}" ]]; then
            family_count=$(echo "${phyla_families[Unknown]}" | wc -l)
            echo "  FAMILIES ($family_count):"
            echo "  --------------------------------------------------"
            echo "${phyla_families[Unknown]}" | sort | sed 's/^/    /'
            echo
        fi
        
        if [[ -n "${phyla_genera[Unknown]}" ]]; then
            genus_count=$(echo "${phyla_genera[Unknown]}" | wc -l)
            echo "  GENERA ($genus_count):"
            echo "  --------------------------------------------------"
            echo "${phyla_genera[Unknown]}" | sort | sed 's/^/    /'
            echo
        fi
        
        echo "========================================================================"
        echo
    fi
    
} > "$output_file"

echo "Done! Results saved to: $output_file"
echo
echo "Summary:"
echo "  Total phyla processed: ${#phyla_names[@]}"
echo "  Families organized: $(echo "${phyla_families[@]}" | wc -w)"
echo "  Genera organized: $(echo "${phyla_genera[@]}" | wc -w)"
