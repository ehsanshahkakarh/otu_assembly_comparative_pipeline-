#!/bin/bash
#
# Generate All Sankey Plots with Colored Flows
# ============================================
#
# This script generates all Sankey diagrams:
# - 16S Bacteria
# - 16S Archaea  
# - 18S Eukaryota
#
# Author: Enhanced Sankey Team
# Date: 2026-01-28
#

echo "========================================================================"
echo "🧬 Generating All Sankey Plots with Colored Flows"
echo "========================================================================"
echo ""

# Change to script directory
cd "$(dirname "$0")"

# 1. Generate Bacteria Sankey
echo "📊 [1/3] Generating Bacteria Sankey diagram..."
echo "----------------------------------------------------------------------"
python colored_sankey_16s.py
echo ""

# 2. Generate Archaea Sankey
echo "📊 [2/3] Generating Archaea Sankey diagram..."
echo "----------------------------------------------------------------------"
# Create temporary script with Archaea configuration
sed 's/PROCESS_DOMAIN = "Bacteria"/PROCESS_DOMAIN = "Archaea"/' colored_sankey_16s.py > temp_archaea_sankey.py
python temp_archaea_sankey.py
rm temp_archaea_sankey.py
echo ""

# 3. Generate Eukaryota Sankey
echo "📊 [3/3] Generating Eukaryota Sankey diagram..."
echo "----------------------------------------------------------------------"
python colored_sankey_18s.py
echo ""

# Summary
echo "========================================================================"
echo "✅ All Sankey Diagrams Generated Successfully!"
echo "========================================================================"
echo ""
echo "📁 Output Files:"
echo "   • sankey_16s_bacteria_colored_flows.html"
echo "   • sankey_16s_archaea_colored_flows.html"
echo "   • sankey_18s_eukaryota_colored_flows.html"
echo ""
echo "💡 Open these HTML files in a web browser to view interactive diagrams"
echo "🎨 All flows are colored to match their taxonomic colors"
echo ""
echo "========================================================================"

