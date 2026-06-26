#!/usr/bin/env python3
"""
18S Tree Annotation Builder
===========================

Emit all iTOL annotation files for the EukCensus-derived 18S tree:

    visuals/trees/18s/18s_nf_branches.txt          (branch coloring by NF bin)
    visuals/trees/18s/18s_division_colorstrip.txt  (genus -> division color)
    visuals/trees/18s/18s_coverage_colorstrip.txt  (NCBI+EukCensus vs EukCensus-only)
    visuals/trees/18s/18s_coverage_branches.txt    (branch highlight where NCBI covered)

Inputs:
    metadata/eukcensus_18S/eukcensus_18S.clusters.97.tsv
    final_merger/outputs/18s_ncbi_merged_{division,family,genus}.csv
    visuals/shared_config/taxonomic_color_mapping.yaml
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "trees" / "py"))
from build_18s_eukcensus_tree import resolve_collisions  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
META = REPO_ROOT / "metadata" / "eukcensus_18S" / "eukcensus_18S.clusters.97.tsv"
MERGED = REPO_ROOT / "final_merger" / "outputs"
PALETTE_YAML = REPO_ROOT / "visuals" / "shared_config" / "taxonomic_color_mapping.yaml"
OUT_DIR = REPO_ROOT / "visuals" / "trees" / "18s"

# NF bins: (lower_exclusive, upper_inclusive, hex_color, legend_label)
NF_BINS = [
    (1.0, 2.0, "#feb24c", "1 < NF <= 2 (near parity)"),
    (2.0, 10.0, "#fd8d3c", "2 < NF <= 10 (moderate)"),
    (10.0, np.inf, "#bd0026", "10 < NF < inf (severe)"),
]
NF_INF_COLOR = "#08306b"
NF_INF_LABEL = "NF = inf (NCBI-absent)"

# Muted fallback palette for divisions not in shared YAML (rare ones in metadata)
FALLBACK_DIV_COLORS = [
    "#aac6d6", "#c8a99f", "#cabad2", "#a8b9a3", "#d6c3a0",
    "#b8a7c0", "#9fb6a8", "#c4a0a8", "#a09fb8", "#b9b2a0",
]


def bin_nf(value):
    """Return (hex_color, legend_label) for an NF value, or None if NF<=1."""
    if pd.isna(value):
        return None
    if np.isinf(value):
        return (NF_INF_COLOR, NF_INF_LABEL)
    if value <= 1.0:
        return None
    for lo, hi, color, label in NF_BINS:
        if lo < value <= hi:
            return (color, label)
    return None


def load_label_maps():
    """Reproduce the tree builder's label maps so annotations align."""
    df = pd.read_csv(META, sep="\t", usecols=["division", "family", "genus"]).dropna()
    return resolve_collisions(
        df["division"].unique(), df["family"].unique(), df["genus"].unique())


def load_division_palette(divisions):
    """Return dict division_name -> hex from shared YAML + fallback."""
    with open(PALETTE_YAML) as f:
        cfg = yaml.safe_load(f)
    pal = dict(cfg.get("eukaryota_colors", {}))
    pal.update({k: v for k, v in cfg.get("special_colors", {}).items()
                if isinstance(k, str)})
    fallback = iter(FALLBACK_DIV_COLORS)
    final = {}
    for d in sorted(divisions):
        final[d] = pal.get(d, next(fallback, "#c0c0c0"))
    return final


def write_header(lines, kind, label, color, legend_labels, legend_colors,
                 strip_width=25):
    """Append a standard iTOL dataset header to `lines` list."""
    lines.append(f"DATASET_{kind}")
    lines.append("SEPARATOR TAB")
    lines.append(f"DATASET_LABEL\t{label}")
    lines.append(f"COLOR\t{color}")
    if kind == "COLORSTRIP":
        lines.append(f"STRIP_WIDTH\t{strip_width}")
        lines.append("MARGIN\t0")
        lines.append("BORDER_WIDTH\t0")
        lines.append("SHOW_INTERNAL\t0")
    lines.append("")
    lines.append(f"LEGEND_TITLE\t{label}")
    lines.append("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_labels)))
    lines.append("LEGEND_COLORS\t" + "\t".join(legend_colors))
    lines.append("LEGEND_LABELS\t" + "\t".join(legend_labels))
    lines.append("")
    lines.append("DATA")


def write_tree_colors_header(lines, label, color, legend_labels, legend_colors):
    lines.append("TREE_COLORS")
    lines.append("SEPARATOR TAB")
    lines.append(f"DATASET_LABEL\t{label}")
    lines.append(f"COLOR\t{color}")
    lines.append("")
    lines.append(f"LEGEND_TITLE\t{label}")
    lines.append("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_labels)))
    lines.append("LEGEND_COLORS\t" + "\t".join(legend_colors))
    lines.append("LEGEND_LABELS\t" + "\t".join(legend_labels))
    lines.append("")
    lines.append("DATA")


def merged_df(rank):
    """Load 18s_ncbi_merged_<rank>.csv with NF coerced to float (inf preserved)."""
    df = pd.read_csv(MERGED / f"18s_ncbi_merged_{rank}.csv")
    df["novelty_factor"] = pd.to_numeric(df["novelty_factor"], errors="coerce")
    return df


def main():
    div_map, fam_map, gen_map = load_label_maps()
    rank_maps = {"division": div_map, "family": fam_map, "genus": gen_map}

    build_nf_branches(rank_maps)
    build_nf_colorstrip(gen_map)
    build_division_colorstrip(gen_map, div_map)
    build_coverage_files(gen_map, fam_map, div_map)
    print("Done.")


def build_nf_colorstrip(gen_map):
    """Per-genus COLORSTRIP companion to 18s_nf_branches.txt — iTOL renders
    legends for COLORSTRIP datasets but not for TREE_COLORS, so this track
    exists primarily to display the NF bin legend on the side panel.
    Genera with NF <= 1 or missing NF are skipped (left blank on the strip)."""
    g = merged_df("genus")
    legend_colors = [c for _, _, c, _ in NF_BINS] + [NF_INF_COLOR]
    legend_labels = [lab for _, _, _, lab in NF_BINS] + [NF_INF_LABEL]
    out = []
    write_header(out, "COLORSTRIP", "NF bin (genus)", "#bd0026",
                 legend_labels=legend_labels,
                 legend_colors=legend_colors,
                 strip_width=12)
    n = 0
    for _, row in g.iterrows():
        taxon = row["genus"]
        if taxon not in gen_map:
            continue
        res = bin_nf(row["novelty_factor"])
        if res is None:
            continue
        color, lab = res
        out.append(f"{gen_map[taxon]}\t{color}\t{lab}")
        n += 1
    (OUT_DIR / "18s_nf_colorstrip.txt").write_text("\n".join(out) + "\n")
    print(f"18s_nf_colorstrip.txt: {n} genera assigned a bin")


def build_nf_branches(rank_maps):
    """One TREE_COLORS file covering NF bins across division/family/genus."""
    legend_colors = [c for _, _, c, _ in NF_BINS] + [NF_INF_COLOR]
    legend_labels = [lab for _, _, _, lab in NF_BINS] + [NF_INF_LABEL]
    out = []
    write_tree_colors_header(out, "Eukaryota NF bins (all ranks)", "#bd0026",
                             legend_labels, legend_colors)
    counts = {lab: 0 for lab in legend_labels}
    counts["NF <= 1 (skipped)"] = 0
    for rank in ("division", "family", "genus"):
        df = merged_df(rank)
        label_map = rank_maps[rank]
        for _, row in df.iterrows():
            taxon = row[rank]
            if taxon not in label_map:
                continue
            res = bin_nf(row["novelty_factor"])
            if res is None:
                counts["NF <= 1 (skipped)"] += 1
                continue
            color, lab = res
            out.append(f"{label_map[taxon]}\tbranch\t{color}\tnormal\t4")
            counts[lab] += 1
    (OUT_DIR / "18s_nf_branches.txt").write_text("\n".join(out) + "\n")
    print("18s_nf_branches.txt:", counts)


def build_division_colorstrip(gen_map, div_map):
    """COLORSTRIP: each genus leaf colored by its division."""
    df = pd.read_csv(META, sep="\t",
                     usecols=["division", "genus"]).dropna()
    pal = load_division_palette(df["division"].unique())
    divs_present = sorted(df["division"].unique())
    out = []
    write_header(out, "COLORSTRIP", "Eukaryota Division (EukCensus)",
                 "#CCCCCC",
                 legend_labels=divs_present,
                 legend_colors=[pal[d] for d in divs_present])
    seen = set()
    for _, r in df.iterrows():
        g, d = r["genus"], r["division"]
        if g in seen:
            continue
        seen.add(g)
        out.append(f"{gen_map[g]}\t{pal[d]}\t{d}")
    (OUT_DIR / "18s_division_colorstrip.txt").write_text("\n".join(out) + "\n")
    print(f"18s_division_colorstrip.txt: {len(seen)} leaves, {len(divs_present)} divisions")


def build_coverage_files(gen_map, fam_map, div_map):
    """COLORSTRIP (per genus leaf) + TREE_COLORS (branches across all ranks)."""
    g = merged_df("genus")
    in_ncbi = g["ncbi_species_count"].fillna(0) > 0
    covered = set(g.loc[in_ncbi, "genus"])
    only = set(g.loc[~in_ncbi, "genus"])

    col_lines = []
    write_header(col_lines, "COLORSTRIP", "Genus genomic coverage",
                 "#cc5500",
                 legend_labels=["NCBI+EukCensus", "EukCensus only"],
                 legend_colors=["#cc5500", "#000000"])
    for genus_name, label in gen_map.items():
        if genus_name in covered:
            col_lines.append(f"{label}\t#cc5500\tNCBI+EukCensus")
        elif genus_name in only:
            col_lines.append(f"{label}\t#000000\tEukCensus only")
    (OUT_DIR / "18s_coverage_colorstrip.txt").write_text("\n".join(col_lines) + "\n")
    print(f"18s_coverage_colorstrip.txt: covered={len(covered)} only={len(only)}")

    branch_lines = []
    write_tree_colors_header(branch_lines,
                             "Eukaryota branch coverage (NCBI+EukCensus)",
                             "#cc5500",
                             legend_labels=["NCBI+EukCensus clade"],
                             legend_colors=["#cc5500"])
    rank_maps = {"genus": gen_map, "family": fam_map, "division": div_map}
    n = 0
    for rank, label_map in rank_maps.items():
        df = merged_df(rank)
        covered_rank = df.loc[df["ncbi_species_count"].fillna(0) > 0, rank]
        for taxon in covered_rank:
            if taxon in label_map:
                branch_lines.append(f"{label_map[taxon]}\tbranch\t#cc5500\tnormal\t4")
                n += 1
    (OUT_DIR / "18s_coverage_branches.txt").write_text("\n".join(branch_lines) + "\n")
    print(f"18s_coverage_branches.txt: {n} covered branches")


if __name__ == "__main__":
    main()
