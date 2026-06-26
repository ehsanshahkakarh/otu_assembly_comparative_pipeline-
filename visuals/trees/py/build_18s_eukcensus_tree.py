#!/usr/bin/env python3
"""
18S EukCensus Tree Builder
==========================

Build a 3-level Newick tree (division -> family -> genus) directly from the
canonical EukCensus 18S metadata:

    metadata/eukcensus_18S/eukcensus_18S.clusters.97.tsv

Every distinct division, family, and genus value (including all *.U.<rank>
placeholders) becomes a node in the tree. Naming collisions across ranks are
resolved by appending _div / _fam / _gen to the offending labels.

Outputs:
    visuals/trees/18s/18s_genus_tree.nwk     (Newick tree, overwrites old)
    visuals/trees/18s/18s_taxon_size.tsv     (taxon, rank, clusters, total_size)
    visuals/trees/18s/archive/18s_genus_tree.OLD.nwk  (archive of previous)
"""

import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
METADATA = REPO_ROOT / "metadata" / "eukcensus_18S" / "eukcensus_18S.clusters.97.tsv"
OUT_DIR = REPO_ROOT / "visuals" / "trees" / "18s"
ARCHIVE_DIR = OUT_DIR / "archive"
TREE_OUT = OUT_DIR / "18s_genus_tree.nwk"
SIZE_OUT = OUT_DIR / "18s_taxon_size.tsv"

_NWK_SAFE = re.compile(r"[(),:;\s]")


def _safe(label: str) -> str:
    """Strip Newick-reserved characters from a node label."""
    return _NWK_SAFE.sub("_", str(label))


def resolve_collisions(divisions, families, genera):
    """
    Return three dict[str -> str] maps (div/fam/gen -> tree label) where
    cross-rank collisions are resolved with _div / _fam / _gen suffixes.
    """
    div_set, fam_set, gen_set = set(divisions), set(families), set(genera)
    collisions = (div_set & fam_set) | (div_set & gen_set) | (fam_set & gen_set)

    def _map(names, rank_suffix):
        return {n: (_safe(n) + rank_suffix if n in collisions else _safe(n))
                for n in names}

    return (_map(div_set, "_div"),
            _map(fam_set, "_fam"),
            _map(gen_set, "_gen"))


def build_newick(df, div_map, fam_map, gen_map):
    """Build a (((leaves)family)division)Eukaryota Newick string."""
    fam_to_genera = defaultdict(list)
    div_to_families = defaultdict(list)
    seen_fam, seen_gen = set(), set()

    grouped = (df.groupby(["division", "family", "genus"], dropna=False)
                 .size().reset_index(name="cluster_count"))

    for _, row in grouped.iterrows():
        d, f, g = row["division"], row["family"], row["genus"]
        fg_key = (d, f, g)
        if fg_key not in seen_gen:
            seen_gen.add(fg_key)
            fam_to_genera[(d, f)].append(gen_map[g])

        df_key = (d, f)
        if df_key not in seen_fam:
            seen_fam.add(df_key)
            div_to_families[d].append((f, fam_map[f]))

    div_blocks = []
    for d, fams in div_to_families.items():
        fam_blocks = []
        for fam_raw, fam_label in fams:
            genera = fam_to_genera[(d, fam_raw)]
            fam_blocks.append(f"({','.join(genera)}){fam_label}")
        div_blocks.append(f"({','.join(fam_blocks)}){div_map[d]}")

    return f"({','.join(div_blocks)})Eukaryota;"


def write_size_table(df, div_map, fam_map, gen_map, out_path):
    """Per-taxon cluster count and summed cluster size (for future tracks)."""
    rows = []
    for rank, col, label_map in (("division", "division", div_map),
                                 ("family", "family", fam_map),
                                 ("genus", "genus", gen_map)):
        agg = df.groupby(col, dropna=False)["size"].agg(["count", "sum"])
        for raw_name, r in agg.iterrows():
            rows.append({"taxon_raw": raw_name,
                         "tree_label": label_map[raw_name],
                         "rank": rank,
                         "cluster_count": int(r["count"]),
                         "total_size": int(r["sum"])})
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def main():
    if not METADATA.exists():
        sys.exit(f"Missing metadata: {METADATA}")

    print(f"Reading {METADATA.name} ...")
    df = pd.read_csv(METADATA, sep="\t",
                     usecols=["size", "division", "family", "genus"])
    df = df.dropna(subset=["division", "family", "genus"])
    print(f"  {len(df):,} clusters | "
          f"{df['division'].nunique()} divisions | "
          f"{df['family'].nunique()} families | "
          f"{df['genus'].nunique()} genera")

    div_map, fam_map, gen_map = resolve_collisions(
        df["division"].unique(), df["family"].unique(), df["genus"].unique())
    suffixed = sum(1 for k, v in {**div_map, **fam_map, **gen_map}.items()
                   if v != _safe(k))
    print(f"  {suffixed} labels suffixed to resolve cross-rank collisions")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    if TREE_OUT.exists():
        archive_path = ARCHIVE_DIR / "18s_genus_tree.OLD.nwk"
        shutil.move(str(TREE_OUT), str(archive_path))
        print(f"  archived previous tree -> {archive_path.relative_to(REPO_ROOT)}")

    print("Building Newick ...")
    newick = build_newick(df, div_map, fam_map, gen_map)
    TREE_OUT.write_text(newick + "\n")
    print(f"  wrote {TREE_OUT.relative_to(REPO_ROOT)}  ({len(newick):,} chars)")

    print("Writing taxon size table ...")
    write_size_table(df, div_map, fam_map, gen_map, SIZE_OUT)
    print(f"  wrote {SIZE_OUT.relative_to(REPO_ROOT)}")

    print("Done.")


if __name__ == "__main__":
    main()
