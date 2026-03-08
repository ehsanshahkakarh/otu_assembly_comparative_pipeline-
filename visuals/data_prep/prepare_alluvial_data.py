#!/usr/bin/env python3
"""
Prepare Alluvial / Sankey Plot Source Data
==========================================
Reads final_merger division-level outputs and produces alluvial-ready CSVs
with column names matching what the R and Python alluvial/sankey scripts expect.

Key transformation: rename 'division' → 'phylum' (the column name every
alluvial and sankey script references).

Inputs:
  - final_merger/outputs/16s_ncbi_merged_division.csv
  - final_merger/outputs/18s_ncbi_merged_division.csv

Outputs:
  - alluvial_plots/source_data/16s_alluvial_division.csv
  - alluvial_plots/source_data/18s_alluvial_division.csv
"""

import pandas as pd
from pathlib import Path

# ── Paths (relative to this script: visuals/data_prep/) ─────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent  # 00parse_database/

MERGER_DIR = BASE_DIR / "final_merger" / "outputs"
OUT_DIR = SCRIPT_DIR.parent / "alluvial_plots" / "source_data"

# Merger files to process
FILES = {
    "16s": "16s_ncbi_merged_division.csv",
    "18s": "18s_ncbi_merged_division.csv",
}


def prepare_alluvial_csv(gene: str, fname: str) -> None:
    """Read a merger division file, rename division→phylum, and write out."""
    src = MERGER_DIR / fname
    if not src.exists():
        print(f"  ⚠ Missing: {src}")
        return

    df = pd.read_csv(src)
    print(f"\n── {fname} ({len(df)} rows) ──")

    # Core rename: division → phylum
    df = df.rename(columns={"division": "phylum"})

    # Ensure output directory exists
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    out_path = OUT_DIR / f"{gene}_alluvial_division.csv"
    df.to_csv(out_path, index=False)
    print(f"  ✓ {out_path.relative_to(SCRIPT_DIR.parent)}")

    # Summary
    domains = df["domain"].unique() if "domain" in df.columns else ["?"]
    matched = (df["match_status"] == "matched").sum() if "match_status" in df.columns else "?"
    print(f"    Domains: {', '.join(str(d) for d in domains)}")
    print(f"    Matched rows: {matched} / {len(df)}")


def main():
    print("=" * 60)
    print("Prepare Alluvial / Sankey Source Data")
    print("=" * 60)

    for gene, fname in FILES.items():
        prepare_alluvial_csv(gene, fname)

    print("\n✅ Alluvial source data generation complete!")
    print(f"   Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()

