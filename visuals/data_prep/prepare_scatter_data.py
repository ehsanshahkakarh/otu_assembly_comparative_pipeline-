#!/usr/bin/env python3
"""
Prepare Scatter Plot Source Data
================================
Reads final_merger outputs and census lineage data to produce
plot-ready CSVs for domain_scatter_panels.R.

Inputs:
  - final_merger/outputs/{16s,18s}_ncbi_merged_{division,family}.csv
  - eukcensus_parse/{16S,18S}_censusparse/output/ (for lineage → phylum/division mapping)

Outputs:
  - scatter_plots/mega_scrip/source_data/{Domain}_{level}_source_data.csv  (6 files)
  - scatter_plots/source_data/{Domain}_{level}_source_data.csv             (mirror)
"""

import pandas as pd
import numpy as np
from pathlib import Path

# ── Paths (relative to this script's location: visuals/data_prep/) ──────────
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent.parent  # 00parse_database/

MERGER_DIR = BASE_DIR / "final_merger" / "outputs"
CENSUS_16S_DIR = BASE_DIR / "eukcensus_parse" / "16S_censusparse" / "output"
CENSUS_18S_DIR = BASE_DIR / "eukcensus_parse" / "18S_censusparse" / "output"

OUT_MEGA = SCRIPT_DIR.parent / "scatter_plots" / "mega_scrip" / "source_data"
OUT_MAIN = SCRIPT_DIR.parent / "scatter_plots" / "source_data"

TOP_N = 10
THRESHOLD = 1.0


# ── Helpers ─────────────────────────────────────────────────────────────────
def calculate_circle_size(isolate_count: pd.Series, genome_count: pd.Series) -> pd.Series:
    """Circle size based on isolate percentage (matches R logic exactly)."""
    iso_pct = np.where(genome_count > 0, (isolate_count / genome_count) * 100, 0)
    return np.where(iso_pct == 0, 25,
           np.where(iso_pct < 10, 20,
           np.where(iso_pct < 50, 15, 10)))


def build_family_to_phylum_map(census_family_file: Path) -> dict:
    """
    Parse census family file to map each family name → its parent phylum/division.
    Uses the lineage + lineage_ranks columns.
    """
    mapping = {}
    if not census_family_file.exists():
        print(f"  ⚠ Census family file not found: {census_family_file}")
        return mapping

    df = pd.read_csv(census_family_file, dtype=str)
    for _, row in df.iterrows():
        name = row.get("Name_to_use", "")
        lineage = str(row.get("lineage", ""))
        ranks = str(row.get("lineage_ranks", ""))

        parts = lineage.split(";")
        rank_parts = ranks.split(";")

        # Find 'phylum' or 'division' in ranks and grab corresponding name
        parent = None
        for target_rank in ("phylum", "division"):
            if target_rank in rank_parts:
                idx = rank_parts.index(target_rank)
                if idx < len(parts):
                    parent = parts[idx]
                    break

        mapping[name] = parent if parent else name
    return mapping


def flag_top_taxa(df: pd.DataFrame) -> pd.DataFrame:
    """Add Is_Top_Novelty and Is_Top_Coverage boolean columns."""
    df = df.copy()

    # Novelty: top N taxa with factor > threshold
    nf_candidates = df.loc[df["novelty_factor"] > THRESHOLD].nlargest(TOP_N, "novelty_factor")
    df["Is_Top_Novelty"] = df.index.isin(nf_candidates.index)

    # Coverage: top N taxa with factor > threshold
    of_candidates = df.loc[df["overrepresentation_factor"] > THRESHOLD].nlargest(TOP_N, "overrepresentation_factor")
    df["Is_Top_Coverage"] = df.index.isin(of_candidates.index)

    return df


def standardise_columns(df: pd.DataFrame, level: str, domain: str,
                         phylum_map: dict | None = None) -> pd.DataFrame:
    """Rename / add columns to match what domain_scatter_panels.R expects."""
    df = df.copy()

    # Rename first column (division or family) → Taxon
    df = df.rename(columns={df.columns[0]: "Taxon"})

    # Circle size
    df["Circle_Size"] = calculate_circle_size(
        df["isolate_count"].fillna(0), df["ncbi_genome_count"].fillna(0)
    )

    # Coverage percentage (not computed in merger; set 0)
    if "coverage_percentage" not in df.columns:
        df["coverage_percentage"] = 0

    # Phylum / Division column
    if domain == "Eukaryota":
        if level == "division":
            df["Division"] = df["Taxon"]
        else:
            if phylum_map:
                df["Division"] = df["Taxon"].map(phylum_map).fillna(df["Taxon"])
            else:
                df["Division"] = df["Taxon"]
    else:
        if level == "division":
            df["Phylum"] = df["Taxon"]
        else:
            if phylum_map:
                df["Phylum"] = df["Taxon"].map(phylum_map).fillna(df["Taxon"])
            else:
                df["Phylum"] = df["Taxon"]

    # Rename match columns for backward compatibility with old source data schema
    rename_map = {
        "direct_name_matches": "direct_matches",
        "direct_taxid_matches": "taxid_matches",
        "hierarchical_taxid_matches": "lineage_matches",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    # Drop columns the scatter script doesn't need
    drop_cols = ["census_taxid", "lineage_name_matches"]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns], errors="ignore")

    return df


def save_csv(df: pd.DataFrame, domain: str, level: str):
    """Write to both output directories."""
    # Map internal 'division' level name to 'phylum' for filenames
    file_level = "phylum" if level == "division" else level
    fname = f"{domain}_{file_level}_source_data.csv"

    for out_dir in (OUT_MEGA, OUT_MAIN):
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / fname
        df.to_csv(path, index=False)
        print(f"  ✓ {path.relative_to(SCRIPT_DIR.parent)}")


# ── Main ────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("Prepare Scatter Plot Source Data")
    print("=" * 60)

    # Pre-build phylum/division lookup from census family files
    print("\nBuilding family → phylum/division lookup from census lineage...")
    phylum_map_16s = build_family_to_phylum_map(
        CENSUS_16S_DIR / "eukcensus16S_by_family.csv"
    )
    print(f"  16S families mapped: {len(phylum_map_16s)}")

    phylum_map_18s = build_family_to_phylum_map(
        CENSUS_18S_DIR / "eukcensus_18S_by_family.csv"
    )
    print(f"  18S families mapped: {len(phylum_map_18s)}")

    # ── Process each merger file ────────────────────────────────────────────
    merger_files = {
        "16s": {"division": "16s_ncbi_merged_division.csv",
                "family":   "16s_ncbi_merged_family.csv"},
        "18s": {"division": "18s_ncbi_merged_division.csv",
                "family":   "18s_ncbi_merged_family.csv"},
    }

    for gene, levels in merger_files.items():
        for level, fname in levels.items():
            fpath = MERGER_DIR / fname
            if not fpath.exists():
                print(f"\n⚠ Missing: {fpath}")
                continue

            raw = pd.read_csv(fpath)
            print(f"\n── {fname} ({len(raw)} rows) ──")

            if gene == "16s":
                # Split by domain → Bacteria & Archaea
                for domain in ("Bacteria", "Archaea"):
                    subset = raw[raw["domain"] == domain].copy()
                    subset = subset[
                        (subset["census_otu_count"] > 0) &
                        (subset["ncbi_species_count"] > 0)
                    ].copy()
                    if subset.empty:
                        print(f"  {domain} {level}: no data after filtering")
                        continue

                    pmap = phylum_map_16s if level == "family" else None
                    subset = standardise_columns(subset, level, domain, pmap)
                    subset = flag_top_taxa(subset)

                    n_nov = subset["Is_Top_Novelty"].sum()
                    n_cov = subset["Is_Top_Coverage"].sum()
                    print(f"  {domain} {level}: {len(subset)} taxa "
                          f"({n_nov} novelty, {n_cov} coverage)")
                    save_csv(subset, domain, level)

            else:  # 18s → Eukaryota
                domain = "Eukaryota"
                subset = raw[
                    (raw["census_otu_count"] > 0) &
                    (raw["ncbi_species_count"] > 0)
                ].copy()
                if subset.empty:
                    print(f"  {domain} {level}: no data after filtering")
                    continue

                pmap = phylum_map_18s if level == "family" else None
                subset = standardise_columns(subset, level, domain, pmap)
                subset = flag_top_taxa(subset)

                n_nov = subset["Is_Top_Novelty"].sum()
                n_cov = subset["Is_Top_Coverage"].sum()
                print(f"  {domain} {level}: {len(subset)} taxa "
                      f"({n_nov} novelty, {n_cov} coverage)")
                save_csv(subset, domain, level)

    print("\n✅ Scatter source data generation complete!")


if __name__ == "__main__":
    main()

