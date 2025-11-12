import pandas as pd

# Load the NCBI assembly summary
df = pd.read_csv("00assembly_summary_genbank.txt", sep="\t", skiprows=1, low_memory=False)
df = df.rename(columns={"#assembly_accession": "assembly_accession"})

# Load taxid → phylum mapping
taxmap = pd.read_csv("taxid_to_phylum.csv")  # Must contain: taxid, phylum

# Merge taxonomy
df = df.merge(taxmap, on="taxid", how="left")
df = df.dropna(subset=["phylum"])

# Optional: strip version number from accession (for versionless match), but here we keep it
df["accession_clean"] = df["assembly_accession"]

# Save accession-to-phylum mapping for merge
df[["accession_clean", "phylum"]].to_csv("ncbi_accession_phylum_table.csv", index=False)
print("✅ Cleaned NCBI accession-to-phylum mapping saved.")

