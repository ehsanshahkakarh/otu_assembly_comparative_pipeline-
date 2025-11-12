import pandas as pd

# Load and concatenate the GTDB taxonomy files
dfs = [pd.read_csv(f, sep="\t", header=None, names=["accession", "taxonomy"])
       for f in ["00bac120_taxonomy.tsv", "00ar53_taxonomy.tsv"]]
df = pd.concat(dfs, ignore_index=True)

# Clean accession and extract phylum
df = (
    df.assign(
        accession_clean=lambda x: x["accession"].str.replace(r"^[A-Z]{2}_", "", regex=True),  # remove RS_/GB_
        phylum=lambda x: x["taxonomy"].str.extract(r"p__([A-Za-z0-9_\-]+)")
    )
    .dropna(subset=["phylum"])
)

# Save detailed mapping for merging
df[["accession_clean", "phylum"]].to_csv("gtdb_accession_phylum_table.csv", index=False)

