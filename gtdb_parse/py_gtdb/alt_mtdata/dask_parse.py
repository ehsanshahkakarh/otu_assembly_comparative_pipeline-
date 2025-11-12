import dask.dataframe as dd

# Load files with Dask
df_bac = dd.read_csv("00bac120_taxonomy.tsv", sep="\t", header=None, names=["accession", "taxonomy"])
df_arc = dd.read_csv("00ar53_taxonomy.tsv", sep="\t", header=None, names=["accession", "taxonomy"])

# Combine
df = dd.concat([df_bac, df_arc])

# Remove prefix
df["accession"] = df["accession"].str.replace(r"^[A-Z]{2}_", "", regex=True)

# Extract phylum
df["phylum"] = df["taxonomy"].str.extract(r"p__([A-Za-z0-9_\-]+)")[0]

# Drop missing phylum rows
df = df.dropna(subset=["phylum"])

# Count genomes per phylum (need .compute() to run)
counts = df["phylum"].value_counts().compute()
counts = counts.reset_index()
counts.columns = ["phylum", "gtdb_genome_count"]
counts.to_csv("gtdb_phylum_counts_dask.csv", index=False)

