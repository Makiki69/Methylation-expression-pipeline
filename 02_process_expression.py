
import pandas as pd
import os

# Path to raw expression matrix
EXPR_FILE = "RNASeq_Expression"

# Path to expression probe annotation (maps expression IDs → gene names)
PROBE_FILE = f"expression_probe"

# Output directory for processed files
OUT_DIR = f"processed"

# Creates processed/ directory if it does not already exist
os.makedirs(OUT_DIR, exist_ok=True)

# Output file for gene-level expression
OUT_EXPR = f"{OUT_DIR}/expression_gene_level.csv"

print("Loading expression matrix...")

# First in RNASeq_Expression is the expression feature ID (e.g., RP11-433M22.1)
expr = pd.read_csv(EXPR_FILE, sep="\t")

# Rename the first column to 'id' so it matches the probe annotation
expr = expr.rename(columns={"sample": "id"})

print("Loading expression probe annotation...")
probe = pd.read_csv(PROBE_FILE, sep="\t")

# Keep only the columns needed for mapping: id → gene
probe = probe[["id", "gene"]]

print("Merging expression matrix with gene names...")
expr_annot = probe.merge(expr, on="id", how="inner")

# Now we have: gene, id, sample1, sample2, ...
# Drop the expression ID and aggregate to gene level
print("Aggregating expression to gene level...")
expr_gene = expr_annot.drop(columns=["id"]).groupby("gene").mean()

# Save final gene-level expression matrix
expr_gene.to_csv(OUT_EXPR)

print("Done. Saved to:", OUT_EXPR)
