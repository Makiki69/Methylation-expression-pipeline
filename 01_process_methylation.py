import pandas as pd
import glob
import os

#Path to methylation matrix
METH_FILE = "DNA_Methylation"

#Path to methylation probe
PROBE_FILE=f"methylation_probe"
# sets output directory for all processed files
OUT_DIR=f"processed"
#Creates the processed/ directory if does not exist already, exist_ok prevents errors if it already exists
os.makedirs(OUT_DIR, exist_ok=True)

OUT_METH = f"{OUT_DIR}/methylation_gene_level.csv"

print("Loading methylation matrix...")

#First column in the methylation_probe is CpG ID, with remaining columns as samples
meth = pd.read_csv(METH_FILE, sep="\t")
meth = meth.rename(columns={"sample": "id"})

print("Loading methylation probe annotation...")    

#Loading and Keeping only id and genecolumns
probe = pd.read_csv(PROBE_FILE, sep="\t", comment="#", header=None)
probe.columns = ["id", "gene", "chrom", "chromStart", "chromEnd", "strand"]

probe = probe[["id", "gene"]]
probe = probe.assign(gene=probe["gene"].str.split(",")).explode("gene")

#Expansion of multi-gene probes
probe = probe.assign(gene=probe["gene"].str.split(",")).explode("gene")

print("Merging CpG methylation with gene names...")
meth_annot = probe.merge(meth, on="id", how="inner")

#should hve gene, id, sample1, sample2,...
#This will allow usto drop CpG id and aggregate to gene level
print("Aggregating methylation to gene level.")
meth_gene = meth_annot.drop(columns=["id"]).groupby("gene").mean()

meth_gene.to_csv(OUT_METH)
print("Done. Saved to:", OUT_METH)
