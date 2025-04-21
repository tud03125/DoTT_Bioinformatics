import os
import pandas as pd

def extract_coordinates(sig_genes_csv, saf_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_coords = os.path.join(output_dir, "unique_extended_genes_coordinates.txt")
    sig_genes_df = pd.read_csv(sig_genes_csv, index_col=0)
    gene_ids = sig_genes_df.index.dropna().unique()
    saf_columns = ["GeneID", "Chr", "Start", "End", "Strand"]
    saf_df = pd.read_csv(saf_file, sep="\t", header=0, names=saf_columns)
    found_coords = saf_df[saf_df["GeneID"].isin(gene_ids)].drop_duplicates()
    found_coords.to_csv(output_coords, sep="\t", index=False)
    print(f"Unique extended gene coordinates saved to {output_coords}")
    return output_coords
