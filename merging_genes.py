import pandas as pd

def merge_intervals(intervals, max_gap=0):
    if not intervals:
        return []
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current[0] <= last[1] + max_gap:
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            merged.append(current)
    return merged

def merge_gene_regions(sig_genes_csv, saf_file, output_file, max_gap=0):
    sig_df = pd.read_csv(sig_genes_csv, index_col=0)
    genes_of_interest = set(sig_df.index.dropna())
    saf_df = pd.read_csv(saf_file, sep="\t")
    saf_df = saf_df[saf_df["GeneID"].isin(genes_of_interest)]
    merged_regions = []
    for gene, group in saf_df.groupby("GeneID"):
        group = group.sort_values("Start")
        intervals = list(zip(group["Start"], group["End"]))
        merged = merge_intervals(intervals, max_gap)
        chrom = group["Chr"].iloc[0]
        strand = group["Strand"].iloc[0]
        for start, end in merged:
            merged_regions.append([gene, chrom, start, end, strand])
    merged_df = pd.DataFrame(merged_regions, columns=["GeneID", "Chr", "Start", "End", "Strand"])
    merged_df.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Merged annotation saved to {output_file}")
    return output_file
