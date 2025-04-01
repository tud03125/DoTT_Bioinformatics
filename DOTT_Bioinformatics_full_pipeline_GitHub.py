#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
import pysam
import matplotlib.pyplot as plt
import seaborn as sns
import contextlib

# rpy2 and mygene imports
import rpy2.robjects as ro
import rpy2.robjects.pandas2ri as pandas2ri
from rpy2.robjects.packages import importr
import mygene

############################################################
# Function: Dynamic Region Extension
############################################################

def dynamic_region_extension(bam_file, chrom, ref_point, strand, max_extension=10000, step=100, count_threshold=5):
    """
    Dynamically extend a region based on read density.
    For plus strand, ref_point is the transcript end.
    For minus strand, ref_point is the transcript start.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        last_valid = ref_point
        for offset in range(step, max_extension + step, step):
            if strand == "+":
                region_start = ref_point
                region_end = ref_point + offset
            else:
                region_start = max(1, ref_point - offset)
                region_end = ref_point
            count = bam.count(contig=chrom, start=region_start, end=region_end)
            if count < count_threshold:
                break
            last_valid = region_end if strand == "+" else region_start
    if strand == "+":
        return ref_point, last_valid
    else:
        return last_valid, ref_point

############################################################
# Function: Custom GTF Loader for hg19
############################################################

def load_gtf_custom(gtf_file):
    """
    Custom loader for a GTF file.
    Splits each line into 9 fields so that the Attributes field captures all remaining text.
    Used primarily for hg19 GTFs.
    """
    records = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if "\t" in line:
                parts = line.split("\t", 8)
            else:
                parts = line.split(None, 8)
            if len(parts) < 9:
                continue
            records.append(parts)
    df = pd.DataFrame(records, columns=["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
    return df

############################################################
# Function: Generate SAF File
############################################################

def generate_saf(gtf_file, extension, output_dir, dynamic, bam_file_for_dynamic, species, id_type, kgx_file="kgXref.txt.gz"):
    """
    Creates an SAF (Simplified Annotation Format) file from a GTF.
    Depending on the species:
      - For mm39: extracts gene_id from the Attributes.
      - For hg38/hg19: extracts gene_id and transcript_id, strips version numbers,
        and maps to official gene symbols using a kgXref mapping table and mygene.
    Also supports dynamic region extension if a BAM file is provided.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_saf = os.path.join(output_dir, f"3utr_{extension}bp_extended_regions.saf")
    
    if species == "hg19":
        print("Loading GTF file using custom loader for hg19...")
        gtf = load_gtf_custom(gtf_file)
    else:
        print("Loading GTF file using pandas read_csv...")
        gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None,
                          names=["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
    
    print("Filtering for transcript features...")
    transcript = gtf[gtf["Feature"] == "transcript"].copy()
    
    if species in ["hg38", "hg19"]:
        print("Extracting gene_id and transcript_id from Attributes for human GTF...")
        transcript["GeneID"] = transcript["Attributes"].str.extract(r'gene_id\s+"([^"]+)"')
        transcript["TranscriptID"] = transcript["Attributes"].str.extract(r'transcript_id\s+"([^"]+)"')
        transcript["GeneID"] = transcript["GeneID"].apply(lambda x: x.split('.')[0] if pd.notna(x) else x)
        transcript["TranscriptID"] = transcript["TranscriptID"].apply(lambda x: x.split('.')[0] if pd.notna(x) else x)
        
        print("Mapping gene IDs to official gene symbols...")
        kgxref = pd.read_csv(kgx_file, sep="\t", compression="gzip", header=None,
                             names=["kgID", "mRNA", "spID", "spDisplayID", "geneSymbol", 
                                    "refseq", "protAcc", "description", "rfamAcc", "tRnaName"])
        kgxref_dict = pd.Series(kgxref.geneSymbol.values, index=kgxref.kgID).to_dict()
        mg = mygene.MyGeneInfo()
        def map_to_official(row):
            gene_id = row["GeneID"]
            transcript_id = row["TranscriptID"]
            if gene_id == transcript_id:
                try:
                    with contextlib.redirect_stderr(open(os.devnull, "w")):
                        result = mg.query(transcript_id, scopes="ensembl.transcript", fields="symbol", species="human", verbose=False)
                    if result.get("hits"):
                        return result["hits"][0].get("symbol", gene_id)
                    else:
                        return gene_id
                except Exception as e:
                    print(f"Error querying mygene for {transcript_id}: {e}")
                    return gene_id
            else:
                return kgxref_dict.get(gene_id, gene_id)
        transcript["OfficialGene"] = transcript.apply(map_to_official, axis=1)
        gene_col = "OfficialGene"
    else:
        print("Extracting gene_id from Attributes for mouse GTF...")
        transcript["GeneID"] = transcript["Attributes"].str.extract(r'gene_id "([^"]+)"')
        transcript["GeneID"] = transcript["GeneID"].apply(lambda x: x.split('.')[0] if pd.notna(x) else x)
        gene_col = "GeneID"
    
    # Apply region extension (dynamic or fixed)
    if dynamic and bam_file_for_dynamic:
        print("Using dynamic region extension...")
        def calc_dynamic(row):
            ref_point = int(row["End"]) if row["Strand"].strip() == "+" else int(row["Start"])
            start_ext, end_ext = dynamic_region_extension(bam_file_for_dynamic, row["Chr"], ref_point, row["Strand"].strip(), max_extension=extension)
            return pd.Series([start_ext, end_ext])
        transcript[["Start_ext", "End_ext"]] = transcript.apply(calc_dynamic, axis=1)
    else:
        print(f"Extending regions by a fixed {extension} bases...")
        def calculate_fixed(row, extension):
            if row["Strand"].strip() == "+":
                return int(row["End"]), int(row["End"]) + extension
            else:
                return max(int(row["Start"]) - extension, 1), int(row["Start"])
        transcript[["Start_ext", "End_ext"]] = transcript.apply(lambda row: pd.Series(calculate_fixed(row, extension)), axis=1)
    
    print("Creating SAF format...")
    saf_df = transcript[[gene_col, "Chr", "Start_ext", "End_ext", "Strand"]].dropna()
    saf_df.columns = ["GeneID", "Chr", "Start", "End", "Strand"]
    saf_df.to_csv(output_saf, sep="\t", index=False)
    print(f"Extended SAF file saved to {output_saf}")
    
    print("First 5 lines of the SAF file:")
    print(saf_df.head())
    
    return output_saf

############################################################
# Function: Check if BAM File is Paired-End
############################################################

def is_bam_paired_end(bam_file, max_reads_check=100000):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for i, read in enumerate(bam):
            if read.is_paired:
                return True
            if i >= max_reads_check:
                break
    return False

############################################################
# Function: Run featureCounts
############################################################

def run_featurecounts(saf_file, extension, bam_files, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_counts = os.path.join(output_dir, f"3utr_{extension}bp_counts.txt")
    first_bam = bam_files[0]
    paired_end = is_bam_paired_end(first_bam)
    print(f"Auto-detect: {first_bam} is {'paired-end' if paired_end else 'single-end'}")
    cmd = [
        "featureCounts",
        "-a", saf_file,
        "-F", "SAF",
        "-o", output_counts,
        "-g", "GeneID"
    ]
    if paired_end:
        cmd.insert(1, "-p")
    cmd.extend(bam_files)
    print(f"Running featureCounts with command:\n{' '.join(cmd)}\n")
    subprocess.run(cmd, check=True)
    print(f"featureCounts output saved to {output_counts}")
    return output_counts

############################################################
# Function: Filter Counts for DESeq2
############################################################

def filter_counts_for_deseq2(counts_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    print(f"Loading counts file: {counts_file}")
    counts_df = pd.read_csv(counts_file, sep="\t", comment="#", low_memory=False)
    print("Dropping unnecessary columns...")
    columns_to_drop = ["Chr", "Start", "End", "Strand", "Length"]
    counts_df_cleaned = counts_df.drop(columns=columns_to_drop, errors='ignore')
    output_cleaned = os.path.join(output_dir, os.path.basename(counts_file).replace(".txt", "_cleaned.txt"))
    counts_df_cleaned.to_csv(output_cleaned, sep="\t", index=False)
    print(f"Cleaned counts matrix saved to: {output_cleaned}")
    return output_cleaned

############################################################
# Function: Run DESeq2 Analysis in Python (using rpy2)
############################################################

def run_deseq2_in_python(cleaned_counts_file, extension, output_dir, condition_list):
    os.makedirs(output_dir, exist_ok=True)
    pandas2ri.activate()
    print("Importing DESeq2 from R...")
    DESeq2 = importr('DESeq2')
    EnhancedVolcano = importr('EnhancedVolcano')
    grdevices = importr('grDevices')
    print(f"Reading cleaned counts from {cleaned_counts_file}...")
    counts_df = pd.read_csv(cleaned_counts_file, sep="\t", index_col=0)
    sample_columns = counts_df.columns
    if len(condition_list) != len(sample_columns):
        raise ValueError("Mismatch between provided condition list and sample columns.")
    counts_matrix = counts_df.astype(int)
    print("counts_matrix.shape =", counts_matrix.shape)
    counts_matrix_r = pandas2ri.py2rpy(counts_matrix)
    coldata_df = pd.DataFrame({"condition": condition_list}, index=counts_matrix.columns)
    coldata_r = pandas2ri.py2rpy(coldata_df)
    print("Creating DESeq2 dataset in R...")
    ro.r('library(DESeq2)')
    dds = DESeq2.DESeqDataSetFromMatrix(countData=counts_matrix_r,
                                        colData=coldata_r,
                                        design=ro.Formula("~ condition"))
    print("Running DESeq2 differential expression analysis...")
    dds = DESeq2.DESeq(dds)
    res = DESeq2.results(dds)
    res_df_r = ro.r('as.data.frame')(res)
    res_df = pandas2ri.rpy2py(res_df_r)
    deseq_results_file = os.path.join(output_dir, f"3_UTR_{extension}bp_deseq2_results.csv")
    res_df.to_csv(deseq_results_file)
    print(f"DESeq2 results saved to {deseq_results_file}")
    print("Creating MA plot...")
    ma_png = os.path.join(output_dir, f"3_UTR_{extension}bp_MA_plot.png")
    grdevices.png(filename=ma_png, width=800, height=600)
    ro.r.plotMA(res, main="MA Plot", ylim=ro.FloatVector([-2, 2]))
    grdevices.dev_off()
    print("Creating Volcano plot...")
    vol_png = os.path.join(output_dir, f"3_UTR_{extension}bp_Volcano_plot.png")
    grdevices.png(filename=vol_png, width=800, height=600)
    volcano_plot = EnhancedVolcano.EnhancedVolcano(
        res,
        lab=ro.r('rownames')(res_df_r),
        x='log2FoldChange',
        y='pvalue',
        title='Differential 3UTR Expression',
        pCutoff=0.05
    )
    ro.r("plot")(volcano_plot)
    grdevices.dev_off()
    print("Filtering for significant genes...")
    res_df.dropna(subset=["padj"], inplace=True)
    sig_genes = res_df[(res_df["padj"] < 0.05) & (res_df["log2FoldChange"] > 1)]
    significant_file = os.path.join(output_dir, "significant_extended_genes.csv")
    sig_genes.to_csv(significant_file, index=True)
    print(f"Significant genes saved to {significant_file}")
    return significant_file

############################################################
# Function: Merge Intervals and Gene Regions
############################################################

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

############################################################
# Function: Extract Coordinates for Significant Genes
############################################################

def extract_coordinates(sig_genes_csv, saf_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_coords = os.path.join(output_dir, "unique_extended_genes_coordinates.txt")
    print(f"Extracting coordinates from {saf_file} for genes in {sig_genes_csv}...")
    sig_genes_df = pd.read_csv(sig_genes_csv, index_col=0)
    gene_ids = sig_genes_df.index.dropna().unique()
    saf_columns = ["GeneID", "Chr", "Start", "End", "Strand"]
    saf_df = pd.read_csv(saf_file, sep="\t", header=0, names=saf_columns)
    found_coords = saf_df[saf_df["GeneID"].isin(gene_ids)].drop_duplicates()
    found_coords.to_csv(output_coords, sep="\t", index=False)
    print(f"Unique extended gene coordinates saved to {output_coords}")
    return output_coords

############################################################
# Helper: Parse Condition Labels
############################################################

def parse_conditions(bam_files, conditions):
    """
    Parses a comma-separated list of condition labels.
    If not provided, splits the BAM files equally into 'Control' and 'Experimental' groups.
    """
    if conditions:
        cond_list = [cond.strip() for cond in conditions.split(",")]
        if len(cond_list) != len(bam_files):
            raise ValueError("Number of conditions provided does not match number of BAM files")
        return cond_list
    else:
        if len(bam_files) % 2 != 0:
            raise ValueError("Number of BAM files must be even if conditions are not provided")
        half = len(bam_files) // 2
        return ["Control"] * half + ["Experimental"] * half

############################################################
# Main Pipeline Entry Point
############################################################

def main():
    parser = argparse.ArgumentParser(description="DoTT Analysis Pipeline")
    parser.add_argument("--gtf-file", required=True, help="Path to the GTF annotation file")
    parser.add_argument("--bam-files", nargs="+", required=True, help="List of BAM file paths")
    parser.add_argument("--species", choices=["mm39", "hg38", "hg19"], required=True,
                        help="Species option (mm39 for mouse, hg38/hg19 for human)")
    parser.add_argument("--id-type", choices=["Ensembl_ID", "Symbol"], default="Symbol",
                        help="Gene identifier type for human annotations")
    parser.add_argument("--extension", type=int, default=10000, help="Fixed extension length in bases")
    parser.add_argument("--dynamic", action="store_true", help="Use dynamic region extension")
    parser.add_argument("--output-dir", required=True, help="Directory to store output files")
    parser.add_argument("--conditions", help="Comma-separated list of condition labels for each BAM file (optional)")
    parser.add_argument("--kgx-file", default="kgXref.txt.gz", help="Path to the kgXref mapping file (for human GTFs)")
    
    args = parser.parse_args()
    
    # Determine condition labels
    condition_list = parse_conditions(args.bam_files, args.conditions)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Generate the SAF file based on species and parameters
    saf_file = generate_saf(
        gtf_file=args.gtf_file,
        extension=args.extension,
        output_dir=args.output_dir,
        dynamic=args.dynamic,
        bam_file_for_dynamic=args.bam_files[0] if args.dynamic else None,
        species=args.species,
        id_type=args.id_type,
        kgx_file=args.kgx_file
    )
    
    # Run featureCounts
    counts_file = run_featurecounts(saf_file, args.extension, args.bam_files, args.output_dir)
    
    # Filter counts for DESeq2
    cleaned_counts = filter_counts_for_deseq2(counts_file, args.output_dir)
    
    # Run DESeq2 analysis
    significant_genes = run_deseq2_in_python(cleaned_counts, args.extension, args.output_dir, condition_list)
    
    # Merge gene regions from significant genes
    merged_annotation_file = os.path.join(args.output_dir, "merged_DoTT_annotation.bed")
    merge_gene_regions(significant_genes, saf_file, merged_annotation_file)
    
    # Extract coordinates for significant genes
    extract_coordinates(significant_genes, saf_file, args.output_dir)
    
    print("\n=== Pipeline Finished Successfully! ===")
    print(f"SAF file: {saf_file}")
    print(f"FeatureCounts output: {counts_file}")
    print(f"Cleaned counts: {cleaned_counts}")
    print(f"DESeq2 significant genes: {significant_genes}")
    print(f"Merged annotation file: {merged_annotation_file}")

if __name__ == "__main__":
    main()
