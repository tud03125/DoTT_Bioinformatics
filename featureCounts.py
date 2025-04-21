import os
import subprocess
import pandas as pd
import pysam

def is_bam_paired_end(bam_file, max_reads_check=100000):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for i, read in enumerate(bam):
            if read.is_paired:
                return True
            if i >= max_reads_check:
                break
    return False

def run_featurecounts(saf_file, extension, bam_files, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_counts = os.path.join(output_dir, f"3utr_{extension}bp_counts.txt")
    first_bam = bam_files[0]
    paired_end = is_bam_paired_end(first_bam)
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
    subprocess.run(cmd, check=True)
    print(f"featureCounts output saved to {output_counts}")
    return output_counts

def filter_counts_for_deseq2(counts_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    counts_df = pd.read_csv(counts_file, sep="\t", comment="#", low_memory=False)
    columns_to_drop = ["Chr", "Start", "End", "Strand", "Length"]
    counts_df_cleaned = counts_df.drop(columns=columns_to_drop, errors='ignore')
    output_cleaned = os.path.join(output_dir, counts_file.split("/")[-1].replace(".txt", "_cleaned.txt"))
    counts_df_cleaned.to_csv(output_cleaned, sep="\t", index=False)
    print(f"Cleaned counts matrix saved to: {output_cleaned}")
    return output_cleaned
