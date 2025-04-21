import os
import pandas as pd
import pysam
import contextlib
import mygene

def load_gtf_custom(gtf_file):
    records = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            parts = line.split("\t", 8)
            if len(parts) < 9:
                continue
            records.append(parts)
    df = pd.DataFrame(records, columns=["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
    return df

def dynamic_region_extension(bam_file, chrom, ref_point, strand, max_extension=10000, step=100, count_threshold=5):
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
    return (ref_point, last_valid) if strand == "+" else (last_valid, ref_point)

def generate_saf(gtf_file, extension, output_dir, dynamic, bam_file_for_dynamic, species, id_type, kgx_file="kgXref.txt.gz"):
    os.makedirs(output_dir, exist_ok=True)
    output_saf = os.path.join(output_dir, f"3utr_{extension}bp_extended_regions.saf")
    
    if species == "hg19":
        gtf = load_gtf_custom(gtf_file)
    else:
        gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None,
                          names=["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
    
    transcript = gtf[gtf["Feature"] == "transcript"].copy()
    
    if species in ["hg38", "hg19"]:
        transcript["GeneID"] = transcript["Attributes"].str.extract(r'gene_id\s+"([^"]+)"')
        transcript["TranscriptID"] = transcript["Attributes"].str.extract(r'transcript_id\s+"([^"]+)"')
        transcript["GeneID"] = transcript["GeneID"].apply(lambda x: x.split('.')[0] if pd.notna(x) else x)
        transcript["TranscriptID"] = transcript["TranscriptID"].apply(lambda x: x.split('.')[0] if pd.notna(x) else x)
        
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
                except Exception:
                    return gene_id
            else:
                return kgxref_dict.get(gene_id, gene_id)
        transcript["OfficialGene"] = transcript.apply(map_to_official, axis=1)
        gene_col = "OfficialGene"
    else:
        transcript["GeneID"] = transcript["Attributes"].str.extract(r'gene_id "([^"]+)"')
        transcript["GeneID"] = transcript["GeneID"].apply(lambda x: x.split('.')[0] if pd.notna(x) else x)
        gene_col = "GeneID"
    
    if dynamic and bam_file_for_dynamic:
        def calc_dynamic(row):
            ref_point = int(row["End"]) if row["Strand"].strip() == "+" else int(row["Start"])
            start_ext, end_ext = dynamic_region_extension(bam_file_for_dynamic, row["Chr"], ref_point, row["Strand"].strip(), max_extension=extension)
            return pd.Series([start_ext, end_ext])
        transcript[["Start_ext", "End_ext"]] = transcript.apply(calc_dynamic, axis=1)
    else:
        def calculate_fixed(row, extension):
            if row["Strand"].strip() == "+":
                return int(row["End"]), int(row["End"]) + extension
            else:
                return max(int(row["Start"]) - extension, 1), int(row["Start"])
        transcript[["Start_ext", "End_ext"]] = transcript.apply(lambda row: pd.Series(calculate_fixed(row, extension)), axis=1)
    
    saf_df = transcript[[gene_col, "Chr", "Start_ext", "End_ext", "Strand"]].dropna()
    saf_df.columns = ["GeneID", "Chr", "Start", "End", "Strand"]
    saf_df.to_csv(output_saf, sep="\t", index=False)
    
    print(f"Extended SAF file saved to {output_saf}")
    return output_saf
