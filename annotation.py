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
    # This function processes ONE bam file path (string) at a time
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # 1. Get chromosome length from BAM header to prevent overflow errors
        try:
            chrom_len = bam.header.get_reference_length(chrom)
        except ValueError:
            # If chrom is not in BAM (e.g., "chr1" vs "1"), try stripping/adding "chr"
            if chrom.startswith("chr"):
                alt_chrom = chrom[3:]
            else:
                alt_chrom = "chr" + chrom
            
            try:
                chrom_len = bam.header.get_reference_length(alt_chrom)
                chrom = alt_chrom # Update chrom to match BAM
            except ValueError:
                # Chromosome missing from BAM entirely; return original point
                return ref_point, ref_point

        last_valid = ref_point
        
        for offset in range(step, max_extension + step, step):
            if strand == "+":
                region_start = ref_point
                # 2. CLAMP: Ensure end does not exceed chromosome length
                region_end = min(ref_point + offset, chrom_len)
            else:
                # 2. CLAMP: Ensure start is not less than 0
                region_start = max(0, ref_point - offset)
                region_end = ref_point
            
            # If our start/end are invalid after clamping (e.g. start > end), stop.
            if region_start >= region_end:
                break

            # 3. Safe count with try-except as final safety net
            try:
                count = bam.count(contig=chrom, start=region_start, end=region_end)
            except ValueError:
                break

            if count < count_threshold:
                break
            
            last_valid = region_end if strand == "+" else region_start
            
            # If we hit the end of the chromosome, stop the loop to avoid next iteration error
            if (strand == "+" and region_end == chrom_len) or (strand == "-" and region_start == 0):
                break

    return (ref_point, last_valid) if strand == "+" else (last_valid, ref_point)

def generate_saf(gtf_file, extension, output_dir, dynamic, bam_file_for_dynamic, species, id_type, kgx_file="kgXref.txt.gz", gap=0):
    os.makedirs(output_dir, exist_ok=True)
    # Include gap in filename for clarity
    output_saf = os.path.join(output_dir, f"3utr_ext{extension}bp_gap{gap}bp_extended_regions.saf")
    
    if species == "hg19":
        gtf = load_gtf_custom(gtf_file)
    else:
        gtf = pd.read_csv(gtf_file, sep="\t", comment="#", header=None,
                          names=["Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"])
    
    transcript = gtf[gtf["Feature"] == "transcript"].copy()
    
    if species in ["hg38", "hg19"]:
        transcript["GeneID"] = transcript["Attributes"].str.extract(r'gene_id\s+"([^"]+)"')
        transcript["TranscriptID"] = transcript["Attributes"].str.extract(r'transcript_id\s+"([^"]+)"')
        
        # --- ROBUST ID FIX ---
        transcript["GeneID"] = transcript["GeneID"].astype(str).str.split('.').str.get(0)
        transcript["TranscriptID"] = transcript["TranscriptID"].astype(str).str.split('.').str.get(0)
        # ---------------------
        
        # Added low_memory=False to prevent DtypeWarning
        kgxref = pd.read_csv(kgx_file, sep="\t", compression="gzip", header=None, low_memory=False,
                             names=["kgID", "mRNA", "spID", "spDisplayID", "geneSymbol", 
                                    "refseq", "protAcc", "description", "rfamAcc", "tRnaName"])
        
        # Ensure kgID is also cleaned of versions
        kgxref["kgID"] = kgxref["kgID"].astype(str).str.split('.').str.get(0)
        kgxref_dict = pd.Series(kgxref.geneSymbol.values, index=kgxref.kgID).to_dict()
        
        mg = mygene.MyGeneInfo()
        
        def map_to_official(row):
            gene_id = str(row["GeneID"])
            transcript_id = str(row["TranscriptID"])
            
            # Check dictionary first
            if gene_id in kgxref_dict:
                return kgxref_dict[gene_id]
                
            if gene_id == transcript_id:
                try:
                    with contextlib.redirect_stderr(open(os.devnull, "w")):
                        result = mg.query(transcript_id, scopes="ensembl.transcript", fields="symbol", species="human", verbose=False)
                    if result and result.get("hits"):
                        return result["hits"].get("symbol", gene_id)
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
        # Apply same robust fix for mouse/other species
        transcript["GeneID"] = transcript["GeneID"].astype(str).str.split('.').str.get(0)
        gene_col = "GeneID"

    if dynamic and bam_file_for_dynamic:
        def calc_dynamic(row):
            ref_point = int(row["End"]) if row["Strand"].strip() == "+" else int(row["Start"])
            strand = row["Strand"].strip()
            
            # 1. Normalize input to always be a list
            if isinstance(bam_file_for_dynamic, list):
                bams = bam_file_for_dynamic
            else:
                bams = [bam_file_for_dynamic]
            
            # 2. Initialize final coordinates to the reference point (i.e., extension = 0)
            final_start = ref_point
            final_end = ref_point
            
            # 3. Iterate through EVERY bam file in the list.
            # This loop structure prevents the "expected str, not list" error 
            # because 'bam_path' is guaranteed to be a single item from the list.
            for bam_path in bams:
                # Extra safety: cast to string to ensure it's a path
                safe_path = str(bam_path)
                
                # Calculate extension for this specific BAM
                s, e = dynamic_region_extension(safe_path, row["Chr"], ref_point, strand, max_extension=extension)
                
                # Update maximum extension logic
                if strand == "+":
                    # For positive strand, we want the largest End coordinate
                    if e > final_end:
                        final_end = e
                else:
                    # For negative strand, we want the smallest Start coordinate (furthest upstream)
                    # Note: We initialize final_start to ref_point. If s is smaller, it means we extended.
                    if s < final_start:
                        final_start = s
                        
            return pd.Series([final_start, final_end])
        
        print("Calculating dynamic extensions across all BAM files. This may take some time...")
        transcript[["Start_ext", "End_ext"]] = transcript.apply(calc_dynamic, axis=1)
    else:
        # Gap logic for fixed extension
        def calculate_fixed(row, extension, gap):
            if row["Strand"].strip() == "+":
                # Start = End + Gap
                return int(row["End"]) + gap, int(row["End"]) + gap + extension
            else:
                # Start = Start - Gap - Extension
                return max(int(row["Start"]) - gap - extension, 1), max(int(row["Start"]) - gap, 1)
        
        transcript[["Start_ext", "End_ext"]] = transcript.apply(lambda row: pd.Series(calculate_fixed(row, extension, gap)), axis=1)

    saf_df = transcript[[gene_col, "Chr", "Start_ext", "End_ext", "Strand"]].dropna()
    saf_df.columns = ["GeneID", "Chr", "Start", "End", "Strand"]
    saf_df.to_csv(output_saf, sep="\t", index=False)
    print(f"Extended SAF file saved to {output_saf}")
    
    return output_saf
