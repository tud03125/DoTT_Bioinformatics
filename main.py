#!/usr/bin/env python3
import argparse
import os
import subprocess

from annotation import generate_saf
from featureCounts import run_featurecounts, filter_counts_for_deseq2
from merging_genes import merge_gene_regions
from coordinate_extraction import extract_coordinates
from utils import parse_conditions

def main():
    parser = argparse.ArgumentParser(description="DoTT Bioinformatics Pipeline")
    parser.add_argument("--gtf-file", required=True, help="Path to the GTF annotation file")
    parser.add_argument("--bam-files", nargs="+", required=True, help="List of BAM file paths")
    parser.add_argument("--species", choices=["mm39", "hg38", "hg19"], required=True,
                        help="Species option (e.g., mm39 for mouse)")
    parser.add_argument("--id-type", choices=["Ensembl_ID", "Symbol"], default="Symbol",
                        help="Gene identifier type for annotations")
    parser.add_argument("--extension", type=int, default=10000, help="Fixed extension length in bases")
    parser.add_argument("--dynamic", action="store_true", help="Use dynamic region extension")
    parser.add_argument("--output-dir", required=True, help="Directory to store output files")
    parser.add_argument("--conditions", required=True,
                        help="Comma-separated list of condition labels for each BAM file")
    parser.add_argument("--kgx-file", default="kgXref.txt.gz",
                        help="Path to the kgXref mapping file (for human GTFs)")
    # Optional modules:
    parser.add_argument("--run_gsea", action="store_true",
                        help="Run GSEA pre-ranked list module")
    parser.add_argument("--unsupervised_ml", action="store_true",
                        help="Run unsupervised ML analysis module")
    parser.add_argument("--supervised_ml", action="store_true",
                        help="Run supervised ML analysis module")
    parser.add_argument("--experimental_condition",
                        help="Label for the experimental condition")
    # Updated ground truth parameter (formerly sim_tx_info)
    parser.add_argument("--ground_truth", help="Path to the ground truth file (CSV or TXT format)")
    # Optional bootstrapping arguments (if applicable)
    parser.add_argument("--bootstrap", action="store_true", help="Perform bootstrapping in DESeq2 analysis")
    parser.add_argument("--n_boot", type=int, default=100, help="Number of bootstrap iterations")
    parser.add_argument("--consensus_threshold", type=float, default=0.5, help="Consensus threshold for bootstrapping")
    
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    
    # Parse conditions (if needed for internal purposes)
    condition_list = parse_conditions(args.bam_files, args.conditions)
    
    # Step 1: Annotation - Generate SAF file
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
    
    # Step 2: FeatureCounts - Run featureCounts and filter counts for DESeq2
    counts_file = run_featurecounts(saf_file, args.extension, args.bam_files, args.output_dir)
    cleaned_counts = filter_counts_for_deseq2(counts_file, args.output_dir)
    
    # Step 3: Differential Expression Analysis - Call the DESeq2 R script
    r_cmd = [
        "Rscript",
        "dott_DESeq2.R",
        "--counts_file", cleaned_counts,
        "--output_dir", args.output_dir,
        "--conditions", args.conditions
    ]
    # Optional Bootstrapping option for DESeq2
    if args.bootstrap:
        r_cmd.extend([
            "--bootstrap", "TRUE",
            "--n_boot", str(args.n_boot),
            "--consensus_threshold", str(args.consensus_threshold)
        ])
    subprocess.run(r_cmd, check=True)
    
    # Define file paths produced by the DESeq2 script:
    deseq2_results_file = os.path.join(args.output_dir, "3_UTR_extended_differential_analysis_results.csv")
    sig_results_file = os.path.join(args.output_dir, "significant_extended_genes.csv")
    absolute_sig_results_file = os.path.join(args.output_dir, "absolute_significant_extended_genes_with_individual_means.csv")
    norm_counts_file = os.path.join(args.output_dir, "normalized_counts.csv")
    
    # Step 4: Merge gene regions and extract coordinates
    merged_annotation_file = os.path.join(args.output_dir, "merged_DoTT_annotation.bed")
    merge_gene_regions(absolute_sig_results_file, saf_file, merged_annotation_file)
    extract_coordinates(absolute_sig_results_file, saf_file, args.output_dir)
    
    # Optional Module: GSEA Pre-ranked List
    if args.run_gsea:
        r_gsea_cmd = [
            "Rscript",
            "dott_gsea_preranked_list.R",
            "--input_file",  absolute_sig_results_file,
            "--output_dir", args.output_dir
        ]
        subprocess.run(r_gsea_cmd, check=True)
        print("GSEA pre-ranked list created. Follow the printed instructions for GenePattern upload.")
    
    # Optional Module: Unsupervised ML Analysis
    if args.unsupervised_ml:
        if not args.experimental_condition:
            raise ValueError("For unsupervised ML, please provide --experimental_condition.")
        from Unsupervised_ML import run_unsupervised_ml
        run_unsupervised_ml(absolute_sig_results_file, norm_counts_file, args.output_dir, 
                            args.experimental_condition, conditions_list=condition_list)
    
    # Optional Module: Supervised ML Analysis
    if args.supervised_ml:
        if not args.ground_truth:
            raise ValueError("For supervised ML, please provide --ground_truth (ground truth file path).")
        if not args.experimental_condition:
            raise ValueError("For supervised ML, please provide --experimental_condition.")
        from Supervised_ML import train_ml_classifier_cv, evaluate_ml_performance, compare_to_ground_truth
        # Compare DESeq2 results to ground truth using the significant results file
        #gt_merged, gt_cm, gt_roc_auc = compare_to_ground_truth(args.ground_truth, sig_results_file, args.gtf_file, args.output_dir)
        gt_merged, gt_cm, gt_roc_auc = compare_to_ground_truth(args.ground_truth, deseq2_results_file, args.gtf_file, args.output_dir)
        print("Ground truth comparison complete. ROC AUC from ground truth comparison before ML classifier:", gt_roc_auc)

        # Train ML classifier using cross-validation on the DESeq2 results file
        #cv_scores, final_model, pred_df, holdout_results = train_ml_classifier_cv(args.ground_truth, sig_results_file, args.gtf_file, args.output_dir, clf_cutoff=0.5, cv=5, test_size=0.2)
        cv_scores, final_model, pred_df, holdout_results = train_ml_classifier_cv(args.ground_truth, deseq2_results_file, args.gtf_file, args.output_dir, clf_cutoff=0.5, cv=5, test_size=0.2)
        # Evaluate ML performance:
        merged_ml, cm_df, metrics_df = evaluate_ml_performance(pred_df, holdout_results, args.ground_truth, args.gtf_file, deseq2_results_file, args.output_dir)
        print("ML classifier ROC AUC:", metrics_df["ROC_AUC"].values[0])

    print("Pipeline finished successfully.")

if __name__ == "__main__":
    main()
