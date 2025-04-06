#!/usr/bin/env python3
import os
import pandas as pd
from scipy.stats import fisher_exact

def compare_enrichment_across_conditions(predictions_csv, output_dir, experimental_condition):
    """
    Compare relative enrichment across conditions using Fisher's exact test.
    This example compares the number of genes predicted as DE (predicted_DE == 1) in the experimental 
    group versus control.
    """
    df = pd.read_csv(predictions_csv)
    # Check that required columns exist
    if 'condition' not in df.columns or 'predicted_DE' not in df.columns:
        print("Missing 'condition' or 'predicted_DE' column; cannot perform enrichment analysis.")
        return None

    # Create a contingency table: rows are conditions, columns are predicted_DE (0 or 1)
    contingency = pd.crosstab(df['condition'], df['predicted_DE'])
    try:
        exp_row = contingency.loc[experimental_condition].values
        # Get the other condition (assumes exactly 2 unique conditions)
        control_condition = [cond for cond in contingency.index if cond != experimental_condition][0]
        ctrl_row = contingency.loc[control_condition].values
    except Exception as e:
        print("Error extracting contingency table:", e)
        return None

    # Assume column 1 corresponds to predicted_DE == 1 and column 0 to predicted_DE == 0.
    table = [[exp_row[1], exp_row[0]], [ctrl_row[1], ctrl_row[0]]]
    oddsratio, pvalue = fisher_exact(table)

    # Save the contingency table and test results to CSV
    enrichment_comparison_file = os.path.join(output_dir, "enrichment_comparison.csv")
    result_df = pd.DataFrame({
        'experimental_condition': [experimental_condition],
        'control_condition': [control_condition],
        'exp_DE': [exp_row[1]],
        'exp_nonDE': [exp_row[0]],
        'ctrl_DE': [ctrl_row[1]],
        'ctrl_nonDE': [ctrl_row[0]],
        'odds_ratio': [oddsratio],
        'p_value': [pvalue]
    })
    result_df.to_csv(enrichment_comparison_file, index=False)

    # Write a summary text file
    summary_file = os.path.join(output_dir, "enrichment_summary.txt")
    with open(summary_file, "w") as f:
        f.write(f"Fisher's exact test p-value: {pvalue}\n")
        f.write(f"Odds Ratio: {oddsratio}\n")
    print(f"Enrichment summary saved to: {summary_file}")
    return enrichment_comparison_file

def assess_replicate_consistency(norm_counts_df, predicted_genes, experimental_samples, output_dir):
    """
    Assess replicate consistency using the coefficient of variation (CV).
    """
    # Select only the rows corresponding to predicted genes and the columns of experimental samples
    sub_df = norm_counts_df.loc[predicted_genes, experimental_samples]
    cv_series = sub_df.std(axis=1) / sub_df.mean(axis=1)
    cv_df = pd.DataFrame({"Gene": cv_series.index, "CV": cv_series.values})
    output_file = os.path.join(output_dir, "replicate_consistency_CV.csv")
    cv_df.to_csv(output_file, index=False)
    print(f"Replicate consistency (CV) saved to: {output_file}")
    return output_file

def run_unsupervised_ml(predictions_file, norm_counts_file, output_dir, experimental_condition, conditions_list=None):
    """
    Run unsupervised ML analyses:
      1. If the predictions file does not include a 'condition' column, create it by comparing
         the 'exp_mean' and 'ctrl_mean' columns.
      2. Update predictions: add a 'predicted_DE' column (1 if the gene's condition equals the experimental condition, else 0).
      3. Compare enrichment across conditions.
      4. Assess replicate consistency.
      
    The optional parameter `conditions_list` (a list of conditions corresponding to the normalized_counts file columns)
    is used to identify experimental sample columns if the predictions file does not include a 'sample' column.
    """
    # Read predictions using the gene names as index.
    dott_df = pd.read_csv(predictions_file, index_col=0)
    
    # If "condition" is missing, create it based on comparing exp_mean and ctrl_mean.
    if "condition" not in dott_df.columns:
        print("Column 'condition' not found in predictions file. Creating it based on exp_mean and ctrl_mean...")
        dott_df["condition"] = dott_df.apply(
            lambda row: experimental_condition if row["exp_mean"] > row["ctrl_mean"] else "Control", axis=1
        )
    
    # Create/Update the predicted_DE column based on the condition.
    dott_df["predicted_DE"] = dott_df.apply(
        lambda row: 1 if row["condition"] == experimental_condition else 0, axis=1
    )
    updated_csv = os.path.join(output_dir, "dott_genes_with_conditions.csv")
    dott_df.to_csv(updated_csv)
    print("Updated DOTT predictions saved to:", updated_csv)
    
    # Run enrichment analysis.
    compare_enrichment_across_conditions(updated_csv, output_dir, experimental_condition)
    
    # Load normalized counts (ensuring gene names are used as the index)
    norm_counts_df = pd.read_csv(norm_counts_file, index_col=0)
    
    # Identify experimental sample columns.
    if "sample" in dott_df.columns:
        experimental_samples = dott_df[dott_df["condition"] == experimental_condition]["sample"].tolist()
    elif conditions_list is not None:
        # Use conditions_list (assumed to be in order as provided via --conditions) to select experimental samples.
        experimental_indices = [i for i, cond in enumerate(conditions_list) if cond == experimental_condition]
        experimental_samples = norm_counts_df.columns[experimental_indices].tolist()
    else:
        # Fallback: assume experimental samples are the second half of the columns.
        experimental_samples = norm_counts_df.columns[len(norm_counts_df.columns)//2:].tolist()
    print("Experimental samples identified:", experimental_samples)
    
    # Filter predictions for the experimental condition.
    dott_df_filtered = dott_df[dott_df["condition"] == experimental_condition]
    filtered_csv = os.path.join(output_dir, "dott_genes_with_conditions_filtered.csv")
    dott_df_filtered.to_csv(filtered_csv)
    predicted_genes = dott_df_filtered.index.tolist()
    
    assess_replicate_consistency(norm_counts_df, predicted_genes, experimental_samples, output_dir)
    
    return updated_csv, filtered_csv
