#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import re
from sklearn.metrics import (confusion_matrix, roc_curve, auc, precision_recall_curve,
                             average_precision_score, f1_score, accuracy_score, roc_auc_score)
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split, cross_val_score
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline
import matplotlib.pyplot as plt

def load_tx2gene_mapping(gtf_file):
    """
    Loads a transcript-to-gene mapping from a GTF file.
    Only considers lines where the feature type is 'transcript'.
    Returns a dictionary: {transcript_id: gene_name}
    """
    mapping = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            # In a GTF file, the 9th field contains the attributes.
            attr = fields[8]
            transcript_match = re.search(r'transcript_id "([^"]+)"', attr)
            gene_match = re.search(r'gene_name "([^"]+)"', attr)
            if transcript_match and gene_match:
                transcript_id = transcript_match.group(1)
                gene_name = gene_match.group(1)
                mapping[transcript_id] = gene_name
    return mapping

def load_ground_truth(ground_truth_file, mapping_gtf_file):
    """
    Loads ground truth data from a file that can be either a tab-delimited TXT file
    (with transcript IDs) or a CSV (with gene names already provided). If the 'gene'
    column is missing, the function will attempt to map 'transcriptid' to gene names
    using the supplied GTF file.
    
    The function also standardizes the ground truth status by looking for a column
    named 'DEstatus.2' (or 'ground_truth') and converting its values to 1 (TRUE) or 0 (FALSE).

    Returns:
        DataFrame with at least the following columns:
            - gene : UCSC-based gene name
            - ground_truth : 1 for TRUE (DoTT) and 0 for FALSE
    """
    ext = os.path.splitext(ground_truth_file)[1].lower()
    if ext == ".csv":
        gt_df = pd.read_csv(ground_truth_file)
    else:
        gt_df = pd.read_csv(ground_truth_file, sep="\t")
    
    # If a 'gene' column is not present, then map 'transcriptid' to gene names.
    if "gene" not in gt_df.columns:
        if "transcriptid" not in gt_df.columns:
            raise ValueError("Ground truth file must contain either a 'gene' or 'transcriptid' column.")
        tx2gene = load_tx2gene_mapping(mapping_gtf_file)
        gt_df["gene"] = gt_df["transcriptid"].map(tx2gene)
        gt_df = gt_df.dropna(subset=["gene"])
    
    # Standardize the ground truth column.
    # The original column (from the simulation) is 'DEstatus.2', but if you've preprocessed your data,
    # it might already be named 'ground_truth'.
    # Check for 'DEstatus.2', 'DEstatus', or 'ground_truth'
    if "DEstatus.2" in gt_df.columns:
        gt_df["ground_truth"] = gt_df["DEstatus.2"].apply(lambda x: 1 if str(x).strip().upper() == "TRUE" else 0)
    elif "DEstatus" in gt_df.columns:
        gt_df["ground_truth"] = gt_df["DEstatus"].apply(lambda x: 1 if str(x).strip().upper() == "TRUE" else 0)
    elif "ground_truth" in gt_df.columns:
        gt_df["ground_truth"] = gt_df["ground_truth"].apply(lambda x: 1 if str(x).strip().upper() == "TRUE" else 0)
    else:
        raise ValueError("Ground truth file must contain 'DEstatus.2', 'DEstatus', or 'ground_truth' column.")
    
    return gt_df

def compare_to_ground_truth(ground_truth_file, deseq_results_file, mapping_gtf_file, output_dir):
    """
    Compares DESeq2 results with ground truth.
    Merges the ground truth data with DESeq2 results (merged on gene names),
    computes a confusion matrix, and saves performance metrics.
    """
    gt_df = load_ground_truth(ground_truth_file, mapping_gtf_file)
    deseq = pd.read_csv(deseq_results_file, index_col=0)
    merged = gt_df.merge(deseq[['padj', 'log2FoldChange']], left_on='gene', right_index=True, how='inner')
    merged['predicted'] = ((merged['padj'] < 0.05) & (abs(merged['log2FoldChange']) > 1)).astype(int)
    merged_file = os.path.join(output_dir, "merged_ground_truth_and_deseq_results.csv")
    merged.to_csv(merged_file, index=False)
    print("Merged ground truth and DESeq2 results saved to:", merged_file)
    
    y_true = merged['ground_truth']
    y_pred = merged['predicted']
    cm = confusion_matrix(y_true, y_pred)
    cm_df = pd.DataFrame(cm, index=["Actual_NonDE", "Actual_DE"], columns=["Pred_NonDE", "Pred_DE"])
    cm_file = os.path.join(output_dir, "confusion_matrix.csv")
    cm_df.to_csv(cm_file, index=False)
    print("Confusion Matrix saved to:", cm_file)
    
    fpr, tpr, _ = roc_curve(y_true, merged['padj'])
    roc_auc_val = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot([0,1], [0,1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve: Ground Truth vs. DESeq2')
    plt.legend(loc="lower right")
    roc_plot_file = os.path.join(output_dir, "ROC_curve.svg")
    plt.savefig(roc_plot_file)
    plt.close()
    print("ROC curve saved to:", roc_plot_file)
    
    return merged, cm_df, roc_auc_val

def train_ml_classifier_cv(ground_truth_file, deseq_results_file, mapping_gtf_file, output_dir, clf_cutoff=0.5, cv=5, test_size=0.2):
    """
    Trains a supervised ML classifier using DESeq2 features and ground truth.
    Merges ground truth (loaded via load_ground_truth) with DESeq2 result data
    (including baseMean, log2FoldChange, pvalue, and padj), and then computes a score.
    
    Returns cross-validation scores, the final model, and prediction DataFrames.
    """
    gt_df = load_ground_truth(ground_truth_file, mapping_gtf_file)
    deseq = pd.read_csv(deseq_results_file, index_col=0)
    data = gt_df.merge(deseq[['baseMean', 'log2FoldChange', 'pvalue', 'padj']], left_on='gene', right_index=True, how='inner')
    data.dropna(subset=['padj'], inplace=True)
    data['padj_adj'] = data['padj'].replace(0, 1e-300)
    data['score'] = -np.log10(data['padj_adj'])
    
    gene_names = data['gene']
    numeric_cols = ['baseMean', 'log2FoldChange', 'pvalue', 'padj', 'score']
    features = data[numeric_cols].copy()
    valid_indices = features.dropna().index
    features = features.loc[valid_indices]
    labels = data.loc[valid_indices, 'ground_truth']
    gene_names = gene_names.loc[valid_indices]
    
    pipeline = Pipeline([
        ('smote', SMOTE(random_state=42)),
        ('clf', RandomForestClassifier(n_estimators=100, random_state=42))
    ])
    
    skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=42)
    cv_scores = cross_val_score(pipeline, features, labels, cv=skf, scoring='roc_auc')
    print("Cross-validation ROC AUC scores:", cv_scores)
    print("Mean ROC AUC: {:.3f}".format(cv_scores.mean()))
    
    X_train, X_holdout, y_train, y_holdout, gene_train, gene_holdout = train_test_split(
        features, labels, gene_names, test_size=test_size, random_state=42, stratify=labels
    )
    pipeline.fit(X_train, y_train)
    pred_prob_train = pipeline.predict_proba(X_train)[:, 1]
    
    pred_df = pd.DataFrame({
        'gene': gene_train,
        'baseMean': X_train['baseMean'],
        'log2FoldChange': X_train['log2FoldChange'],
        'pvalue': X_train['pvalue'],
        'padj': X_train['padj'],
        'predicted_score': pred_prob_train,
        'true_label': y_train
    })
    pred_df['predicted_DE'] = (pred_df['predicted_score'] >= clf_cutoff).astype(int)
    pred_file = os.path.join(output_dir, "ml_predictions_with_genes.csv")
    pred_df.to_csv(pred_file, index=False)
    print("ML predictions saved to:", pred_file)
    
    holdout_pred_prob = pipeline.predict_proba(X_holdout)[:, 1]
    holdout_pred = (holdout_pred_prob >= clf_cutoff).astype(int)
    holdout_results = {
        'gene_holdout': gene_holdout,
        'predicted_prob': holdout_pred_prob,
        'predicted_label': holdout_pred,
        'y_holdout': y_holdout
    }
    
    return cv_scores, pipeline, pred_df, holdout_results

def evaluate_ml_performance(pred_df, holdout_results, ground_truth_file, mapping_gtf_file, output_dir):
    """
    Evaluates the ML classifier performance by merging ML predictions with ground truth,
    and generating confusion matrices, ROC curves, and precision-recall curves.
    Saves all performance metrics to CSV files.
    """
    gt_df = load_ground_truth(ground_truth_file, mapping_gtf_file)
    merged_ml = gt_df[['gene', 'ground_truth']].merge(pred_df, on='gene', how='inner')
    merged_ml_file = os.path.join(output_dir, "merged_ground_truth_and_ml_predictions.csv")
    merged_ml.to_csv(merged_ml_file, index=False)
    print("Merged ML predictions with ground truth saved to:", merged_ml_file)
    
    y_true = merged_ml['ground_truth']
    y_pred = merged_ml['predicted_DE']
    cm = confusion_matrix(y_true, y_pred)
    cm_df = pd.DataFrame(cm, index=["Actual_NonDE", "Actual_DE"], columns=["Pred_NonDE", "Pred_DE"])
    cm_file = os.path.join(output_dir, "ml_confusion_matrix.csv")
    cm_df.to_csv(cm_file, index=False)
    print("ML Confusion Matrix saved to:", cm_file)
    
    fpr, tpr, _ = roc_curve(y_true, merged_ml['predicted_score'])
    roc_auc_val = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot([0,1], [0,1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ML ROC Curve')
    plt.legend(loc="lower right")
    roc_plot_file = os.path.join(output_dir, "ml_ROC_curve.svg")
    plt.savefig(roc_plot_file)
    plt.close()
    print("ML ROC curve saved to:", roc_plot_file)
    
    pr_precision, pr_recall, _ = precision_recall_curve(y_true, merged_ml['predicted_score'])
    avg_precision = average_precision_score(y_true, merged_ml['predicted_score'])
    plt.figure()
    plt.plot(pr_recall, pr_precision, label=f'PR curve (AP = {avg_precision:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('ML Precision-Recall Curve')
    plt.legend(loc="upper right")
    pr_plot_file = os.path.join(output_dir, "ml_PR_curve.svg")
    plt.savefig(pr_plot_file)
    plt.close()
    print("ML PR curve saved to:", pr_plot_file)
    
    try:
        tn, fp, fn, tp = cm.ravel()
    except Exception:
        tn, fp, fn, tp = 0, 0, 0, 0
    sensitivity = tp/(tp+fn) if (tp+fn) > 0 else float('nan')
    specificity = tn/(tn+fp) if (tn+fp) > 0 else float('nan')
    precision_val = tp/(tp+fp) if (tp+fp) > 0 else float('nan')
    accuracy = (tp+tn)/(tp+tn+fp+fn) if (tp+tn+fp+fn) > 0 else float('nan')
    f1 = 2 * precision_val * sensitivity / (precision_val + sensitivity) if (precision_val + sensitivity) > 0 else float('nan')
    
    metrics_df = pd.DataFrame([{
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "Precision": precision_val,
        "Accuracy": accuracy,
        "F1_Score": f1,
        "ROC_AUC": roc_auc_val,
        "Average_Precision": avg_precision
    }])
    metrics_file = os.path.join(output_dir, "ml_performance_metrics.csv")
    metrics_df.to_csv(metrics_file, index=False)
    print("ML performance metrics saved to:", metrics_file)
    
    return merged_ml, cm_df, metrics_df
