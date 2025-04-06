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
            # In a GTF file, the 9th field (index 8) contains the attributes.
            attr = fields[8]
            transcript_match = re.search(r'transcript_id "([^"]+)"', attr)
            gene_match = re.search(r'gene_name "([^"]+)"', attr)
            if transcript_match and gene_match:
                transcript_id = transcript_match.group(1)
                gene_name = gene_match.group(1)
                mapping[transcript_id] = gene_name
    return mapping

def compare_to_ground_truth(sim_tx_info_file, deseq_results_file, mapping_gtf_file, output_dir):
    """
    Compares DESeq2 results with simulation ground truth.
    Merges the ground truth with DESeq2 results on gene names,
    computes a confusion matrix, and saves performance metrics.
    """
    sim_tx = pd.read_csv(sim_tx_info_file, sep="\t")
    # Build transcript-to-gene mapping.
    tx2gene = load_tx2gene_mapping(mapping_gtf_file)
    sim_tx['gene'] = sim_tx['transcriptid'].map(tx2gene)
    sim_tx = sim_tx.dropna(subset=['gene'])
    sim_tx['ground_truth'] = sim_tx['DEstatus.2'].apply(lambda x: 1 if str(x).strip().upper() == "TRUE" else 0)
    
    deseq = pd.read_csv(deseq_results_file, index_col=0)
    merged = sim_tx.merge(deseq[['padj', 'log2FoldChange']], left_on='gene', right_index=True, how='inner')
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

def train_ml_classifier_cv(sim_tx_info_file, deseq_results_file, mapping_gtf_file, output_dir, clf_cutoff=0.5, cv=5, test_size=0.2):
    """
    Trains a supervised ML classifier using DESeq2 features and simulation ground truth.
    Returns cross-validation scores, the final model, and predictions.
    """
    # Load simulation ground truth and mapping.
    sim_tx = pd.read_csv(sim_tx_info_file, sep="\t")
    tx2gene = load_tx2gene_mapping(mapping_gtf_file)
    sim_tx['gene'] = sim_tx['transcriptid'].map(tx2gene)
    sim_tx = sim_tx.dropna(subset=['gene'])
    sim_tx['ground_truth'] = sim_tx['DEstatus.2'].apply(lambda x: 1 if str(x).strip().upper() == "TRUE" else 0)
    
    deseq = pd.read_csv(deseq_results_file, index_col=0)
    data = sim_tx.merge(deseq[['baseMean', 'log2FoldChange', 'pvalue', 'padj']], left_on='gene', right_index=True, how='inner')
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

def evaluate_ml_performance(pred_df, holdout_results, sim_tx_info_file, mapping_gtf_file, output_dir):
    """
    Evaluates ML classifier performance by merging predictions with ground truth,
    computing confusion matrices, ROC and PR curves, and saving performance metrics.
    """
    sim_tx = pd.read_csv(sim_tx_info_file, sep="\t")
    tx2gene = load_tx2gene_mapping(mapping_gtf_file)
    sim_tx['gene'] = sim_tx['transcriptid'].map(tx2gene)
    sim_tx = sim_tx.dropna(subset=['gene'])
    sim_tx['ground_truth'] = sim_tx['DEstatus.2'].apply(lambda x: 1 if str(x).strip().upper() == "TRUE" else 0)
    
    merged_ml = sim_tx[['gene', 'ground_truth']].merge(pred_df, on='gene', how='inner')
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
