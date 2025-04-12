#!/usr/bin/env python3
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (confusion_matrix, roc_curve, auc, precision_recall_curve,
                             average_precision_score, f1_score, accuracy_score, precision_score)
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split, cross_val_score, cross_val_predict
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline


def load_tx2gene_mapping(gtf_file):
    """
    Loads a transcript-to-gene mapping from a GTF file.
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
    Loads ground truth data from a file that can be either a CSV or tab-delimited TXT.
    If a 'gene' column is missing, the function maps 'transcriptid' to gene names using the provided GTF file.
    It also standardizes the ground truth status.
    Returns a DataFrame with at least columns: 'gene' and 'ground_truth'.
    """
    ext = os.path.splitext(ground_truth_file)[1].lower()
    if ext == ".csv":
        gt_df = pd.read_csv(ground_truth_file)
    else:
        gt_df = pd.read_csv(ground_truth_file, sep="\t")
    
    if "gene" not in gt_df.columns:
        if "transcriptid" not in gt_df.columns:
            raise ValueError("Ground truth file must contain either a 'gene' or 'transcriptid' column.")
        tx2gene = load_tx2gene_mapping(mapping_gtf_file)
        gt_df["gene"] = gt_df["transcriptid"].map(tx2gene)
        gt_df = gt_df.dropna(subset=["gene"])
    
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
    Compares DESeq2 results to ground truth and generates performance curves and metrics.
    """
    gt = load_ground_truth(ground_truth_file, mapping_gtf_file)
    gt["gene"] = gt["gene"].astype(str)
    
    deseq = pd.read_csv(deseq_results_file, index_col=0)
    merged = gt.merge(deseq[['padj', 'log2FoldChange']], left_on='gene', right_index=True, how='inner')
    merged['predicted'] = ((merged['padj'] < 0.05) & (merged['log2FoldChange'] > 1)).astype(int)
    
    merged_file = os.path.join(output_dir, "merged_ground_truth_and_deseq_results.csv")
    merged.to_csv(merged_file, index=False)
    print(f"Merged ground truth and DESeq2 results saved to: {merged_file}")
    
    y_true = merged['ground_truth']
    y_pred = merged['predicted']
    cm = confusion_matrix(y_true, y_pred)
    cm_df = pd.DataFrame(cm, index=["Actual_NonDE", "Actual_DE"], columns=["Pred_NonDE", "Pred_DE"])
    cm_file = os.path.join(output_dir, "Original_confusion_matrix.csv")
    cm_df.to_csv(cm_file)
    print("Confusion Matrix:")
    print(cm_df)
    
    deseq['padj_adj'] = deseq['padj'].replace(0, 1e-300)
    deseq['score'] = -np.log10(deseq['padj_adj'])
    merged = gt.merge(deseq[['score']], left_on='gene', right_index=True, how='inner')
    merged = merged.dropna(subset=['score'])
    
    fpr, tpr, _ = roc_curve(merged['ground_truth'], merged['score'])
    roc_auc_val = auc(fpr, tpr)
    print(f"ROC AUC: {roc_auc_val:.2f}")
    
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve: Ground Truth vs. DESeq2 Predictions')
    plt.legend(loc="lower right")
    roc_plot_file = os.path.join(output_dir, "Original_ROC_curve.svg")
    plt.savefig(roc_plot_file)
    plt.close()
    print(f"ROC curve plot saved to: {roc_plot_file}")
    
    pr_precision, pr_recall, _ = precision_recall_curve(merged['ground_truth'], merged['score'])
    avg_precision = average_precision_score(merged['ground_truth'], merged['score'])
    print(f"Average Precision: {avg_precision:.2f}")
    
    plt.figure()
    plt.plot(pr_recall, pr_precision, label=f'PR curve (AP = {avg_precision:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve: Ground Truth vs. DESeq2 Predictions')
    plt.legend(loc="upper right")
    pr_plot_file = os.path.join(output_dir, "Original_PR_curve.svg")
    plt.savefig(pr_plot_file)
    plt.close()
    print(f"Precision-Recall curve plot saved to: {pr_plot_file}")
    
    try:
        tn, fp, fn, tp = cm.ravel()
    except ValueError:
        tn, fp, fn, tp = 0, 0, 0, 0
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float('nan')
    specificity = tn / (tn + fp) if (tn + fp) > 0 else float('nan')
    precision_val = tp / (tp + fp) if (tp + fp) > 0 else float('nan')
    accuracy = (tp + tn) / (tn + fp + fn + tp) if (tn + fp + fn + tp) > 0 else float('nan')
    f1 = 2 * precision_val * sensitivity / (precision_val + sensitivity) if (precision_val + sensitivity) > 0 else float('nan')
    
    metrics = {
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "Precision": precision_val,
        "Accuracy": accuracy,
        "F1 Score": f1,
        "ROC AUC": roc_auc_val,
        "Average Precision": avg_precision
    }
    
    print("Performance Metrics:")
    for metric, value in metrics.items():
        print(f"{metric}: {value:.3f}")
    
    metrics_df = pd.DataFrame([metrics])
    metrics_file = os.path.join(output_dir, "Original_performance_metrics.csv")
    metrics_df.to_csv(metrics_file, index=False)
    print(f"Original performance metrics saved to: {metrics_file}")
    
    return merged_file, cm_file, roc_auc_val


def train_ml_classifier_cv(ground_truth_file, deseq2_results_file, mapping_gtf_file, output_dir, clf_cutoff=0.5, cv=5, test_size=0.2):
    """
    Trains a supervised ML classifier using DESeq2 features and ground truth,
    and returns cross-validation scores, the final trained model, training predictions,
    holdout results, and cross-validation predictions.
    """
    gt_df = load_ground_truth(ground_truth_file, mapping_gtf_file)
    deseq = pd.read_csv(deseq2_results_file, index_col=0)
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
    
    # Get cross-validation predictions for additional evaluation
    cv_pred_prob = cross_val_predict(pipeline, features, labels, cv=skf, method='predict_proba')[:, 1]
    cv_pred = (cv_pred_prob >= clf_cutoff).astype(int)
    cv_results = {
        'cv_pred_prob': cv_pred_prob,
        'cv_pred': cv_pred,
        'labels': labels,
        'gene': gene_names  # Include gene information for merging
    }
    
    cv_df = pd.DataFrame({
        'gene': cv_results['gene'],  
        'cv_pred_prob': cv_results['cv_pred_prob'],
        'cv_pred': cv_results['cv_pred'],
        'true_label': cv_results['labels']
    })
    
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
    
    holdout_df = pd.DataFrame({
        'gene': holdout_results['gene_holdout'],   # Ensure holdout_results includes gene names
        'predicted_prob': holdout_results['predicted_prob'],
        'predicted_label': holdout_results['predicted_label'],
        'true_label': holdout_results['y_holdout']
    })
    
    return cv_scores, pipeline, pred_df, holdout_results, cv_results


def evaluate_ml_performance(pred_df, holdout_results, ground_truth_file, mapping_gtf_file, deseq2_results_file, output_dir, cv_results=None):
    """
    Evaluates the ML classifier performance by merging ML predictions with the ground truth and then generating 
    confusion matrices, ROC curves, and precision-recall curves for the training set, holdout set, and cross-validation (if provided).
    """
    ## Evaluate on the training (merged) predictions
    gt_df = load_ground_truth(ground_truth_file, mapping_gtf_file)
    gt_df["gene"] = gt_df["gene"].astype(str)
    pred_df["gene"] = pred_df["gene"].astype(str)
    
    ##############################
    ## 1. Merge Training Set (pred_df)
    ##############################
    merged_ml = gt_df[['gene', 'ground_truth']].merge(pred_df, on='gene', how='right')
    merged_ml['ground_truth'] = merged_ml['ground_truth'].fillna(0)
    merged_ml['padj'] = merged_ml['padj'].fillna(1)
    merged_ml = merged_ml.dropna(subset=['log2FoldChange'])
    merged_ml_file = os.path.join(output_dir, "merged_ground_truth_and_ml_predictions.csv")
    merged_ml.to_csv(merged_ml_file, index=False)
    print("Merged ML predictions with ground truth saved to:", merged_ml_file)
    
    ##############################
    ## 2. Merge Holdout Results
    ##############################
    # Create a DataFrame for holdout predictions. Make sure that holdout_results contains gene names.
    holdout_df = pd.DataFrame({
        'gene': holdout_results['gene_holdout'],   # Ensure your train_split returns gene_holdout properly.
        'predicted_prob': holdout_results['predicted_prob'],
        'predicted_label': holdout_results['predicted_label'],
        'true_label': holdout_results['y_holdout']
    })
    # Merge with ground truth using the gene column.
    merged_holdout = gt_df[['gene', 'ground_truth']].merge(holdout_df, on='gene', how='right')
    merged_holdout['ground_truth'] = merged_holdout['ground_truth'].fillna(0)
    holdout_merged_file = os.path.join(output_dir, "merged_ground_truth_and_holdout_predictions.csv")
    merged_holdout.to_csv(holdout_merged_file, index=False)
    print("Merged holdout predictions with ground truth saved to:", holdout_merged_file)
    
    ##############################
    ## 3. Merge Cross-Validation Results (Optional)
    ##############################
    if cv_results is not None:
        # Here we assume that cv_results was updated in train_ml_classifier_cv() to include gene names.
        cv_df = pd.DataFrame({
            'gene': cv_results['gene'],             # Requires that cv_results includes the key 'gene'
            'cv_pred_prob': cv_results['cv_pred_prob'],
            'cv_pred': cv_results['cv_pred'],
            'true_label': cv_results['labels']
        })
        merged_cv = gt_df[['gene', 'ground_truth']].merge(cv_df, on='gene', how='right')
        merged_cv['ground_truth'] = merged_cv['ground_truth'].fillna(0)
        cv_merged_file = os.path.join(output_dir, "merged_ground_truth_and_cv_predictions.csv")
        merged_cv.to_csv(cv_merged_file, index=False)
        print("Merged CV predictions with ground truth saved to:", cv_merged_file)
    
    # Training set confusion matrix
    y_true_train = merged_ml['ground_truth']
    y_pred_train = merged_ml['predicted_DE']
    cm = confusion_matrix(y_true_train, y_pred_train)
    cm_df = pd.DataFrame(cm, index=["Actual_NonDE", "Actual_DE"], columns=["Pred_NonDE", "Pred_DE"])
    cm_file = os.path.join(output_dir, "ml_confusion_matrix.csv")
    cm_df.to_csv(cm_file)
    print("ML Confusion Matrix saved to:", cm_file)
    
    # ROC and PR curves for training
    fpr, tpr, _ = roc_curve(y_true_train, merged_ml['predicted_score'])
    roc_auc_val = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, label=f'ROC curve (AUC = {roc_auc_val:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ML ROC Curve (Training Set)')
    plt.legend(loc="lower right")
    roc_plot_file = os.path.join(output_dir, "ml_ROC_curve.svg")
    plt.savefig(roc_plot_file)
    plt.close()
    print("ML ROC curve saved to:", roc_plot_file)
    
    pr_precision, pr_recall, _ = precision_recall_curve(y_true_train, merged_ml['predicted_score'])
    avg_precision = average_precision_score(y_true_train, merged_ml['predicted_score'])
    plt.figure()
    plt.plot(pr_recall, pr_precision, label=f'PR curve (AP = {avg_precision:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('ML Precision-Recall Curve (Training Set)')
    plt.legend(loc="upper right")
    pr_plot_file = os.path.join(output_dir, "ml_PR_curve.svg")
    plt.savefig(pr_plot_file)
    plt.close()
    print("ML PR curve saved to:", pr_plot_file)
    
    try:
        tn, fp, fn, tp = cm.ravel()
    except:
        tn, fp, fn, tp = 0, 0, 0, 0
    sensitivity = tp/(tp+fn) if (tp+fn) > 0 else float('nan')
    specificity = tn/(tn+fp) if (tn+fp) > 0 else float('nan')
    accuracy = accuracy_score(y_true_train, y_pred_train)
    precision_val = precision_score(y_true_train, y_pred_train)
    f1 = f1_score(y_true_train, y_pred_train)
    
    training_metrics = {
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "Accuracy": accuracy,
        "Precision": precision_val,
        "F1_Score": f1,
        "ROC_AUC": roc_auc_val,
        "Average_Precision": avg_precision
    }
    training_metrics_df = pd.DataFrame([training_metrics])
    training_metrics_file = os.path.join(output_dir, "ml_training_performance_metrics.csv")
    training_metrics_df.to_csv(training_metrics_file, index=False)
    print("ML training performance metrics saved to:", training_metrics_file)
    
    ## Evaluate on the holdout set
    y_holdout = holdout_results['y_holdout']
    holdout_pred_prob = holdout_results['predicted_prob']
    holdout_pred = holdout_results['predicted_label']
    cm_holdout = confusion_matrix(y_holdout, holdout_pred)
    cm_holdout_df = pd.DataFrame(cm_holdout, index=["Actual_NonDE", "Actual_DE"], columns=["Pred_NonDE", "Pred_DE"])
    cm_holdout_file = os.path.join(output_dir, "holdout_confusion_matrix.csv")
    cm_holdout_df.to_csv(cm_holdout_file)
    print("Holdout Confusion Matrix saved to:", cm_holdout_file)
    
    fpr_h, tpr_h, _ = roc_curve(y_holdout, holdout_pred_prob)
    roc_auc_h = auc(fpr_h, tpr_h)
    plt.figure()
    plt.plot(fpr_h, tpr_h, label=f'Holdout ROC (AUC = {roc_auc_h:.2f})')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Holdout ROC Curve')
    plt.legend(loc="lower right")
    roc_holdout_file = os.path.join(output_dir, "holdout_ROC_curve.svg")
    plt.savefig(roc_holdout_file)
    plt.close()
    print("Holdout ROC curve saved to:", roc_holdout_file)
    
    pr_precision_h, pr_recall_h, _ = precision_recall_curve(y_holdout, holdout_pred_prob)
    avg_precision_h = average_precision_score(y_holdout, holdout_pred_prob)
    plt.figure()
    plt.plot(pr_recall_h, pr_precision_h, label=f'Holdout PR (AP = {avg_precision_h:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Holdout Precision-Recall Curve')
    plt.legend(loc="upper right")
    pr_holdout_file = os.path.join(output_dir, "holdout_PR_curve.svg")
    plt.savefig(pr_holdout_file)
    plt.close()
    print("Holdout PR curve saved to:", pr_holdout_file)
    
    try:
        tn_h, fp_h, fn_h, tp_h = cm_holdout.ravel()
    except:
        tn_h, fp_h, fn_h, tp_h = 0, 0, 0, 0
    holdout_sensitivity = tp_h/(tp_h+fn_h) if (tp_h+fn_h) > 0 else float('nan')
    holdout_specificity = tn_h/(tn_h+fp_h) if (tn_h+fp_h) > 0 else float('nan')
    holdout_accuracy = accuracy_score(y_holdout, holdout_pred)
    holdout_precision = precision_score(y_holdout, holdout_pred)
    holdout_f1 = f1_score(y_holdout, holdout_pred)
    
    holdout_metrics = {
        "Sensitivity": holdout_sensitivity,
        "Specificity": holdout_specificity,
        "Accuracy": holdout_accuracy,
        "Precision": holdout_precision,
        "F1_Score": holdout_f1,
        "ROC_AUC": roc_auc_h,
        "Average_Precision": avg_precision_h
    }
    holdout_metrics_df = pd.DataFrame([holdout_metrics])
    holdout_metrics_file = os.path.join(output_dir, "holdout_performance_metrics.csv")
    holdout_metrics_df.to_csv(holdout_metrics_file, index=False)
    print("Holdout performance metrics saved to:", holdout_metrics_file)
    
    ## Optionally, evaluate cross-validation performance if cv_results is provided
    if cv_results is not None:
        cv_labels = cv_results['labels']
        cv_pred_prob = cv_results['cv_pred_prob']
        cv_pred = cv_results['cv_pred']
        cv_cm = confusion_matrix(cv_labels, cv_pred)
        cv_cm_df = pd.DataFrame(cv_cm, index=["Actual_NonDE", "Actual_DE"], columns=["Pred_NonDE", "Pred_DE"])
        cv_cm_file = os.path.join(output_dir, "cv_confusion_matrix.csv")
        cv_cm_df.to_csv(cv_cm_file)
        print("CV Confusion Matrix saved to:", cv_cm_file)
        
        fpr_cv, tpr_cv, _ = roc_curve(cv_labels, cv_pred_prob)
        roc_auc_cv = auc(fpr_cv, tpr_cv)
        plt.figure()
        plt.plot(fpr_cv, tpr_cv, label=f'CV ROC (AUC = {roc_auc_cv:.2f})')
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('CV ROC Curve')
        plt.legend(loc="lower right")
        roc_cv_file = os.path.join(output_dir, "cv_ROC_curve.svg")
        plt.savefig(roc_cv_file)
        plt.close()
        print("CV ROC curve saved to:", roc_cv_file)
        
        pr_precision_cv, pr_recall_cv, _ = precision_recall_curve(cv_labels, cv_pred_prob)
        avg_precision_cv = average_precision_score(cv_labels, cv_pred_prob)
        plt.figure()
        plt.plot(pr_recall_cv, pr_precision_cv, label=f'CV PR (AP = {avg_precision_cv:.2f})')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('CV Precision-Recall Curve')
        plt.legend(loc="upper right")
        pr_cv_file = os.path.join(output_dir, "cv_PR_curve.svg")
        plt.savefig(pr_cv_file)
        plt.close()
        print("CV PR curve saved to:", pr_cv_file)
        
        try:
            tn_cv, fp_cv, fn_cv, tp_cv = cv_cm.ravel()
        except:
            tn_cv, fp_cv, fn_cv, tp_cv = 0, 0, 0, 0
        cv_sensitivity = tp_cv/(tp_cv+fn_cv) if (tp_cv+fn_cv) > 0 else float('nan')
        cv_specificity = tn_cv/(tn_cv+fp_cv) if (tn_cv+fp_cv) > 0 else float('nan')
        cv_accuracy = accuracy_score(cv_labels, cv_pred)
        cv_precision_val = precision_score(cv_labels, cv_pred)
        cv_f1 = f1_score(cv_labels, cv_pred)
        cv_metrics = {
            "Sensitivity": cv_sensitivity,
            "Specificity": cv_specificity,
            "Accuracy": cv_accuracy,
            "Precision": cv_precision_val,
            "F1_Score": cv_f1,
            "ROC_AUC": roc_auc_cv,
            "Average_Precision": avg_precision_cv
        }
        cv_metrics_df = pd.DataFrame([cv_metrics])
        cv_metrics_file = os.path.join(output_dir, "cv_performance_metrics.csv")
        cv_metrics_df.to_csv(cv_metrics_file, index=False)
        print("CV performance metrics saved to:", cv_metrics_file)
    
    return merged_ml, merged_holdout, merged_cv, cm_df, training_metrics_df
