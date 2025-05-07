# === METHYLATION FEATURE COMPARISON - MULTI-MOTIF VERSION ===

import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Set Working Directory ===
os.chdir("C:/set/your/path")
print(f"Working directory set to: {os.getcwd()}")

# === Feature Parsing ===
def parse_features_biopython(gbk_file, feature_types=["CDS", "tRNA", "rRNA"]):
    print(f"[INFO] Parsing GenBank file: {gbk_file}")
    features = []
    for record in SeqIO.parse(gbk_file, "genbank"):
        for f in record.features:
            if f.type in feature_types:
                label = f.qualifiers.get("gene", ["NA"])[0] if "gene" in f.qualifiers else f.qualifiers.get("product", ["NA"])[0]
                features.append({
                    "contig": record.id,
                    "start": int(f.location.start),
                    "end": int(f.location.end),
                    "strand": f.location.strand,
                    "type": f.type,
                    "label": label
                })
    print(f"[INFO] Extracted {len(features)} features from {gbk_file}")
    return pd.DataFrame(features)

# === Annotate Methylation Sites ===
def annotate_methylation(meth_df, features_df):
    meth_df['feature_type'] = 'intergenic'
    for _, feature in features_df.iterrows():
        in_feature = (
            (meth_df['Position'] >= feature['start']) &
            (meth_df['Position'] <= feature['end'])
        )
        meth_df.loc[in_feature, 'feature_type'] = feature['type']
    return meth_df

# === Feature Summary ===
def summarize_by_feature(df):
    return df.groupby('feature_type').agg(
        total_sites=('Position', 'count'),
        mean_percent_modified=('Percent_modified', lambda x: round(x.mean(), 3)),
        forward_strand=('Strand', lambda x: (x == '+').sum()),
        reverse_strand=('Strand', lambda x: (x == '-').sum())
    ).reset_index()

# === Annotation + Summary Workflow ===
def run_annotation_and_summary(meth_file, gbk_file, tag):
    meth_df = pd.read_csv(meth_file, sep='\t')
    meth_df = meth_df[meth_df['Total_coverage'] >= 5].copy()
    features_df = parse_features_biopython(gbk_file)
    annotated = annotate_methylation(meth_df.copy(), features_df)
    summary = summarize_by_feature(annotated)
    
    annotated.to_csv(f"outputs/annotated_sites/filtered_annotated_methylation_{tag}.csv", index=False)
    summary.to_csv(f"outputs/feature_summaries/feature_summary_{tag}.csv", index=False)
    
    return annotated, summary

# === Comparative Feature Summary ===
def compare_feature_summaries(hyb_summary, ont_summary):
    merged = pd.merge(hyb_summary, ont_summary, on='feature_type', how='outer', suffixes=('_Hyb', '_ONT'))
    merged.fillna(0, inplace=True)
    merged['motif_ratio_ONT_Hyb'] = merged['total_sites_ONT'] / merged['total_sites_Hyb'].replace(0, pd.NA)
    merged.to_csv("outputs/comparisons/comparative_feature_summary.csv", index=False)
    return merged

# === Multi-Motif Confusion Matrix Computation ===
def compute_confusion_matrix(ann_hyb, ann_ont, threshold=20):
    features = ['CDS', 'rRNA', 'tRNA']
    summary = []

    hyb_filtered = ann_hyb[ann_hyb['Percent_modified'] >= threshold]
    ont_filtered = ann_ont[ann_ont['Percent_modified'] >= threshold]

    all_motifs = set(hyb_filtered['motif']).union(set(ont_filtered['motif']))

    for feature in features:
        for motif in sorted(all_motifs):
            hyb_sites = set(hyb_filtered[(hyb_filtered['feature_type'] == feature) & (hyb_filtered['motif'] == motif)]['Position'])
            ont_sites = set(ont_filtered[(ont_filtered['feature_type'] == feature) & (ont_filtered['motif'] == motif)]['Position'])

            tp = len(hyb_sites & ont_sites)
            fp = len(ont_sites - hyb_sites)
            fn = len(hyb_sites - ont_sites)

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

            summary.append({
                'feature_type': feature,
                'motif': motif,
                'TP': tp,
                'FP': fp,
                'FN': fn,
                'Precision': round(precision, 3),
                'Recall': round(recall, 3),
                'F1_score': round(f1, 3)
            })

    confusion_df = pd.DataFrame(summary)
    confusion_df.to_csv("outputs/comparisons/confusion_matrix_by_motif_and_feature.csv", index=False)
    return confusion_df

# === Plot Comparative Summaries ===
def plot_comparison(merged):
    features = merged['feature_type']
    x = range(len(features))

    plt.figure(figsize=(12, 6))
    plt.bar([i - 0.2 for i in x], merged['total_sites_Hyb'], width=0.4, label='Hybrid', color='#1f77b4')
    plt.bar([i + 0.2 for i in x], merged['total_sites_ONT'], width=0.4, label='ONT', color='#ff7f0e')
    plt.xticks(x, features, rotation=45)
    plt.ylabel("GATC Site Count")
    plt.title("Motif Count by Feature Type")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/plots/motif_count_comparison.png", dpi=300)
    plt.close()

    plt.figure(figsize=(12, 6))
    plt.bar([i - 0.2 for i in x], merged['mean_percent_modified_Hyb'], width=0.4, label='Hybrid', color='#1f77b4')
    plt.bar([i + 0.2 for i in x], merged['mean_percent_modified_ONT'], width=0.4, label='ONT', color='#ff7f0e')
    plt.xticks(x, features, rotation=45)
    plt.ylabel("Mean % Modification")
    plt.title("Mean Modification by Feature Type")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/plots/mean_modification_comparison.png", dpi=300)
    plt.close()

# === Boxplot of % Modification ===
def plot_boxplot_by_feature(annotated_df, tag):
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=annotated_df, x='feature_type', y='Percent_modified', palette='Set3')
    plt.title(f"Distribution of Percent Modified by Feature Type ({tag})")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"outputs/plots/boxplot_percent_modified_by_feature_{tag}.png", dpi=300)
    plt.close()

# === RUN HERE ===
if __name__ == "__main__":
    # Create output folders
    os.makedirs("outputs/annotated_sites", exist_ok=True)
    os.makedirs("outputs/feature_summaries", exist_ok=True)
    os.makedirs("outputs/comparisons", exist_ok=True)
    os.makedirs("outputs/plots", exist_ok=True)

    # Input files
    ont_csv = "ONT_methylation_sites.tsv"
    ont_gbk = "ONT_genbank.gbk"
    hyb_csv = "Hybrid_methylation_sites.tsv"
    hyb_gbk = "Hybrid_genbank.gbk"

    # Annotation and Summaries
    ann_hyb, sum_hyb = run_annotation_and_summary(hyb_csv, hyb_gbk, "Hyb")
    ann_ont, sum_ont = run_annotation_and_summary(ont_csv, ont_gbk, "ONT")

    # Boxplots
    plot_boxplot_by_feature(ann_hyb, "Hyb")
    plot_boxplot_by_feature(ann_ont, "ONT")

    # Feature Summary Comparison
    merged_summary = compare_feature_summaries(sum_hyb, sum_ont)

    # Comparative Plots
    plot_comparison(merged_summary)

    # Confusion Matrix per Feature × Motif
    confusion_df = compute_confusion_matrix(ann_hyb, ann_ont, threshold=20)

    print("✅ Methylation annotation, multi-motif comparison, and confusion matrix complete. All outputs saved.")
