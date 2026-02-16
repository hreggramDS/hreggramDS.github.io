---
layout: pipeline
title: "scRNA-seq Analysis Pipeline"
date: 2025-01-20
version: "2.0"
domain: "Bioinformatics"
technologies: ["Python", "Scanpy", "AnnData", "Cell Ranger", "NumPy", "Matplotlib"]
repository: "https://github.com/hreggramDS/scrna-pipeline"
description: "End-to-end single-cell RNA sequencing analysis pipeline from raw FASTQ files to biological interpretation"
---

## Overview

A comprehensive, reproducible pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data. This workflow covers every step from raw sequencing reads to cell type annotation, differential expression, and trajectory inference — designed for 10x Genomics Chromium data but adaptable to other platforms.

---

## Pipeline Architecture

```
FASTQ files
    │
    ▼
┌─────────────────┐
│  Cell Ranger     │  Alignment, barcode counting, UMI deduplication
│  (or STARsolo)   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Quality Control │  Cell/gene filtering, doublet detection, QC metrics
│                  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Normalization   │  Library size normalization, log transform, HVG selection
│  & Feature Sel.  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Dim. Reduction  │  PCA → Batch correction (optional) → UMAP/t-SNE
│                  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Clustering      │  KNN graph → Leiden/Louvain community detection
│                  │
└────────┬────────┘
         │
         ├──────────────────────────┐
         ▼                          ▼
┌─────────────────┐     ┌───────────────────────┐
│  Marker Genes &  │     │  Trajectory Inference  │
│  Cell Type Anno. │     │  (PAGA, Diffusion Map) │
└────────┬────────┘     └───────────┬───────────┘
         │                          │
         └──────────┬───────────────┘
                    ▼
            ┌──────────────┐
            │  Reporting &  │  Figures, tables, annotated AnnData object
            │  Export       │
            └──────────────┘
```

---

## Step 1: Raw Data Processing (Cell Ranger)

### Input Requirements

| Input | Format | Description |
|---|---|---|
| Sequencing reads | FASTQ | R1 (barcode+UMI), R2 (cDNA), I1 (sample index) |
| Reference genome | FASTA + GTF | Pre-built or custom genome reference |
| Chemistry | Auto-detect | 10x Chromium v2, v3, or v3.1 |

### Running Cell Ranger

```bash
cellranger count \
    --id=sample_01 \
    --transcriptome=/ref/refdata-gex-GRCh38-2024-A \
    --fastqs=/data/fastqs/sample_01/ \
    --sample=Sample_01 \
    --localcores=16 \
    --localmem=64 \
    --expect-cells=10000
```

### Key Outputs

- `filtered_feature_bc_matrix.h5` — Count matrix (cells × genes)
- `web_summary.html` — QC report with estimated cell count, sequencing saturation, mapping rates
- `molecule_info.h5` — Per-molecule data for secondary analyses

---

## Step 2: Quality Control

### Loading and Initial QC

```python
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

# Load count matrix
adata = sc.read_10x_h5('sample_01/outs/filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()

# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt', 'ribo'], 
    percent_top=None, log1p=False, inplace=True
)

# Plot QC distributions
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
sc.pl.violin(adata, 'pct_counts_ribo', ax=axes[3], show=False)
plt.tight_layout()
plt.savefig('qc_violins.png', dpi=150)
```

### Filtering Criteria

```python
# Apply QC filters
print(f"Cells before filtering: {adata.n_obs}")

# Cell-level filters
sc.pp.filter_cells(adata, min_genes=200)       # Remove empty droplets
sc.pp.filter_cells(adata, max_genes=8000)       # Remove potential doublets
adata = adata[adata.obs.pct_counts_mt < 20, :]  # Remove dying cells

# Gene-level filters
sc.pp.filter_genes(adata, min_cells=3)           # Remove rarely detected genes

print(f"Cells after filtering: {adata.n_obs}")
print(f"Genes after filtering: {adata.n_vars}")
```

### Doublet Detection

```python
import scrublet as scr

# Run Scrublet for doublet detection
scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2, min_cells=3,
    min_gene_variability_pctl=85,
    n_prin_comps=30
)

adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = predicted_doublets

# Remove predicted doublets
adata = adata[~adata.obs['predicted_doublet'], :]
print(f"Cells after doublet removal: {adata.n_obs}")
```

### QC Thresholds Summary

| Metric | Threshold | Rationale |
|---|---|---|
| Min genes/cell | 200 | Remove empty droplets |
| Max genes/cell | 8,000 | Remove likely doublets |
| Max % mitochondrial | 20% | Remove dying/stressed cells |
| Min cells/gene | 3 | Remove noise genes |
| Doublet score | Scrublet auto | Remove computational doublets |

---

## Step 3: Normalization & Feature Selection

```python
# Save raw counts for DE analysis later
adata.layers['counts'] = adata.X.copy()

# Library size normalization
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform
sc.pp.log1p(adata)

# Store normalized data
adata.raw = adata.copy()

# Highly variable gene selection
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    n_top_genes=3000  # Alternative: select top N
)

print(f"Highly variable genes: {adata.var['highly_variable'].sum()}")

# Scale to unit variance (clipped)
sc.pp.scale(adata, max_value=10)
```

---

## Step 4: Dimensionality Reduction

### PCA

```python
# PCA on highly variable genes
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# Determine number of PCs to use
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
# Typically use 30-50 PCs based on elbow
n_pcs = 40
```

### Batch Correction (if applicable)

```python
# If samples have batch effects, apply Harmony or scVI
# Option 1: Harmony
import harmonypy as hm

ho = hm.run_harmony(
    adata.obsm['X_pca'][:, :n_pcs],
    adata.obs,
    'batch'
)
adata.obsm['X_pca_harmony'] = ho.Z_corr
use_rep = 'X_pca_harmony'

# Option 2: scVI (for more complex batch structures)
# import scvi
# scvi.model.SCVI.setup_anndata(adata, batch_key='batch', layer='counts')
# model = scvi.model.SCVI(adata)
# model.train()
# adata.obsm['X_scVI'] = model.get_latent_representation()
# use_rep = 'X_scVI'
```

### Neighborhood Graph & UMAP

```python
# Build KNN graph
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, use_rep='X_pca')

# UMAP embedding for visualization
sc.tl.umap(adata, min_dist=0.3)

# t-SNE as alternative
sc.tl.tsne(adata, n_pcs=n_pcs)
```

---

## Step 5: Clustering

```python
# Leiden clustering (preferred over Louvain)
# Try multiple resolutions
for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')

# Visualize clusters at different resolutions
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, res in zip(axes, [0.3, 0.8, 1.5]):
    sc.pl.umap(adata, color=f'leiden_{res}', ax=ax, 
               title=f'Resolution {res}', show=False)
plt.tight_layout()
plt.savefig('clustering_resolutions.png', dpi=150)

# Select optimal resolution (typically 0.5-1.0)
adata.obs['leiden'] = adata.obs['leiden_0.8']
```

### Clustering Quality Assessment

```python
from sklearn.metrics import silhouette_score

# Silhouette score
sil = silhouette_score(
    adata.obsm['X_pca'][:, :n_pcs],
    adata.obs['leiden'],
    sample_size=5000
)
print(f"Silhouette score: {sil:.3f}")
```

---

## Step 6: Marker Genes & Cell Type Annotation

### Differential Expression

```python
# Find marker genes per cluster
sc.tl.rank_genes_groups(
    adata, 'leiden', 
    method='wilcoxon',
    use_raw=True,
    pts=True  # Include percentage expressing
)

# View top markers
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

# Extract results as DataFrame
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers_filtered = markers[
    (markers['pvals_adj'] < 0.05) & 
    (markers['logfoldchanges'] > 1)
]
```

### Manual Annotation with Known Markers

```python
# Define canonical marker genes by cell type
marker_dict = {
    'T cells': ['CD3D', 'CD3E', 'CD3G'],
    'CD4+ T': ['CD4', 'IL7R', 'CCR7'],
    'CD8+ T': ['CD8A', 'CD8B', 'GZMK'],
    'B cells': ['CD79A', 'MS4A1', 'CD19'],
    'NK cells': ['NKG7', 'GNLY', 'NCAM1'],
    'Monocytes': ['CD14', 'LYZ', 'S100A9'],
    'Dendritic': ['FCER1A', 'CST3', 'ITGAX'],
    'Platelets': ['PPBP', 'PF4'],
}

# Dot plot of markers
sc.pl.dotplot(
    adata, marker_dict, groupby='leiden',
    dendrogram=True, standard_scale='var'
)

# Assign cell type labels
cluster_to_celltype = {
    '0': 'CD4+ T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    # ... extend for all clusters
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)
```

### Automated Annotation (Optional)

```python
# Using celltypist for automated annotation
import celltypist
from celltypist import models

# Download and use pre-trained model
models.download_models(model='Immune_All_Low.pkl')
model = models.Model.load(model='Immune_All_Low.pkl')

predictions = celltypist.annotate(
    adata, model=model, majority_voting=True
)
adata.obs['celltypist_label'] = predictions.predicted_labels.majority_voting
```

---

## Step 7: Trajectory Inference (Optional)

```python
# PAGA for trajectory analysis
sc.tl.paga(adata, groups='cell_type')
sc.pl.paga(adata, threshold=0.03, show=True)

# Diffusion pseudotime
sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_dcs=10)

# Visualize trajectories
sc.pl.umap(adata, color=['cell_type', 'dpt_pseudotime'], 
           wspace=0.4)
```

---

## Step 8: Reporting & Export

### Generate Summary Figures

```python
# Final UMAP with cell type annotations
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color='cell_type', ax=axes[0], 
           title='Cell Types', show=False, frameon=False)
sc.pl.umap(adata, color='leiden', ax=axes[1], 
           title='Clusters', show=False, frameon=False)
plt.tight_layout()
plt.savefig('final_umap.png', dpi=300, bbox_inches='tight')

# Heatmap of top markers
sc.pl.rank_genes_groups_heatmap(
    adata, n_genes=5, groupby='cell_type',
    show_gene_labels=True, vmin=-3, vmax=3
)
```

### Export Results

```python
# Save annotated AnnData object
adata.write('results/annotated_data.h5ad')

# Export cluster markers to CSV
markers.to_csv('results/cluster_markers.csv', index=False)

# Export cell metadata
adata.obs.to_csv('results/cell_metadata.csv')

# Export UMAP coordinates
umap_df = pd.DataFrame(
    adata.obsm['X_umap'], 
    columns=['UMAP1', 'UMAP2'],
    index=adata.obs_names
)
umap_df.to_csv('results/umap_coordinates.csv')
```

---

## Configuration & Parameters

### Recommended Defaults

```yaml
# pipeline_config.yml
qc:
  min_genes: 200
  max_genes: 8000
  max_pct_mt: 20
  min_cells_per_gene: 3

normalization:
  target_sum: 10000
  n_top_genes: 3000

reduction:
  n_pcs: 40
  n_neighbors: 15

clustering:
  algorithm: leiden
  resolution: 0.8

de:
  method: wilcoxon
  min_logfc: 1.0
  max_pval: 0.05
```

---

## Dependencies

```
# requirements.txt
scanpy>=1.9.0
anndata>=0.9.0
numpy>=1.21
scipy>=1.7
pandas>=1.3
matplotlib>=3.5
seaborn>=0.12
scrublet>=0.2.3
harmonypy>=0.0.9
celltypist>=1.3.0
leidenalg>=0.9.0
```

---

## Troubleshooting

| Issue | Likely Cause | Solution |
|---|---|---|
| Too few cells after QC | Aggressive filtering | Relax thresholds, check distribution plots |
| Clusters don't separate in UMAP | Insufficient HVGs or PCs | Increase n_top_genes, try more PCs |
| Batch effects dominate UMAP | Uncorrected technical variation | Apply Harmony or scVI batch correction |
| No clear marker genes | Resolution too low/high | Try multiple resolutions, check DE method |
| Memory errors | Large dataset | Use backed mode, subsample for QC |

---

## References

- Wolf, F. A., Angerer, P. & Theis, F. J. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol.* **19**, 15 (2018).
- Wolock, S. L., Lopez, R. & Klein, A. M. Scrublet: Computational Identification of Cell Doublets in Single-Cell Transcriptomic Data. *Cell Systems* **8**, 281–291 (2019).
- Korsunsky, I. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat. Methods* **16**, 1289–1296 (2019).
