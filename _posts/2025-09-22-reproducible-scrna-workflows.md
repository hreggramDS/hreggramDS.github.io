---
layout: post
title: "Building Reproducible scRNA-seq Workflows"
date: 2025-09-22
categories: ["bioinformatics", "tutorial"]
tags: ["scRNA-seq", "scanpy", "reproducibility", "python", "single-cell"]
author: "Greg Hamilton"
excerpt: "Lessons learned from building production-grade single-cell RNA sequencing pipelines — from reproducibility patterns to common pitfalls."
---

Single-cell RNA sequencing has become the workhorse of modern genomics, but building analysis workflows that are reproducible, scalable, and maintainable is harder than the tutorials make it look. After running dozens of scRNA-seq analyses across different datasets and biological contexts, here are the patterns and practices that have served me well.

## The Reproducibility Problem

Most scRNA-seq tutorials show you the *happy path*: load data, filter, normalize, cluster, plot UMAP, done. In practice, you're often:

- Re-running analyses with different QC thresholds after discussing with biologists
- Integrating multiple samples with batch effects
- Trying different clustering resolutions and comparing results
- Going back to earlier steps when downstream results look wrong
- Sharing notebooks with collaborators who have different Python environments

Without intentional structure, this becomes a mess of Jupyter notebooks with run-out-of-order cells, hardcoded paths, and results that can't be regenerated.

## Pattern 1: Separate Configuration from Code

The single biggest improvement I've made is extracting all parameters into a configuration file:

```yaml
# config.yml
sample:
  name: "experiment_2025"
  input_path: "data/raw/filtered_feature_bc_matrix.h5"
  output_dir: "results/experiment_2025"

qc:
  min_genes: 200
  max_genes: 6000
  max_pct_mt: 15
  min_cells_per_gene: 3
  doublet_rate: 0.06

normalization:
  target_sum: 10000
  n_top_genes: 3000
  hvg_flavor: "seurat_v3"  # or "seurat", "cell_ranger"

dimensionality_reduction:
  n_pcs: 40
  n_neighbors: 15
  umap_min_dist: 0.3

clustering:
  algorithm: "leiden"  # or "louvain"
  resolutions: [0.3, 0.5, 0.8, 1.0, 1.5]
  selected_resolution: 0.8

batch_correction:
  enabled: false
  method: "harmony"  # or "scvi", "bbknn"
  batch_key: "sample_id"
```

```python
import yaml

def load_config(path='config.yml'):
    with open(path) as f:
        return yaml.safe_load(f)

config = load_config()
```

Now when a collaborator asks "what QC thresholds did you use?", the answer is a version-controlled YAML file, not "let me find the notebook cell where I set that."

## Pattern 2: Checkpoint Your AnnData Objects

scRNA-seq workflows are iterative. You don't want to re-run Cell Ranger every time you change your clustering resolution. Save intermediate states:

```python
import scanpy as sc
from pathlib import Path

def save_checkpoint(adata, step_name, output_dir):
    """Save AnnData checkpoint with metadata."""
    path = Path(output_dir) / f"{step_name}.h5ad"
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write(path)
    print(f"Checkpoint saved: {path} ({adata.n_obs} cells, {adata.n_vars} genes)")

def load_checkpoint(step_name, output_dir):
    """Load AnnData from checkpoint."""
    path = Path(output_dir) / f"{step_name}.h5ad"
    if path.exists():
        adata = sc.read(path)
        print(f"Loaded checkpoint: {path} ({adata.n_obs} cells, {adata.n_vars} genes)")
        return adata
    raise FileNotFoundError(f"No checkpoint found: {path}")

# Usage
save_checkpoint(adata, '01_raw_loaded', config['sample']['output_dir'])
save_checkpoint(adata, '02_qc_filtered', config['sample']['output_dir'])
save_checkpoint(adata, '03_normalized', config['sample']['output_dir'])
save_checkpoint(adata, '04_clustered', config['sample']['output_dir'])
```

This lets you restart from any step without re-running everything upstream.

## Pattern 3: Structured Logging over Print Statements

Replace ad-hoc print statements with structured logging that creates an audit trail:

```python
import logging
from datetime import datetime

def setup_logger(output_dir, name='scrna_pipeline'):
    """Set up logging to both file and console."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # File handler
    log_path = Path(output_dir) / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)s | %(message)s'
    ))
    logger.addHandler(fh)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(ch)
    
    return logger

logger = setup_logger(config['sample']['output_dir'])

# Now use throughout the pipeline
logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
logger.info(f"QC filter: min_genes={config['qc']['min_genes']}, max_pct_mt={config['qc']['max_pct_mt']}")
logger.info(f"After QC: {adata.n_obs} cells retained ({pct_retained:.1f}%)")
```

## Pattern 4: Defensive QC with Visual Checks

Don't just apply filters — generate plots that let you *verify* your thresholds are appropriate for each dataset:

```python
def qc_report(adata, config, output_dir):
    """Generate comprehensive QC report with filter visualization."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    # Distribution plots with filter lines
    metrics = [
        ('n_genes_by_counts', config['qc']['min_genes'], config['qc']['max_genes']),
        ('total_counts', None, None),
        ('pct_counts_mt', None, config['qc']['max_pct_mt']),
    ]
    
    for ax, (metric, lo, hi) in zip(axes[0], metrics):
        ax.hist(adata.obs[metric], bins=100, edgecolor='black', alpha=0.7)
        if lo is not None:
            ax.axvline(lo, color='red', linestyle='--', label=f'min={lo}')
        if hi is not None:
            ax.axvline(hi, color='red', linestyle='--', label=f'max={hi}')
        ax.set_xlabel(metric)
        ax.set_ylabel('Cells')
        ax.legend()
    
    # Joint distributions
    axes[1, 0].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], 
                        s=1, alpha=0.3)
    axes[1, 0].set_xlabel('Total Counts')
    axes[1, 0].set_ylabel('N Genes')
    
    axes[1, 1].scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'], 
                        s=1, alpha=0.3)
    axes[1, 1].set_xlabel('Total Counts')
    axes[1, 1].set_ylabel('% Mitochondrial')
    
    # Summary statistics
    summary_text = (
        f"Total cells: {adata.n_obs}\n"
        f"Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}\n"
        f"Median UMIs/cell: {adata.obs['total_counts'].median():.0f}\n"
        f"Median % MT: {adata.obs['pct_counts_mt'].median():.1f}%"
    )
    axes[1, 2].text(0.1, 0.5, summary_text, fontsize=14, 
                     transform=axes[1, 2].transAxes, verticalalignment='center',
                     fontfamily='monospace')
    axes[1, 2].axis('off')
    
    plt.suptitle('Quality Control Report', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'qc_report.png', dpi=150, bbox_inches='tight')
    plt.close()
```

## Pattern 5: Pin Your Environment

This seems obvious but is constantly overlooked. Every analysis should have a pinned environment:

```bash
# Create environment from scratch
conda create -n scrna python=3.10
conda activate scrna
pip install scanpy==1.9.6 anndata==0.10.3 numpy==1.24.4 \
    scipy==1.11.4 matplotlib==3.8.2 scrublet==0.2.3

# Lock it
pip freeze > requirements.txt

# Or better yet, use conda
conda env export > environment.yml
```

I've seen analyses break months later because a `scanpy` update changed default parameters for `sc.pp.highly_variable_genes`. Pinning prevents this.

## Common Pitfalls

### 1. Running Cells Out of Order

The most common source of irreproducible notebooks. Mitigations:
- Use `nbconvert --execute` to verify notebooks run top-to-bottom
- Or convert critical analysis steps to `.py` scripts that run linearly

### 2. Forgetting to Store Raw Counts

Once you normalize and scale, you need the original counts for differential expression:

```python
# Do this BEFORE normalizing
adata.layers['counts'] = adata.X.copy()

# Later, for DE analysis
sc.tl.rank_genes_groups(adata, 'leiden', use_raw=True)
```

### 3. Over-clustering or Under-clustering

No single resolution is "correct." Always try multiple and compare with known biology:

```python
for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
```

### 4. Ignoring Batch Effects

If you have multiple samples, check for batch effects before clustering:

```python
sc.pl.umap(adata, color='sample_id')  # Do samples mix or separate?
```

If samples cluster by batch rather than biology, apply correction (Harmony, scVI, BBKNN) before proceeding.

## Putting It Together

A well-structured scRNA-seq project looks like this:

```
project/
├── config.yml                 # All parameters
├── environment.yml            # Pinned dependencies
├── scripts/
│   ├── 01_qc.py              # Quality control
│   ├── 02_normalize.py       # Normalization & HVG
│   ├── 03_reduce.py          # PCA, batch correction, UMAP
│   ├── 04_cluster.py         # Clustering & markers
│   └── 05_annotate.py        # Cell type annotation
├── notebooks/
│   └── exploration.ipynb     # Interactive analysis & figures
├── data/
│   ├── raw/                  # Input data (not version controlled)
│   └── processed/            # Intermediate checkpoints
├── results/
│   ├── figures/              # Publication-ready plots
│   ├── tables/               # Marker genes, cell counts
│   └── annotated.h5ad        # Final annotated object
└── README.md                 # How to run the analysis
```

Each numbered script reads a checkpoint, does one thing, and saves the next checkpoint. The notebook is for exploration and figure generation, not for the core pipeline logic.

## Want the Full Pipeline?

I've documented my complete scRNA-seq pipeline in the [Pipelines section](/pages/pipelines/) of this site, including Cell Ranger configuration, QC thresholds, and downstream analysis steps.

---

*What reproducibility patterns have worked for you? I'm always looking to improve my workflows. Reach out on [GitHub](https://github.com/hreggramDS) or the [contact page](/pages/contact/).*
