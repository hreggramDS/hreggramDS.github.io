---
layout: post
title: "A Practical Guide to Topological Data Analysis"
date: 2025-06-10
categories: ["data-science", "tutorial"]
tags: ["topological-data-analysis", "persistent-homology", "python", "mathematics"]
author: "Greg Hamilton"
excerpt: "An accessible introduction to topological data analysis (TDA) — what it is, why it matters, and how to use it in Python with real data."
---

Topological Data Analysis (TDA) is one of those fields that sounds intimidating but is built on surprisingly intuitive ideas. At its core, TDA asks: **what is the *shape* of your data?** Not in the "is it normally distributed?" sense, but in a deeper geometric sense — does your data have clusters? Loops? Voids? And are those features real or just noise?

In this post, I'll walk through the key ideas behind TDA and show you how to apply them in Python.

## Why Care About Shape?

Consider a dataset sampled from a circle versus a blob. Standard summary statistics (mean, variance) won't distinguish them. PCA will project them both down to line segments. But topologically they're fundamentally different — the circle has a **hole** and the blob doesn't.

That hole is a topological feature, and persistent homology gives us a principled way to detect it.

## The Big Idea: Persistent Homology

Here's the intuition in four steps:

1. **Start with your point cloud** — a set of data points in some metric space
2. **Grow balls** around each point, gradually increasing the radius $\epsilon$
3. **Track topological features** as balls overlap and form connections:
   - **$H_0$**: Connected components (clusters merge as $\epsilon$ grows)
   - **$H_1$**: Loops (holes form and then fill in)
   - **$H_2$**: Voids (cavities in 3D structures)
4. **Record when features are born and die** — this gives you a *persistence diagram*

Features that persist across a wide range of $\epsilon$ values are likely real structure in your data. Short-lived features are likely noise. This is the fundamental principle of persistent homology.

## Getting Started in Python

### Setup

```bash
pip install ripser persim scikit-tda matplotlib numpy
```

### Your First Persistence Diagram

```python
import numpy as np
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt

# Generate a noisy circle
np.random.seed(42)
n = 100
theta = 2 * np.pi * np.random.rand(n)
noise = 0.1 * np.random.randn(n, 2)
circle = np.column_stack([np.cos(theta), np.sin(theta)]) + noise

# Compute persistent homology
result = ripser(circle, maxdim=1)

# Plot persistence diagram
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: original data
axes[0].scatter(circle[:, 0], circle[:, 1], s=10)
axes[0].set_title('Point Cloud (Noisy Circle)')
axes[0].set_aspect('equal')

# Right: persistence diagram
plot_diagrams(result['dgms'], ax=axes[1], show=False)
axes[1].set_title('Persistence Diagram')

plt.tight_layout()
plt.savefig('tda_intro.png', dpi=150)
plt.show()
```

### Reading the Persistence Diagram

In the resulting diagram:
- **$H_0$ points** (blue) correspond to connected components. Most die quickly (noise), but one lives forever (the single connected component of the circle).
- **$H_1$ points** (orange) correspond to loops. **One point stands far from the diagonal** — that's the circle's hole! Its high persistence (death - birth) tells us it's a robust feature.

The diagonal line represents features with zero lifetime (born and immediately die). Points far from the diagonal represent significant topological features.

## A More Interesting Example: Detecting Structure in High Dimensions

TDA really shines when your data is high-dimensional and you can't just eyeball it.

```python
from sklearn.datasets import make_moons, make_circles
from ripser import ripser
from persim import plot_diagrams

# Create three different datasets
datasets = {
    'Blobs': np.random.randn(200, 2) * 0.3 + np.array([[0, 0]]),
    'Two Moons': make_moons(200, noise=0.05)[0],
    'Nested Circles': make_circles(200, noise=0.05, factor=0.5)[0]
}

fig, axes = plt.subplots(3, 2, figsize=(12, 15))

for idx, (name, data) in enumerate(datasets.items()):
    # Scatter plot
    axes[idx, 0].scatter(data[:, 0], data[:, 1], s=10)
    axes[idx, 0].set_title(f'{name} — Point Cloud')
    axes[idx, 0].set_aspect('equal')
    
    # Persistence diagram
    result = ripser(data, maxdim=1)
    plot_diagrams(result['dgms'], ax=axes[idx, 1], show=False)
    axes[idx, 1].set_title(f'{name} — Persistence Diagram')

plt.tight_layout()
plt.show()
```

Each dataset produces a characteristically different persistence diagram:
- **Blobs**: Short-lived $H_0$ features (noise) and no significant $H_1$ (no loops)
- **Two Moons**: Two persistent $H_0$ components (two clusters) and mild $H_1$ activity
- **Nested Circles**: Two very persistent $H_1$ features (two separate loops!)

## Using TDA Features for Machine Learning

Persistence diagrams aren't directly usable as feature vectors (they're multi-sets of varying size). But several vectorization methods exist:

### Persistence Images

```python
from persim import PersistenceImager

pimgr = PersistenceImager(pixel_size=0.1, birth_range=(0, 1.5), pers_range=(0, 1.5))

# Convert persistence diagram to a fixed-size image
dgm_h1 = result['dgms'][1]
pimg = pimgr.transform(dgm_h1)

# pimg is now a 2D array that can be used as features for any ML model
```

### Persistence Landscapes

```python
def persistence_landscape(dgm, k=1, resolution=100):
    """Compute the k-th persistence landscape."""
    finite = dgm[np.isfinite(dgm[:, 1])]
    births, deaths = finite[:, 0], finite[:, 1]
    
    t_max = np.max(deaths) if len(deaths) > 0 else 1
    t = np.linspace(0, t_max, resolution)
    
    tent_functions = np.maximum(
        np.minimum(
            t[None, :] - births[:, None],
            deaths[:, None] - t[None, :]
        ), 0
    )
    
    sorted_tents = np.sort(tent_functions, axis=0)[::-1]
    
    if k <= len(sorted_tents):
        return t, sorted_tents[k-1]
    return t, np.zeros(resolution)
```

## When to Use TDA

TDA is particularly useful when:

- Your data has **complex geometric structure** that simple statistics miss
- You're working in **high dimensions** where visualization fails
- You need features that are **robust to noise** and **coordinate-free**
- You want to detect **phase transitions** or **structural changes** in parameterized data
- Standard clustering struggles because your data has **non-convex shapes**

## When *Not* to Use TDA

Be honest about its limitations:

- **Computational cost**: Rips complex construction is $O(n^2)$ in memory and up to $O(2^n)$ in theory (though persistent homology algorithms are much more efficient in practice)
- **Interpretability**: Explaining what $H_1$ features mean to a business stakeholder requires more effort than explaining a regression coefficient
- **Small datasets**: With few points, the persistence diagrams are noisy and unstable
- **Tabular data**: If your data is purely tabular with no inherent geometry, TDA may not add value over standard ML

## Further Reading

- Edelsbrunner & Harer, *Computational Topology: An Introduction* — The theoretical foundation
- Carlsson, *Topology and Data* (2009) — The seminal survey paper
- [Ripser documentation](https://ripser.scikit-tda.org/) — The fastest persistent homology library
- [GUDHI library](https://gudhi.inria.fr/) — Comprehensive TDA toolkit with excellent tutorials
- My project on [TDA for Quantum Systems](/projects/tda-quantum-systems/) — How I applied these ideas in my research

---

*This post is part of a planned series on applied mathematics for data science. Next up: using the Mapper algorithm for unsupervised data exploration.*
