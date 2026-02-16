---
layout: project
title: "Topological Data Analysis for Quantum Systems"
date: 2022-05-15
team: ["Greg Hamilton"]
status: "Completed"
domain: ["Physics", "Topological Data Analysis"]
technologies: ["Python", "NumPy", "SciPy", "Ripser", "GUDHI", "Matplotlib"]
github: "https://github.com/hreggramDS/tda-quantum"
description: "Applied persistent homology and topological invariants to detect phase transitions in many-body localized quantum systems"
featured: true
---

## Overview

This project — the core of my doctoral research — applies techniques from topological data analysis (TDA) to quantum many-body physics. By computing persistent homology on quantum state data, we can detect phase transitions, characterize entanglement structure, and extract topological features that are invisible to conventional order parameters.

## Problem Statement

Identifying phase transitions in disordered quantum systems is notoriously difficult:

- **Many-body localization (MBL)** transitions lack a simple local order parameter
- Traditional spectral and entanglement-based diagnostics can be noisy and system-size dependent
- High-dimensional Hilbert spaces make visualization and intuitive understanding challenging

We needed a tool that could capture the *shape* of quantum data — and TDA provides exactly that.

## Methodology

### Persistent Homology Pipeline

1. **State Preparation** — Generate ground and excited states for disordered Hamiltonians across disorder strengths
2. **Point Cloud Construction** — Embed quantum states (or their reduced density matrices) as point clouds in metric spaces
3. **Filtration** — Build Vietoris-Rips or alpha complexes across a range of scales
4. **Persistence Computation** — Track birth and death of topological features (connected components, loops, voids)
5. **Feature Extraction** — Compute persistence diagrams, Betti curves, and persistence landscapes
6. **Classification** — Use topological features to classify phases and detect transitions

### Mathematical Framework

The key insight is that persistent homology captures multi-scale topological structure. For a point cloud $X$ and distance parameter $\epsilon$:

- $H_0(\epsilon)$ tracks connected components (clustering)
- $H_1(\epsilon)$ tracks loops (cyclical structure)
- $H_2(\epsilon)$ tracks voids (higher-order cavities)

The persistence diagram $\text{Dgm}_k(X)$ records the birth-death pairs of $k$-dimensional features.

## Implementation

### Core Computation

```python
import numpy as np
from ripser import ripser
from persim import plot_diagrams
import gudhi

def compute_persistence(states, max_dim=2):
    """Compute persistent homology of quantum state point cloud."""
    # Build distance matrix from state overlaps
    n = len(states)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            overlap = np.abs(np.vdot(states[i], states[j]))**2
            dist_matrix[i, j] = np.sqrt(1 - overlap)
            dist_matrix[j, i] = dist_matrix[i, j]
    
    # Compute persistent homology
    result = ripser(dist_matrix, maxdim=max_dim, distance_matrix=True)
    return result['dgms']

def persistence_landscape(dgm, num_landscapes=5, resolution=1000):
    """Compute persistence landscapes for statistical analysis."""
    births = dgm[:, 0]
    deaths = dgm[:, 1]
    t_range = np.linspace(0, np.max(deaths[np.isfinite(deaths)]), resolution)
    
    landscapes = np.zeros((num_landscapes, resolution))
    for t_idx, t in enumerate(t_range):
        tent_values = np.maximum(
            np.minimum(t - births, deaths - t), 0
        )
        sorted_vals = np.sort(tent_values)[::-1]
        for k in range(min(num_landscapes, len(sorted_vals))):
            landscapes[k, t_idx] = sorted_vals[k]
    
    return landscapes, t_range
```

### Phase Transition Detection

```python
def detect_transition(disorder_values, persistence_features):
    """Detect MBL transition from topological features."""
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(persistence_features)
    
    # Two-phase clustering
    kmeans = KMeans(n_clusters=2, random_state=42)
    labels = kmeans.fit_predict(features_scaled)
    
    # Find transition point
    transitions = np.where(np.diff(labels) != 0)[0]
    critical_disorder = disorder_values[transitions]
    
    return critical_disorder, labels
```

## Results

### Key Findings

- **Phase Detection**: Persistent homology successfully distinguishes thermal and MBL phases using $H_1$ features
- **Critical Point**: Topological features identify the MBL transition at disorder strengths consistent with exact diagonalization studies
- **Scaling Behavior**: Persistence landscapes show characteristic scaling near the critical point
- **Robustness**: TDA-based diagnostics are more robust to finite-size effects than traditional entanglement entropy measures

### Performance Comparison

| Method | System Size (L=12) | System Size (L=16) | Finite-Size Stability |
|---|---|---|---|
| Entanglement Entropy | Critical W ≈ 3.5 | Critical W ≈ 3.7 | Moderate drift |
| Level Statistics | Critical W ≈ 3.6 | Critical W ≈ 3.7 | Some drift |
| **TDA (This Work)** | **Critical W ≈ 3.6** | **Critical W ≈ 3.65** | **Minimal drift** |

## Deliverables

- **Ph.D. Dissertation**: Full treatment published via ProQuest
- **Computational Pipeline**: Reusable Python codebase for TDA on quantum data
- **Methodology Guide**: Step-by-step protocol for applying persistent homology to physics problems

## Learnings & Future Work

- TDA provides a genuinely complementary perspective to spectral and entanglement-based methods
- Persistence landscapes enable rigorous statistical comparison across parameter regimes
- Future directions include applying similar methods to quantum error correction codes and topological quantum computing
- Opportunities to extend to time-series topological analysis (crocker plots) for dynamical systems
