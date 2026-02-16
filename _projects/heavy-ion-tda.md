---
layout: project
title: "TDA for Relativistic Heavy-Ion Collisions"
date: 2021-08-01
team: ["Greg Hamilton"]
status: "Completed"
domain: ["Physics", "Topological Data Analysis"]
technologies: ["Python", "C++", "ROOT", "NumPy", "Ripser", "Matplotlib"]
description: "Topological analysis of particle production patterns in relativistic heavy-ion collision data"
featured: true
---

## Overview

Applied topological data analysis methods to study the spatial and momentum distributions of particles produced in relativistic heavy-ion collisions. By treating particle production events as point clouds and computing their persistent homology, we extracted geometric and topological signatures of the quark-gluon plasma (QGP) — a state of matter where quarks and gluons are deconfined at extreme temperatures and densities.

## Problem Statement

Relativistic heavy-ion collisions at facilities like RHIC and the LHC produce thousands of particles per event, creating complex multi-particle final states. Key challenges include:

- **Identifying collective flow patterns** that signal QGP formation
- **Characterizing event-by-event fluctuations** in particle distributions
- **Distinguishing** signal (QGP signatures) from background (hadronic rescattering, detector effects)
- Traditional observables (flow harmonics $v_n$, correlations) capture only partial geometric information

TDA provides a complementary lens by capturing the full topological structure of particle distributions.

## Methodology

### Analysis Pipeline

1. **Event Selection** — Apply centrality and kinematic cuts to collision events
2. **Point Cloud Construction** — Map each event's particles to $(\eta, \phi)$ pseudorapidity-azimuthal angle space
3. **Persistent Homology** — Compute Rips complex filtration for each event
4. **Feature Extraction** — Extract Betti numbers, persistence diagrams, and persistence images
5. **Statistical Analysis** — Compare topological features across centrality classes and collision energies

### Topological Observables

For each collision event with $N$ particles at positions $\{(\eta_i, \phi_i)\}_{i=1}^N$:

- **$\beta_0$ persistence**: Tracks merging of particle clusters — sensitive to jet structure and flow
- **$\beta_1$ persistence**: Tracks formation and filling of holes — sensitive to ring-like flow patterns
- **Persistence entropy**: $S = -\sum_i p_i \log p_i$ where $p_i = \ell_i / L$ measures topological complexity

## Implementation

### Event Processing

```python
import numpy as np
from ripser import ripser
from persim import PersistenceImager

def process_event(particles, max_dim=1):
    """Compute persistent homology for a single collision event.
    
    Args:
        particles: array of shape (N, 2) with (eta, phi) coordinates
        max_dim: maximum homological dimension
    
    Returns:
        persistence diagrams and derived features
    """
    # Handle periodic boundary in phi
    # Use product metric: d = sqrt(d_eta^2 + d_phi^2)
    # where d_phi accounts for 2pi periodicity
    n = len(particles)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            d_eta = particles[i, 0] - particles[j, 0]
            d_phi = particles[i, 1] - particles[j, 1]
            # Periodic boundary condition in azimuth
            d_phi = min(abs(d_phi), 2*np.pi - abs(d_phi))
            dist[i, j] = dist[j, i] = np.sqrt(d_eta**2 + d_phi**2)
    
    result = ripser(dist, maxdim=max_dim, distance_matrix=True)
    return result['dgms']

def extract_features(dgms):
    """Extract summary statistics from persistence diagrams."""
    features = {}
    for dim, dgm in enumerate(dgms):
        finite = dgm[np.isfinite(dgm[:, 1])]
        lifetimes = finite[:, 1] - finite[:, 0]
        
        features[f'H{dim}_count'] = len(finite)
        features[f'H{dim}_max_persistence'] = np.max(lifetimes) if len(lifetimes) > 0 else 0
        features[f'H{dim}_mean_persistence'] = np.mean(lifetimes) if len(lifetimes) > 0 else 0
        features[f'H{dim}_total_persistence'] = np.sum(lifetimes)
        features[f'H{dim}_entropy'] = persistence_entropy(lifetimes)
    
    return features

def persistence_entropy(lifetimes):
    """Compute persistence entropy."""
    if len(lifetimes) == 0 or np.sum(lifetimes) == 0:
        return 0.0
    probs = lifetimes / np.sum(lifetimes)
    return -np.sum(probs * np.log(probs + 1e-10))
```

### Centrality Dependence Analysis

```python
def analyze_centrality_dependence(events_by_centrality):
    """Compare topological features across centrality classes."""
    results = {}
    for cent_class, events in events_by_centrality.items():
        features_list = []
        for particles in events:
            dgms = process_event(particles)
            features_list.append(extract_features(dgms))
        
        df = pd.DataFrame(features_list)
        results[cent_class] = {
            'mean': df.mean(),
            'std': df.std(),
            'n_events': len(events)
        }
    
    return results
```

## Results

### Key Findings

- **Flow Sensitivity**: $H_1$ persistence is strongly correlated with elliptic flow coefficient $v_2$, confirming that TDA captures collective behavior
- **Centrality Scaling**: Total persistence scales monotonically with centrality, reflecting increasing particle multiplicity and collectivity
- **Event Classification**: Topological features can distinguish central from peripheral collisions with ~90% accuracy using simple classifiers
- **Complementary Information**: TDA features contain geometric information beyond what is captured by standard flow harmonics alone

### Centrality Dependence of Topological Features

| Centrality | Avg $H_0$ Persistence | Avg $H_1$ Persistence | Persistence Entropy |
|---|---|---|---|
| 0–10% (central) | 12.4 ± 0.3 | 3.8 ± 0.2 | 4.2 ± 0.1 |
| 10–30% | 9.7 ± 0.2 | 3.1 ± 0.2 | 3.8 ± 0.1 |
| 30–50% | 6.8 ± 0.2 | 2.1 ± 0.1 | 3.3 ± 0.1 |
| 50–80% (peripheral) | 3.2 ± 0.1 | 0.9 ± 0.1 | 2.5 ± 0.1 |

## Deliverables

- **Analysis Code**: Python/C++ codebase for TDA analysis of collision data
- **Publication**: Results included in Ph.D. dissertation
- **Methodology**: Generalizable framework for TDA on point process data

## Learnings & Future Work

- Periodic boundary conditions require careful handling in persistent homology computations
- Persistence images provide an effective vectorization for machine learning on topological features
- Natural extensions include time-dependent TDA for hydrodynamic evolution and application to event generators
