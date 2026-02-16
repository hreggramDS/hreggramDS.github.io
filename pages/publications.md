---
layout: page
title: Publications
description: Academic papers, technical reports, and conference presentations
permalink: /pages/publications/
last_modified_at: 2026-02-16
---

# Publications & Presentations

{% if site.data.publications %}

## Peer-Reviewed Publications

{% assign peer_reviewed = site.data.publications.papers | where: "type", "journal" %}
{% for pub in peer_reviewed %}
<div class="publication-entry">

**{{ pub.title }}**
- **Authors**: {{ pub.authors }}
- **{{ pub.venue }}**, {{ pub.year }}
{% if pub.abstract %}- {{ pub.abstract }}{% endif %}
- {% if pub.pdf %}[PDF]({{ pub.pdf }}){% endif %}{% if pub.doi %} \| [DOI]({{ pub.doi }}){% endif %}{% if pub.arxiv %} \| [arXiv]({{ pub.arxiv }}){% endif %}{% if pub.github %} \| [Code]({{ pub.github }}){% endif %}

</div>
---
{% endfor %}

{% else %}

## Academic Papers

### Ph.D. Dissertation

**Topological Data Analytics in Quantum Systems**
- **Author**: Greg Hamilton
- *University of Illinois at Urbana-Champaign*, 2022
- Applied persistent homology and topological invariants to many-body localization, relativistic heavy-ion collisions, and quantum information theory.
- [View on ProQuest](#)

---

### Peer-Reviewed Publications

**Persistent Homology of Many-Body Localization**
- **Authors**: G. Hamilton et al.
- *Physical Review B* (expected)
- Applied topological data analysis to detect phase transitions in disordered quantum spin chains. Demonstrated that persistent homology features provide robust, finite-size-stable diagnostics for the MBL transition.
- [arXiv](#) \| [Code](https://github.com/hreggramDS/tda-quantum)

---

**Topological Signatures in Heavy-Ion Collisions**
- **Authors**: G. Hamilton et al.
- *Working paper*
- Used persistent homology to analyze the spatial structure of particle production in relativistic heavy-ion collisions, extracting topological observables sensitive to collective flow and jet activity.
- [arXiv](#)

---

{% endif %}

## Conference Presentations

**Topological Data Analysis for Quantum Phase Transitions**
- **Conference**: APS March Meeting, Las Vegas, NV, 2023
- **Type**: Contributed Talk
- Presented TDA-based methodologies for detecting many-body localization transitions.

**Persistent Homology in Relativistic Heavy-Ion Collisions**
- **Conference**: APS Division of Nuclear Physics, New Orleans, LA, 2021
- **Type**: Contributed Talk
- Introduced topological observables for characterizing particle production geometry.

**Applications of TDA to Disordered Quantum Systems**
- **Conference**: Topology, Algebra, and Geometry in Data Science (TAGDS), Columbus, OH, 2022
- **Type**: Poster
- Demonstrated connections between persistent homology features and entanglement measures.

---

## Research Highlights

### Topological Data Analysis
Research applying persistent homology and topological methods to complex data structures in physics and quantum systems. Developed novel computational pipelines for computing and analyzing persistence diagrams of quantum state spaces.

### Many-Body Localization
Investigation of localization phenomena in quantum many-body systems. Used TDA to identify the MBL phase transition with improved finite-size stability compared to traditional entanglement entropy diagnostics.

### Single-Cell Genomics
Development of computational methods for single-cell RNA sequencing data analysis, including quality control frameworks, clustering validation, and reproducible pipeline architecture.

---

## External Profiles

- [Google Scholar](#) — *Add your profile link*
- [ORCID](#) — *Add your ORCID*
- [arXiv](#) — *Add your author page*

---

*This page will be updated as new publications become available. Publication details are placeholders — update with actual citation information.*
