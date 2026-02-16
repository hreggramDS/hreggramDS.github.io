---
layout: project
title: "Advanced Analytics for Management Consulting"
date: 2023-03-01
team: ["Greg Hamilton"]
status: "Completed"
domain: ["Business Intelligence", "Machine Learning"]
technologies: ["Python", "scikit-learn", "SQL", "Tableau", "Excel"]
description: "Data-driven frameworks for strategic decision-making across Fortune 500 engagements at BCG"
featured: true
---

## Overview

During my time at Boston Consulting Group, I developed and deployed analytical frameworks that translated complex data into actionable business strategy. These projects spanned industries including manufacturing, healthcare, and retail — each requiring rapid iteration, clear stakeholder communication, and measurable business outcomes.

## Problem Context

Management consulting engagements have unique constraints:

- **Tight timelines** — Deliver actionable insights in 4–8 week sprints
- **Diverse domains** — Each engagement may involve a new industry and dataset
- **Stakeholder translation** — Results must be interpretable by C-suite executives, not just analysts
- **Impact orientation** — Models and analyses must connect directly to dollar-value decisions

## Selected Engagements

### Customer Segmentation & Targeting (Retail)

**Challenge**: A major retailer needed to optimize marketing spend across customer segments with declining response rates.

**Approach**:
1. Integrated transaction, demographic, and behavioral data across 10M+ customer records
2. Engineered RFM (Recency, Frequency, Monetary) features and augmented with behavioral signals
3. Applied K-means and Gaussian Mixture Models for multi-dimensional segmentation
4. Built propensity-to-purchase models per segment using logistic regression and gradient boosting

```python
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

# Multi-dimensional customer segmentation
features = ['recency', 'frequency', 'monetary', 
            'channel_diversity', 'category_breadth',
            'avg_basket_size', 'promo_sensitivity']

scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_customers[features])

# Optimal cluster selection via BIC
bic_scores = []
for k in range(2, 12):
    gmm = GaussianMixture(n_components=k, random_state=42)
    gmm.fit(X_scaled)
    bic_scores.append(gmm.bic(X_scaled))

optimal_k = np.argmin(bic_scores) + 2
gmm_final = GaussianMixture(n_components=optimal_k, random_state=42)
df_customers['segment'] = gmm_final.fit_predict(X_scaled)
```

**Impact**: Identified 3 high-value micro-segments comprising 8% of the customer base but driving 35% of margin. Recommended targeted campaigns projected to increase ROI on marketing spend by 20%+.

---

### Supply Chain Optimization (Manufacturing)

**Challenge**: A global manufacturer experienced frequent stockouts in some product lines while carrying excess inventory in others.

**Approach**:
1. Built demand forecasting models at the SKU-location level using historical sales, seasonality, and leading indicators
2. Implemented safety stock optimization using service-level-aware newsvendor models
3. Created a simulation framework to stress-test inventory policies under demand uncertainty

```python
from scipy.optimize import minimize_scalar
from scipy.stats import norm

def optimal_safety_stock(demand_mean, demand_std, lead_time, 
                          service_level=0.95, holding_cost=1.0, 
                          stockout_cost=10.0):
    """Compute optimal safety stock balancing holding and stockout costs."""
    z = norm.ppf(service_level)
    lt_demand_std = demand_std * np.sqrt(lead_time)
    safety_stock = z * lt_demand_std
    return safety_stock

def simulate_policy(demand_dist, reorder_point, order_quantity, 
                     lead_time, n_days=365, n_sims=1000):
    """Monte Carlo simulation of inventory policy performance."""
    results = []
    for _ in range(n_sims):
        inventory = order_quantity
        stockout_days = 0
        total_holding = 0
        
        for day in range(n_days):
            demand = demand_dist.rvs()
            inventory -= demand
            if inventory < 0:
                stockout_days += 1
                inventory = 0
            if inventory <= reorder_point:
                inventory += order_quantity  # simplified
            total_holding += max(inventory, 0)
        
        results.append({
            'fill_rate': 1 - stockout_days/n_days,
            'avg_inventory': total_holding/n_days
        })
    return pd.DataFrame(results)
```

**Impact**: Recommended inventory policy changes projected to reduce stockouts by 40% while decreasing average inventory carrying costs by 15%.

---

### Operational Efficiency (Healthcare)

**Challenge**: A hospital network needed to reduce patient wait times and optimize staff scheduling across departments.

**Approach**:
1. Analyzed 2 years of patient flow data (arrivals, service times, discharge patterns)
2. Built queuing models to identify bottleneck departments and peak demand periods
3. Developed staff scheduling optimization using mixed-integer programming
4. Created interactive dashboards for real-time monitoring

**Impact**: Recommended scheduling changes projected to reduce average ED wait times by 25% without increasing staffing costs.

## Methodology

### Consulting Analytics Framework

Every engagement followed a structured analytical approach:

| Phase | Activities | Deliverables |
|---|---|---|
| **Discovery** (Week 1) | Data audit, stakeholder interviews, hypothesis generation | Data landscape map, initial hypotheses |
| **Analysis** (Weeks 2–4) | Feature engineering, modeling, validation | Models, key findings, sensitivity analysis |
| **Synthesis** (Weeks 5–6) | Insight extraction, business translation, recommendations | Executive presentation, implementation roadmap |
| **Handoff** | Documentation, knowledge transfer, tool deployment | Technical docs, trained client teams |

## Key Takeaways

- **Speed over perfection**: In consulting, a good-enough model delivered Tuesday beats a perfect model delivered next month
- **Storytelling matters**: The best analysis fails if it can't be communicated in a 20-minute exec presentation
- **Domain immersion**: Spending time on shop floors, in hospitals, and in stores is as important as time in Jupyter notebooks
- **Reusable frameworks**: Building modular analytical templates dramatically accelerates subsequent engagements

## Deliverables

- **Analytical Frameworks**: Reusable templates for segmentation, demand forecasting, and optimization
- **Client Presentations**: C-suite ready decks translating analysis into strategic recommendations
- **Implementation Guides**: Step-by-step deployment plans for recommended changes
- **Knowledge Transfer**: Training materials for client data teams
