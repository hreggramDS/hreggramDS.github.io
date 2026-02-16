---
layout: project
title: "Predictive Maintenance for Heavy Equipment"
date: 2024-06-01
team: ["Greg Hamilton"]
status: "In Progress"
domain: ["Machine Learning", "Industrial IoT"]
technologies: ["Python", "XGBoost", "LightGBM", "Spark", "AWS", "Docker"]
description: "Machine learning system for predicting component failures in construction and mining equipment using sensor telemetry data"
featured: true
---

## Overview

Developed and deployed a machine learning system that predicts component failures in heavy construction and mining equipment, enabling proactive maintenance scheduling. The system ingests real-time sensor telemetry, extracts degradation signals, and generates failure probability forecasts at the component level.

## Problem Statement

Unplanned equipment downtime in construction and mining operations is extremely costly:

- A single day of downtime for a large mining truck can cost **$10,000–$50,000+** in lost productivity
- Reactive maintenance leads to cascading failures that compound repair costs
- Over-maintenance on fixed schedules wastes parts and labor while still missing unexpected failures

The goal: **predict failures before they happen**, targeting the optimal maintenance window that minimizes both downtime and unnecessary servicing.

## Methodology

### System Architecture

```
Sensor Data → Ingestion Layer → Feature Engineering → ML Models → Alert System
                  │                     │                  │            │
              AWS S3/Kinesis      Time-series          XGBoost      Dashboard
                                 aggregations         LightGBM     + API
                                 rolling stats        Survival
                                 spectral features    Analysis
```

### Feature Engineering

Key engineered features from raw sensor telemetry:

1. **Rolling Statistics** — Mean, std, skew, kurtosis over 1h/6h/24h/7d windows
2. **Rate-of-Change Features** — First and second derivatives of key signals
3. **Spectral Features** — FFT-based frequency domain features for vibration data
4. **Operational Context** — Load factor, duty cycle, environmental conditions
5. **Cumulative Degradation** — Exponentially-weighted moving averages capturing long-term drift

### Model Approach

Two complementary modeling strategies:

**Binary Classification** — Will a failure occur in the next N days?
```python
import xgboost as xgb
from sklearn.model_selection import TimeSeriesSplit

model = xgb.XGBClassifier(
    n_estimators=500,
    max_depth=6,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    scale_pos_weight=ratio_neg_to_pos,
    eval_metric='aucpr'
)

tscv = TimeSeriesSplit(n_splits=5)
for train_idx, val_idx in tscv.split(X):
    model.fit(
        X.iloc[train_idx], y.iloc[train_idx],
        eval_set=[(X.iloc[val_idx], y.iloc[val_idx])],
        early_stopping_rounds=50,
        verbose=False
    )
```

**Survival Analysis** — Estimate remaining useful life (RUL)
```python
from lifelines import CoxPHFitter, WeibullAFTFitter

# Weibull Accelerated Failure Time model
aft = WeibullAFTFitter()
aft.fit(
    df_survival,
    duration_col='time_to_event',
    event_col='failure_observed',
    formula="temperature + vibration + pressure + hours_since_service"
)

# Predict median remaining life
rul_predictions = aft.predict_median(new_data)
```

## Implementation

### Data Pipeline

```python
from pyspark.sql import SparkSession
import pyspark.sql.functions as F
from pyspark.sql.window import Window

spark = SparkSession.builder.appName("PredMaint").getOrCreate()

# Define rolling window specifications
windows = {
    '1h': Window.partitionBy('equipment_id')
               .orderBy('timestamp')
               .rangeBetween(-3600, 0),
    '24h': Window.partitionBy('equipment_id')
                .orderBy('timestamp')
                .rangeBetween(-86400, 0),
}

# Compute rolling features
df_features = df_raw.select(
    'equipment_id', 'timestamp',
    *[F.avg(col).over(windows['1h']).alias(f'{col}_avg_1h')
      for col in sensor_cols],
    *[F.stddev(col).over(windows['24h']).alias(f'{col}_std_24h')
      for col in sensor_cols],
)
```

### Model Serving

- Models containerized with Docker and served via REST API
- Batch predictions run nightly on full fleet; real-time scoring for critical alerts
- Model registry tracks versions, metrics, and A/B test assignments
- Automated retraining triggered on data drift detection

## Results

### Key Metrics

| Metric | Value |
|---|---|
| Precision @ 7-day horizon | 82% |
| Recall @ 7-day horizon | 74% |
| PR-AUC | 0.79 |
| False alarm rate | < 5% |
| Avg. lead time before failure | 4.2 days |

### Business Impact

- **30%+ reduction** in unplanned downtime events for monitored components
- Maintenance scheduling improved from reactive to proactive for top failure modes
- Framework extensible to additional equipment types and component categories

## Deliverables

- **ML Pipeline**: End-to-end training and serving infrastructure
- **Feature Store**: Reusable feature engineering library for telemetry data
- **Dashboard**: Real-time equipment health monitoring and alert management
- **Documentation**: Model cards, data dictionaries, and operational runbooks

## Learnings & Future Work

- Time-series feature engineering is more impactful than model complexity for this domain
- Survival analysis models provide more actionable output than binary classifiers for maintenance planning
- Future work includes incorporating computer vision on inspection images and expanding to fleet-level optimization
