# ULM vs MLM: Understanding the Difference

## Overview

This document explains the difference between **Univariate Linear Model (ULM)** and **Multivariate Linear Model (MLM)** as implemented in decoupler, and why we use MLM in this pipeline.

## Key Differences

### Univariate Linear Model (ULM)

- **Approach**: Fits a separate linear model for each regulator/source independently
- **Model**: For each sample and each regulator, fits: `genes = β × regulator_weights + ε`
- **Assumption**: Each regulator's activity is assessed independently, without considering interactions with other regulators
- **Use case**: Good for identifying which individual regulators are active, but may miss interactions

### Multivariate Linear Model (MLM)

- **Approach**: Fits a single multivariate linear model incorporating all regulators simultaneously
- **Model**: For each sample, fits: `genes = β₁ × regulator₁_weights + β₂ × regulator₂_weights + ... + ε`
- **Assumption**: Accounts for the combined effect of multiple regulators and their interactions
- **Use case**: Better captures regulatory interactions and provides more accurate activity estimates when regulators are correlated

## Why MLM for This Pipeline?

We use MLM because:

1. **Better handles correlated regulators**: Transcription factors and pathways often regulate overlapping sets of genes. MLM accounts for these correlations.
2. **More accurate activity estimates**: By considering all regulators simultaneously, MLM reduces false positives from independent analysis.
3. **Captures interactions**: MLM can reveal how multiple regulators work together to control gene expression.

## Toy Example: 3 Genes and 2 Transcription Factors

Let's consider a simple example with:
- **3 genes**: G1, G2, G3
- **2 transcription factors**: TF1, TF2

### Regulatory Network

The regulatory network defines how TFs regulate genes:

```
TF1 → G1 (weight: 0.8)
TF1 → G2 (weight: 0.5)
TF2 → G2 (weight: 0.6)
TF2 → G3 (weight: 0.7)
```

Notice that both TF1 and TF2 regulate G2 (co-regulation).

### Observed Gene Expression

For a given cell, we observe:
- G1 expression: 2.0
- G2 expression: 3.5
- G3 expression: 1.8

### ULM Approach

ULM fits separate models for each TF:

**For TF1:**
- Model: `G1, G2 = β₁ × TF1_weights + ε`
- Using weights [0.8, 0.5] for G1 and G2
- Estimates TF1 activity independently

**For TF2:**
- Model: `G2, G3 = β₂ × TF2_weights + ε`
- Using weights [0.6, 0.7] for G2 and G3
- Estimates TF2 activity independently

**Problem**: Both TF1 and TF2 regulate G2, but ULM doesn't account for this overlap. Each TF's activity estimate for G2 is based on its own model, ignoring the other TF.

### MLM Approach

MLM fits a single multivariate model:

**Model**: `G1, G2, G3 = β₁ × TF1_weights + β₂ × TF2_weights + ε`

Where:
- TF1_weights = [0.8, 0.5, 0.0] (for G1, G2, G3)
- TF2_weights = [0.0, 0.6, 0.7] (for G1, G2, G3)

**Advantage**: MLM simultaneously estimates both TF1 and TF2 activities, accounting for the fact that G2 is regulated by both factors. This provides:

1. **More accurate estimates**: The model knows that G2's expression is influenced by both TFs
2. **Proper attribution**: If G2 is highly expressed, MLM can determine how much is due to TF1 vs TF2
3. **Interaction detection**: If TF1 and TF2 typically work together, MLM captures this relationship

### Mathematical Intuition

In the ULM case:
- TF1 activity estimate based on G1 and G2 expression, but G2 is also influenced by TF2
- This "contamination" can lead to inaccurate TF1 activity estimates
- Similar issue for TF2

In the MLM case:
- The model accounts for both TFs simultaneously
- It can "subtract out" the effect of TF2 when estimating TF1 activity
- Results in cleaner, more accurate activity estimates

## When to Use Each?

### Use ULM when:
- Regulators are known to be independent
- You want a quick, simple analysis
- You're exploring which regulators might be involved (initial screening)

### Use MLM when:
- Regulators may be correlated or interact
- You want accurate activity estimates
- You're doing downstream analysis requiring precise scores
- **This is the default for our pipeline**

## Implementation in decoupler

In decoupler:
- **ULM**: `dc.mt.ulm(data=adata, net=network)` → stores in `adata.obsm['score_ulm']`
- **MLM**: `dc.mt.mlm(data=adata, net=network)` → stores in `adata.obsm['score_mlm']`

Our pipeline uses MLM (`dc.mt.mlm`) for all enrichment scoring to ensure accurate and robust activity estimates.

## References

- [decoupler documentation](https://decoupler.readthedocs.io/)
- [decoupler MLM tutorial](https://decoupler.readthedocs.io/en/latest/notebooks/scell/rna_sc.html)

