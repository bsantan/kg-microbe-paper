# Comparison to Publication Standard (Figure 6)

**Date:** 2026-01-30

This report compares both the current and NERSC repository results against the publication standard values extracted from Figure 6.

## Executive Summary

### Inflammatory Bowel Disease

- **Current Repository:** DIFFERS
- **NERSC Repository:** MISSING

### Parkinson's Disease

- **Current Repository:** DIFFERS
- **NERSC Repository:** MISSING

---

## Detailed Comparison

### Inflammatory Bowel Disease

#### Publication Standard (Figure 6)

| Metric | Value |
|--------|-------|
| Increased likelihood (total) | 15,398 |
| Increased likelihood (butyrate producers) | 514 |
| Decreased likelihood (total) | 39,554 |
| Decreased likelihood (butyrate producers) | 221 |
| Total microbes | 54,952 |
| Chi-square | 647 |
| P-value | 1.00e-142 |

#### Current Repository

| Metric | Standard | Repository | Difference | Match |
|--------|----------|------------|------------|-------|
| Increased Total | 15,398 | 7,764 | -7,634 | ❌ |
| Increased Butyrate | 514 | 186 | -328 | ❌ |
| Decreased Total | 39,554 | 3,036 | -36,518 | ❌ |
| Decreased Butyrate | 221 | 478 | +257 | ❌ |
| Chi-square | 647 | 671.68 | 24.68 | ❌ |
| P-value | 1.00e-142 | 4.30e-148 | - | ❌ |

**Source:** `data/Intermediate_Files/IBD_Classification_butyrate_producers_summary.csv`

#### NERSC Repository

❌ **No data found in repository**

---

### Parkinson's Disease

#### Publication Standard (Figure 6)

| Metric | Value |
|--------|-------|
| Increased likelihood (total) | 2,990 |
| Increased likelihood (butyrate producers) | 299 |
| Decreased likelihood (total) | 22,531 |
| Decreased likelihood (butyrate producers) | 152 |
| Total microbes | 25,521 |
| Chi-square | 1317 |
| P-value | 2.00e-288 |

#### Current Repository

| Metric | Standard | Repository | Difference | Match |
|--------|----------|------------|------------|-------|
| Increased Total | 2,990 | 7,401 | +4,411 | ❌ |
| Increased Butyrate | 299 | 150 | -149 | ❌ |
| Decreased Total | 22,531 | 1,240 | -21,291 | ❌ |
| Decreased Butyrate | 152 | 302 | +150 | ❌ |
| Chi-square | 1317 | 1063.60 | -253.40 | ❌ |
| P-value | 2.00e-288 | 2.70e-233 | - | ❌ |

**Source:** `data/Intermediate_Files/PD_Classification_butyrate_producers_summary.csv`

#### NERSC Repository

❌ **No data found in repository**

---

## Interpretation

### Key Findings

**Note:** The repository data uses different metrics than the publication figure.

The publication figure shows:
- Microbes with **increased/decreased likelihood of disease** (epidemiological association)
- Total counts include all strain-level representations

The repository CSV files show:
- Microbes with **increased/decreased abundance in disease state** (differential abundance)
- Potentially different aggregation or filtering

This explains why the numerical values differ significantly between the publication figure and the repository files.

### Recommendations

1. **Clarify metric definitions:** Document whether the analysis measures:
   - Disease association (likelihood)
   - Differential abundance (increased/decreased in disease)

2. **Verify data source:** Confirm which data files were used to generate the publication figure

3. **Document analysis pipeline:** Clarify any filtering, aggregation, or transformation steps
