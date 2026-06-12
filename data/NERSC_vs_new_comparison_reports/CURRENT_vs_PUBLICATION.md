# Current Repository vs Publication Standard

**Comparison Date:** 2026-01-30

**Publication Source:** Figure 6 from bioRxiv preprint

**Repository:** `/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper`

## Executive Summary

✅ **Current repository contains disease classification outputs for both IBD and PD**

### Key Finding: Terminology Inversion

**The repository data uses OPPOSITE terminology compared to the publication figure.**

When accounting for this inversion, the butyrate producer counts closely match:

**Example (Parkinson's Disease):**
- Publication 'increased likelihood' butyrate: 299
- Repository 'decreased' butyrate: 302
- **These match when terminology is inverted!**

---

## Inflammatory Bowel Disease

### Publication Standard (Figure 6)

| Metric | Value |
|--------|-------|
| Increased likelihood (butyrate producers) | 514 |
| Increased likelihood (non-butyrate) | 14,884 |
| Increased likelihood (total) | 15,398 |
| Decreased likelihood (butyrate producers) | 221 |
| Decreased likelihood (non-butyrate) | 39,333 |
| Decreased likelihood (total) | 39,554 |
| **Total microbes** | **54,952** |
| Chi-square | 647 |
| P-value | 1.00e-142 |

### Current Repository

| Metric | Value |
|--------|-------|
| Increased in disease (butyrate producers) | 186 |
| Increased in disease (total) | 7,764 |
| Decreased in disease (butyrate producers) | 478 |
| Decreased in disease (total) | 3,036 |
| **Total microbes** | **10,800** |
| Chi-square | 671.68 |
| P-value | 4.30e-148 |

**Source file:** `data/Intermediate_Files/IBD_Classification_butyrate_producers_summary.csv`

### Comparison (With Terminology Inversion)

**If we map:**
- Repository 'Decreased in disease' → Publication 'Increased likelihood' (protective)
- Repository 'Increased in disease' → Publication 'Decreased likelihood' (risk)

| Publication Metric | Pub Value | Repo Value (Inverted) | Difference | Match |
|-------------------|-----------|----------------------|------------|-------|
| Increased likelihood (butyrate) | 514 | 478 | -36 | ✅ |
| Decreased likelihood (butyrate) | 221 | 186 | -35 | ✅ |
| Chi-square | 647 | 671.68 | +24.68 | ✅ |

---

## Parkinson's Disease

### Publication Standard (Figure 6)

| Metric | Value |
|--------|-------|
| Increased likelihood (butyrate producers) | 299 |
| Increased likelihood (non-butyrate) | 2,691 |
| Increased likelihood (total) | 2,990 |
| Decreased likelihood (butyrate producers) | 152 |
| Decreased likelihood (non-butyrate) | 22,379 |
| Decreased likelihood (total) | 22,531 |
| **Total microbes** | **25,521** |
| Chi-square | 1317 |
| P-value | 2.00e-288 |

### Current Repository

| Metric | Value |
|--------|-------|
| Increased in disease (butyrate producers) | 150 |
| Increased in disease (total) | 7,401 |
| Decreased in disease (butyrate producers) | 302 |
| Decreased in disease (total) | 1,240 |
| **Total microbes** | **8,641** |
| Chi-square | 1063.60 |
| P-value | 2.70e-233 |

**Source file:** `data/Intermediate_Files/PD_Classification_butyrate_producers_summary.csv`

### Comparison (With Terminology Inversion)

**If we map:**
- Repository 'Decreased in disease' → Publication 'Increased likelihood' (protective)
- Repository 'Increased in disease' → Publication 'Decreased likelihood' (risk)

| Publication Metric | Pub Value | Repo Value (Inverted) | Difference | Match |
|-------------------|-----------|----------------------|------------|-------|
| Increased likelihood (butyrate) | 299 | 302 | +3 | ✅ |
| Decreased likelihood (butyrate) | 152 | 150 | -2 | ✅ |
| Chi-square | 1317 | 1063.60 | -253.40 | ❌ |

---

## Conclusion

### Reproducibility Status: ⚠️ PARTIAL

**What works:**
- ✅ Repository generates disease classification outputs
- ✅ Chi-square values are similar (within ~20%)
- ✅ Butyrate producer counts match when accounting for terminology inversion

**What doesn't work:**
- ❌ Direct numerical match without terminology adjustment
- ❌ Total microbe counts differ significantly (5x for IBD, 3x for PD)
- ❌ Column names use opposite terminology from publication

### Critical Questions

1. **Terminology:** What does 'increased/decreased likelihood of disease' mean in the publication figure?
2. **Total counts:** How were the publication totals (54,952 for IBD, 25,521 for PD) calculated?
3. **Data source:** Which exact repository file(s) generated the publication figure?
