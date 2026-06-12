# NERSC Repository vs Publication Standard

**Comparison Date:** 2026-01-30

**Publication Source:** Figure 6 from bioRxiv preprint

**Repository:** `/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper_NERSC/kg-microbe-paper`

## Executive Summary

### ❌ CRITICAL FINDING: NERSC Repository Missing Disease Classification Outputs

**The NERSC repository snapshot does NOT contain any disease classification output files.**

**Missing files for both IBD and PD:**
- Summary CSV files
- Outcome TSV files
- Ranks CSV files
- Classification JSON files

**Implication:** The NERSC snapshot was taken BEFORE the disease classification analysis (Figure 6) was performed.

---

## Inflammatory Bowel Disease

### Publication Standard (Figure 6)

| Metric | Value |
|--------|-------|
| Increased likelihood (butyrate producers) | 514 |
| Increased likelihood (total) | 15,398 |
| Decreased likelihood (butyrate producers) | 221 |
| Decreased likelihood (total) | 39,554 |
| **Total microbes** | **54,952** |
| Chi-square | 647 |
| P-value | 1.00e-142 |

### NERSC Repository

❌ **No data found in NERSC repository**

**Searched locations:**
- `src/Intermediate_Files/IBD_Classification_butyrate_producers_summary.csv`
- `src/Intermediate_Files_but/IBD_Classification_butyrate_producers_summary.csv`

**Status:** File does not exist in NERSC snapshot

---

## Parkinson's Disease

### Publication Standard (Figure 6)

| Metric | Value |
|--------|-------|
| Increased likelihood (butyrate producers) | 299 |
| Increased likelihood (total) | 2,990 |
| Decreased likelihood (butyrate producers) | 152 |
| Decreased likelihood (total) | 22,531 |
| **Total microbes** | **25,521** |
| Chi-square | 1317 |
| P-value | 2.00e-288 |

### NERSC Repository

❌ **No data found in NERSC repository**

**Searched locations:**
- `src/Intermediate_Files/PD_Classification_butyrate_producers_summary.csv`
- `src/Intermediate_Files_but/PD_Classification_butyrate_producers_summary.csv`

**Status:** File does not exist in NERSC snapshot

---

## Conclusion

### Status: ❌ CANNOT COMPARE

**The NERSC repository snapshot does not contain the disease classification analysis outputs.**

This means:
1. NERSC snapshot was taken BEFORE Figure 6 analysis was performed
2. Current repository is the ONLY source of disease classification outputs
3. NERSC cannot be used to validate publication results
4. The current repository IS the source of the publication data

### Timeline Inference

Based on file presence:

```
NERSC Snapshot
     |
     v
[Contains: butyrate competency analysis]
[Missing: disease classification]
     |
     | Disease classification
     | analysis performed
     v
Current Repository
[Contains: all analyses including Figure 6]
     |
     v
Publication Figure 6 generated
```

### Recommendation

**Use the current repository** as the authoritative source for:
- Disease classification analysis
- Figure 6 generation
- All bioRxiv preprint outputs

**NERSC snapshot can be archived or removed** as it does not contain the publication-critical analysis.
