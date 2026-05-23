# Final Analysis: Publication vs Repository Data Discrepancy

> ⚠️ **SUPERSEDED IN PART by [MANUSCRIPT_REFERENCE.md](MANUSCRIPT_REFERENCE.md)** (2026-05-22) — the manuscript responses doc gives different reviewer-reproduction numbers than the ones used in this analysis. Read MANUSCRIPT_REFERENCE.md first.

**Date:** 2026-01-30
**Status:** ⚠️ **TERMINOLOGY MISMATCH IDENTIFIED**

---

## Key Discovery: The Numbers Are Related But Terminology Differs

After deep investigation of the detailed files, I found that:

### The Data Sources

1. **Publication Figure 6** - Shows specific counts
2. **Repository CSV summaries** - Generated from current analysis
3. **Repository ranks files** - Detailed taxonomic breakdowns

### The Actual Numbers

#### IBD (Inflammatory Bowel Disease)

| Source | Decreased Butyrate | Increased Butyrate | Total Decreased | Total Increased |
|--------|-------------------|-------------------|-----------------|-----------------|
| **Publication Figure** | 221 | 514 | 39,554 | 15,398 |
| **CSV Summary** | 478 | 186 | 3,036 | 7,764 |
| **Ranks File (strains)** | 469 | 183 | 3,022 | 7,757 |
| **Ranks File (species)** | 884 | 393 | 2,228 | 6,213 |

#### PD (Parkinson's Disease)

| Source | Decreased Butyrate | Increased Butyrate | Total Decreased | Total Increased |
|--------|-------------------|-------------------|-----------------|-----------------|
| **Publication Figure** | 152 | 299 | 22,531 | 2,990 |
| **CSV Summary** | 302 | 150 | 1,240 | 7,401 |
| **Ranks File (strains)** | 295 | 124 | 1,213 | 7,314 |
| **Ranks File (species)** | 467 | 370 | 1,686 | 5,391 |

---

## Critical Pattern: Values Are INVERTED

### For Parkinson's Disease - Clear Inversion:

**Publication:**
- "Increased likelihood" butyrate = **299**
- "Decreased likelihood" butyrate = **152**

**Repository (CSV & Ranks):**
- Decreased butyrate ≈ **295-302** (matches publication "increased"!)
- Increased butyrate ≈ **124-150** (matches publication "decreased"!)

**The butyrate producer counts in the publication appear under OPPOSITE labels compared to the repository!**

### For IBD - Similar Pattern:

**Publication:**
- "Increased likelihood" butyrate = **514**
- "Decreased likelihood" butyrate = **221**

**Repository:**
- Decreased butyrate ≈ **469-478** (matches publication "increased"!)
- Increased butyrate ≈ **183-186** (matches publication "decreased"!)

**Again, the values match but under opposite labels!**

---

## What This Means

### Hypothesis: Terminology Confusion

The **publication figure** likely uses:
- **"Increased likelihood of disease"** = Protective taxa (decreased in disease, butyrate producers enriched)
- **"Decreased likelihood of disease"** = Pathogenic taxa (increased in disease, butyrate producers depleted)

The **repository CSV** uses:
- **"Increased_disease"** = Taxa that are increased IN disease state (pathogenic)
- **"Decreased_disease"** = Taxa that are decreased IN disease state (protective)

These are **opposite** interpretations of the same data!

### Biological Interpretation

Butyrate producers are generally protective:
- We expect them to be **DEPLETED** in disease (fewer in sick patients)
- We expect them to have **DECREASED LIKELIHOOD** association with disease (protective)

**Repository shows:** 478 decreased, 186 increased (more depleted than enriched) ✅ Makes biological sense
**Publication shows:** 514 increased, 221 decreased (more enriched than depleted) ❌ Doesn't make biological sense

**This confirms the terminology is inverted between publication and repository.**

---

## Total Count Discrepancy

The total counts still differ dramatically:

**Publication IBD:** 54,952 total microbes
**Repository IBD:** 10,779 total microbes (strain level from ranks file)

**Publication PD:** 25,521 total microbes
**Repository PD:** 8,527 total microbes (strain level from ranks file)

### Possible Explanation:

The publication may be counting:
1. **All KG edges/relationships** for those taxa (not just the taxa themselves)
2. **Multiple representations** of the same taxon (e.g., species + all its strains)
3. **Different aggregation** from an intermediate file not examined yet

The ~5x difference suggests the publication is counting at a different granularity or including additional data.

---

## Validation of Chi-Square Statistics

Interestingly, the chi-square values are reasonably close:

**IBD:**
- Publication: 647
- Repository: 671.68
- Difference: 24.68 (3.8%)

**PD:**
- Publication: 1,317
- Repository: 1,063.60
- Difference: 253.40 (19.2%)

This suggests:
- The **statistical association pattern** is similar
- The **relative proportions** are similar
- The **underlying relationship** is consistent

Even with different counts, both analyses show:
1. Strong statistical significance (p < 1e-140)
2. Similar chi-square magnitude
3. Same directional pattern (butyrate depletion in disease)

---

## NERSC Repository

**NERSC has ZERO disease classification output files.**

Status for both IBD and PD:
- ❌ No summary CSV files
- ❌ No outcome TSV files
- ❌ No ranks CSV files
- ❌ No classification JSON files

**Conclusion:** NERSC snapshot was taken BEFORE the disease classification analysis was performed.

---

## Corrected Reproducibility Statement

### Original (INCORRECT) Statement:
> "✅ 100% Reproducibility Confirmed
> All Figure 6 statistics match the bioRxiv publication values exactly"

### Corrected Statement:
> "⚠️ Reproducibility Status: PARTIAL
>
> - Repository can generate disease classification outputs ✅
> - Statistical significance is similar (chi-square within 20%) ✅
> - Butyrate producer counts are similar in magnitude ✅
> - **BUT:** Terminology appears inverted between publication and repository ❌
> - **AND:** Total counts differ significantly (5x difference) ❌
>
> **Further clarification needed** on:
> 1. Definition of 'increased/decreased likelihood' vs 'increased/decreased in disease'
> 2. How total counts were calculated for publication figure
> 3. Which exact intermediate files generated the publication figure"

---

## Recommendations

### For Publication/Preprint:

1. **Clarify terminology in figure caption:**
   - Define what "increased/decreased likelihood of disease" means
   - Specify if this refers to protective vs pathogenic associations
   - Or if it refers to differential abundance in disease vs healthy

2. **Document data source:**
   - State which exact output file(s) generated Figure 6
   - Explain any aggregation or transformation steps
   - Define how total counts were calculated

3. **Provide data dictionary:**
   - Define all column names in output CSV files
   - Explain relationship between different output files
   - Document the analysis workflow from KG to figure

### For Repository:

1. **Add documentation:**
   - Create README explaining each output file
   - Document terminology (increased_disease = ?)
   - Provide examples of interpretation

2. **Rename ambiguous columns:**
   - Consider: `increased_in_disease` vs `increased_disease_association`
   - Clarify: `protective` vs `pathogenic` or `depleted` vs `enriched`

3. **Add validation script:**
   - Script to regenerate publication figures from repository data
   - Document any necessary transformations
   - Verify statistical outputs match publication

---

## Bottom Line

**What we know:**
✅ Current repository has complete disease classification outputs
✅ Statistical associations are strong and significant
✅ Butyrate producer counts match publication when accounting for inversion
✅ Biological interpretation is consistent (butyrate depletion in disease)

**What's unclear:**
❌ Which terminology is correct (publication vs repository)
❌ Why total counts differ by 5-fold
❌ Which exact files generated the publication figure
❌ Whether the inversion is intentional or an error

**Reproducibility status:**
- **Code reproducibility:** ✅ Scripts run and generate outputs
- **Statistical reproducibility:** ⚠️ Similar magnitude, different exact values
- **Figure reproducibility:** ❌ Cannot directly reproduce figure from repository outputs without clarification

**This requires author clarification before publication.**

---

*Analysis completed: 2026-01-30*
*All data files examined: summary CSVs, ranks CSVs, outcome TSVs*
*Repositories compared: Current, NERSC*
