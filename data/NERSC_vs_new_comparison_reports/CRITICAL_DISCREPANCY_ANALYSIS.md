# CRITICAL DISCREPANCY: Repository Data Does NOT Match Publication Figure

> ⚠️ **PARTIALLY RESOLVED — See [MANUSCRIPT_REFERENCE.md](MANUSCRIPT_REFERENCE.md)** (2026-05-22). After fixes `716066f` (strain-dict contamination) and `82912ac` (per-disease results file), the current repo run matches the manuscript-acknowledged PD reproduction (1337 / 8e-293) and is within ~4% of the paper IBD value. The "5x total" discrepancy this file describes was caused by a contamination short-circuit in `create_species_strains_dictionary` that has since been fixed.

**Date:** 2026-01-30
**Status:** ❌ **MAJOR REPRODUCIBILITY ISSUE IDENTIFIED**

---

## CORRECTION TO EARLIER ANALYSIS

**My earlier statement was WRONG:**
> "✅ Critical Finding: 100% Reproducibility Confirmed
> All Figure 6 statistics match the bioRxiv publication values exactly"

This was completely incorrect. The repository data and publication figure show **dramatically different values**.

---

## The Actual Discrepancies

### Inflammatory Bowel Disease (IBD)

| Metric | Publication Figure | Current Repo | Difference |
|--------|-------------------|--------------|------------|
| **Increased (total)** | 15,398 | 7,764 | **-7,634** (-50%) |
| **Increased (butyrate)** | 514 | 186 | **-328** (-64%) |
| **Decreased (total)** | 39,554 | 3,036 | **-36,518** (-92%) |
| **Decreased (butyrate)** | 221 | 478 | **+257** (+116%) |
| **Total microbes** | 54,952 | 10,800 | **-44,152** (-80%) |
| **Chi-square** | 647 | 671.68 | +24.68 (+4%) |
| **P-value** | 1.0e-142 | 4.3e-148 | Similar magnitude |

### Parkinson's Disease (PD)

| Metric | Publication Figure | Current Repo | Difference |
|--------|-------------------|--------------|------------|
| **Increased (total)** | 2,990 | 7,401 | **+4,411** (+148%) |
| **Increased (butyrate)** | 299 | 150 | **-149** (-50%) |
| **Decreased (total)** | 22,531 | 1,240 | **-21,291** (-94%) |
| **Decreased (butyrate)** | 152 | 302 | **+150** (+99%) |
| **Total microbes** | 25,521 | 8,641 | **-16,880** (-66%) |
| **Chi-square** | 1,317 | 1,063.60 | -253.40 (-19%) |
| **P-value** | 2.0e-288 | 2.7e-233 | Different magnitude |

---

## The Numbers Are Inverted!

### Critical Pattern Identified:

**IBD:**
- Publication: Decreased butyrate = **221**
- Repository: Decreased butyrate = **478**
- Publication: Increased butyrate = **514**
- Repository: Increased butyrate = **186**

**The butyrate producer counts are OPPOSITE in magnitude!**

**PD:**
- Publication: Decreased butyrate = **152**
- Repository: Decreased butyrate = **302**
- Publication: Increased butyrate = **299**
- Repository: Increased butyrate = **150**

**Again, the counts are inverted!**

---

## Possible Explanations

### Hypothesis 1: Different Metric Definitions

**Publication Figure** may show:
- "Increased/Decreased **likelihood of disease**" (epidemiological association from literature)
- Counts based on strain-level KG representation
- Total: 54,952 for IBD, 25,521 for PD

**Repository CSV** may show:
- "Increased/Decreased **in disease state**" (differential abundance analysis)
- Counts based on species-level or filtered data
- Total: 10,800 for IBD, 8,641 for PD

### Hypothesis 2: Data Transformation or Filtering

The repository may apply:
- Taxonomic rank aggregation (species vs strain)
- Confidence thresholds
- Statistical filtering
- Different edge types from KG

### Hypothesis 3: Terminology Confusion

Possible that:
- "Decreased_disease" in CSV = decreased abundance in disease (depleted, protective)
- "Decreased likelihood" in figure = protective association from literature
- These could be measuring **opposite** directions!

### Hypothesis 4: Different Analysis Runs

The publication figure may represent:
- An earlier analysis version
- Different KG version
- Different data sources
- Different statistical approach

---

## Evidence for Inverted Terminology

Looking at the pattern:

**If we SWAP the repository's "increased" and "decreased" labels:**

IBD (swapped):
- Decreased butyrate: 186 (vs publication 221) - closer!
- Increased butyrate: 478 (vs publication 514) - closer!

PD (swapped):
- Decreased butyrate: 150 (vs publication 152) - very close!
- Increased butyrate: 302 (vs publication 299) - very close!

**This suggests the repository might use OPPOSITE terminology for the direction!**

But the total counts still don't match (10,800 vs 54,952), suggesting additional differences in:
- Aggregation level (species vs strain)
- Filtering criteria
- Data source

---

## What Chi-Square Values Tell Us

Interestingly, chi-square values are relatively close:
- IBD: 647 (publication) vs 671.68 (repository) = 4% difference
- PD: 1,317 (publication) vs 1,063.60 (repository) = 19% difference

This suggests:
- Both analyses show **similar statistical strength**
- The **association pattern** is similar
- But the **underlying data** is different

---

## NERSC Repository Status

**NERSC has NO disease classification output files** (IBD or PD summaries)

This confirms:
- NERSC snapshot was taken BEFORE final disease analysis
- Current repository is the only source of these specific CSV files
- We cannot use NERSC to validate against publication

---

## Critical Questions That Need Answering

1. **Which data generated the publication figure?**
   - Was it from a different analysis run?
   - Different scripts not in the repository?
   - Different intermediate files?

2. **What do the CSV column names actually mean?**
   - "Num_Decreased_disease_Total" = decreased in disease? Or protective?
   - "Num_Increased_disease_Total" = increased in disease? Or risk factor?

3. **Why are total counts so different?**
   - Publication: 54,952 (IBD), 25,521 (PD)
   - Repository: 10,800 (IBD), 8,641 (PD)
   - This is a **5x difference** for IBD!

4. **Is there another output file we should be comparing to?**
   - Maybe the publication used `outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_{disease}.csv`?
   - Or different aggregation?

---

## Recommended Next Steps

### Immediate:
1. ✅ **Document this discrepancy** (this file)
2. ⚠️ **Search for alternative output files** that might match publication
3. ⚠️ **Check if figure was generated from different scripts**
4. ⚠️ **Verify column name interpretations** with original authors

### Investigation:
1. Read the `outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_*.csv` files
2. Check if they contain the publication values
3. Examine the `classification_butyrate_produces_{disease}_microbes_*.json` files
4. Review the classification scripts for terminology mapping

### Documentation:
1. Clarify in paper methods which exact files produced the figure
2. Define what "increased/decreased" means in context
3. Document any filtering or aggregation steps
4. Provide data dictionary for all output files

---

## Current Status Summary

| Aspect | Status |
|--------|--------|
| **Reproducibility** | ❌ **FAILED** - Numbers don't match |
| **Data Availability** | ✅ Current repo has data (NERSC missing) |
| **Statistical Significance** | ✅ Both show strong associations |
| **Biological Conclusion** | ⚠️ **Unclear** - depends on metric interpretation |
| **Code Availability** | ✅ Scripts present and runnable |
| **Transparency** | ❌ **Gap** - unclear which outputs = publication |

---

## Bottom Line

**We CANNOT confirm reproducibility** because:
1. Repository values ≠ Publication figure values
2. Differences are too large to be rounding errors
3. Even the directionality is questionable (inverted counts)
4. Total counts differ by 5-fold

**Further investigation required** to:
- Locate the actual data source for the publication figure
- Clarify metric definitions and terminology
- Reconcile the dramatic differences in counts

**This is a significant reproducibility concern** that needs to be resolved before publication.

---

*Analysis by: Claude Code comprehensive comparison*
*Date: 2026-01-30*
*Repository: kg-microbe-paper*
