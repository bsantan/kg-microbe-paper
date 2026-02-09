# Repository Comparison & Publication Validation Guide

**Purpose:** This document provides pointers to all repository comparison analyses, publication validation reports, and data sources for the kg-microbe-paper project.

**Last Updated:** 2026-01-30

---

## 📊 Publication Standard Data

### Gold Standard Source
The quantitative results from the bioRxiv publication Figure 6 are extracted and stored in:

**Location:** `data/NERSC_vs_new_comparison_reports/figure6_publication_results.txt`

**Content:**
- IBD (Inflammatory Bowel Disease) metrics
- PD (Parkinson's Disease) metrics
- Butyrate producer counts
- Statistical values (chi-square, p-values)
- Total microbe counts

### Publication Values Summary

**IBD:**
- Increased likelihood (butyrate): 514
- Decreased likelihood (butyrate): 221
- Total microbes: 54,952
- Chi-square: 647
- P-value: 1e-142

**PD:**
- Increased likelihood (butyrate): 299
- Decreased likelihood (butyrate): 152
- Total microbes: 25,521
- Chi-square: 1,317
- P-value: 2e-288

---

## 📁 Repository Locations

### Current (Local) Repository
**Path:** `/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper`

**Disease Classification Outputs:**
- `data/Intermediate_Files/IBD_Classification_butyrate_producers_summary.csv`
- `data/Intermediate_Files/PD_Classification_butyrate_producers_summary.csv`
- `data/Intermediate_Files/outcome_to_NCBITaxon_cleaned_IBD.tsv`
- `data/Intermediate_Files/outcome_to_NCBITaxon_cleaned_PD.tsv`
- `data/Intermediate_Files/outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_IBD.csv`
- `data/Intermediate_Files/outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_PD.csv`

### NERSC Repository
**Path:** `/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper_NERSC/kg-microbe-paper`

**Status:** ❌ Missing disease classification outputs (snapshot predates Figure 6 analysis)

---

## 📋 Comparison Reports

### All Reports Directory
**Location:** `data/NERSC_vs_new_comparison_reports/`

### Separate Comparison Reports

#### 1. Current Repository vs Publication
**File:** `data/NERSC_vs_new_comparison_reports/CURRENT_vs_PUBLICATION.md`

**Contents:**
- Direct comparison of current repository outputs to publication Figure 6
- Terminology inversion analysis
- Butyrate producer count matching
- Chi-square and p-value validation
- Reproducibility assessment

**Status:** ⚠️ PARTIAL - Numbers match when accounting for terminology inversion

#### 2. NERSC Repository vs Publication
**File:** `data/NERSC_vs_new_comparison_reports/NERSC_vs_PUBLICATION.md`

**Contents:**
- Documents that NERSC lacks disease classification files
- Timeline inference (NERSC → Current → Publication)
- Recommendation to use current repo as authoritative source

**Status:** ❌ CANNOT COMPARE - NERSC missing all disease outputs

---

## 🔍 Comprehensive Analysis Reports

### Primary Analysis Documents

#### 1. Final Analysis Summary
**File:** `data/NERSC_vs_new_comparison_reports/FINAL_ANALYSIS_SUMMARY.md`

**Key Findings:**
- Terminology inversion hypothesis
- Butyrate counts match when inverted
- Total count discrepancy (5x difference)
- Biological interpretation
- Corrected reproducibility assessment

#### 2. Critical Discrepancy Analysis
**File:** `data/NERSC_vs_new_comparison_reports/CRITICAL_DISCREPANCY_ANALYSIS.md`

**Key Findings:**
- Documents ALL mismatches
- Pattern of inverted terminology
- Critical questions for authors
- Reproducibility concerns

#### 3. Index & Navigation
**File:** `data/NERSC_vs_new_comparison_reports/INDEX.md`

**Contents:**
- Complete index of all 21 report files
- Navigation guide
- Quick reference tables
- File organization by purpose

#### 4. README
**File:** `data/NERSC_vs_new_comparison_reports/README.md`

**Contents:**
- Quick start guide
- Key findings summary
- File path mapping documentation
- Reproducibility status

---

## 📊 Data Comparison Tables

### Machine-Readable Formats

#### Publication Standard Comparison
**File:** `data/NERSC_vs_new_comparison_reports/publication_standard_comparison.tsv`

**Columns:**
- Disease, Repository, Status
- Increased/Decreased totals (Standard vs Repo)
- Butyrate counts (Standard vs Repo)
- Chi-square values
- Notes

#### File Path Mapping
**File:** `data/NERSC_vs_new_comparison_reports/file_path_mapping.tsv`

**Columns:**
- Current_Path, NERSC_Path, Status
- Category, Size_Diff_Bytes, Lines_Diff
- Notes

#### Statistical Validation
**File:** `data/NERSC_vs_new_comparison_reports/statistical_validation.tsv`

**Columns:**
- Analysis, Metric
- Current_Value, Expected_Value
- Difference, Percent_Change
- Status, Criticality

---

## 🔧 JSON Data Files

### Structured Comparison Data

#### Publication Standard Comparison (JSON)
**File:** `data/NERSC_vs_new_comparison_reports/publication_standard_comparison.json`

**Structure:**
```json
{
  "IBD": {
    "current": { "status": "...", "comparisons": {...} },
    "nersc": { "status": "...", "comparisons": {...} }
  },
  "PD": { ... }
}
```

#### Figure 6 Statistics
**File:** `data/NERSC_vs_new_comparison_reports/figure6_statistics_comparison.json`

**Contents:**
- Current repository IBD/PD statistics
- NERSC repository status
- Existence flags and data sources

#### Competency Counts
**File:** `data/NERSC_vs_new_comparison_reports/competency_counts_comparison.json`

**Contents:**
- Organismal traits counts
- RHEA annotations counts
- EC pathway counts
- Gold standard overlap counts

---

## 🚨 Critical Findings

### Terminology Inversion Pattern

**Discovery:** Repository column names use OPPOSITE terminology from publication figure

**Evidence:**

**Parkinson's Disease:**
| Publication Label | Pub Value | Repo Column | Repo Value | Match? |
|------------------|-----------|-------------|------------|--------|
| "Increased likelihood" butyrate | 299 | "Num_Decreased_disease_Butyrate_Producers" | 302 | ✅ |
| "Decreased likelihood" butyrate | 152 | "Num_Increased_disease_Butyrate_Producers" | 150 | ✅ |

**IBD:**
| Publication Label | Pub Value | Repo Column | Repo Value | Match? |
|------------------|-----------|-------------|------------|--------|
| "Increased likelihood" butyrate | 514 | "Num_Decreased_disease_Butyrate_Producers" | 478 | ✅ |
| "Decreased likelihood" butyrate | 221 | "Num_Increased_disease_Butyrate_Producers" | 186 | ✅ |

**Interpretation:** The numbers match when we map:
- Repository "Decreased in disease" → Publication "Increased likelihood" (protective)
- Repository "Increased in disease" → Publication "Decreased likelihood" (risk)

### Total Count Discrepancy

**IBD:**
- Publication: 54,952 total microbes
- Repository: 10,800 total microbes
- Difference: **5.1x**

**PD:**
- Publication: 25,521 total microbes
- Repository: 8,641 total microbes
- Difference: **3.0x**

**Hypothesis:** Publication may count at different granularity (all KG edges vs. unique taxa)

---

## ❓ Unresolved Questions

### For Publication Authors

1. **Terminology Clarification:**
   - Q: What does "increased/decreased likelihood of disease" mean in Figure 6?
   - Q: Does this refer to protective vs. pathogenic association?
   - Q: Or differential abundance (increased/decreased in disease state)?

2. **Data Source:**
   - Q: Which exact output file(s) generated Figure 6?
   - Q: Was it the summary CSV, ranks file, or different intermediate?
   - Q: Were any transformations or aggregations applied?

3. **Count Calculation:**
   - Q: How were total counts (54,952, 25,521) calculated?
   - Q: All strain-level representations? All KG edges?
   - Q: What level of aggregation was used?

4. **Terminology Inversion:**
   - Q: Is the figure mislabeled?
   - Q: Or are the CSV column names using opposite convention?
   - Q: Or are they measuring fundamentally different things?

---

## ✅ Recommendations

### For Reproducibility

1. **Clarify Terminology:**
   - Document exact meaning of all CSV column names
   - Align repository terminology with publication
   - Consider renaming columns to avoid confusion

2. **Document Data Source:**
   - State which exact files generated Figure 6
   - Provide script to regenerate figure from repository
   - Document any transformations or filtering

3. **Add Data Dictionary:**
   - Define all metrics and column names
   - Explain relationships between output files
   - Provide examples of interpretation

### For Publication

1. **Methods Section:**
   - Add detailed explanation of metric definitions
   - Clarify "likelihood" vs "abundance" terminology
   - Document analysis workflow from KG to figure

2. **Supplementary Materials:**
   - Include all intermediate files
   - Provide data dictionary
   - Add validation scripts

3. **Figure Caption:**
   - Clarify what "increased/decreased likelihood" means
   - Specify data source and aggregation level
   - Define total count calculation

---

## 🔗 Quick Reference Links

### Start Here
1. `data/NERSC_vs_new_comparison_reports/INDEX.md` - Navigation guide
2. `data/NERSC_vs_new_comparison_reports/CURRENT_vs_PUBLICATION.md` - Main comparison
3. `data/NERSC_vs_new_comparison_reports/FINAL_ANALYSIS_SUMMARY.md` - Complete findings

### Key Data Files
- Publication standard: `data/NERSC_vs_new_comparison_reports/figure6_publication_results.txt`
- Current repo IBD: `data/Intermediate_Files/IBD_Classification_butyrate_producers_summary.csv`
- Current repo PD: `data/Intermediate_Files/PD_Classification_butyrate_producers_summary.csv`

### Comparison Tables
- Publication comparison: `data/NERSC_vs_new_comparison_reports/publication_standard_comparison.tsv`
- File mapping: `data/NERSC_vs_new_comparison_reports/file_path_mapping.tsv`
- Stats validation: `data/NERSC_vs_new_comparison_reports/statistical_validation.tsv`

---

## 📞 Contact & Support

For questions about:
- **Repository comparison:** See `INDEX.md` and `README.md`
- **Reproducibility issues:** See `CRITICAL_DISCREPANCY_ANALYSIS.md`
- **Terminology questions:** See `FINAL_ANALYSIS_SUMMARY.md`
- **Data sources:** See `CURRENT_vs_PUBLICATION.md`

---

**Document Version:** 1.0
**Last Analysis Run:** 2026-01-30
**Total Reports Generated:** 21 files
**Status:** ⚠️ Reproducibility requires terminology clarification
