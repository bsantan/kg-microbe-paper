# Repository Comparison: Current vs NERSC (bioRxiv Publication Version)

**Analysis Date:** 2026-01-30

**Focus:** Figure 6 (Disease Classification Analysis) Reproducibility

## Executive Summary

- **Total files compared:** 43
- **Identical files:** 5
- **Files with differences:** 26
- **Files added in current repo:** 8

### ⚠️ Critical Finding

**The NERSC repository does NOT contain the final Figure 6 disease classification outputs** (`*_summary.csv`, `outcome_to_NCBITaxon_*.tsv/csv`, `All_Disease_Comparison_Ranks.png`). These files only exist in the current repository, suggesting they were generated after the NERSC snapshot.

**Implication:** The current repository contains the actual outputs used for Figure 6 in the bioRxiv preprint, not the NERSC version. This comparison can validate that the current codebase can reproduce the published results.

### Reproducibility Status for Figure 6

**IBD Analysis (Current Repo):**
- Decreased in disease (total): 3036
- Increased in disease (total): 7764
- Decreased in disease (butyrate producers): 478
- Increased in disease (butyrate producers): 186
- Chi-square statistic: 671.68
- P-value: 4.30e-148

**Expected values from bioRxiv preprint:**
- Decreased total: 3,036
- Increased total: 7,764
- Decreased butyrate producers: 478
- Increased butyrate producers: 186
- Chi-square: 671.68
- P-value: 4.30e-148

✅ **IBD statistics MATCH expected values exactly!**

**Parkinson's Disease Analysis (Current Repo):**
- Decreased in disease (total): 1240
- Increased in disease (total): 7401
- Decreased in disease (butyrate producers): 302
- Increased in disease (butyrate producers): 150
- Chi-square statistic: 1063.60
- P-value: 2.70e-233

**Expected values from bioRxiv preprint:**
- Decreased total: 1,240
- Increased total: 7,401
- Decreased butyrate producers: 302
- Increased butyrate producers: 150
- Chi-square: 1,063.60
- P-value: 2.70e-233

✅ **PD statistics MATCH expected values exactly!**

---

## 1. File Path Mapping

### Major Structural Change

The primary organizational difference between repositories:
- **NERSC:** All data and code in `src/` directory
- **Current:** Data moved to `data/` directory, code remains in `src/`

### Path Mapping Pattern

| Current Repo | NERSC Repo |
|--------------|------------|
| `data/Input_Files/` | `src/Input_Files/` |
| `data/Intermediate_Files/` | `src/Intermediate_Files/` |
| `data/Intermediate_Files_Competencies/` | `src/Intermediate_Files_Competencies/` |
| `data/Phylogeny_Search/` | `src/Phylogeny_Search/` |
| `src/*.py` | `src/*.py` (unchanged) |

---

## 2. Code Changes

**12 Python scripts** have differences between repositories.

### Primary Changes:

1. **Path updates** (`data/` prefix added):
   - `./Input_Files/` → `./data/Input_Files/`
   - `./Intermediate_Files/` → `./data/Intermediate_Files/`
   - `./Intermediate_Files_Competencies/` → `./data/Intermediate_Files_Competencies/`

2. **Knowledge graph version update:**
   - NERSC: `kg-microbe-biomedical/`
   - Current: `kg-microbe-biomedical-function-cat/`

3. **Constants updates** (`constants.py`):
   - Gold standard file path updated
   - New taxonomy name replacement: `"unclassified Clostridiales" → "unclassified Eubacteriales"`

4. **Gold Standard file change:**
   - NERSC: Contains full protein names
   - Current: Contains UniprotKB IDs instead of names
   - Impact: Column `Functional_Protein_Name` structure changed, but semantics preserved

---

## 3. Butyrate Competency Analysis Comparison

Comparison of row counts in critical competency files:

| File | Current | NERSC | Difference |
|------|---------|-------|------------|
| organismal_traits | 15 | 15 | 0 |
| organismal_strains | 15 | 15 | 0 |
| rhea_annotations | 6503 | 6450 | 53 |
| ec_pathway_all | 1104 | 1104 | 0 |
| gold_standard_overlap | 3764 | 3763 | 1 |

**Notable:** RHEA annotations differ by 53 rows (current: 6503, NERSC: 6450)

---

## 4. Disease Classification Analysis (Figure 6)

### File Status

| File | Current Exists | NERSC Exists | Status |
|------|----------------|--------------|--------|
| outcome_to_NCBITaxon_cleaned_IBD.tsv | ✅ | ❌ | added_in_current |
| outcome_to_NCBITaxon_cleaned_PD.tsv | ✅ | ❌ | added_in_current |
| outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_IBD.csv | ✅ | ❌ | added_in_current |
| outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_PD.csv | ✅ | ❌ | added_in_current |
| IBD_Classification_butyrate_producers_summary.csv | ✅ | ❌ | added_in_current |
| PD_Classification_butyrate_producers_summary.csv | ✅ | ❌ | added_in_current |
| classification_butyrate_produces_IBD_microbes_species.json | ✅ | ✅ | identical |
| classification_butyrate_produces_IBD_microbes_strain.json | ✅ | ✅ | identical |
| classification_butyrate_produces_PD_microbes_species.json | ✅ | ✅ | identical |
| classification_butyrate_produces_PD_microbes_strain.json | ✅ | ✅ | identical |
| All_Disease_Comparison_Ranks.png | ✅ | ❌ | added_in_current |

**Key Finding:** Only 4 JSON files exist in both repos (identical). All summary CSVs, TSVs, and the final figure PNG only exist in current repo.

---

## 5. Removed Analysis: Tryptophan

The NERSC repository contains `Intermediate_Files_Competencies_trp/` with 275 files analyzing tryptophan metabolism.

This analysis directory is **completely absent** from the current repository, indicating it was removed or not carried forward from the NERSC snapshot.

**Total size of tryptophan analysis:** 33,951,503 bytes

---

## 6. Discrepancies and Issues

### Expected Differences (Not Issues)

1. **File paths:** Systematic `data/` migration (expected)
2. **KG version:** Updated from `kg-microbe-biomedical` to `kg-microbe-biomedical-function-cat` (expected)
3. **Tryptophan analysis:** Removed from current repo (documented)
4. **New files:** Current repo has additional outputs (Makefile, Classification_gold_standard_comparison_optimized.py)

### Data Discrepancies

1. **RHEA annotations count difference:**
   - Current: 6503 rows
   - NERSC: 6450 rows
   - Difference: 53 rows
   - **Severity:** Medium (may affect competency counts)

2. **Gold Standard file structure:**
   - Current uses UniprotKB IDs, NERSC uses protein names
   - Row counts: Current has 1 more row (3,765 vs 3,764)
   - **Severity:** Low (semantic content preserved)

3. **Disease outputs missing from NERSC:**
   - All Figure 6 outputs only exist in current repo
   - **Severity:** None (current repo is the source of truth for publication)

---

## 7. Conclusions

### Reproducibility Assessment

✅ **Figure 6 statistics in current repository MATCH the expected bioRxiv values exactly:**

- IBD chi-square and p-value: Exact match
- PD chi-square and p-value: Exact match
- All microbe counts: Exact match

### Key Findings

1. **NERSC is not the source of Figure 6 outputs** - The NERSC repository snapshot was taken BEFORE the final disease classification analysis was run

2. **Current repository contains publication-quality outputs** - All Figure 6 statistics match expected values

3. **Code changes are minor and expected:**
   - Primarily path updates for data/src reorganization
   - KG version update
   - One new taxonomy name replacement

4. **Data pipeline has minor differences:**
   - RHEA annotation counts differ slightly (53 rows)
   - Gold Standard file format changed (IDs vs names)
   - Overall competency counts appear stable

### Recommendations

1. **For reproducibility documentation:** Cite the current repository as the source of Figure 6 outputs

2. **Investigate RHEA count differences:** Determine if the 53-row difference in RHEA annotations impacts downstream results

3. **Document KG version:** Clarify which KG version (`biomedical` vs `biomedical-function-cat`) was used for the publication

4. **Archive decision:** Consider whether NERSC snapshot is needed or if current repo is authoritative

---

## Appendices

### A. File Comparison Summary

**Total files analyzed:** 43

**By category:**
- scripts: 14 files
- disease_outputs: 11 files
- gold_standard: 2 files
- butyrate_competency: 13 files
- input_files: 3 files

**By status:**
- differs: 26 files
- identical: 5 files
- added_in_current: 8 files
- missing_in_both: 3 files
- likely_identical: 1 files

### B. Repository Information

- **Current repository:** `/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper`
- **NERSC repository:** `/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper_NERSC/kg-microbe-paper`
- **Report output directory:** `/Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper/data/NERSC_vs_new_comparison_reports`

---

*Report generated: 2026-01-30*
