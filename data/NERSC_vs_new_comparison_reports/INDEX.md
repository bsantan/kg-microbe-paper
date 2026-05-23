# Complete Index of Repository Comparison Reports

**Generated:** 2026-01-30
**Total Reports:** 19 files

> ⚠️ **READ FIRST: [MANUSCRIPT_REFERENCE.md](MANUSCRIPT_REFERENCE.md)** (added 2026-05-22)
>
> The reports below predate the second reproducibility review's conclusion and use second-hand transcribed Figure 6C numbers. The authoritative source is `revisions2/KG-Microbe_Responses2_mpj2.docx`, which states PD reproduction χ² = 1337 / p = 8e-293 (vs paper 1317 / 2e-288) and describes IBD as "consistent" between paper and reviewer reproduction. The "Reviewer #4 = IBD 672 / PD 1064" claim that appears in several files in this directory is **not consistent** with the manuscript's stated reviewer reproduction value.

---

## 🚨 START HERE - Critical Findings

**Read these first to understand the major issues:**

1. **[FINAL_ANALYSIS_SUMMARY.md](FINAL_ANALYSIS_SUMMARY.md)** ⭐ **MOST IMPORTANT**
   - Complete investigation results
   - Terminology inversion hypothesis
   - Biological interpretation
   - Corrected reproducibility assessment

2. **[CRITICAL_DISCREPANCY_ANALYSIS.md](CRITICAL_DISCREPANCY_ANALYSIS.md)** ⚠️ **KEY ISSUE**
   - Documents that repository ≠ publication
   - Shows dramatic count differences
   - Identifies inverted terminology pattern
   - Lists critical questions for authors

3. **[README.md](README.md)** 📖 **OVERVIEW**
   - Quick start guide
   - File organization
   - Navigation help

---

## 📊 Publication Standard Comparison (New Analysis)

These files compare repository data to the publication figure values:

4. **[figure6_publication_results.txt](figure6_publication_results.txt)**
   - Extracted quantitative values from publication Figure 6
   - Gold standard reference data
   - Used for all comparisons

5. **[publication_standard_comparison_report.md](publication_standard_comparison_report.md)**
   - Formal comparison report
   - Side-by-side comparison tables
   - Current vs NERSC vs Publication
   - Interpretation and recommendations

6. **[publication_standard_comparison.tsv](publication_standard_comparison.tsv)**
   - Machine-readable comparison table
   - All metrics in structured format
   - Status flags for each comparison

7. **[publication_standard_comparison.json](publication_standard_comparison.json)**
   - Complete comparison data with metadata
   - Includes all calculations
   - Detailed comparison objects

---

## 📝 Repository-to-Repository Comparison (Original Analysis)

These files compare current repository to NERSC snapshot:

8. **[comparison_report.md](comparison_report.md)**
   - Original comprehensive narrative report
   - File path mapping documentation
   - Code changes analysis
   - ⚠️ Contains incorrect "100% reproducibility" claim (corrected in newer files)

9. **[comparison_table.tsv](comparison_table.tsv)**
   - Master comparison table (43 files)
   - Current vs NERSC status for each file
   - Organized by analysis pipeline step

10. **[statistical_validation.tsv](statistical_validation.tsv)**
    - Validates Figure 6 metrics
    - Compares current repo to expected values
    - Shows ALL metrics MATCH (IBD and PD)
    - ⚠️ But these are DIFFERENT from publication figure!

---

## 🔧 Code and Path Analysis

11. **[file_path_mapping.tsv](file_path_mapping.tsv)**
    - Complete path mapping: current ↔ NERSC
    - Documents data/ → src/ migration
    - Status for each file pair
    - Size and line count differences

12. **[code_diff_summary.tsv](code_diff_summary.tsv)**
    - Summary of code changes by file
    - Categorization: path changes, logic changes
    - Impact assessment (high/medium/low)
    - Key changes description

13. **[code_comparison.tsv](code_comparison.tsv)**
    - Full diff output for all changed code files
    - Line-by-line changes
    - 441 lines of detailed diffs

14. **[file_comparison.tsv](file_comparison.tsv)**
    - Raw file-by-file comparison data
    - Size, line count, MD5 hashes
    - Match/differ status
    - Foundation for other reports

15. **[discrepancy_summary.tsv](discrepancy_summary.tsv)**
    - All identified discrepancies
    - Severity ratings (low/medium/high)
    - Impact descriptions
    - Categorized by type

---

## 📈 Data Analysis

16. **[figure6_statistics_comparison.json](figure6_statistics_comparison.json)**
    - IBD and PD disease statistics
    - Current repo values
    - NERSC repo status (missing)

17. **[competency_counts_comparison.json](competency_counts_comparison.json)**
    - Butyrate competency analysis row counts
    - Organismal traits, RHEA, EC pathway counts
    - Current vs NERSC differences

18. **[tryptophan_analysis.json](tryptophan_analysis.json)**
    - Documents tryptophan analysis in NERSC (removed from current)
    - 275 files documented
    - ~1,105 lines of file listings

19. **[critical_statistics.json](critical_statistics.json)**
    - Empty file (artifact from initial analysis)

---

## 📁 File Organization by Purpose

### For Understanding the Major Issue:
1. FINAL_ANALYSIS_SUMMARY.md
2. CRITICAL_DISCREPANCY_ANALYSIS.md
3. figure6_publication_results.txt
4. publication_standard_comparison_report.md

### For Detailed Data:
- publication_standard_comparison.tsv (publication vs repos)
- comparison_table.tsv (current vs NERSC)
- statistical_validation.tsv (metrics validation)

### For Code Changes:
- code_diff_summary.tsv (summary)
- code_comparison.tsv (full diffs)
- file_path_mapping.tsv (path changes)

### For JSON/Structured Data:
- publication_standard_comparison.json
- figure6_statistics_comparison.json
- competency_counts_comparison.json
- tryptophan_analysis.json

---

## 🎯 Key Conclusions Across All Reports

### ✅ What Works:
- Current repository has complete analysis outputs
- NERSC repository missing disease classification files
- Code can be run and generates outputs
- Statistical associations are strong and significant
- Chi-square values within 20% of each other

### ❌ What's Broken:
- **Repository data ≠ Publication figure data**
- Butyrate producer counts appear under OPPOSITE labels
- Total counts differ by 3-5x
- Unclear which files generated the publication figure
- Cannot directly reproduce publication figure from repository

### ⚠️ What Needs Clarification:
1. Definition of "increased/decreased likelihood of disease"
2. Definition of "increased/decreased in disease state"
3. Which exact output files generated Figure 6
4. How were total counts (54,952, 25,521) calculated
5. Is the terminology inversion intentional or an error

---

## 📊 The Numbers At A Glance

### IBD (Inflammatory Bowel Disease)

| Source | Dec Butyrate | Inc Butyrate | Total |
|--------|--------------|--------------|-------|
| Publication Figure | 221 | 514 | 54,952 |
| Current CSV | 478 | 186 | 10,800 |
| Ranks File (strains) | 469 | 183 | 10,779 |

**Pattern:** Repository counts match publication when SWAPPED!

### PD (Parkinson's Disease)

| Source | Dec Butyrate | Inc Butyrate | Total |
|--------|--------------|--------------|-------|
| Publication Figure | 152 | 299 | 25,521 |
| Current CSV | 302 | 150 | 8,641 |
| Ranks File (strains) | 295 | 124 | 8,527 |

**Pattern:** Repository counts match publication when SWAPPED!

---

## 🔍 How to Use These Reports

**For quick overview:**
→ Start with README.md

**To understand the reproducibility issue:**
→ Read FINAL_ANALYSIS_SUMMARY.md
→ Then CRITICAL_DISCREPANCY_ANALYSIS.md

**For detailed comparisons:**
→ Use the TSV files in Excel/spreadsheet software
→ publication_standard_comparison.tsv (main comparison)
→ comparison_table.tsv (file-by-file)

**For programmatic access:**
→ Use JSON files
→ publication_standard_comparison.json (complete data)

**For specific analyses:**
→ Code changes: code_diff_summary.tsv
→ File paths: file_path_mapping.tsv
→ Statistics: statistical_validation.tsv

---

## ⚡ Quick Links

**Most Important File:** [FINAL_ANALYSIS_SUMMARY.md](FINAL_ANALYSIS_SUMMARY.md)

**Quick Data Check:** [publication_standard_comparison.tsv](publication_standard_comparison.tsv)

**Original Comparison:** [comparison_report.md](comparison_report.md)

**Publication Values:** [figure6_publication_results.txt](figure6_publication_results.txt)

---

## 📞 Questions to Answer Before Publication

These questions emerged from the comprehensive analysis:

1. **Terminology:**
   - What does "increased/decreased likelihood of disease" measure?
   - What does "Num_Increased_disease_Total" in CSV mean?
   - Are these measuring the same thing?

2. **Data Source:**
   - Which exact file(s) generated Figure 6 in the publication?
   - Was it the summary CSV, ranks file, or something else?
   - Were there any transformations applied?

3. **Counts:**
   - How were the total counts (54,952 for IBD, 25,521 for PD) calculated?
   - Why do they differ 5x from the repository totals?
   - What granularity level was used (species, strain, edges)?

4. **Validation:**
   - Can the publication figure be regenerated from current repository?
   - If so, what script/workflow does this?
   - If not, where is the missing data/code?

---

*This index covers all 19 files in the comparison analysis*
*Last updated: 2026-01-30*
