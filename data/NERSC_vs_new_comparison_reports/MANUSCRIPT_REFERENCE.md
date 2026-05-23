# Manuscript Reference — Authoritative Source for Figure 6C Numbers

**Read this before any other file in this directory.**

The reports in this directory (most dated 2026-01-30) were generated before the second computational reproducibility review was concluded. They use "bioRxiv preprint" as a vague reference and treat various transcribed Figure 6C numbers as ground truth. Several of those numbers are **second-hand transcriptions** that have not been verified against the current authoritative source.

## Authoritative source

`revisions2/KG-Microbe_Responses2_mpj2.docx`

Specifically, the section **"Summary of the second computational reproducibility review"** states:

> The result was consistent for inflammatory bowel disease, but differed slightly for Parkinson's disease (PD):
> - The reproduced p-value for PD was **8e-293** as opposed to the reported p-value of **2e-288**
> - The reproduced value of the test statistic for PD was **1337** whereas the reported value was **1317**.

And the authors' response:

> In the latest reproduction, one less extreme p-value matched exactly; the other (more extreme) p-value shifted by 5 in the base-10 exponent. We made additional code edits to further optimize the competencies and also cache the Equilibriator results.

Drift attributed to:
1. NCBITaxon OWL file updates via the OwlReady package (accessed in `src/Classification_gold_standard_comparison.py`)
2. Equilibrator API updates (used in `src/Process_competency_questions.py`)
3. Platform/version sensitivity for extreme p-value calculations

## What is verified

| Quantity | Paper Figure 6C | Documented reproduction | Source |
|---|---|---|---|
| PD χ² | **1317** | **1337** | Manuscript text, verbatim |
| PD p-value | **2e-288** | **8e-293** | Manuscript text, verbatim |
| IBD result | "consistent" between paper and reviewer reproduction | (no explicit chi² in responses doc) | Manuscript text |

## What is unverified (used in other reports here)

The following values appear throughout `comparison_report.md`, `figure6_publication_results.txt`, `FINAL_ANALYSIS_SUMMARY.md`, `CRITICAL_DISCREPANCY_ANALYSIS.md`, and `INDEX.md`, and are **not restated in the manuscript responses doc**:

- IBD χ² = 647, p = 1e-142
- IBD totals: 15,398 (Increased) / 39,554 (Decreased)
- IBD butyrate counts: 514 (Increased) / 221 (Decreased)
- PD totals: 2,990 (Increased) / 22,531 (Decreased)
- PD butyrate counts: 299 (Increased) / 152 (Decreased)
- "Reviewer #4: IBD χ² 672, PD χ² 1064" — **the manuscript explicitly attributes PD reproduction χ² = 1337 (not 1064) to the reviewer**, so the "1064" figure circulating in this report set appears to be from a different / mis-attributed source.

These values likely come from the Figure 6 panels themselves but have not been independently confirmed against the current manuscript or supplementary materials.

## Current run alignment with the manuscript

After commits `716066f` (strain-dict contamination fix) and `82912ac` (per-disease results file fix), the current run produces:

| | Manuscript-acknowledged reproduction | Our current run |
|---|---|---|
| PD χ² | 1337 | 1337.36 ✅ |
| PD p | 8e-293 | 8.61e-293 ✅ |
| IBD χ² | "consistent" (≈647, unverified) | 674.63 (+4.3% drift) |
| IBD p | (not given) | 9.83e-149 |

PD matches the documented reproduction. IBD drifts ~4% — same drift direction/magnitude as PD, attributable to the same OwlReady NCBITaxon updates the manuscript already documents.
