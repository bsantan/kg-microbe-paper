# Testing Strategy for Data Reorganization

This document outlines the testing strategy to verify that the data directory reorganization works correctly.

## Pre-Test Setup

### 1. Clean Environment
Remove any existing downloaded/generated files to start fresh:

```bash
# Remove downloaded files
rm -f data/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz
rm -f data/Input_Files/ontologies.tar.gz
rm -f data/Input_Files/ncbitaxon_nodes.tsv
rm -rf data/Input_Files/kg-microbe-biomedical-function-cat/

# Remove generated intermediate files (keep git-tracked reference files)
rm -rf data/Intermediate_Files/*.tsv
rm -rf data/Intermediate_Files/*.csv
rm -rf data/Intermediate_Files_Competencies/butyrate_produces/*.tsv
rm -rf data/Phylogeny_Search/*.owl
rm -rf data/Phylogeny_Search/*.sqlite3
```

### 2. Verify Git-Tracked Files Exist
These files should be present (tracked in git):

```bash
# Should exist (373 KB)
ls -lh data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv

# Should exist (4 files, ~123 KB each)
ls -lh data/ML_model_shap/

# Should exist (12 KB each)
ls -lh data/Intermediate_Files_Competencies/vital_ids*.csv
```

**Expected output:**
```
data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv: 373K
data/ML_model_shap/shap_feature_importance_class_temperature_*: 123K each (4 files)
data/Intermediate_Files_Competencies/vital_ids.csv: 12K
data/Intermediate_Files_Competencies/vital_ids_manual.csv: 12K
```

## Test 1: Makefile Targets

### Test 1.1: Download Knowledge Graph
```bash
cd /Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper
make download_kg
```

**Expected outputs:**
- `data/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz` (~4.3 GB)
- `data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv` (~29 GB)
- `data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv` (~15 GB)

**Verification:**
```bash
ls -lh data/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz
ls -lh data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_*.tsv
```

### Test 1.2: Download Ontologies
```bash
make download_ontologies
```

**Expected outputs:**
- `data/Input_Files/ontologies.tar.gz` (~75 MB)
- `data/Input_Files/ncbitaxon_nodes.tsv` (~130 MB)

**Verification:**
```bash
ls -lh data/Input_Files/ontologies.tar.gz
ls -lh data/Input_Files/ncbitaxon_nodes.tsv
```

### Test 1.3: Create Edge Subfiles
```bash
make rhea_chebi_competencies
make ec_competencies
make taxonomy_competencies
```

**Expected outputs:**
- `data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv` (~9.5 GB)
- `data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_competency_specific_ec.tsv` (~9.2 GB)
- `data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv` (~45 MB)

**Verification:**
```bash
ls -lh data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_*.tsv
head -n 2 data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv
head -n 2 data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_competency_specific_ec.tsv
head -n 2 data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv
```

All three files should have `subject	predicate	object` header.

### Test 1.4: Setup Gold Standard
```bash
make setup_gold_standard
```

**Expected output:**
- `data/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv` (copied from Input_Files)

**Verification:**
```bash
ls -lh data/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv
diff data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv \
     data/Intermediate_Files_Competencies/butyrate_produces/Gold_Standard_Species_Overlap_butyrate_produces.csv
```

Should show no differences.

### Test 1.5: Run All Makefile Targets
```bash
make all
```

Should complete all above steps without errors.

## Test 2: Python Analysis Scripts

Run from the `src/` directory in the correct order:

```bash
cd src
```

### Test 2.1: HMP Gut Microbiome Analysis
```bash
uv run python gut_microbes_competencies.py
```

**Expected outputs in `data/Intermediate_Files/`:**
- `HMP_species.tsv`
- `HMP_organismal_traits.tsv`
- `HMP_functional_annotations.tsv`
- Related visualization files (PNG)

**Verification:**
```bash
ls -lh ../data/Intermediate_Files/HMP_*.tsv
wc -l ../data/Intermediate_Files/HMP_*.tsv
```

Should show multiple TSV files with reasonable line counts.

### Test 2.2: Metabolite Competency Analysis
```bash
uv run python Process_competency_questions.py
```

**Expected outputs in `data/Intermediate_Files_Competencies/butyrate_produces/`:**
- `butyrate_produces_organismal_traits.tsv`
- `butyrate_produces_functional_annotations.tsv`
- `butyrate_produces_ec_competencies.tsv`
- `butyrate_produces_all_methods.tsv`
- Visualization files (Venn diagrams, treemaps)

**Expected outputs in `data/Phylogeny_Search/`:**
- `ncbitaxon_rank.tsv` (~57 MB - cached taxonomy ranks)

**Verification:**
```bash
ls -lh ../data/Intermediate_Files_Competencies/butyrate_produces/butyrate_produces_*.tsv
ls -lh ../data/Phylogeny_Search/ncbitaxon_rank.tsv
head -n 5 ../data/Intermediate_Files_Competencies/butyrate_produces/butyrate_produces_all_methods.tsv
```

### Test 2.3: Gold Standard Comparison
```bash
uv run python Gold_standard_Competency_analysis.py
```

**Expected outputs in `data/Intermediate_Files_Competencies/butyrate_produces/`:**
- `gold_standard_ids.tsv`
- `gold_standard_ids_manual.tsv`
- `gold_standard_overlap_*.tsv` (multiple files for different methods)
- Monte Carlo simulation results
- Comparison plots (PNG)

**Verification:**
```bash
ls -lh ../data/Intermediate_Files_Competencies/butyrate_produces/gold_standard_*.tsv
wc -l ../data/Intermediate_Files_Competencies/butyrate_produces/gold_standard_overlap_*.tsv
```

### Test 2.4: Disease Classification Analysis

**Option A: Full version (requires more memory):**
```bash
uv run python classification.py
```

**Option B: Optimized version (recommended for memory-constrained systems):**
```bash
uv run python Classification_gold_standard_comparison_optimized.py
```

**Expected outputs in `data/Intermediate_Files/`:**
- `outcome_to_NCBITaxon_cleaned.tsv`
- `NCBITaxon_to_disease_associations.tsv`
- Classification model outputs
- SHAP plots and feature importance files

**Verification:**
```bash
ls -lh ../data/Intermediate_Files/outcome_to_NCBITaxon_*.tsv
ls -lh ../data/Intermediate_Files/*_shap_*.png
```

## Test 3: File Location Verification

### Test 3.1: All Files in Correct Locations

Run this comprehensive check:

```bash
cd /Users/marcin/Documents/VIMSS/ontology/KG-Hub/KG-Microbe/paper_KGM/kg-microbe-paper

echo "=== Checking data/Input_Files/ ==="
ls -lh data/Input_Files/ | grep -v "^d" | awk '{print $9, $5}'

echo -e "\n=== Checking data/Input_Files/kg-microbe-biomedical-function-cat/ ==="
ls -lh data/Input_Files/kg-microbe-biomedical-function-cat/*.tsv | awk '{print $9, $5}'

echo -e "\n=== Checking data/Phylogeny_Search/ ==="
ls -lh data/Phylogeny_Search/*.tsv 2>/dev/null | awk '{print $9, $5}'

echo -e "\n=== Checking data/Intermediate_Files/ ==="
ls -lh data/Intermediate_Files/*.tsv 2>/dev/null | awk '{print $9, $5}'

echo -e "\n=== Checking data/Intermediate_Files_Competencies/butyrate_produces/ ==="
ls -lh data/Intermediate_Files_Competencies/butyrate_produces/*.tsv 2>/dev/null | head -20 | awk '{print $9, $5}'

echo -e "\n=== Checking data/ML_model_shap/ ==="
ls -lh data/ML_model_shap/*.tsv 2>/dev/null | awk '{print $9, $5}'
```

### Test 3.2: No Files in Old Locations

Verify that no data files remain in old `src/` locations:

```bash
echo "=== Checking for files in old src/ locations (should be empty) ==="
find src/Input_Files src/Intermediate_Files src/Intermediate_Files_Competencies src/Phylogeny_Search -type f 2>/dev/null
```

**Expected:** No output (directories may exist but should be empty or not exist).

## Test 4: Git Status Check

Verify that generated files are properly git-ignored:

```bash
git status
```

**Expected output:**
```
On branch main
Your branch is up to date with 'origin/main'.

Untracked files:
  (use "git add <file>..." to include in what will be committed)
	data/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz
	data/Input_Files/kg-microbe-biomedical-function-cat/
	data/Input_Files/ncbitaxon_nodes.tsv
	data/Input_Files/ontologies.tar.gz

nothing added to commit but untracked files present (use "git add" to track)
```

**Important:** Generated intermediate files in `data/Intermediate_Files/` should NOT appear as untracked (they're tracked in git).

## Test 5: NERSC-Specific Tests

If testing on NERSC, also verify:

### Test 5.1: SLURM Job Submission
```bash
cd /global/cfs/cdirs/m4689/kg-microbe-paper
sbatch run_classification.sl
```

**Verification:**
```bash
squeue -u $USER
tail -f logs/classification_*.out
```

### Test 5.2: UV Cache Location
```bash
echo $UV_CACHE_DIR
ls -lh $UV_CACHE_DIR
```

Should show local scratch directory, not NFS location.

## Success Criteria

✅ **Pass**: All tests complete without errors and all expected files are created in correct locations

| Test | Expected Files | Location | Size Check |
|------|---------------|----------|------------|
| 1.1 | KG edges/nodes | data/Input_Files/kg-microbe-biomedical-function-cat/ | ~44 GB total |
| 1.2 | Ontology files | data/Input_Files/ | ~205 MB total |
| 1.3 | Edge subfiles | data/Input_Files/kg-microbe-biomedical-function-cat/ | ~18 GB total |
| 1.4 | Gold Standard | data/Intermediate_Files_Competencies/butyrate_produces/ | 373 KB |
| 2.1 | HMP outputs | data/Intermediate_Files/ | Multiple TSV files |
| 2.2 | Competency outputs | data/Intermediate_Files_Competencies/butyrate_produces/ | Multiple TSV files |
| 2.2 | Taxonomy cache | data/Phylogeny_Search/ | ~57 MB |
| 2.3 | Gold standard comparison | data/Intermediate_Files_Competencies/butyrate_produces/ | Multiple TSV files |
| 2.4 | Disease analysis | data/Intermediate_Files/ | Multiple TSV/PNG files |

## Troubleshooting

### Issue: "File not found" errors
**Solution:** Check that all paths reference `data/` not `src/`

### Issue: OOM errors on classification scripts
**Solution:** Use `Classification_gold_standard_comparison_optimized.py` instead of `classification.py`

### Issue: Git showing large files as untracked
**Solution:** Verify .gitignore patterns with `git check-ignore -v <file>`

### Issue: NFS locking errors on NERSC
**Solution:** Ensure `UV_CACHE_DIR` is set to local scratch in SLURM script

## Quick Test Script

Save this as `test_reorganization.sh`:

```bash
#!/bin/bash

set -e

echo "=== Testing Data Reorganization ==="

# Test 1: Check git-tracked files exist
echo "1. Checking git-tracked reference files..."
test -f data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv
test -f data/ML_model_shap/shap_feature_importance_class_temperature_hyperthermophilic_temperature__EC_RHEA_v5.tsv
echo "✓ Git-tracked files present"

# Test 2: Run Makefile
echo "2. Running Makefile..."
make all
echo "✓ Makefile completed"

# Test 3: Check KG files created
echo "3. Verifying KG downloads..."
test -f data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv
test -f data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv
test -f data/Input_Files/ncbitaxon_nodes.tsv
echo "✓ KG files downloaded and processed"

# Test 4: Run Python scripts
echo "4. Running Python analysis scripts..."
cd src
uv run python gut_microbes_competencies.py
echo "✓ HMP analysis complete"

uv run python Process_competency_questions.py
echo "✓ Competency analysis complete"

uv run python Gold_standard_Competency_analysis.py
echo "✓ Gold standard comparison complete"

uv run python Classification_gold_standard_comparison_optimized.py
echo "✓ Disease classification complete"

cd ..

# Test 5: Verify no old locations
echo "5. Checking for files in old src/ locations..."
if find src/Input_Files src/Intermediate_Files -type f 2>/dev/null | grep -q .; then
    echo "✗ Files found in old src/ locations!"
    exit 1
fi
echo "✓ No files in old locations"

echo -e "\n=== All Tests Passed ==="
```

Run with:
```bash
chmod +x test_reorganization.sh
./test_reorganization.sh
```
