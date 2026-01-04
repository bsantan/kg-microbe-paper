# Testing Quick Reference

## Automated Test (Recommended)

Run the full automated test suite:

```bash
./test_reorganization.sh
```

This will:
1. Verify git-tracked reference files exist
2. Run all Makefile targets
3. Execute all 4 main Python analysis scripts
4. Check files are in correct locations
5. Verify no files remain in old `src/` paths

**Runtime:** ~30-60 minutes (depending on download speed and compute)

---

## Manual Testing

### Step 1: Run Makefile
```bash
make all
```

**Expected outputs:**
- `data/Input_Files/KGMicrobe-biomedical-function-20250222.tar.gz` (4.3 GB)
- `data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges*.tsv` (3 files, ~44 GB total)
- `data/Input_Files/ncbitaxon_nodes.tsv` (130 MB)

### Step 2: Run Python Analysis Scripts (in order)

```bash
cd src

# 1. HMP gut microbiome analysis
python gut_microbes_competencies.py

# 2. Metabolite competency analysis
python Process_competency_questions.py

# 3. Gold standard comparison
python Gold_standard_Competency_analysis.py

# 4. Disease classification (optimized version)
python Classification_gold_standard_comparison_optimized.py
```

### Step 3: Verify File Locations

```bash
cd ..

# Check all expected output files exist
ls -lh data/Input_Files/Gold_Standard_*.csv
ls -lh data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges*.tsv
ls -lh data/Phylogeny_Search/ncbitaxon_rank.tsv
ls -lh data/Intermediate_Files/*.tsv | head -10
ls -lh data/Intermediate_Files_Competencies/butyrate_produces/*.tsv | head -10
ls -lh data/ML_model_shap/*.tsv
```

### Step 4: Check Git Status

```bash
git status
```

**Should show ONLY these untracked files:**
- `data/KGMicrobe-biomedical-function-20250222.tar.gz`
- `data/kg-microbe-biomedical-function-cat/`
- `data/ncbitaxon_nodes.tsv`
- `data/ontologies.tar.gz`

**Should NOT show:**
- Any files in `data/Intermediate_Files/`
- Any files in `data/Intermediate_Files_Competencies/`
- Any files in `data/Phylogeny_Search/`
- Any files in `data/ML_model_shap/`

(These are tracked in git)

---

## Expected File Counts

| Directory | Files | Total Size |
|-----------|-------|------------|
| `data/Input_Files/` | 5 items | ~4.5 GB |
| `data/Input_Files/kg-microbe-biomedical-function-cat/` | 5 TSV files | ~62 GB |
| `data/Phylogeny_Search/` | 1 TSV file | ~57 MB |
| `data/Intermediate_Files/` | 10-20 TSV/CSV files | ~100-500 MB |
| `data/Intermediate_Files_Competencies/butyrate_produces/` | 15-25 TSV files | ~10-50 MB |
| `data/ML_model_shap/` | 4 TSV files | ~492 KB |

---

## Quick Validation Commands

### Verify all main outputs exist:
```bash
test -f data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv && \
test -f data/Input_Files/ncbitaxon_nodes.tsv && \
test -f data/Phylogeny_Search/ncbitaxon_rank.tsv && \
test -f data/Intermediate_Files_Competencies/butyrate_produces/butyrate_produces_all_methods.tsv && \
echo "✓ All key files present" || echo "✗ Missing files"
```

### Check no old locations have data:
```bash
find src/Input_Files src/Intermediate_Files -type f 2>/dev/null | wc -l
```
**Expected:** 0

### Verify file headers:
```bash
# All edge files should have this header
head -n 1 data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv
head -n 1 data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv
```
**Expected:** `subject	predicate	object`

---

## Troubleshooting

### "FileNotFoundError" in Python scripts
**Cause:** Path references still pointing to old `src/` locations
**Fix:** Check which file is missing and verify path updates

### "DuckDB error: Could not open file"
**Cause:** KG files not downloaded or in wrong location
**Fix:** Run `make all` to download files

### Git showing large files as changes
**Cause:** .gitignore not configured correctly
**Fix:** Run `git check-ignore -v <file>` to debug

### OOM errors
**Cause:** Insufficient memory for full classification script
**Fix:** Use `Classification_gold_standard_comparison_optimized.py` instead of `classification.py`

---

## NERSC-Specific Testing

On NERSC, test with SLURM:

```bash
cd /global/cfs/cdirs/m4689/kg-microbe-paper
sbatch run_classification.sl
```

Monitor:
```bash
squeue -u $USER
tail -f logs/classification_*.out
```

Verify UV cache is on local scratch:
```bash
grep UV_CACHE_DIR run_classification.sl
```
Should show `${TMPDIR:-/tmp}/uv-cache-$$`, NOT an NFS path.
