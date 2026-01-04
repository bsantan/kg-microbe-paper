#!/bin/bash

set -e

echo "==================================================================="
echo "Testing Data Reorganization for kg-microbe-paper"
echo "==================================================================="

# Test 1: Check git-tracked files exist
echo -e "\n[1/5] Checking git-tracked reference files..."
if [ -f data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv ]; then
    echo "  ✓ Gold_Standard file exists ($(ls -lh data/Input_Files/Gold_Standard_Species_Overlap_butyrate_produces.csv | awk '{print $5}'))"
else
    echo "  ✗ Gold_Standard file missing!"
    exit 1
fi

shap_count=$(ls data/ML_model_shap/*.tsv 2>/dev/null | wc -l)
if [ "$shap_count" -eq 4 ]; then
    echo "  ✓ All 4 SHAP files present in data/ML_model_shap/"
else
    echo "  ✗ Expected 4 SHAP files, found $shap_count"
    exit 1
fi

if [ -f data/Intermediate_Files_Competencies/vital_ids.csv ]; then
    echo "  ✓ vital_ids files present"
else
    echo "  ✗ vital_ids files missing!"
    exit 1
fi

# Test 2: Run Makefile
echo -e "\n[2/5] Running Makefile targets..."
echo "  → make download_kg (this will download ~4.3 GB)..."
make download_kg
echo "  ✓ Knowledge graph downloaded"

echo "  → make download_ontologies..."
make download_ontologies
echo "  ✓ Ontologies downloaded"

echo "  → make rhea_chebi_competencies..."
make rhea_chebi_competencies
echo "  ✓ RHEA-CHEBI edge file created"

echo "  → make ec_competencies..."
make ec_competencies
echo "  ✓ EC edge file created"

echo "  → make taxonomy_competencies..."
make taxonomy_competencies
echo "  ✓ Taxonomy edge file created"

echo "  → make setup_gold_standard..."
make setup_gold_standard
echo "  ✓ Gold standard copied to Intermediate_Files_Competencies"

# Test 3: Check KG files created
echo -e "\n[3/5] Verifying downloaded and generated files..."
if [ -f data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv ]; then
    echo "  ✓ merged-kg_edges.tsv exists ($(ls -lh data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv | awk '{print $5}'))"
else
    echo "  ✗ merged-kg_edges.tsv missing!"
    exit 1
fi

if [ -f data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv ]; then
    echo "  ✓ merged-kg_edges_noEC.tsv exists ($(ls -lh data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv | awk '{print $5}'))"
else
    echo "  ✗ merged-kg_edges_noEC.tsv missing!"
    exit 1
fi

if [ -f data/Input_Files/ncbitaxon_nodes.tsv ]; then
    echo "  ✓ ncbitaxon_nodes.tsv exists ($(ls -lh data/Input_Files/ncbitaxon_nodes.tsv | awk '{print $5}'))"
else
    echo "  ✗ ncbitaxon_nodes.tsv missing!"
    exit 1
fi

# Test 4: Run Python scripts
echo -e "\n[4/5] Running Python analysis pipeline..."
cd src

echo "  → python gut_microbes_competencies.py..."
if python gut_microbes_competencies.py > ../logs/test_hmp.log 2>&1; then
    echo "  ✓ HMP analysis complete"
    if [ -f ../data/Intermediate_Files/HMP_species.tsv ]; then
        echo "    - HMP_species.tsv created ($(wc -l < ../data/Intermediate_Files/HMP_species.tsv) lines)"
    fi
else
    echo "  ✗ HMP analysis failed! Check logs/test_hmp.log"
    exit 1
fi

echo "  → python Process_competency_questions.py..."
if python Process_competency_questions.py > ../logs/test_competency.log 2>&1; then
    echo "  ✓ Competency analysis complete"
    if [ -f ../data/Intermediate_Files_Competencies/butyrate_produces/butyrate_produces_all_methods.tsv ]; then
        echo "    - butyrate_produces_all_methods.tsv created ($(wc -l < ../data/Intermediate_Files_Competencies/butyrate_produces/butyrate_produces_all_methods.tsv) lines)"
    fi
    if [ -f ../data/Phylogeny_Search/ncbitaxon_rank.tsv ]; then
        echo "    - ncbitaxon_rank.tsv cached ($(ls -lh ../data/Phylogeny_Search/ncbitaxon_rank.tsv | awk '{print $5}'))"
    fi
else
    echo "  ✗ Competency analysis failed! Check logs/test_competency.log"
    exit 1
fi

echo "  → python Gold_standard_Competency_analysis.py..."
if python Gold_standard_Competency_analysis.py > ../logs/test_gold_standard.log 2>&1; then
    echo "  ✓ Gold standard comparison complete"
    gs_files=$(ls ../data/Intermediate_Files_Competencies/butyrate_produces/gold_standard*.tsv 2>/dev/null | wc -l)
    echo "    - $gs_files gold standard output files created"
else
    echo "  ✗ Gold standard comparison failed! Check logs/test_gold_standard.log"
    exit 1
fi

echo "  → python Classification_gold_standard_comparison_optimized.py..."
if python Classification_gold_standard_comparison_optimized.py > ../logs/test_classification.log 2>&1; then
    echo "  ✓ Disease classification complete"
    if [ -f ../data/Intermediate_Files/outcome_to_NCBITaxon_cleaned.tsv ]; then
        echo "    - outcome_to_NCBITaxon_cleaned.tsv created ($(wc -l < ../data/Intermediate_Files/outcome_to_NCBITaxon_cleaned.tsv) lines)"
    fi
else
    echo "  ✗ Disease classification failed! Check logs/test_classification.log"
    exit 1
fi

cd ..

# Test 5: Verify no old locations
echo -e "\n[5/5] Verifying no data files in old src/ locations..."
old_files=$(find src/Input_Files src/Intermediate_Files src/Intermediate_Files_Competencies src/Phylogeny_Search -type f 2>/dev/null || true)
if [ -n "$old_files" ]; then
    echo "  ✗ Files found in old src/ locations:"
    echo "$old_files"
    exit 1
else
    echo "  ✓ No files in old src/ locations"
fi

# Summary
echo -e "\n==================================================================="
echo "✅ All Tests Passed!"
echo "==================================================================="
echo ""
echo "Summary of created files:"
echo ""
echo "data/Input_Files/:"
ls -lh data/Input_Files/ | grep -v "^d" | grep -v "^total" | awk '{printf "  - %-60s %8s\n", $9, $5}'
echo ""
echo "data/Input_Files/kg-microbe-biomedical-function-cat/:"
ls -lh data/Input_Files/kg-microbe-biomedical-function-cat/*.tsv 2>/dev/null | awk '{printf "  - %-60s %8s\n", $9, $5}' | head -10
echo ""
echo "data/Phylogeny_Search/:"
ls -lh data/Phylogeny_Search/*.tsv 2>/dev/null | awk '{printf "  - %-60s %8s\n", $9, $5}'
echo ""
echo "data/Intermediate_Files/:"
ls -lh data/Intermediate_Files/*.tsv 2>/dev/null | awk '{printf "  - %-60s %8s\n", $9, $5}' | head -10
echo ""
echo "data/Intermediate_Files_Competencies/butyrate_produces/:"
ls -lh data/Intermediate_Files_Competencies/butyrate_produces/*.tsv 2>/dev/null | awk '{printf "  - %-60s %8s\n", $9, $5}' | head -10
echo ""
echo "==================================================================="
