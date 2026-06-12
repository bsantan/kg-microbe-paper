#!/bin/bash
# Pre-filter the large edges file to create a smaller disease-focused subset
# This can be run if the Python script still runs out of memory

INPUT_FILE="./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv"
OUTPUT_FILE="./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_disease_subset.tsv"

echo "Pre-filtering edges file to disease-relevant subset..."
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

# Extract header
head -1 "$INPUT_FILE" > "$OUTPUT_FILE"

# Filter for edges involving NCBITaxon or MONDO (disease) entities
# This reduces the file from ~150M rows to ~few thousand rows
echo "Filtering for NCBITaxon and MONDO entities..."
grep -E "NCBITaxon:|MONDO:" "$INPUT_FILE" >> "$OUTPUT_FILE"

echo ""
echo "Filtering complete!"
echo "Original file lines: $(wc -l < "$INPUT_FILE")"
echo "Filtered file lines: $(wc -l < "$OUTPUT_FILE")"
echo "Reduction: $(echo "scale=1; 100 * (1 - $(wc -l < "$OUTPUT_FILE") / $(wc -l < "$INPUT_FILE"))" | bc)%"
echo ""
echo "To use the filtered file, update Classification_gold_standard_comparison.py:"
echo "  Change: merged-kg_edges.tsv"
echo "  To:     merged-kg_edges_disease_subset.tsv"
