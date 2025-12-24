# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains competency and biomedical analyses of the kg-microbe-biomedical-function knowledge graph. The analyses focus on microbial metabolite production (particularly butyrate), gut microbiome competencies, and disease associations (IBD, Parkinson's Disease).

## System Requirements

- Unix-based OS (not tested on Windows)
- Python >= 3.11, < 3.13 (catboost v1.2.7 is not compatible with Python >= 3.13)
- [uv](https://docs.astral.sh/uv/) package manager

## Initial Setup

### 1. Install Dependencies

Install project dependencies using uv:

```bash
uv sync
```

### 2. Download All Data

Run the Makefile from the repository root to download all required data:

```bash
uv run make all
```

This single command will:
- Download the KG-Microbe biomedical function knowledge graph from NERSC
- Extract and create three essential edge files in `src/Input_Files/kg-microbe-biomedical-function-cat/`:
  - `merged-kg_edges_noEC.tsv` - NCBI Taxon, CHEBI, UniprotKB, and RHEA edges
  - `merged-kg_edges_competency_specific_ec.tsv` - EC enzyme classification edges
  - `merged-kg_edges_ncbitaxon.tsv` - NCBI Taxonomy edges
- Download NCBI Taxonomy ontologies from kg-microbe releases
- Extract `ncbitaxon_nodes.tsv` to `src/Input_Files/`

## Running Analysis Scripts

All Python scripts should be run from the repository root using `uv run`:

```bash
uv run python src/<script_name>.py
```

### Main Analysis Pipeline (in order)

Run these scripts sequentially from the repository root:

1. **Gut Microbiome Competencies (HMP)**
   ```bash
   uv run python src/gut_microbes_competencies.py
   ```
   Analyzes Human Microbiome Project taxa for organismal traits and functional annotations.

2. **Metabolite Competencies**
   ```bash
   uv run python src/Process_competency_questions.py
   ```
   Identifies taxa with metabolic traits (e.g., butyrate production) using multiple semantic representations.

3. **Gold Standard Comparison**
   ```bash
   uv run python src/Gold_standard_Competency_analysis.py
   ```
   Compares KG results to literature (Vital et al.) for butyrate producers.

4. **Biomedical Analysis**
   ```bash
   uv run python src/Classification_gold_standard_comparison.py
   ```
   Analyzes microbial metabolism in disease contexts (IBD, PD).

### Machine Learning Analysis

The Jupyter notebooks in `src/` provide additional analyses:
- `kg_microbe_train_taxa_to_temperature__EC_RHEA.ipynb` - Train models to predict temperature preferences
- `trait_bubble_plot.ipynb` - Visualize traits

## Code Architecture

### Core Modules

**`constants.py`**
- Central configuration file containing:
  - Metabolite mappings (CHEBI, GO term IDs)
  - File naming conventions
  - Disease mappings (IBD, PD)
  - Model hyperparameters
  - Gut microbiome reference lists

**`duckdb_utils.py`**
- DuckDB query utilities for knowledge graph operations
- Key functions: `duckdb_load_table()`, `output_table_to_file()`, `get_node_label()`
- All KG queries use DuckDB in-memory databases

**`ncbi_phylogeny_search.py`**
- NCBI Taxonomy phylogeny traversal using owlready2
- Functions for finding parent/child taxa at different ranks
- Downloads and caches NCBI Taxonomy OWL ontology

**`Competencies.py`**
- Main competency analysis logic (~95KB, largest module)
- Implements three semantic representations of metabolite competencies:
  1. Organismal traits (direct biolink:produces edges)
  2. Functional annotations via RHEA reactions
  3. EC pathway-based annotations
- Uses eQuilibrator API for reaction direction prediction

### Analysis Flow

1. **Knowledge Graph Loading**: DuckDB loads TSV edge/node files into in-memory tables
2. **Competency Queries**: Multiple DuckDB queries identify taxa with metabolite-related annotations
3. **Phylogeny Resolution**: owlready2 resolves taxonomic hierarchies and ranks
4. **Gold Standard Mapping**: String matching + manual curation maps literature taxa to NCBI Taxonomy
5. **Statistical Analysis**: Monte Carlo simulations test significance of overlaps
6. **Visualization**: matplotlib/seaborn create Venn diagrams, treemaps, bar plots

### Key Design Patterns

- **Directory Structure**: Each metabolite/direction creates subdirectories under `Intermediate_Files_Competencies/`
- **File Naming**: Consistent patterns like `{metabolite}_{direction}` (e.g., `butyrate_produces`)
- **Graph Patterns**: Edge patterns like `NCBITaxon -> produces -> CHEBI` or `NCBITaxon -> derives_from -> UniprotKB -> participates_in -> RHEA`
- **Caching**: Intermediate files avoid re-downloading/re-computing (check existence before downloading)

### Important Constants

Update `constants.py` to modify:
- `ALL_METABOLITES`: List of metabolites to analyze (default: `["butyrate"]`)
- `ALL_DIRECTIONS`: Analysis directions (default: `["produces"]`)
- `RANDOM_SEED`: For reproducibility (default: 12)
- `MODEL_PARAMETERS`: CatBoost hyperparameters

## Output Files

Analysis outputs are organized into:
- `src/Intermediate_Files/` - HMP analysis, disease analysis
- `src/Intermediate_Files_Competencies/{metabolite}_{direction}/` - Competency results per metabolite
- `src/Phylogeny_Search/` - NCBI Taxonomy rank files, phylogeny searches
- `data/` - SHAP feature importance from ML models

## Common Workflows

### Adding a New Metabolite

1. Add CHEBI/GO mappings to `constants.py` under `CHEBI_mappings` and `GO_mappings`
2. Add relevant term synonyms to `METABOLITES_RELEVANT_TERMS`
3. Add metabolite name to `ALL_METABOLITES` list
4. Run `uv run python Process_competency_questions.py`

### Modifying EC Pathways

EC pathway definitions are hardcoded in `Competencies.py` in the `genomic_ec_competency()` function. Each pathway is a list of EC numbers representing a biosynthetic pathway.

### Changing Disease Contexts

Update `COMPETENCY_DISEASE_MAP` in `constants.py` with MONDO disease IDs.

## Dependencies

The codebase relies heavily on:
- **DuckDB**: All KG queries (in-memory SQL database)
- **owlready2**: NCBI Taxonomy OWL ontology loading and traversal
- **equilibrator_api**: Thermodynamic reaction direction prediction
- **CatBoost**: Machine learning for trait prediction
- **pandas**: Data manipulation throughout
- **matplotlib/seaborn**: All visualizations

## Known Quirks

- Scripts assume execution from repository root directory
- The Makefile is located at the repository root (not in `src/`)
- Large files like `ncbitaxon_rank.tsv` (60MB) are created on first run and cached
- The `.venv` directory (created by uv) should not be committed (already in `.gitignore`)
- NCBI Taxonomy names are sometimes outdated (see `REPLACED_TAXA_NAMES` in `constants.py`)
- All commands should be run with `uv run` to ensure the correct Python environment is used
