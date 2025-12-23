# kg-microbe-paper

Repository for competency and biomedical analyses of kg-microbe-biomedical-function.

## Setup

### Prerequisites
- Unix-based OS (not tested on Windows)
- Python >= 3.11, < 3.13 (catboost v1.2.7 is not compatible with Python >= 3.13)
- [uv](https://docs.astral.sh/uv/) package manager

### Installation

1. Install uv if you haven't already:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

2. Clone the repository and install dependencies:
```bash
git clone https://github.com/bsantan/kg-microbe-paper.git
cd kg-microbe-paper
uv sync
```

### Dependencies
All dependencies are managed via `pyproject.toml`:
- Python >= 3.11, < 3.13
- tqdm==4.67.1
- matplotlib==3.10.0
- matplotlib_venn==1.1.1
- statsmodels==0.14.4
- equilibrator_api==0.6.0
- owlready2==0.47
- catboost==1.2.7
- scikit-learn==1.6.1
- shap==0.46.0
- seaborn==0.13.2
- kaleido==0.2.1
- numpy==1.26.4
- pandas==2.2.3
- duckdb==1.0.0

## Running the Analysis Scripts

All commands should be run with `uv run` to ensure the correct Python environment is used.

In order to create subfiles of the necessary edges in the graph, first run the following command:

```bash
cd src
uv run make all
```

This will create 3 files that are used in this analysis.

```
merged-kg_edges_noEC.tsv
merged-kg_edges_competency_specific_ec.tsv
merged-kg_edges_ncbitaxon.tsv
```

Next, the following file will need to be downloaded into /src/Input_Files directory:

https://github.com/Knowledge-Graph-Hub/kg-microbe/releases/download/2025-03-07/ontologies.tar.gz

Then extract it:

```
cd src/Input_Files
tar -xvzf ontologies.tar.gz ncbitaxon_nodes.tsv
```

### Gut Microbiome Competencies (Human Microbiome Project)

The first script will perform a series of DuckDB queries to evaluate the existance of taxa from the HMP, and examine the presence of organismal traits or functional annotations for those taxa.

```bash
cd src
uv run python gut_microbes_competencies.py
```

Note: in order for this to run, the 'ncbitaxon_nodes.tsv' file must be present in the Input_Files directory, which can be accessed at '/data/transformed/ontologies/ncbitaxon_nodes.tsv' by running the kg-microbe transform step: https://github.com/Knowledge-Graph-Hub/kg-microbe.git.

#### Expected Outputs

Step 1: identify which taxa have organismal traits or functional annotations. The organismal traits are identified by a DuckDB query internally, and the functional annotations are output to './Phylogeny_Search':

```
unique_ncbitaxon_uniprot_ids.txt
```

Step 2: download HMP supplementary file, map all HMP taxa to NCBI Taxonomy IDs, and evaluate which have organismal traits or functional annotations. These summary files are output to './Intermediate_Files':

```
HMP_Microbes_Mapped_Summary.csv
HMP_Microbes_Mapped.csv
```

### Metabolite Competencies

The second script will perform a series of DuckDB queries to identify taxa with with a given metabolic trait. Currently, this supports finding taxa with one of 4 semantic representations of butyrate production.

```bash
cd src
uv run python Process_competency_questions.py
```

#### Expected Outputs

Step 1: perform competencies in KG that involve microbe-metabolite annotation or microbe-protein-biochemical_reaction-metabolite annotation. This supports any metabolite or any direction, and also produces a 'term_mappings.csv' file for cases where multiple mappings terms are available for a given metabolite. For butyrate production, these are output to './Intermediate_Files_Competencies/butyrate_produces':

```
NCBI_organismal_traits.tsv: edges with pattern 'NCBITaxon/strain,biolink:produces,butyrate'
NCBI_organismal_traits_strains.tsv: strain rolled up to NCBITaxon node, edges with pattern 'NCBITaxon,biolink:produces,butyrate'
NCBI_organismal_traits_strains_proteomes.tsv: taxa that have any functional annotations from UniProt in the KG (i.e. have a proteome), edges with pattern 'NCBITaxon,biolink:produces,butyrate'
NCBI_genomic_traits_GO.tsv: taxa that have at least one protein with a GO annotation that has 'butyrate' in the label
NCBI_genomic_traits_RHEA.tsv: taxa that have at least one protein with a RHEA annotation that has 'butyrate' as a participant
```

Step 2: Predict reaction direction for microbe-protein-biochemical_reaction-metabolite annotations. This step uses the eQuilibrator API. For butyrate production, these are output to './Intermediate_Files_Competencies/butyrate_produces':

```
reaction_direction_dict.json: A dictionary of all RHEA terms from NCBI_genomic_traits_RHEA.tsv and the predicted predicate representing reaction direction for butyrate (input or output)
NCBI_genomic_traits_RHEA.tsv: UPDATED version, taxa that have at least one protein with a RHEA annotation that has 'butyrate' as an output
```

Step 3: perform competencies in KG that involve a set of ECs that are defined as being part of a pathway involved in the metabolite/direction (butyrate production in this case). 4 butyrate production pathways are represented in this case with a unique set of ECs. For butyrate production, these are output to './Intermediate_Files_Competencies/butyrate_produces':

```
NCBI_genomic_traits_EC_pathway_N.tsv: taxa that have at least one functional protein with at least all but one EC in the corresponding pathway (N)
NCBI_genomic_traits_EC_pathway_all.tsv: taxa that have at least one functional protein any pathway annotation from the 4 defined EC sets
```

Step 4: Create summary files. These are output to './Intermediate_Files_Competencies':

```
All_Competencies_produces_raw.tsv: for each metabolite given (butyrate in this case), the following values are recorded:
- Total_Traits_Annotations: total edges with pattern 'NCBITaxon,biolink:produces,butyrate'
- Total_Traits_Proteome_Overlap: total edges with pattern 'NCBITaxon,biolink:produces,butyrate' that also have any functional annotation from Uniprot (i.e. have a proteome)
- Total_Rhea-Chebi_Annotations: total edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate'
- Total_Rhea-Chebi_Traits_Overlap: taxa in edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' that also have edges with pattern 'NCBITaxon,biolink:produces,butyrate'
- EC_Annotations,Total_EC_Traits_Overlap: taxa in edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
- Total_EC_Rhea-Chebi_Overlap: taxa in edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' that also have edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate'
- Total_EC_and_Rhea-Chebi: taxa in edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' or that have edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate'
- Total_Rhea-Chebi_and_EC_Traits_Overlap: taxa in edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' that also have edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' that also have edges with pattern 'NCBITaxon,biolink:produces,butyrate'
- Total_Traits_and_Rhea-Chebi_and_EC_Annotations: taxa in edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' or that have edges with pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' or that have edges with pattern 'NCBITaxon,biolink:produces,butyrate'
All_Competencies_produces_summary.tsv: for each metabolite given (butyrate in this case), the following values are summarized:
- traits_only: total edges with only the pattern 'NCBITaxon,biolink:produces,butyrate'
- ec_only: total edges with only the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
- rhea_chebi_only: total edges with only the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate'
- traits_rhea_chebi_only: total edges with only both the pattern 'NCBITaxon,biolink:produces,butyrate' and 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate'
- traits_ec_only: total edges with only both the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' and 'NCBITaxon,biolink:produces,butyrate'
- ec_rhea_chebi_only: total edges with only both the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' and 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate'
- traits_rhea_chebi_ec_only: total edges with only all patterns 'NCBITaxon,biolink:produces,butyrate' and 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' and 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' 
```

Step 5: Create visualizations of competencies. These are output to './Intermediate_Files_Competencies':

```
All_Sources_Competencies_Comparison_butyrate_produces.png: A visualization of all semantic representations of butyrate production in the KG (from All_Competencies_produces_summary.tsv)
Monte_Carlos_Distribution.png: Monte Carlo simulation of 1,000 iterations to test significance of overlap between taxa in edges with the pattern 'NCBITaxon,biolink:produces,butyrate' vs. taxa in edges with the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
```
  
Step 6: The second part of this script compares the butyrate metabolite competency results in the KG to those of another analysis in Vital et al. Ensure that the ncbitaxon nodes file exists as 'Input_Files/ncbitaxon_nodes.tsv'. 

Note: For systems with 24GB of RAM or less, this step may not be able to execute. The Gold_Standard_Species_Overlap_butyrate_produces.csv file can therefore be accessed in the '/data' directory of this repository in order to proceed to the next step (move this to the './Intermediate_Files' directory).

```
gold_standard_ids.tsv: Taxa from Input_Files/Vital_etal_butyrate+producing_microbes.csv mapped to NCBITaxon IDs automatically using string matching.
gold_standard_ids_manual.tsv: For taxa that could not be mapped, manual mappings added to NCBITaxon IDs.
- Family/Family_ID: taxon at the family level.
- Name/Name_ID: taxon at the species level.
Gold_Standard_Species_Overlap_butyrate_produces.csv: Summary of presence of taxa at the species level in each semantic representation of the KG or Vital et al. set.
- Value: Species NCBITaxon ID
- Proteome: 1 if taxon has a proteome, 0 if not
- Gold_Standard: 1 if taxon in Vital et al., 0 if not
- Organismal: 1 if taxon in edge of the pattern 'NCBITaxon,biolink:produces,butyrate', 0 if not
- Functional: 1 if taxon in edge of the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate', 0 if not
- Functional_EC: 1 if taxon in edge of the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>', 0 if not
- Functional_Protein_Name: name of involved protein if taxon in edge of the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
- Functional_EC_Pathway: pathway that taxon is involved in if taxon in edge of the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' (note: butyrate_pathway_1 = acetyl_coa_pathway_1, butyrate_pathway_2 = acetyl_coa_pathway_2, butyrate_pathway_3 = glutarate_pathway, butyrate_pathway_4 = 4_aminobutyrate_pathway, butyrate_pathway_5 = lysine_pathway)
Gold_Standard_Species_Venn_Diagrams_EC_RHEA_Org.png: Visualization of number of taxa at species level in KG vs. Vital et al.
Gold_Standard_Families_Overlap_butyrate_produces.csv: Summary of presence of taxa at the family level in each semantic representation of the KG or Vital et al. set.
- Value: Family NCBITaxon ID
- Gold_Standard: 1 if taxon in Vital et al., 0 if not
- Organismal: 1 if taxon in edge of the pattern 'NCBITaxon,biolink:produces,butyrate', 0 if not
- Functional: 1 if taxon in edge of the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate', 0 if not
Gold_Standard_Families_Venn_Diagrams.png: Visualization of number of taxa at the family level in KG vs. Vital et al. (excluding pathway EC sets)
```

### Comparison to literature set

The third script will use the outputs from the Metabolite Competencies to compare to the Vital et al.

```bash
cd src
uv run python Gold_standard_Competency_analysis.py
```

#### Expected Outputs

Step 1: identify ranks of relevant taxa. These are output to './Phylogeny_Search':

```
ncbitaxon_rank.tsv: The corresponding rank of each taxa in NCBI Taxonomy, found using owlready2. 
- NCBITaxon_ID: taxon ID
- Rank: taxon rank
Gold_Standard_Species_Overlap_all_butyrate_produces_microbes_phylum.json: phylum of each taxa from either the semantic representations in the KG or Vital et al. that can produce butyrate (using Gold_Standard_Species_Overlap_butyrate_produces.csv), phylum: taxon pairs
Gold_Standard_Species_Overlap_all_butyrate_produces_microbes_family.json: family of each taxa from either the semantic representations in the KG or Vital et al. that can produce butyrate (using Gold_Standard_Species_Overlap_butyrate_produces.csv), family: taxon pairs
Gold_Standard_Species_Overlap_all_butyrate_produces_microbes_genus.json: genus of each taxa from either the semantic representations in the KG or Vital et al. that can produce butyrate (using Gold_Standard_Species_Overlap_butyrate_produces.csv), genus: taxon pairs
```

The children of each taxa at the species or strain level is found next. The following files include each taxa at the family level from either the semantic representations in the KG or Vital et al. that can produce butyrate (using Gold_Standard_Species_Overlap_butyrate_produces.csv). This is output to './Intermediate_Files':

```
competencies_all_microbes_families_butyrate_produces_microbes_species.json: species of each taxon, taxon: species pairs
competencies_all_microbes_families_butyrate_produces_microbes_strain.json: strains of each taxon, taxon: strain pairs
competencies_all_microbes_families_butyrate_produces_microbes_strains_and_species.json: all children (strains if they exist or species if not) of each taxa, family: taxon pairs
```

Note: Steps 2 and 3 are done for the following sets of taxa, and the end of the filename is listed for each step:

```
proteomes_no_gs: taxa with proteome but no butyrate production semantic representation in KG
organismal_no_gs: taxa with edges of the pattern 'NCBITaxon,biolink:produces,butyrate' but not in Vital et al. set
gs_proteomes_total: taxa with proteome and in Vital et al. set
gs_proteomes_no_annotation: taxa with proteome and in Vital et al. set but no butyrate production semantic representation in KG
gs_no_kg: taxa in the Vital et al. set but no butyrate production semantic representation in KG
gs_and_kg: taxa in the Vital et al. set and any butyrate production semantic representation in KG
ec_no_gs: taxa with edges of the pattern 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>' but not in Vital et al. set
all_no_gs: taxa with any butyrate production semantic representation in KG but not in Vital et al. set
all_kg: taxa with any butyrate production semantic representation in KG
```

Step 2: summarize the ranks of relevant taxa, and assess the importance in the human gut. Each of the following files consists of a different set of taxa defined below with the following columns:
- Value: taxon NCBITaxon ID
- Rank: taxon rank
- Genus: genus that the taxon belongs to
- Family: family that the taxon belongs to
- Phylum: phylum that the taxon belongs to
- Location: 'human' if phylum is in the list [Bacteroidota, Bacillota, Actinomycetota, Pseudomonadota, Fusobacteriota], 'other' if not
- Impt_Family: 'impt_fam' if phylum is in the list [Ruminococcaceae, Lachnospiraceae, Bacteroidaceae, Prevotellaceae, Clostridaceae, Oscillospiraceae, Eubacteriaceae], 'other' if not

The first file is output to './Phylogeny_Search':

```
all_microbes_butyrate_produces_families.tsv: taxa from KG or Vital et al. set
```

The other files are output to './Intermediate_Files_Competencies/butyrate_produces/Competency_Analysis'. This is done for the set of taxa listed above and creates files with the substring:

```
_butyrate_produces_families.tsv
```

A treemap of all families organized by the phylum level is also output for each set of taxa list above and creates files with the substring:

```
_ncbi_taxonomy_treemap.png
```

Step 3: Analyze the number of butyrate producers per genus. This is done for the set of taxa listed above and creates files with the substring:

```
_Genus_Traits_Comparison_species_and_strain_all.tsv
```

Next this is done for only those with the 'impt_fam' annotation. This is done with a threshold for including only the top X percent of families by the number of butyrate producers per family (100% or 10%). These are output to './Intermediate_Files_Competencies/butyrate_produces/Competency_Analysis/Thresold_<1.0 or 0.1>':

```
_Genus_Traits_Comparison_species_and_strain_all.tsv
```

### Biomedical Analysis

The fourth script will use the outputs from the Metabolite Competencies to analyze microbial metabolism in the context of disease.

```bash
cd src
uv run python Classification_gold_standard_comparison.py
```

#### Expected Outputs

Step 1: identify all taxa involved in diseases according to KG (for PD and IBD) and identify all children (strains if they exist otherwise species) for relevant taxa. These are output to './Intermdiate_Files':

```
classification_butyrate_produces_<IBD/PD_microbes_strain.json: species of each taxon, taxon: species pairs
classification_butyrate_produces_<IBD/PD_microbes_strain.json: strains of each taxon, taxon: strain pairs
```

Step 2: summarize the number of species or strains per taxa in the disease set. These are output to './Intermdiate_Files':

```
outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_<IBD/PD>.csv: each taxa annotated as being involved with disease and the number of butyrate producers by their children
- Name: taxon NCBITaxon ID
- Disease_Relationship: direction of disease involvement (increase or decrease likelihood of disease)
- Rank: taxon rank
- Num_Species: number of children that are species level
- Num_Species_Butyrate_Producers: number of children that are species level and are butyrate producers according to the KG, where taxa are in edges of the patterns 'NCBITaxon,biolink:produces,butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
- Num_Strains: number of children that are strain level
- Num_Strains_Butyrate_Producers: number of children that are strain level and are butyrate producers according to the KG, where taxa are in edges of the patterns 'NCBITaxon,biolink:produces,butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
```

Step 3: determine significance of difference between butyrate producers in set of taxa that increase vs. decrease disease. These are output to './Intermediate_Files':

```
<IBD/PD>_Classification_butyrate_producers_summary.csv
Num_Decreased_disease_Total: total taxa that decrease liklihood of disease
Num_Increased_disease_Total: total taxa that increase liklihood of disease
Num_Decreased_disease_Butyrate_Producers: total taxa that decrease liklihood of disease and are butyrate producers according to the KG, where taxa are in edges of the patterns 'NCBITaxon,biolink:produces,butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
Num_Increased_disease_Butyrate_Producers: total taxa that increase liklihood of disease and are butyrate producers according to the KG, where taxa are in edges of the patterns 'NCBITaxon,biolink:produces,butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:participates_in, RHEA, <predicted to produce>, butyrate' or 'NCBITaxon, biolink:derives_from, UniprotKB,biolink:enables, <set of EC in defined pathway>'
chi2: chi2 value of difference between butyrate producers in increased vs decreased group after a chi-squared test
P_Val: p value of difference between butyrate producers in increased vs decreased group after a chi-squared test
```

