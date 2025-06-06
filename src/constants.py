CHEBI_mappings = {
    "butyrate" : ['CHEBI:17968', 'CHEBI:85867', 'CHEBI:88764', 'CHEBI:87684', 'CHEBI:64103', 'CHEBI:88806', 'CHEBI:90150', 'CHEBI:87683', 'CHEBI:89719', 'CHEBI:87318', 'CHEBI:453', 'CHEBI:32097', 'CHEBI:50477', 'CHEBI:87422', 'CHEBI:31415']
}

GO_mappings = {
    "butyrate" : ['GO:0046358', 'GO:0046359', 'GO:0047761', 'GO:1903544', 'GO:0019605', 'GO:1903545', 'GO:0047371', 'GO:0030645', 'GO:0033508', 'GO:0044581', 'GO:0052638', 'GO:0019672', 'GO:1900500', 'GO:1900502', 'GO:1900501', 'GO:0034338', 'GO:1900751'],
    "tryptophan" : ['GO:0009034', 'GO:0006569', 'GO:0006568', 'GO:0036469', 'GO:0000162', 'GO:0004830', 'GO:0006436']
}

# ALL_METABOLITES = ["indole"]
ALL_METABOLITES = ["butyrate"]#, "acetate", "propionate"] #, "tryptophan", "indole"] # "tryptamine", "3-(1H-indol-3-yl)propanoic acid"
ALL_DIRECTIONS = ["produces"]#, "consumes", "produces_consumes"]

METABOLITES_RELEVANT_TERMS = {
    "butyrate": ["butyrate"],
    "acetate": ["acetate"],
    "propionate": ["propionate"],
    "tryptophan": ["tryptophan", "L-tryptophan","L-tryptophan zwitterion"],
    "indole": ["indole", "1H-indole","3H-indole"]
}

# Traits relevant files
ORGANISMAL_TRAITS_ANNOTATIONS_FILE = "NCBI_organismal_traits"
ORGANISMAL_TRAITS_STRAINS_ANNOTATIONS_FILE = "NCBI_organismal_traits_strains"
ORGANISMAL_TRAITS_STRAINS_PROTEOMES_ANNOTATIONS_FILE = "NCBI_organismal_traits_strains_proteomes"

# UniProt relevant files
RHEA_CHEBI_ANNOTATIONS_FILE = "NCBI_genomic_traits_RHEA"
ORGANISMAL_RHEA_CHEBI_OVERLAP_SPECIES_FILE = "NCBI_organismal_genomic_rhea_comparison_species"
ORGANISMAL_RHEA_CHEBI_OVERLAP_STRAINS_FILE = "NCBI_organismal_genomic_rhea_comparison_strains"
RHEA_ALL_PARTICIPANTS = "rhea_all_participants"
EC_ANNOTATIONS_FILE_SUBSTRING = "NCBI_genomic_traits_EC_pathway_"

ALL_COMPETENCIES_DF_FILE = "All_Competencies"

RANDOM_SEED = 12

MODEL_PARAMETERS = {
        'iterations': 10000,
        'learning_rate': 0.05,
        'l2_leaf_reg': 4,
        'bagging_temperature': 1,
        'random_strength': 6,
        'loss_function': 'Logloss',
        'random_seed': RANDOM_SEED,
        'verbose': 100,
        'early_stopping_rounds': 50,
        'use_best_model': True,
    }

GOLD_STANDARD_FILES = {
    "butyrate_produces": "Input_Files/Vital_etal_butyrate+producing_microbes.csv"
}

COMPETENCY_DISEASE_MAP = {
    "IBD": ["MONDO:0005011", "MONDO:0005265", "MONDO:0005101"],
    "PD": ["MONDO:0005180"]
}

GUT_PHYLA_LIST = ["Bacteroidota", "Bacillota", "Actinomycetota", "Pseudomonadota", "Fusobacteriota"]

GUT_FAMILIES_LIST = ["Ruminococcaceae", "Lachnospiraceae", "Bacteroidaceae", "Prevotellaceae", "Clostridaceae", "Oscillospiraceae", "Eubacteriaceae"]