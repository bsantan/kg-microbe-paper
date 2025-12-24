from collections import Counter, defaultdict
import duckdb
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import chi2_contingency
from duckdb_utils import duckdb_load_table
from classification_utils import differentiate_edge_direction, remove_conflicting_directionality, subset_by_features
from constants import COMPETENCY_DISEASE_MAP
from ncbi_phylogeny_search import find_microbes_strain, get_all_ranks, search_lower_subclass_phylogeny

def get_disease_pairs(output_dir, data_edges, data_nodes, disease_id,disease_name): #MONDO:0005180, #"MONDO:0005011", "MONDO:0005265", "MONDO:0005101"

    # Remove NaN rows
    data_edges = data_edges[~data_edges['predicate'].apply(lambda x: isinstance(x, float))]
    data_edges = data_edges[~data_edges['subject'].apply(lambda x: isinstance(x, float))]
    data_edges = data_edges[~data_edges['object'].apply(lambda x: isinstance(x, float))]

    data_edges.replace(to_replace="NCBITaxon:Actinomyces sp. ICM47", value="NCBITaxon:936548", inplace=True)

    data = differentiate_edge_direction(data_edges)

    data_pairs = data[['subject','object']].drop_duplicates() 
    
    # Target Class
    # data_subset = subset_by_features(generic_input_dir, data_pairs, "NCBITaxon:", True, "MONDO:0005011", "MONDO:0005265", "MONDO:0005101")
    data_subset = subset_by_features(output_dir, data_pairs, "NCBITaxon:", disease_id, True)

    print("len original")
    print(len(data_subset))
    data_pairs_cleaned = remove_conflicting_directionality(data_subset)
    # 239 bugs leftover after this for IBD

    print(len(data_pairs_cleaned))
    data_pairs_cleaned.to_csv(output_dir + "/outcome_to_NCBITaxon_cleaned_" + disease_name + ".tsv", sep="\t", header=True, index=False)

    return data_pairs_cleaned

def combine_labels(taxonomic_labels_counts, keep_label, new_label):

    try:
        taxonomic_labels_counts[keep_label] = taxonomic_labels_counts[keep_label] + taxonomic_labels_counts[new_label]
        del taxonomic_labels_counts[new_label]
    except KeyError:
        return taxonomic_labels_counts

    return taxonomic_labels_counts

def main():

    # feature_table = pd.read_csv("./src/Intermediate_Files_func/feature_table.csv", index_col="subject")
    # disease_microbes = feature_table.index.unique().to_list()

    data_edges = pd.read_csv("./src/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", header=0, sep="\t")
    data_nodes = pd.read_csv("./src/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", header=0, sep="\t")

    output_dir = "./src/Intermediate_Files"

    original_taxon_labels_counts = []
    original_taxon_labels = []

    representative_taxon_labels_counts = []
    representative_taxon_labels = []

    for disease_id in COMPETENCY_DISEASE_MAP.values():

        disease_name = [k for k,v in COMPETENCY_DISEASE_MAP.items() if v == disease_id][0]

        print(disease_name)

        disease_microbes_df = get_disease_pairs(output_dir, data_edges, data_nodes, disease_id, disease_name)
        disease_microbes = disease_microbes_df["subject"].unique().tolist()
        
        # Get all bugs that are in disbiome disease set and gold standard analysis (butyrate produces)
        butyrate_production_output_dir = "./src/Intermediate_Files_Competencies/butyrate_produces"

        # Get all microbes annotated to butyrate production in KG using Vital et al analysis file
        gs_analysis_microbes_file = butyrate_production_output_dir + '/Gold_Standard_Species_Overlap_butyrate_produces.csv'
        gs_analysis_microbes_df = pd.read_csv(gs_analysis_microbes_file)
        kg_analysis_microbes = gs_analysis_microbes_df[(gs_analysis_microbes_df['Organismal'] == 1) | (gs_analysis_microbes_df['Functional'] == 1) | (gs_analysis_microbes_df['Functional_EC'] == 1)]["Value"].tolist()

        gs_analysis_disease_overlap = set(kg_analysis_microbes) & set(disease_microbes)
        print(len(gs_analysis_disease_overlap))

        phylogeny_output_dir = "./src/Phylogeny_Search"
        ncbi_taxa_ranks_df = get_all_ranks(phylogeny_output_dir)

        disease_microbes_ranks = []
        disease_microbes_disease_relationship = []
        # disease_microbes_children_numbers = []
        # disease_microbes_butyrate_producers = []
        disease_microbes_species_numbers = []
        disease_microbes_species_butyrate_producers = []
        disease_microbes_strains_numbers = []
        disease_microbes_strains_butyrate_producers = []

        conn = duckdb.connect(":memory:")
        duckdb_load_table(conn, "./src/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv", "ncbitaxon_edges", ["subject", "predicate", "object"])
        # duckdb_load_table(conn, "./src/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", "ncbitaxon_edges", ["subject", "predicate", "object"])

        microbes_strain_dict, microbes_species_dict = find_microbes_strain(conn, ncbi_taxa_ranks_df, disease_microbes, output_dir, "classification_butyrate_produces_" + disease_name)

        for microbe in disease_microbes: #["NCBITaxon:853"]:#tqdm.tqdm(relevant_ncbitaxa):
            ibd_relationship = disease_microbes_df.loc[disease_microbes_df["subject"] == microbe, "object"].values[0]
            disease_microbes_disease_relationship.append(ibd_relationship)
            microbe_rank = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].get(microbe, "not_found")
            disease_microbes_ranks.append(microbe_rank)
            # # Search traits from only children
            # children = search_lower_subclass_phylogeny(conn, microbe)
            # if children == "not found": children = []
            # # See how many of children are butyrate producers in graph
            # children_with_trait = set(kg_analysis_microbes) & set(children)
            # disease_microbes_children_numbers.append(len(children))
            # disease_microbes_butyrate_producers.append(len(children_with_trait))
            
            if microbe_rank != "species":
                # Search traits from species
                ibd_species = microbes_species_dict.get(microbe, [])
                #ibd_strains.append([key for key, values in microbes_strain_dict.items() if microbe in values][0])
                if ibd_species == ['not found']: ibd_species = []
            elif microbe_rank == "species":
                ibd_species = [microbe]
            # See how many of strains are butyrate producers in graph
            species_with_trait = set(kg_analysis_microbes) & set(ibd_species)
            disease_microbes_species_numbers.append(len(ibd_species))
            disease_microbes_species_butyrate_producers.append(len(species_with_trait))

            if microbe_rank not in ["strain", "subspecies"]:
                # Search traits from strains
                ibd_strains = microbes_strain_dict.get(microbe, [])
                #ibd_strains.append([key for key, values in microbes_strain_dict.items() if microbe in values][0])
                if ibd_strains == ['not found']: ibd_strains = []
            elif microbe_rank in ["strain", "subspecies"]:
                ibd_strains = [microbe]
            # See how many of strains are butyrate producers in graph
            strains_with_trait = set(kg_analysis_microbes) & set(ibd_strains)
            disease_microbes_strains_numbers.append(len(ibd_strains))
            disease_microbes_strains_butyrate_producers.append(len(strains_with_trait))

            # ibd_microbe_children[microbe].extend(strains)

        ranks_df = pd.DataFrame(columns=["Name","Rank","Num_Strains","Num_Strains_Butyrate_Producers"]) #"Num_Children","Num_Butyrate_Producers",
        ranks_df["Name"] = disease_microbes
        ranks_df["Disease_Relationship"] = disease_microbes_disease_relationship
        ranks_df["Rank"] = disease_microbes_ranks
        # ranks_df["Num_Children"] = disease_microbes_children_numbers
        # ranks_df["Num_Butyrate_Producers"] = disease_microbes_butyrate_producers
        ranks_df["Num_Species"] = disease_microbes_species_numbers
        ranks_df["Num_Species_Butyrate_Producers"] = disease_microbes_species_butyrate_producers
        ranks_df["Num_Strains"] = disease_microbes_strains_numbers
        ranks_df["Num_Strains_Butyrate_Producers"] = disease_microbes_strains_butyrate_producers

        ranks_df.to_csv(output_dir + "/outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_" + disease_name + ".csv",sep=",",index=False)


        num_increased_total = 0
        num_decreased_total = 0
        num_increased_disease_butyrate_producers = 0
        num_decreased_disease_butyrate_producers = 0
        total_strains = 0
        total_species = 0

        for i in range(len(ranks_df)):
            if ranks_df.iloc[i].loc["Num_Strains"] > 0:
                num_total = ranks_df.iloc[i].loc["Num_Strains"]
                num_producers = ranks_df.iloc[i].loc["Num_Strains_Butyrate_Producers"]
                total_strains += num_total
            elif ranks_df.iloc[i].loc["Num_Strains"] == 0:
                num_total = ranks_df.iloc[i].loc["Num_Species"]
                num_producers = ranks_df.iloc[i].loc["Num_Species_Butyrate_Producers"]
                total_species += num_total
            if "increased" in ranks_df.iloc[i].loc["Disease_Relationship"]:
                num_increased_total += num_total
                num_increased_disease_butyrate_producers += num_producers
            elif "decreased" in ranks_df.iloc[i].loc["Disease_Relationship"]:
                num_decreased_total += num_total
                num_decreased_disease_butyrate_producers += num_producers

        # Chi Square test 1s then 0s
        contingency_table = [
            [num_increased_disease_butyrate_producers, num_increased_total-num_increased_disease_butyrate_producers],
            [num_decreased_disease_butyrate_producers, num_decreased_total-num_decreased_disease_butyrate_producers]
        ]
        chi2, pval, dof, expected = chi2_contingency(contingency_table)

        classification_summary_df = pd.DataFrame()
        classification_summary_df["Num_Decreased_disease_Total"] = [num_decreased_total]
        classification_summary_df["Num_Increased_disease_Total"] = [num_increased_total]
        classification_summary_df["Num_Decreased_disease_Butyrate_Producers"] = [num_decreased_disease_butyrate_producers]
        classification_summary_df["Num_Increased_disease_Butyrate_Producers"] = [num_increased_disease_butyrate_producers]
        classification_summary_df["chi2"] = [chi2]
        classification_summary_df["P_Val"] = [pval]

        classification_summary_file = output_dir + '/' + disease_name + '_Classification_butyrate_producers_summary.csv'
        
        classification_summary_df.to_csv(classification_summary_file,sep=",",index=False)

        # For testing
        # ranks_df = pd.read_csv(output_dir + "/outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_" + disease_name + ".csv",sep=',')
        # disease_microbes_ranks = ranks_df["Rank"].tolist()
        # print(len(disease_microbes_ranks))

        # Keep track of original ranks
        taxonomic_labels_counts = Counter(disease_microbes_ranks)
        taxonomic_labels_counts = combine_labels(taxonomic_labels_counts, "strain", "subspecies")
        taxonomic_labels_counts = combine_labels(taxonomic_labels_counts, "genus", "NCBITaxon#_species_group")
        taxonomic_labels = list(taxonomic_labels_counts.keys())
        original_taxon_labels_counts.append(taxonomic_labels_counts)
        original_taxon_labels.append(taxonomic_labels)

        # Keep track of representative species or strains ranks
        new_labels_counts = {'species' : total_species, 'strain' : total_strains}
        # For testing
        # new_labels_counts = {'species' : 3567, 'strain' : 9543}
        representative_taxon_labels_counts.append(new_labels_counts)
        representative_taxon_labels.append(list(new_labels_counts.keys()))

    # Create a piechart
    # Note NCBITaxon:335058 doesn't exist but perhaps is NCBITaxon:93681
    combined_original_taxon_labels = list(set([item for sublist in original_taxon_labels for item in sublist]))
    order_list = ['phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'not_found']
    combined_original_taxon_labels = sorted(combined_original_taxon_labels, key=lambda x: (x not in order_list, order_list.index(x) if x in order_list else float('inf')))
    representative_taxon_labels = sorted(representative_taxon_labels, key=lambda x: (x not in order_list, order_list.index(x) if x in order_list else float('inf')))
    # Set colors (same color for the same label across pie charts)
    colors = plt.cm.Paired(range(len(combined_original_taxon_labels)))
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    for i,disease_name in enumerate(COMPETENCY_DISEASE_MAP.keys()):
        # For original taxa
        wedges, texts = ax[0, i].pie([original_taxon_labels_counts[i].get(label, 0) for label in combined_original_taxon_labels], labels=None, autopct=None, startangle=90, colors = colors, wedgeprops={'edgecolor': 'black'}) #autopct='%1.0f%%'
        ax[0, i].axis('equal')
        ax[0, i].set_title(disease_name + ": " + str(sum(original_taxon_labels_counts[i].values())) + " total", fontsize=18)
        
        # For representative taxa
        wedges, texts = ax[1, i].pie([representative_taxon_labels_counts[i].get(label, 0) for label in combined_original_taxon_labels], labels=None, autopct=None, startangle=90, colors = colors, wedgeprops={'edgecolor': 'black'}) #autopct='%1.0f%%'
        ax[1, i].axis('equal')
        ax[1, i].set_title(disease_name + ": " + str(sum(representative_taxon_labels_counts[i].values())) + " total", fontsize=18)


    fig.legend(wedges, combined_original_taxon_labels, title="Rank", loc="upper left", fontsize=12, bbox_to_anchor=(0.8, 0.8))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    fig.savefig(output_dir + '/All_Disease_Comparison_Ranks.png')
    plt.close()

if __name__ == '__main__':
    main()    

