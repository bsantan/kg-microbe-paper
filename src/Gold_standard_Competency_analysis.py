




from collections import defaultdict
import json
import math
import os
import duckdb
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import tqdm

from duckdb_utils import duckdb_load_table, get_node_label
from Competencies import convert_to_species
from constants import ALL_DIRECTIONS, ALL_METABOLITES, GUT_FAMILIES_LIST, GUT_PHYLA_LIST
from ncbi_phylogeny_search import find_microbes_family, find_microbes_phylum, find_microbes_rank, find_microbes_strain, get_all_ranks, get_ncbitaxon_with_uniprot

import plotly.express as px
import plotly.io as pio
from scipy.stats import binomtest 

def create_treemap(conn, df, filename, metab_direction_output_dir):

    df['Family'] = df['Family'].apply(lambda x: get_node_label(conn, x))
    df['Phylum'] = df['Phylum'].apply(lambda x: get_node_label(conn, x))

    phyla_colors = {'human': 'rgba(255,0,0,0.2)', 'other': 'rgba(0,0,255,0.2)'}

    # Create the edges DataFrame with levels and values
    edges = (
        df.groupby(['Phylum', 'Family', 'Location']) #, 'Value'])
        .size()  # Count occurrences
        .reset_index(name='value')  # Add count as 'value'
    )

    # Plot treemap with 3 levels
    fig = px.treemap(
        edges,
        path=['Phylum', 'Family'], #, 'Species_or_Strain'],  # Define the 3 levels of hierarchy
        values='value',                          # Optional: size of the nodes
        title="Distribution of Phyla and Families",
        color='Location',              # Use the Color column for node color
        color_discrete_map=phyla_colors  # Semi-transparent colors
    )

    # Add the numbers inside the boxes and customize outlines
    fig.update_traces(
        textinfo='label+value',
    )

    # Optionally, save as an image (e.g., PNG)
    pio.write_image(fig, metab_direction_output_dir + "/" + filename + "_ncbi_taxonomy_treemap.png", format="png", scale=2)

def combine_microbes_dicts(dict1, dict2, output_dir, new_filename):

    if not os.path.exists(output_dir + "/" + new_filename):
        # Combine the dictionaries
        combined_dict = defaultdict(list)
        for key, value in dict1.items():
            combined_dict[key].extend(value)
        for key, value in dict2.items():
            combined_dict[key].extend(value)
        with open(output_dir + "/" + new_filename, 'w') as json_file:
            json.dump(combined_dict, json_file, indent=4)
    else:
        with open(output_dir + "/" + new_filename, 'r') as json_file:
            data = json.load(json_file)
            combined_dict = defaultdict(list, data)

    return combined_dict

def create_ordered_subset( df, filename, all_microbes_families_species, grouped_rank, group_by_rank, threshold = 1.0):

    print(filename)
    # Remove families with no children
    all_microbes_families_species = {key: value for key, value in all_microbes_families_species.items() if len(value) > 0}

    ranked_value_dict = df.groupby([grouped_rank.capitalize()])['Value'].apply(list).to_dict()
    ranked_mapping = df[[grouped_rank.capitalize(), group_by_rank.capitalize()]].drop_duplicates().set_index(grouped_rank.capitalize())[group_by_rank.capitalize()].to_dict()
    # Now, sort ranked_value_dict by phylum
    ranked_value_dict = {k: ranked_value_dict[k] for k in sorted(ranked_value_dict, key=lambda x: ranked_mapping.get(x, ''))}

    print("len families dict with non proteomes")
    print(len(ranked_value_dict.items()))
    ranked_value_dict = {key: value for key, value in ranked_value_dict.items() if key in all_microbes_families_species.keys()}
    print("len families dict with only proteomes")
    print(len(ranked_value_dict.items()))

    if not ranked_value_dict:  # Check if the input dict is empty
        threshold_ranked_value_dict = {}
    else:
        lengths = {key: len(value) for key, value in ranked_value_dict.items()}
        threshold_index = max(0, int(len(lengths) * threshold) - 1)  # Avoid out-of-range index
        threshold_value = sorted(lengths.values(), reverse=True)[threshold_index]
        threshold_ranked_value_dict = {key: value for key, value in ranked_value_dict.items() if len(value) >= threshold_value}
        
    return threshold_ranked_value_dict, ranked_mapping

def create_families_piechart(conn, metabolite, direction, filename, all_microbes_families_species, child_rank, ncbitaxon_func_ids, output_dir, grouped_rank, group_by_rank, threshold_ranked_value_dict, ranked_mapping):

    if len(threshold_ranked_value_dict.items()) > 0:
        ranked_microbe_labels = []
        ranked_microbe_ids = []
        ranked_labels = []
        fractions_produces_but = []
        total_fam_sizes = []
        traits_data = []

        # Determine grid dimensions (roughly square)
        num_cols = math.ceil(math.sqrt(len(threshold_ranked_value_dict.keys())))
        num_rows = math.ceil(len(threshold_ranked_value_dict.keys()) / num_cols)
        fig, ax = plt.subplots(num_rows, num_cols, figsize=(16, 12))
        # Ensure ax is always an array
        # Flatten the axes array for easy indexing, will be array only if num_cols is >1
        if isinstance(ax, np.ndarray):
            ax = ax.flatten()
        else:
            ax = [ax]

        for i, (ranked_microbe, values) in enumerate(threshold_ranked_value_dict.items()):
            # Subset by only values in proteomes
            values = [j for j in values if j in ncbitaxon_func_ids]
            ranked_microbe_label = get_node_label(conn, ranked_microbe)
            ranked_labels.append(ranked_mapping[ranked_microbe])
            if ranked_microbe == "not found": 
                not_found_children += 1
            else:
                # Remove values that did not have rank
                all_children = all_microbes_families_species[ranked_microbe]
                values = [j for j in values if j in all_children]
                total_children = len(all_microbes_families_species[ranked_microbe])
                with_trait = len(values)
                without_trait = total_children - with_trait
                ranked_microbe_labels.append(ranked_microbe_label)
                ranked_microbe_ids.append(ranked_microbe)
                fractions_produces_but.append(with_trait / total_children)
                total_fam_sizes.append(total_children)
                traits_data.append([with_trait, without_trait])

        colors = ["#56B4E9","#9C27B0"]
        for m in range(len(ranked_microbe_labels)):
            # For original taxa
            wedges, texts = ax[m].pie(traits_data[m], labels=None, autopct=None, startangle=90, wedgeprops={'edgecolor': 'black'}, colors = colors) #autopct='%1.0f%%'
            ax[m].axis('equal')
            ax[m].set_title(ranked_microbe_labels[m] + ": " + str(total_fam_sizes[m]) + " total", fontsize=12)
        # Hide unused subplots
        phrase = metabolite.capitalize() + " " + direction.capitalize() + "r"
        traits_labels = [phrase, "Non " + phrase]
        # traits_labels = ["Butyrate Producer", "Non Butyrate Producer"]
        for j in range(i + 1, len(ax)):
            fig.delaxes(ax[j])
        fig.legend(wedges, traits_labels, title="Microbial Trait", loc="upper left", fontsize=12, bbox_to_anchor=(0.6, 0.2))
        plt.tight_layout(rect=[0, 0, 0.95, 1])
        fig.savefig(output_dir + '/' + filename + '_' + grouped_rank.capitalize() + '_Traits_Comparison_' + child_rank + '.png')
        plt.close()

        # Calculate total number of 1s and total observations
        total_but_producers = int(sum(p * n for p, n in zip(fractions_produces_but, total_fam_sizes)))
        total_observations = sum(total_fam_sizes)

        # Perform binomial test
        p_value = binomtest(total_but_producers, total_observations, p=0.5, alternative='greater')
        families_results_stats_df = pd.DataFrame()
        families_results_stats_df["Total " + phrase + "s"] = [total_but_producers]
        families_results_stats_df["Total_Children"] = [total_observations]
        families_results_stats_df["Binomial_p_value"] = [p_value]
        families_results_stats_df.to_csv(output_dir + '/' + filename + '_' + grouped_rank.capitalize() + '_Traits_Comparison_Statistics_' + child_rank + '.tsv')

        # Output to tsv
        families_results_df = pd.DataFrame()
        families_results_df[grouped_rank.capitalize()] = ranked_microbe_labels
        families_results_df[grouped_rank.capitalize() + "_ID"] = ranked_microbe_ids
        families_results_df["Ranked_Mapping_" + group_by_rank] = ranked_labels
        families_results_df["fractions_produces_but"] = fractions_produces_but
        families_results_df["total_fam_sizes"] = total_fam_sizes
        families_results_df.to_csv(output_dir + '/' + filename + '_' + grouped_rank.capitalize() + '_Traits_Comparison_' + child_rank + '.tsv')

def get_rank(microbe, all_microbes_df, ncbi_taxa_ranks_df):

    if microbe in all_microbes_df["Value"].values:
        rank = all_microbes_df.loc[all_microbes_df["Value"] == microbe, "Rank"].values[0]
    else:
        rank = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].get(microbe, None)
    
    return rank

def post_competency_analysis(conn, metabolite, direction, microbes_family_dict, microbes_phylum_dict, microbes_genus_dict, ncbi_taxa_ranks_df, microbes_list, metab_direction_output_dir, filename, all_microbes_df):

    new_filename = metab_direction_output_dir + "/" + filename + "_" + metabolite + "_" + direction + "_families.tsv"

    if not os.path.exists(new_filename):

        new_df = pd.DataFrame()
        rank_list = []
        microbes_list_genera = []
        microbes_list_families = []
        microbes_list_phyla = []
        microbe_location = []
        microbes_impt_families = []

        for m in tqdm.tqdm(microbes_list):
            rank = get_rank(m, all_microbes_df, ncbi_taxa_ranks_df)
            rank_list.append(rank)
            gen = next((k for k, v in microbes_genus_dict.items() if m in v), None)
            fam = next((k for k, v in microbes_family_dict.items() if m in v), None)
            phy = next((k for k, v in microbes_phylum_dict.items() if m in v), None)
            if gen:
                microbes_list_genera.append(gen)
            else:
                microbes_list_genera.append("none")
                #microbes_impt_families.append("none")
            if fam:
                fam_lab = get_node_label(conn, fam)
                impt_fam = "impt_fam" if fam_lab in GUT_FAMILIES_LIST else "other"
                microbes_list_families.append(fam)
                microbes_impt_families.append(impt_fam)
            else:
                microbes_list_families.append("none")
                microbes_impt_families.append("none")
            if phy:
                lab = get_node_label(conn, phy)
                loc = "human" if lab in GUT_PHYLA_LIST else "other"
                microbes_list_phyla.append(phy)
                microbe_location.append(loc)
            else:
                microbes_list_phyla.append("none")
                microbe_location.append("none")

        new_df["Value"] = microbes_list
        new_df["Rank"] = rank_list
        new_df["Genus"] = microbes_list_genera
        new_df["Family"] = microbes_list_families
        new_df["Phylum"] = microbes_list_phyla
        new_df["Location"] = microbe_location
        new_df["Impt_Family"] = microbes_impt_families

        new_df = new_df.sort_values(by=["Phylum","Location"])
        new_df.to_csv(new_filename,sep='\t',index=False)

    else:
        new_df = pd.read_csv(new_filename,sep='\t')

    return new_df

def get_species_overlap(conn, metabolite, direction, gs_analysis_microbes_df, metab_direction_output_dir):

    # Get all microbes from KG
    all_values = gs_analysis_microbes_df[(gs_analysis_microbes_df['Organismal'] == 1) | (gs_analysis_microbes_df['Functional'] == 1) | (gs_analysis_microbes_df["Functional_EC"] == 1) |  (gs_analysis_microbes_df["Gold_Standard"] == 1)]["Value"].tolist()

    gold_standard_list = gs_analysis_microbes_df[(gs_analysis_microbes_df['Gold_Standard'] == 1)]["Value"].tolist()

    species_overlap_list = []

    # Create dictionary of species: [strains]
    all_values_species_dict = convert_to_species(conn, all_values, "all_kg_gs_values_" + metabolite + "_" + direction)

    new_df = pd.DataFrame()

    for l in all_values:
        overlapping_species = []
        for st in gold_standard_list:
            for key, values in all_values_species_dict.items(): #gold_standard_species_dict
                # If strain from gold standard has species that is also in kg
                if st in values and key in l:
                    overlapping_species.append(l)

        species_overlap_list.append(len(overlapping_species))

    new_df["KG_Microbe"] = all_values
    new_df["Species_Overlap_GS"] = species_overlap_list

    new_df.to_csv(metab_direction_output_dir + "/overlapping_kg_gs_species_" + metabolite + "_" + direction + ".tsv",sep='\t',index=False)

    return new_df

def main():

    for direction in ALL_DIRECTIONS:
        for metabolite in ALL_METABOLITES:
            metab_direction_output_dir = "./Intermediate_Files_Competencies/" + metabolite + "_" + direction

            gs_analysis_microbes_file = metab_direction_output_dir + '/Gold_Standard_Species_Overlap_'+ metabolite + "_" + direction + '.csv'
            gs_analysis_microbes_df = pd.read_csv(gs_analysis_microbes_file)
            
            phylogeny_output_dir = "./Phylogeny_Search"
            ncbi_taxa_ranks_df = get_all_ranks(phylogeny_output_dir)

            conn = duckdb.connect(":memory:")
            duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv", "edges", ["subject", "predicate", "object"])
            duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])

            # First create list of colors that will be used in every treemap to align family boxes
            all_microbes_list = gs_analysis_microbes_df["Value"].unique().tolist()
            # Get genus of each bug
            microbes_genus_dict = find_microbes_rank(conn, ncbi_taxa_ranks_df, all_microbes_list, phylogeny_output_dir, "/Gold_Standard_Species_Overlap_all_"+ metabolite + "_" + direction, "genus")
            # Get family of each bug
            microbes_family_dict = find_microbes_rank(conn, ncbi_taxa_ranks_df, all_microbes_list, phylogeny_output_dir, "/Gold_Standard_Species_Overlap_all_"+ metabolite + "_" + direction, "family")
            # Get phylum of each bug
            microbes_phylum_dict = find_microbes_rank(conn, ncbi_taxa_ranks_df, all_microbes_list, phylogeny_output_dir, "/Gold_Standard_Species_Overlap_all_"+ metabolite + "_" + direction, "phylum")

            all_microbes_df = post_competency_analysis(conn, metabolite, direction, microbes_family_dict, microbes_phylum_dict, microbes_genus_dict, ncbi_taxa_ranks_df, all_microbes_list, metab_direction_output_dir, "all_microbes", pd.DataFrame(columns = ["Value"]))
            families_list = all_microbes_df['Family'].unique().tolist()
            families_list = sorted(families_list)

            genera_list = all_microbes_df['Genus'].unique().tolist()
            genera_list = sorted(genera_list)

            # Get all species & strains from each family found
            microbes_strain_dict, microbes_species_dict = find_microbes_strain(conn, ncbi_taxa_ranks_df, families_list, "./Intermediate_Files", "competencies_all_microbes_families_"+ metabolite + "_" + direction)

            # Only include bugs that have a proteome in the graph
            ncbitaxon_func_ids = get_ncbitaxon_with_uniprot(conn, "./Phylogeny_Search")

            # Combine species and strain dictionaries
            microbes_species_and_strain_dict = combine_microbes_dicts(microbes_species_dict, microbes_strain_dict, "./Intermediate_Files", "competencies_all_microbes_families_"+ metabolite + "_" + direction + "_microbes_strains_and_species.json")
            filtered_microbes_species_and_strain_dict = {
                key: [value for value in values if value in ncbitaxon_func_ids]
                for key, values in microbes_species_and_strain_dict.items()
            }

            microbial_subsets = {
                "gs_no_kg" : (gs_analysis_microbes_df['Organismal'] == 0) & (gs_analysis_microbes_df['Functional'] == 0) & (gs_analysis_microbes_df['Functional_EC'] == 0) & (gs_analysis_microbes_df['Gold_Standard'] == 1),
                "proteomes_no_gs" : (gs_analysis_microbes_df['Proteome'] == 1) & (gs_analysis_microbes_df['Gold_Standard'] == 0),
                "gs_and_kg": ((gs_analysis_microbes_df['Organismal'] == 1) | (gs_analysis_microbes_df['Functional'] == 1) | (gs_analysis_microbes_df["Functional_EC"] == 1)) & (gs_analysis_microbes_df['Gold_Standard'] == 1),
                "organismal_no_gs" : (gs_analysis_microbes_df['Organismal'] == 1) & (gs_analysis_microbes_df['Gold_Standard'] == 0),
                "rhea_chebi_no_gs" : (gs_analysis_microbes_df['Functional'] == 1) & (gs_analysis_microbes_df['Gold_Standard'] == 0),
                "ec_no_gs" : (gs_analysis_microbes_df['Functional_EC'] == 1) & (gs_analysis_microbes_df['Gold_Standard'] == 0),
                "all_no_gs" : gs_analysis_microbes_df['Gold_Standard'] == 0,
                "all_kg" : ((gs_analysis_microbes_df['Organismal'] == 1) | (gs_analysis_microbes_df['Functional'] == 1) | (gs_analysis_microbes_df["Functional_EC"] == 1)),
                "gs_proteomes_no_annotation" : (gs_analysis_microbes_df['Proteome'] == 1) & (gs_analysis_microbes_df['Gold_Standard'] == 1) & (gs_analysis_microbes_df['Functional'] == 0) & (gs_analysis_microbes_df['Functional_EC'] == 0),
                "gs_proteomes_total" : (gs_analysis_microbes_df['Proteome'] == 1) & (gs_analysis_microbes_df['Gold_Standard'] == 1),
                "proteomes_probable_protein_name" : gs_analysis_microbes_df['Functional_Protein_Name'].str.contains(r'\bprobable\b', case=False),
            }

            final_data = {}

            output_dir = metab_direction_output_dir + "/Competency_Analysis"
            os.makedirs(output_dir, exist_ok=True)

            for microbial_subset, constraint in microbial_subsets.items():
                subset_list = gs_analysis_microbes_df[constraint]["Value"].tolist()
                df = post_competency_analysis(conn, metabolite, direction, microbes_family_dict, microbes_phylum_dict, microbes_genus_dict, ncbi_taxa_ranks_df, subset_list, output_dir, microbial_subset, all_microbes_df)
                threshold_ranked_value_dict, family_mapping = create_ordered_subset(df, microbial_subset, filtered_microbes_species_and_strain_dict, "genus", "family", 0.1)
                create_families_piechart(conn, metabolite, direction, microbial_subset, filtered_microbes_species_and_strain_dict, "species_and_strain_all", ncbitaxon_func_ids, output_dir, "genus", "family",threshold_ranked_value_dict, family_mapping)
                create_treemap(conn, df, microbial_subset, output_dir)

                # # Only include families in human gut phyla
                # df_human = df.loc[df["Location"] == "human"]
                # threshold_ranked_value_dict, family_mapping = create_ordered_subset(df_human, microbial_subset, filtered_microbes_species_and_strain_dict, "genus", "family", 0.1)
                # create_families_piechart(conn, microbial_subset, filtered_microbes_species_and_strain_dict, "species_and_strain_human", ncbitaxon_func_ids, output_dir, "genus", threshold_ranked_value_dict, family_mapping)

                # Only include impt families for metabolite trait 
                df_impt = df.loc[df["Impt_Family"] == "impt_fam"]
                for threshold in [0.1, 1.0]:
                    output_subdir = output_dir + "/Threshold_" + str(threshold)
                    os.makedirs(output_subdir, exist_ok=True)
                    threshold_ranked_value_dict, family_mapping = create_ordered_subset(df_impt, microbial_subset, filtered_microbes_species_and_strain_dict, "genus", "family", threshold)
                    create_families_piechart(conn, metabolite, direction, microbial_subset, filtered_microbes_species_and_strain_dict, "species_and_strain_impt_families_" + str(threshold), ncbitaxon_func_ids, output_subdir, "genus", "family", threshold_ranked_value_dict, family_mapping)

                subset_list_human = df.loc[df["Location"] == "human", "Value"].tolist()
                final_data[microbial_subset] = [len(subset_list), len(subset_list_human)]

            # Convert to DataFrame
            final_data_df = pd.DataFrame([(key, values[0], values[1]) for key, values in final_data.items()],
                        columns=["Microbial_Subset", "Num_Non-Human_" + metabolite + "_" + direction, "Num_Human_" + metabolite + "_" + direction])
            final_data_df.to_csv(output_dir + "/competency_summary_" + metabolite + "_" + direction + ".tsv",sep='\t')

if __name__ == '__main__':
    main()    

