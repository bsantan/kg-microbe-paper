

from collections import defaultdict
import csv
import json
import math
import random
import re
import textwrap
import duckdb
from matplotlib.pylab import norm
import numpy as np
import tqdm
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

from duckdb_utils import duckdb_load_table, get_node_label, get_table_count, output_table_to_file
import os
from constants import ALL_COMPETENCIES_DF_FILE, ALL_DIRECTIONS, EC_ANNOTATIONS_FILE_SUBSTRING, GOLD_STANDARD_FILES, METABOLITES_RELEVANT_TERMS, ORGANISMAL_RHEA_CHEBI_OVERLAP_SPECIES_FILE, ORGANISMAL_RHEA_CHEBI_OVERLAP_STRAINS_FILE, ORGANISMAL_TRAITS_ANNOTATIONS_FILE, ORGANISMAL_TRAITS_STRAINS_ANNOTATIONS_FILE, ORGANISMAL_TRAITS_STRAINS_PROTEOMES_ANNOTATIONS_FILE, RHEA_ALL_PARTICIPANTS, RHEA_CHEBI_ANNOTATIONS_FILE, GO_mappings, CHEBI_mappings
import pandas as pd
import sys

from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

from equilibrator_api import ComponentContribution, Q_

from ncbi_phylogeny_search import find_microbes_family, find_microbes_species, get_all_kg_taxa, get_all_ranks, get_ncbitaxon_with_uniprot, get_taxa_per_rank, load_graph
cc = ComponentContribution()

def get_rhea_participants(metabolite):

    directory = "./Intermediate_Files_Competencies/" + metabolite
    df = pd.read_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", sep = "\t")

    rhea_ids = df["rhea"].unique().tolist()
    rhea_duckdb_query =  '(' + ' OR '.join([f"e.subject = '{term}'" for term in rhea_ids]) + ')'

    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])

    query = (
        f"""
        CREATE OR REPLACE TABLE rhea_to_chebi_participants AS
        SELECT *
        FROM edges e
        WHERE {rhea_duckdb_query}
        AND split_part(e.object, ':', 1) = 'CHEBI';
        """
    )

    # print(query)
    conn.execute(query)

    query = (
        f"""
        CREATE OR REPLACE TABLE chebi_atp_participants AS
        SELECT * FROM rhea_to_chebi_participants
        WHERE object = 'CHEBI:30616';
        """
    )

    # print(query)
    conn.execute(query)

    print(metabolite)
    get_table_count(conn, "chebi_atp_participants")
    print(len(rhea_ids))

    query = (
        f"""
        CREATE OR REPLACE TABLE chebi_nadph_participants AS
        SELECT * FROM rhea_to_chebi_participants
        WHERE object = 'CHEBI:16474';
        """
    )

    # print(query)
    conn.execute(query)

    print(metabolite)
    get_table_count(conn, "chebi_nadph_participants")
    print(len(rhea_ids))

    output_table_to_file(conn, "(SELECT DISTINCT * FROM rhea_to_chebi_participants ORDER BY subject DESC)", directory + "/" + RHEA_ALL_PARTICIPANTS + ".tsv")

    # Get chebi labels
    df = pd.read_csv(directory + "/" + RHEA_ALL_PARTICIPANTS + ".tsv", sep = "\t")
    chebi_labels = []
    for chebi in df["object"]:
        lab = get_node_label(conn, chebi)
        chebi_labels.append(lab)
    df["Object_Label"] = chebi_labels

    print("writing to ",directory + "/" + RHEA_ALL_PARTICIPANTS + ".tsv")
    df[["subject", "Object_Label"]].to_csv(directory + "/" + RHEA_ALL_PARTICIPANTS + ".tsv", sep='\t')

def plot_competencies_barplot(df, direction):

    directory = "./Intermediate_Files_Competencies"

    print(df)
    ax = df.plot(kind='bar', x='Metabolite', y=['Rhea-Chebi_Traits_Overlap_Percentage', 'Traits_Rhea-Chebi_Overlap_Percentage'], figsize=(8, 6))

    # Add labels and title
    # ax.set_title('Percentage Coverage Across Sources: Production Only')
    ax.set_title('Percentage Coverage Across Sources: ' + direction.capitalize().replace("_",", "))
    ax.set_xlabel('Metabolite')
    # current_labels = [label.get_text() for label in plt.gca().get_xticklabels()]
    # # new_labels = [label.replace("produces_consumes_", "") for label in current_labels]
    # plt.gca().set_xticklabels(new_labels)
    ax.set_ylabel('Percentage Overlap')
    plt.xticks(rotation=0)
    # ax.legend(['Rhea-Chebi: Traits Overlap (out of all Traits)', 'Traits:Rhea-Chebi Overlap (out of all Rhea-Chebi)'], title='Comparison')
    ax.legend(['Coverage of Traits in Uniprot', 'Coverage of Uniprot in Traits'], title='Comparison')

    # Annotate the bars
    for i in range(len(df)):
        # Annotate the first bar (Coverage of Traits in Uniprot)
        ax.text(i-0.2, df['Rhea-Chebi_Traits_Overlap_Percentage'][i] + 1,  # Adjust the offset (+1) as needed
                f"{int(df['Traits_Annotations_Proteome_Overlap'][i])}",  # Format the value as needed
                ha='center', va='bottom')

        # Annotate the second bar (Coverage of Uniprot in Traits)
        ax.text(i+0.2, df['Traits_Rhea-Chebi_Overlap_Percentage'][i] + 1,  # Adjust the offset (+1) as needed
                f"{int(df['Rhea-Chebi_Annotations'][i])}",  # Format the value as needed
                ha='center', va='bottom')

    plt.savefig(directory + "/All_Competencies_barplot" + "_" + direction + ".png")

def get_total_proteomes_from_graph():

    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", "edges", ["subject", "predicate", "object"])

    # Get only NCBITaxon from organismal traits that have proteomes in genomic traits
    query = (
        f"""
        SELECT DISTINCT e.object
        FROM edges e
        WHERE e.subject LIKE '%UniprotKB:%'
        AND e.object LIKE '%NCBITaxon:%'
        """
    )

    conn.execute(query)
    result = [row[0] for row in conn.execute(query).fetchall()]

    ncbitaxa_with_proteomes = result
    print("ncbitaxa_with_proteomes: ",len(ncbitaxa_with_proteomes))

    return ncbitaxa_with_proteomes

def monte_carlo_simulations(metabolite, direction, comparison_group, actual_overlap, num_inputs,  all_possible_inputs, ax, random_seeds):

    directory = "./Intermediate_Files_Competencies"
    print("all_possible_inputs")
    print(all_possible_inputs)
    print("comparison_group")
    print(comparison_group)
    print("actual_overlap")
    print(actual_overlap)
    print("num_inputs")
    print(num_inputs)
    print("num seeds")
    print(len(random_seeds))

    all_permutations = []
    for seed in random_seeds:
        np.random.seed(seed)
        random_selection = np.random.choice(all_possible_inputs, num_inputs)
        random_selection_overlap = len(list(set([i for i in random_selection if i in comparison_group])))
        all_permutations.append(random_selection_overlap)

    # Compute p-value (one-tailed: greater than or equal to the value)
    p_value = np.sum(np.array(all_permutations) >= actual_overlap) / len(all_permutations)
    # Plot histogram
    # plt.figure()
    ax.hist(all_permutations, density=False, alpha=0.7, color='gray', edgecolor='black', label="Permutations") #bins=30, 

    # Add vertical line for the value of interest
    ax.axvline(actual_overlap, color='red', linestyle='--', linewidth=1.5, label=f"Value = {actual_overlap}")

    # # Annotate p-value on the plot
    ax.text(p_value + 1, plt.gca().get_ylim()[1] * 0.2, f'P-value = {p_value:.4f}', color='red')

    return p_value

def get_microbial_list(filename, col_name):

    df = pd.read_csv(filename, sep="\t")
    l = df[col_name].unique().tolist()

    return l

def plot_competencies_venn_diagrams_with_proteomes(conn):

    directory = "./Intermediate_Files_Competencies"

    # # Get total number of bugs with proteomes
    # total_proteomes = get_total_proteomes_from_graph()
    ncbitaxon_func_ids = get_ncbitaxon_with_uniprot(conn, "./Phylogeny_Search")
    print("Len Taxa with a functional annotation")
    print(len(ncbitaxon_func_ids))
    total_proteomes = ncbitaxon_func_ids
    p_values = {}

    all_dfs = pd.DataFrame()
    for direction in ALL_DIRECTIONS:
        final_df_filename = directory + "/" + ALL_COMPETENCIES_DF_FILE + "_" + direction + "_raw.tsv"
        final_df = pd.read_csv(final_df_filename)
        final_df["Direction"] = direction
        all_dfs = pd.concat([all_dfs,final_df],axis=0)

    num_columns = len(ALL_DIRECTIONS)
    num_rows = math.ceil(len(all_dfs) / num_columns)
    fig_hist, axes_hist = plt.subplots(num_rows, len(ALL_DIRECTIONS), figsize=(12, 8))
    # Ensure axes is an array, even if there's only one row or column
    if not isinstance(axes_hist, (list, np.ndarray)):
        axes_hist = [axes_hist]
    else:
        axes_hist = axes_hist.flatten()
    fig_venn, axes_venn = plt.subplots(num_rows, len(ALL_DIRECTIONS), figsize=(12, 8))
    # Ensure axes is an array, even if there's only one row or column
    if not isinstance(axes_venn, (list, np.ndarray)):
        axes_venn = [axes_venn]
    else:
        axes_venn = axes_venn.flatten()

    num_permutations = 1000
    RANDOM_SEEDS = [random.randint(0, 10000) for _ in range(num_permutations)]

    for i in range(len(all_dfs)):
        metabolite = all_dfs.iloc[i].loc["Metabolite"]
        direction = all_dfs.iloc[i].loc["Direction"]
        organismal_traits = int(all_dfs.iloc[i].loc["Total_Traits_Annotations"])
        organismal_traits_proteomes_list = get_microbial_list(directory + "/" + metabolite + "_" + direction + "/" + ORGANISMAL_TRAITS_STRAINS_PROTEOMES_ANNOTATIONS_FILE + ".tsv", "updated_subject")
        # functional_annotations_proteomes_list = get_microbial_list(directory + "/" + metabolite + "_" + direction + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", "ncbitaxon")
        organismal_traits_proteomes = int(all_dfs.iloc[i].loc["Total_Traits_Proteome_Overlap"])
        if metabolite == "butyrate" and direction == "produces":
            functional_annotations = int(all_dfs.iloc[i].loc["Total_EC_and_Rhea-Chebi"])
            overlap = int(all_dfs.iloc[i].loc["Total_Rhea-Chebi_Traits_Overlap"])
        else:
            functional_annotations = int(all_dfs.iloc[i].loc["Total_Rhea-Chebi_Annotations"])
            overlap = int(all_dfs.iloc[i].loc["Total_Rhea-Chebi_Traits_Overlap"])

        p_value = monte_carlo_simulations(metabolite, direction, organismal_traits_proteomes_list, overlap, functional_annotations, total_proteomes, axes_hist[i], RANDOM_SEEDS)
        p_values[metabolite + "_" + direction] = p_value

        # Create Venn diagram
        venn = venn3(
            subsets=(
                organismal_traits, # A, 100
                functional_annotations, # B, 010
                overlap, # AB, 110
                organismal_traits_proteomes, # C, 001
                organismal_traits, # AC, 101
                overlap, # BC, 011
                overlap # ABC, 111
            ),
            # set_labels=("Organismal Trait", "Functional Trait"),
            set_labels=("", ""),
            ax=axes_venn[i]
        )

        # axes[i].set_title(metabolite.capitalize() + "/ " + direction.capitalize().replace("_",", "))
        # Make labels larger, shift the vertical position of overlap labels for better clarity
        for region_id in ['100', '010', '001']:
            label = venn.get_label_by_id(region_id)
            if label:
                label.set_fontsize(14)
        for region_id in ['111']:
            label = venn.get_label_by_id(region_id)
            if label: 
                label.set_fontsize(14); label.set_y(0.2)
        for region_id in ['101', '011', '110']:
            label = venn.get_label_by_id(region_id)
            if label:
                label.set_text('')


    # Add a legend (key) to explain the colors
    fig_venn.legend(
        loc="upper right",
        labels=["Only Organismal Trait", "Only Functional Trait", "Organismal Trait with Proteome", "Both Traits"],
        handles=[
            venn.get_patch_by_id('100'), # get_patch_by_id
            venn.get_patch_by_id('010'),
            venn.get_patch_by_id('001'),
            venn.get_patch_by_id('111')
        ],
        fontsize=14
    )

    row_labels = [", ".join(word.capitalize() for word in direction.split("_")) for direction in ALL_DIRECTIONS]
    col_labels = [i.capitalize() for i in pd.Index(all_dfs["Metabolite"]).unique().tolist()]
    # Add row labels on the left side of the figure
    for i, label in enumerate(row_labels):
        wrapped_label = "\n".join(textwrap.wrap(label, width=10))
        for f in [fig_hist, fig_venn]:
            f.text(
                0.04, 0.78 - i * 0.33, wrapped_label, 
                va='center', ha='center', 
                rotation='vertical', fontsize=16
            )

    # Add column labels on the top of the figure
    for j, label in enumerate(col_labels):
        wrapped_label = "\n".join(textwrap.wrap(label, width=10))
        fig_hist.text(
            0.17 + j * 0.33, 0.9, wrapped_label, 
            va='center', ha='center', 
            fontsize=16
        )
        fig_venn.text(
            0.17 + j * 0.33, 0.8, wrapped_label, 
            va='center', ha='center', 
            fontsize=16
        )

    # Add p-value labels beneath each Venn diagram
    for i, row in enumerate(ALL_DIRECTIONS):
        for j, col in enumerate(pd.Index(all_dfs["Metabolite"]).unique().tolist()):
            p_val = p_values[col + "_" + row]
            # fig_hist.text(
            #     0.17 + j * 0.33, 0.7 - i * 0.26, f"= {p_val}",  # Adjust y for placement
            #     va='center', ha='center', 
            #     fontsize=12
            # )
            fig_venn.text(
                0.17 + j * 0.33, 0.55 - i * 0.26, f"p = {p_val}",  # Adjust y for placement
                va='center', ha='center', 
                fontsize=12
            )

    # Add title
    fig_venn.suptitle("Coverage of Traits Across Sources",fontsize=18)
    # Adjust spacing manually if needed (optional)
    fig_venn.subplots_adjust(top=0.8, bottom=0.02, left=0.04, right=0.98, wspace=0.1, hspace=0.1)
    fig_venn.savefig(directory + "/All_Competencies_Venn_Diagrams_With_Proteomes.png")
    plt.close(fig_venn)

    # Create histogram figure
    fig_hist.suptitle("Monte Carlo Simulations for Functional Traits")
    fig_hist.savefig(directory + "/Monte_Carlos_Distribution.png")
    plt.close(fig_hist)

def plot_competencies_venn_diagrams():

    directory = "./Intermediate_Files_Competencies"

    all_dfs = pd.DataFrame()
    for direction in ALL_DIRECTIONS:
        final_df_filename = directory + "/" + ALL_COMPETENCIES_DF_FILE + "_" + direction + "_raw.tsv"
        final_df = pd.read_csv(final_df_filename)
        final_df["Direction"] = direction
        all_dfs = pd.concat([all_dfs,final_df],axis=0)

    num_columns = len(ALL_DIRECTIONS)
    num_rows = math.ceil(len(all_dfs) / num_columns)
    fig, axes = plt.subplots(num_rows, len(ALL_DIRECTIONS), figsize=(12, 8))
    # Ensure axes is an array, even if there's only one row or column
    if not isinstance(axes, (list, np.ndarray)):
        axes = [axes]
    else:
        axes = axes.flatten()

    for i in range(len(all_dfs)):
        metabolite = all_dfs.iloc[i].loc["Metabolite"]
        direction = all_dfs.iloc[i].loc["Direction"]
        organismal_traits = int(all_dfs.iloc[i].loc["Total_Traits_Annotations"])
        functional_annotations = int(all_dfs.iloc[i].loc["Total_Rhea-Chebi_Annotations"])
        overlap = int(all_dfs.iloc[i].loc["Total_Rhea-Chebi_Traits_Overlap"])

        # Create Venn diagram
        venn = venn2(
            subsets=(organismal_traits, functional_annotations, overlap),
            # set_labels=("Organismal Trait", "Functional Trait"),
            set_labels=("", ""),
            ax=axes[i]
        )

        # axes[i].set_title(metabolite.capitalize() + "/ " + direction.capitalize().replace("_",", "))
        # Make labels larger, shift the vertical position of overlap labels for better clarity
        for region_id in ['10', '01']:
            label = venn.get_label_by_id(region_id)
            if label:
                label.set_fontsize(14)
        label = venn.get_label_by_id('11')
        if label: label.set_fontsize(14); label.set_y(0.2)

    # Add a legend (key) to explain the colors
    fig.legend(
        loc="upper right",
        labels=["Only Organismal Trait", "Only Functional Trait", "Both Traits"],
        handles=[
            venn.get_patch_by_id('10'), # get_patch_by_id
            venn.get_patch_by_id('01'),
            venn.get_patch_by_id('11')
        ],
        fontsize=14
    )

    # col_labels = [i.capitalize().replace("_",", ") for i in ALL_DIRECTIONS]
    row_labels = [", ".join(word.capitalize() for word in direction.split("_")) for direction in ALL_DIRECTIONS]
    col_labels = [i.capitalize() for i in pd.Index(all_dfs["Metabolite"]).unique().tolist()]
    # Add row labels on the left side of the figure
    for i, label in enumerate(row_labels):
        wrapped_label = "\n".join(textwrap.wrap(label, width=10))
        fig.text(
            0.04, 0.78 - i * 0.33, wrapped_label, 
            va='center', ha='center', 
            rotation='vertical', fontsize=16
        )

    # Add column labels on the top of the figure
    for j, label in enumerate(col_labels):
        wrapped_label = "\n".join(textwrap.wrap(label, width=10))
        fig.text(
            0.17 + j * 0.33, 0.85, wrapped_label, 
            va='center', ha='center', 
            fontsize=16
        )

    # Add title
    fig.suptitle("Coverage of Traits Across Sources",fontsize=18)
    # Adjust spacing manually if needed (optional)
    fig.subplots_adjust(top=0.85, bottom=0.02, left=0.04, right=0.98, wspace=0.1, hspace=0.1)
    fig.savefig(directory + "/All_Competencies_Venn_Diagrams.png")
    plt.close()

# def get_taxonomic_rank(df, column):

#     rank_df = pd.read_csv("./ncbitaxon_rank.tsv")
#     df[column] = df[column].map(rank_df.set_index('Rank')[column])

def round_percentages(x):
    if isinstance(x, (int, float)):
        if x >= 1:  # If the value is 1 or more, round to the nearest whole number
            return int(round(x))
        else:  # If the value is less than 1, round to two decimal places
            return round(x, 2)
    return x

def create_metabolite_competency_df(metabolite, direction, reaction_direction_dict):

    directory = "./Intermediate_Files_Competencies/" + metabolite + "_" + direction

    new_row = {}
    # competency_df = pd.DataFrame(columns = ["Metabolite", "Traits_Annotations", "Rhea-Chebi_Annotations",  "Traits_Annotations_Proteome_Overlap", "Total Rhea-Chebi_Traits_Overlap", "Rhea-Chebi_Traits_Overlap_Percentage", "Traits_Rhea-Chebi_Overlap_Percentage"])

    new_row["Metabolite"] = metabolite
    # competency_df.loc[0,"Metabolite"] = metabolite

    organismal_strains = pd.read_csv(directory + "/" + ORGANISMAL_TRAITS_STRAINS_ANNOTATIONS_FILE + ".tsv",delimiter="\t")
    organismal_strains_list = organismal_strains["updated_subject"].unique().tolist()
    total_traits = len(organismal_strains.drop_duplicates(subset=["updated_subject"]))

    new_row["Total_Traits_Annotations"] = total_traits
    # competency_df.loc[0,"Traits_Annotations"] = len(organismal_strains.drop_duplicates(subset=["updated_subject"]))

    organismal_strains_proteomes = pd.read_csv(directory + "/" + ORGANISMAL_TRAITS_STRAINS_PROTEOMES_ANNOTATIONS_FILE + ".tsv",delimiter="\t")
    total_traits_proteomes_overlap = len(organismal_strains_proteomes)
    new_row["Total_Traits_Proteome_Overlap"] = total_traits_proteomes_overlap
    # competency_df.loc[0,"Traits_Annotations_Proteome_Overlap"] = len(organismal_strains_proteomes)

    rhea_chebi_strains = pd.read_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", delimiter="\t")
    if direction == "produces":
        rhea_chebi_strains = rhea_chebi_strains[~((rhea_chebi_strains['rhea'].isin(reaction_direction_dict.keys())) & (rhea_chebi_strains['rhea'].map(reaction_direction_dict) == "input"))]
        # Overwrite original rhea data with new direction information
        rhea_chebi_strains.to_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", sep="\t")
    if direction == "consumes":
        rhea_chebi_strains = rhea_chebi_strains[~((rhea_chebi_strains['rhea'].isin(reaction_direction_dict.keys())) & (rhea_chebi_strains['rhea'].map(reaction_direction_dict) == "output"))]
        # Overwrite original rhea data with new direction information
        rhea_chebi_strains.to_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", sep="\t")
    total_rhea_chebi = len(rhea_chebi_strains.drop_duplicates(subset=["ncbitaxon"]))
    rhea_chebi_annotations_list = rhea_chebi_strains["ncbitaxon"].unique().tolist()
    new_row["Total_Rhea-Chebi_Annotations"] = total_rhea_chebi
    # competency_df.loc[0,"Rhea-Chebi_Annotations"] = len(rhea_chebi_strains.drop_duplicates(subset=["ncbitaxon"]))

    total_traits_rhea_chebi_overlap = len(list(set(rhea_chebi_annotations_list) & set(organismal_strains_list)))
    new_row["Total_Rhea-Chebi_Traits_Overlap"] = total_traits_rhea_chebi_overlap
    # competency_df["Total Rhea-Chebi_Traits_Overlap"] = organismal_rhea_chebi_strains_overlap

    print(metabolite, direction)
    if metabolite == "butyrate" and direction == "produces":
        print("adding ECs")
        ec_strains_list = pd.read_csv(directory + "/" + EC_ANNOTATIONS_FILE_SUBSTRING + "all.tsv", delimiter="\t").drop_duplicates(subset=["subject"])["subject"].tolist()
    else:
        ec_strains_list = []
    total_ec = len(ec_strains_list)
    new_row["EC_Annotations"] = total_ec

    total_ec_traits_overlap = len(list(set(ec_strains_list) & set(organismal_strains_list)))
    new_row["Total_EC_Traits_Overlap"] = total_ec_traits_overlap
    # competency_df.loc[0,"EC_Annotations"] = len(ec_strains_list)
    # competency_df.loc[0,"Total EC_Traits_Overlap"] = len(list(set(ec_strains_list) & set(organismal_annotations)))
    total_ec_rhea_chebi_overlap = len(list(set(ec_strains_list) & set(rhea_chebi_annotations_list)))
    new_row["Total_EC_Rhea-Chebi_Overlap"] = total_ec_rhea_chebi_overlap

    total_ec_and_rhea_chebi = len(list(set(ec_strains_list) | set(rhea_chebi_annotations_list)))
    new_row["Total_EC_and_Rhea-Chebi"] = total_ec_and_rhea_chebi

    # competency_df.loc[0,"Total EC_Rhea-Chebi_Overlap"] = len(list(set(ec_strains_list) & set(rhea_chebi_annotations)))

    total_ec_traits_rhea_chebi_overlap = len(list((set(ec_strains_list) & set(rhea_chebi_annotations_list)) & set(organismal_strains_list)))
    new_row["Total_Rhea-Chebi_and_EC_Traits_Overlap"] = total_ec_traits_rhea_chebi_overlap

    # competency_df.loc[0,"Total Rhea-Chebi_and_EC_Traits_Overlap"] = len(list((set(ec_strains_list) | set(rhea_chebi_annotations)) & set(organismal_annotations)))

    total = len(list(set(ec_strains_list) | set(rhea_chebi_annotations_list) | set(organismal_strains_list)))
    new_row["Total_Traits_and_Rhea-Chebi_and_EC_Annotations"] = total
    # competency_df.loc[0,"Traits_and_Rhea-Chebi_and_EC_Annotations"] = len(list(set(ec_strains_list) | set(rhea_chebi_annotations) | set(organismal_annotations)))

    # Metrics for paper
    summary_row = {
        "Metabolite" : metabolite,
        "traits_only" : (total_traits - total_ec_traits_overlap - total_traits_rhea_chebi_overlap) + total_ec_traits_rhea_chebi_overlap,
        "ec_only" : (total_ec - total_ec_rhea_chebi_overlap - total_ec_traits_overlap) + total_ec_traits_rhea_chebi_overlap,
        "rhea_chebi_only" : (total_rhea_chebi - total_traits_rhea_chebi_overlap - total_ec_rhea_chebi_overlap) + total_ec_traits_rhea_chebi_overlap,
        "traits_rhea_chebi_only" : total_traits_rhea_chebi_overlap - total_ec_traits_rhea_chebi_overlap,
        "traits_ec_only" : total_ec_traits_overlap - total_ec_traits_rhea_chebi_overlap,
        "ec_rhea_chebi_only" : total_ec_rhea_chebi_overlap - total_ec_traits_rhea_chebi_overlap,
        "traits_rhea_chebi_ec_only" : total_ec_traits_rhea_chebi_overlap,
    }

    competency_df = pd.DataFrame()  # Initialize if not already existing
    competency_df = pd.concat([competency_df, pd.DataFrame([new_row])], ignore_index=True)
    competency_df = competency_df.applymap(round_percentages)

    combined_competency_kg_summary_df = pd.DataFrame()
    combined_competency_kg_summary_df = pd.concat([combined_competency_kg_summary_df, pd.DataFrame([summary_row])], ignore_index=True)
    return competency_df, combined_competency_kg_summary_df

def combine_all_competency_dfs(competency_dfs, direction, filename):

    directory = "./Intermediate_Files_Competencies"

    cols = competency_dfs[0].columns
    final_df = pd.DataFrame(columns = cols)

    for df in competency_dfs:
        final_df = pd.concat([final_df, df], ignore_index=True)

    final_df.to_csv(directory + "/" + ALL_COMPETENCIES_DF_FILE + "_" + direction + "_" + filename + ".tsv", index=False)

    return final_df

def visualize_all_competencies(final_df, direction):

    directory = "./Intermediate_Files_Competencies"
    metabolite = final_df.iloc[0].loc["Metabolite"]

    values = final_df.iloc[0].tolist()
    values.remove(metabolite)
    labels = list(final_df.columns)
    labels.remove("Metabolite")

    plt.figure(figsize=(10,7))
    wedges, texts, autotexts =  plt.pie(
        values, 
        labels = None, 
        autopct=lambda pct: f"{int(round(pct * sum(values) / 100))}" if pct > 0 else "", 
        startangle=90, 
        wedgeprops={'edgecolor': 'black'}
        )
    plt.legend(wedges, labels, title = "Source", loc="center left", bbox_to_anchor=(1, 0.5))
    plt.title("Representation in KG for: " + metabolite.capitalize() + direction.capitalize())

    plt.savefig(directory + '/All_Sources_Competencies_Comparison_' + metabolite + '_' + direction + '.png')
    plt.close()

def process_metabolite_competency_questions(metabolite):

    directory = "./Intermediate_Files_Competencies/" + metabolite

    # Open a text file in write mode
    with open(directory + '/' + metabolite + '_output.txt', 'w') as f:
        # Redirect standard output to the file
        sys.stdout = f

        # Read in species only organismal traits 
        print("Organismal Traits")
        print("Looking at NCBITaxon - CHEBI strains for only proteomes")
        # organismal_species = pd.read_csv(directory + "/NCBI_organismal_traits.tsv",delimiter="\t")
        organismal_species = pd.read_csv(directory + "/NCBI_organismal_traits_strains_proteomes.tsv",delimiter="\t")

        organismal_species_produces = organismal_species[organismal_species["predicate"] == "biolink:produces"]
        organismal_species_consumes = organismal_species[organismal_species["predicate"] == "biolink:consumes"]

        print('length total NCBITaxon IDs with organismal trait:', len(organismal_species))
        print('length total then wNCBITaxon IDs with produces organismal trait:', len(organismal_species_produces))
        print('length total NCBITaxon IDs with consumes organismal trait:', len(organismal_species_consumes))

        # Read in strains only organismal traits
        organismal_strains = pd.read_csv(directory + "/NCBI_organismal_traits_strains.tsv",delimiter="\t")

        duplicate_rows = organismal_strains[organismal_strains["updated_subject"].duplicated(keep=False)]
        duplicate_rows = duplicate_rows.sort_values(by=["updated_subject"])

        organismal_strains = organismal_strains.drop_duplicates(subset=["updated_subject"])

        print('length total strains with organismal trait:', len(organismal_strains))

        # Read in genomic traits through RHEA
        print("Looking at NCBITaxon - related to - RHEA - CHEBI")
        genomic_rhea = pd.read_csv(directory + "/"+ RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv",delimiter="\t")

        print('length total unique NCBITaxon IDs with genomic RHEA-CHEBI traits:', len(genomic_rhea.drop_duplicates(subset=["ncbitaxon"])))

        # Read in overlap of genomic RHEA and organismal species
        overlap_rhea_organismal_species = pd.read_csv(directory + "/NCBI_organismal_genomic_rhea_comparison_species.tsv",delimiter="\t")

        print('length total unique NCBITaxon IDs with both genomic RHEA-CHEBI traits and organismal trait:', len(overlap_rhea_organismal_species.drop_duplicates(subset=["updated_subject"])))

        if len(organismal_species) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_organismal_species.drop_duplicates(subset=["updated_subject"])) / len(organismal_species)  * 100
        print("Percent of NCBITaxon IDs with organismal trait that also have genomic RHEA-CHEBI traits: ",percentage_overlap)

        if len(genomic_rhea) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_organismal_species.drop_duplicates(subset=["updated_subject"])) / len(genomic_rhea.drop_duplicates(subset=["ncbitaxon"]))  * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-CHEBIs that also have organismal trait: ",percentage_overlap)

        # Read in genomic traits through GO
        print("Genomic Traits")
        print("Looking at NCBITaxon - associated with - GO")
        genomic_go = pd.read_csv(directory + "/NCBI_genomic_traits_GO.tsv",delimiter="\t")

        print('length total unique NCBITaxon IDs with genomic GO traits:', len(genomic_go.drop_duplicates(subset=["subject"])))

        # Read in overlap of genomic GO and organsmal species
        overlap_go_organismal_species = pd.read_csv(directory + "/NCBI_organismal_genomic_go_comparison_species.tsv",delimiter="\t")

        print('length total unique NCBITaxon IDs with both genomic GO traits and organismal trait:', len(overlap_go_organismal_species.drop_duplicates(subset=["subject_id"])))

        if len(organismal_species) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_go_organismal_species.drop_duplicates(subset=["subject_id"])) / len(organismal_species)  * 100
        print("Percent of NCBITaxon IDs with organismal trait that also have genomic GO traits: ",percentage_overlap)

        percentage_overlap = len(overlap_go_organismal_species.drop_duplicates(subset=["subject_id"])) / len(genomic_go.drop_duplicates(subset=["subject"]))  * 100
        print("Percent of NCBITaxon IDs with genomic GOs that also have organismal trait: ",percentage_overlap)

        # Look at how many NCBITaxon ID's have both RHEA-CHEBI and GO traits
        overlap_rhea_go = set(genomic_rhea['ncbitaxon']).intersection(set(genomic_go['subject']))
        print("Number of NCBITaxon ID's that have RHEA-CHEBI traits and GO traits: ",len(set(overlap_rhea_go)))

        # Read in overlap of genomic RHEA and genomic go
        if len(genomic_rhea) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go) / len(genomic_rhea.drop_duplicates(subset=["ncbitaxon"])) * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-CHEBI that also have genomic GO traits: ",percentage_overlap)
        if len(genomic_go) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go) / len(genomic_go.drop_duplicates(subset=["subject"])) * 100
        print("Percent of NCBITaxon IDs with genomic GO that also have genomic RHEA-CHEBI traits: ",percentage_overlap)

        # Read in genomic traits through CHEBI
        print("Looking at NCBITaxon - binds - CHEBI")
        genomic_chebi = pd.read_csv(directory + "/NCBI_genomic_traits_CHEBI.tsv",delimiter="\t")

        print('length total unique NCBITaxon IDs with genomic CHEBI traits:', len(genomic_chebi.drop_duplicates(subset=["object"])))

        # Read in overlap of genomic CHEBI and organismal species
        overlap_chebi_organismal_species = pd.read_csv(directory + "/NCBI_organismal_genomic_chebi_comparison_species.tsv",delimiter="\t")

        print('length total unique genomic chebi, organismal overlap species level:', len(overlap_chebi_organismal_species.drop_duplicates(subset=["updated_subject"])))

        if len(organismal_species) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(set(organismal_species['updated_subject']).intersection(set(genomic_chebi['object']))) / len(set(organismal_species['updated_subject'])) * 100
        print("Percent of NCBITaxon IDs with organismal trait that also have genomic CHEBI traits: ",percentage_overlap)

        # Handling cases where there is no chebi overlap, where denominator is 0
        try:
            percentage_overlap = len(set(organismal_species['updated_subject']).intersection(set(genomic_chebi['object']))) / len(set(genomic_chebi['object'])) * 100
        except ZeroDivisionError:
            percentage_overlap = 0
        print("Percent of NCBITaxon IDs with genomic CHEBI traits that also have organismal trait: ",percentage_overlap)

        # Look at how many NCBITaxon ID's have both RHEA-CHEBI and CHEBI traits
        overlap_rhea_chebi = set(genomic_rhea['ncbitaxon']).intersection(set(genomic_chebi['object']))
        print("Number of NCBITaxon ID's that have RHEA-CHEBI traits and CHEBI traits: ",len(set(overlap_rhea_chebi)))

        # Read in overlap of genomic CHEBI and genomic RHEA
        if len(genomic_rhea) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_chebi) / len(set(genomic_rhea['ncbitaxon'])) * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-CHEBI that also have genomic CHEBI traits: ",percentage_overlap)

        # Handling cases where there is no chebi overlap, where denominator is 0
        try:
            percentage_overlap = len(overlap_rhea_chebi) / len(set(genomic_chebi['object'])) * 100
        except ZeroDivisionError:
            percentage_overlap = 0
        print("Percent of NCBITaxon IDs with genomic CHEBI that also have genomic RHEA-CHEBI traits: ",percentage_overlap)

        # Look at how many NCBITaxon ID's have both GO and CHEBI traits
        overlap_go_chebi = set(genomic_go['subject']).intersection(set(genomic_chebi['object']))
        print("Number of NCBITaxon ID's that have GO traits and CHEBI traits: ",len(set(overlap_go_chebi)))

        # Read in overlap of genomic CHEBI and genomic GO
        if len(genomic_go) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_go_chebi) / len(set(genomic_go['subject'])) * 100
        print("Percent of NCBITaxon IDs with genomic GO that also have genomic CHEBI traits: ",percentage_overlap)

        # Handling cases where there is no chebi overlap, where denominator is 0
        try:
            percentage_overlap = len(overlap_go_chebi) / len(set(genomic_chebi['object'])) * 100
        except ZeroDivisionError:
            percentage_overlap = 0
        print("Percent of NCBITaxon IDs with genomic CHEBI that also have genomic GO traits: ",percentage_overlap)

        # Look at genomic Rhea-GO traits
        print("Looking at NCBITaxon - RHEA - GO")
        genomic_rhea_go = pd.read_csv(directory + "/NCBI_genomic_traits_RHEA_GO.tsv",delimiter="\t")

        print('length total unique NCBITaxon IDs with genomic RHEA-GO traits:', len(genomic_rhea_go.drop_duplicates(subset=["ncbitaxon"])))

        # Read in overlap of genomic GO and organsmal species
        overlap_rhea_go_organismal_species = pd.read_csv(directory + "/NCBI_organismal_genomic_rhea_go_comparison_species.tsv",delimiter="\t")

        print('length total unique NCBITaxon IDs with both genomic GO traits and organismal trait:', len(overlap_rhea_go_organismal_species.drop_duplicates(subset=["updated_subject"])))

        if len(organismal_species) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go_organismal_species.drop_duplicates(subset=["updated_subject"])) / len(organismal_species)  * 100
        print("Percent of NCBITaxon IDs with organismal trait that also have genomic RHEA-GO traits: ",percentage_overlap)

        if len(genomic_rhea_go) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go_organismal_species.drop_duplicates(subset=["updated_subject"])) / len(genomic_rhea_go.drop_duplicates(subset=["ncbitaxon"]))  * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-GO traits that also have organismal trait: ",percentage_overlap)

        # Look at how many NCBITaxon ID's have both RHEA-CHEBI and RHEA-GO traits
        overlap_rhea_go_chebi = set(genomic_rhea_go['ncbitaxon']).intersection(set(genomic_rhea['ncbitaxon']))
        print("Number of NCBITaxon ID's that have RHEA-CHEBI traits and RHEA-GO traits: ",len(set(overlap_rhea_go_chebi)))

        # Read in overlap of genomic RHEA and genomic go
        if len(genomic_rhea) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go_chebi) / len(genomic_rhea.drop_duplicates(subset=["ncbitaxon"])) * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-CHEBI that also have genomic RHEA-GO traits: ",percentage_overlap)
        if len(genomic_rhea_go) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go_chebi) / len(genomic_rhea_go.drop_duplicates(subset=["ncbitaxon"])) * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-GO that also have genomic RHEA-CHEBI traits: ",percentage_overlap)

        # Look at how many NCBITaxon ID's have both RHEA-GO and GO traits
        overlap_rhea_go_go = set(genomic_rhea_go['ncbitaxon']).intersection(set(genomic_go['subject']))
        print("Number of NCBITaxon ID's that have RHEA-GO traits and GO traits: ",len(set(overlap_rhea_go_go)))

        # Read in overlap of genomic RHEA and genomic go
        if len(genomic_rhea_go) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go_go) / len(genomic_rhea_go.drop_duplicates(subset=["ncbitaxon"])) * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-GO that also have genomic GO traits: ",percentage_overlap)
        if len(genomic_go) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go_go) / len(genomic_go.drop_duplicates(subset=["subject"])) * 100
        print("Percent of NCBITaxon IDs with genomic GO that also have genomic RHEA-GO traits: ",percentage_overlap)

        # Read in overlap of genomic CHEBI and RHEA-GO
        # Look at how many NCBITaxon ID's have both RHEA-GO and CHEBI traits
        overlap_rhea_go_binds_chebi = set(genomic_rhea['ncbitaxon']).intersection(set(genomic_chebi['object']))
        print("Number of NCBITaxon ID's that have RHEA-GO traits and CHEBI traits: ",len(set(overlap_rhea_go_binds_chebi)))
        # Handling cases where there is no chebi overlap, where denominator is 0
        try:
            percentage_overlap = len(overlap_rhea_go_binds_chebi) / len(set(genomic_chebi['object'])) * 100
        except ZeroDivisionError:
            percentage_overlap = 0
        print("Percent of NCBITaxon IDs with genomic CHEBI that also have genomic RHEA-GO traits: ",percentage_overlap)

        if len(genomic_rhea_go) == 0:
            percentage_overlap = 0  # Handle divide by zero case
        else:
            percentage_overlap = len(overlap_rhea_go_binds_chebi) / len(genomic_rhea_go.drop_duplicates(subset=["ncbitaxon"])) * 100
        print("Percent of NCBITaxon IDs with genomic RHEA-GO that also have genomic CHEBI traits: ",percentage_overlap)


        # Reset standard output back to the console
        sys.stdout = sys.__stdout__

def get_term_labels(conn, operator, term, prefix, output_dir):

    if operator == "LIKE":
        term = "%" + term + "%"
    query = (
        f"""
        CREATE OR REPLACE TABLE term_labels AS
        SELECT * FROM nodes
        WHERE split_part(id, ':', 1) = '{prefix}'
        AND name {operator} '{term}';
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "term_labels", output_dir + "/" + term + "_" + prefix + "_term_labels.tsv")

    terms_list = pd.read_csv(output_dir + "/" + term + "_" + prefix + "_term_labels.tsv", sep = "\t")["id"].tolist()

    return terms_list

# Function to create the dynamic query for each pathway list
def create_query(pathway_list, pathway_index):
    query_parts = []
    for sublist in pathway_list:
        # Create the OR condition for each sublist (terms in the sublist should be ORed)
        or_conditions = " OR ".join([f"e.object = 'EC:{term}'" for term in sublist])
        # If the sublist contains more than one term, wrap it in parentheses
        query_parts.append(f"({or_conditions})")
    
    # Combine all sublist conditions with AND
    combined_conditions = " AND ".join(query_parts)

    
    # Return the full query for this pathway
    return f"SELECT DISTINCT u.ncbi_taxon AS subject FROM edges AS e JOIN uniprot_to_ncbi u ON e.subject = u.uniprotkb WHERE e.subject LIKE '%UniprotKB:%' AND ({combined_conditions})"

# Function to generate the WHERE clause for a pathway list
def generate_exists_conditions(pathway_list):
    conditions = []
    for sublist in pathway_list:
        or_conditions = " OR ".join([f"e.object = 'EC:{ec}'" for ec in sublist])
        conditions.append(f"EXISTS (SELECT 1 FROM edges e  WHERE ({or_conditions} AND e.subject = u.uniprotkb))")#e.subject = u.uniprotkb AND (
    return " AND ".join(conditions)

# Generate the query dynamically
def create_query_for_pathway(pathway_list,idx):
    exists_conditions = generate_exists_conditions(pathway_list)
    query = f"""
    CREATE OR REPLACE TABLE ncbi_butyrate_pathway_{str(idx+1)} AS
    SELECT DISTINCT u.ncbi_taxon AS subject
    FROM edges e
    LEFT JOIN uniprot_to_ncbi u ON e.subject = u.uniprotkb
    WHERE {exists_conditions};
    """
    return query

# Function to generate SQL for a single pathway
def generate_pathway_sql(pathway_name, pathway_list):
    # Create CASE statement for filtering
    case_statements = []
    for idx, group in enumerate(pathway_list, start=1):
        conditions = " OR ".join([f"e.object = 'EC:{ec}'" for ec in group])
        case_statements.append(f"WHEN {conditions} THEN 'group{idx}'")
    
    case_sql = "\n        ".join(case_statements)
    
    # Full SQL for the pathway
    sql = f"""
CREATE OR REPLACE TABLE {pathway_name} AS
WITH filtered_edges AS (
    SELECT
        u.ncbi_taxon,
        e.subject AS uniprotkb,
        CASE
            {case_sql}
        ELSE NULL
        END AS group_name
    FROM edges e
    JOIN uniprot_to_ncbi u ON e.subject = u.uniprotkb
),
grouped_edges AS (
    SELECT
        ncbi_taxon,
        group_name,
        COUNT(DISTINCT uniprotkb) AS unique_uniprot_count
    FROM filtered_edges
    WHERE group_name IS NOT NULL
    GROUP BY ncbi_taxon, group_name
),
final_taxa AS (
    SELECT ncbi_taxon
    FROM grouped_edges
    WHERE unique_uniprot_count > 0
    GROUP BY ncbi_taxon
    HAVING COUNT(DISTINCT group_name) = {len(pathway_list) - 1}
)
SELECT DISTINCT ncbi_taxon AS subject
FROM final_taxa;
"""
    return sql

def genomic_ec_competency(metabolite, direction):

    output_dir = "./Intermediate_Files_Competencies" + "/" + metabolite + "_" + direction

    conn = duckdb.connect(":memory:")

    print("Loading EC, RHEA relevant table.")

    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_competency_specific_ec.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])
    output_dir = "./Intermediate_Files_Competencies" + "/" + metabolite + "_" + direction

    query = (
        f"""
        CREATE OR REPLACE TABLE uniprot_to_ncbi AS
        SELECT subject AS uniprotkb, object AS ncbi_taxon
        FROM edges
        WHERE split_part(subject, ':', 1) = 'UniprotKB'
        AND split_part(object, ':', 1) = 'NCBITaxon'
        AND predicate = 'biolink:derives_from';
        """
    )

    print(query)
    conn.execute(query)

    # Only have EC competencies for produces butyrate
    if metabolite == "butyrate" and direction == "produces":
        acetyl_coa_pathway_list = [['2.3.1.9','1.8.3.5'], ['1.1.1.35', '1.1.1.36', '1.1.1.157'], ['4.2.1.17', '4.2.1.55'], ['1.3.1.86', '1.3.1.44'], ['2.8.3.9']]
        acetyl_coa_pathway_2_list = [['2.3.1.9','1.8.3.5'], ['1.1.1.35', '1.1.1.36', '1.1.1.157'], ['4.2.1.17', '4.2.1.55'], ['1.3.1.86', '1.3.1.44'],['2.7.2.7'],['2.3.1.19']]
        glutarate_pathway_list = [['2.8.3.12'], ['4.2.1.-'], ['1.3.1.86', '1.3.1.44'], ['4.1.1.70']]
        four_aminobutyrate_pathway_list = [['1.3.1.86', '1.3.1.44'], ['1.1.1.61'], ['4.2.1.120'], ['5.3.3.3'],['2.8.3.-']]
        lysine_pathway_list = [['5.4.3.2'],['5.4.3.3'],['2.3.1.247'], ['1.3.1.86', '1.3.1.44'], ['2.8.3.9']]

        for idx, pathway in enumerate([acetyl_coa_pathway_list,
                        acetyl_coa_pathway_2_list,
                        glutarate_pathway_list,
                        four_aminobutyrate_pathway_list,
                        lysine_pathway_list]):

            # Generate the query for acetyl_coa_pathway_list
            #query = create_query_for_pathway(pathway,idx)
            
            #print(query)
            #conn.execute(query)
            
            sql_query = generate_pathway_sql("ncbi_butyrate_pathway_" + str(idx+1), pathway)
            print(sql_query)
            conn.execute(sql_query)
            output_table_to_file(conn, "(SELECT * FROM ncbi_butyrate_pathway_" + str(idx+1) + ")", output_dir + "/" + EC_ANNOTATIONS_FILE_SUBSTRING + str(idx+1) + ".tsv") #ncbi_butyrate_pathway_" + str(idx+1) + "

        all_dfs = pd.DataFrame()
        for idx, pathway in enumerate([acetyl_coa_pathway_list,
            acetyl_coa_pathway_2_list,
            glutarate_pathway_list,
            four_aminobutyrate_pathway_list,
            lysine_pathway_list]):

            df = pd.read_csv(output_dir + "/" + EC_ANNOTATIONS_FILE_SUBSTRING + str(idx+1) + ".tsv", sep='\t')
            df["Pathway"] = "butyrate_pathway_" + str(idx + 1)
            all_dfs = pd.concat([all_dfs,df],axis=0)

        new_ec_filename = EC_ANNOTATIONS_FILE_SUBSTRING + "all"
        all_dfs.to_csv(output_dir + "/" + new_ec_filename + ".tsv",sep='\t')

    conn.close()

def organismal_genomic_competency(metabolite, direction):

    output_dir = "./Intermediate_Files_Competencies" + "/" + metabolite + "_" + direction

    os.makedirs(output_dir, exist_ok=True)

    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_noEC.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])

    # Get corresponding GO and CHEBI mappings
    # GO_terms_list = GO_mappings[metabolite]
    GO_terms_list = get_term_labels(conn, "LIKE", metabolite, "GO", output_dir)
    GO_duckdb_query =  '(' + ' OR '.join([f"e.object = '{term}'" for term in GO_terms_list]) + ')'
    # CHEBI_terms_list = CHEBI_mappings[metabolite]
    CHEBI_exact_terms_list = [] #get_term_labels(conn, "=", metabolite, "CHEBI", output_dir)
    term_labels = METABOLITES_RELEVANT_TERMS[metabolite]
    term_ids = []
    for l in term_labels:
        new_term = get_term_labels(conn, "=", l, "CHEBI", output_dir)
        term_ids.extend(new_term)
        CHEBI_exact_terms_list.extend(new_term)
    df = pd.DataFrame({"id": term_ids, "name": term_labels})
    df.to_csv(output_dir + "/term_mappings.csv")
    if metabolite == "3-(1H-indol-3-yl)propanoic acid":
        GO_terms_list = get_term_labels(conn, "LIKE", "indole", "GO", output_dir)
        GO_duckdb_query =  '(' + ' OR '.join([f"e.object = '{term}'" for term in GO_terms_list]) + ')'
    CHEBI_terms_list = get_term_labels(conn, "LIKE", metabolite, "CHEBI", output_dir)
    CHEBI_object_duckdb_query =  '(' + ' OR '.join([f"e.object = '{term}'" for term in CHEBI_terms_list]) + ')'
    CHEBI_subject_duckdb_query =  '(' + ' OR '.join([f"e.subject = '{term}'" for term in CHEBI_terms_list]) + ')'
    CHEBI_exact_object_duckdb_query =  '(' + ' OR '.join([f"e.object = '{term}'" for term in CHEBI_exact_terms_list]) + ')'
    CHEBI_exact_subject_duckdb_query =  '(' + ' OR '.join([f"e.subject = '{term}'" for term in CHEBI_exact_terms_list]) + ')'
    # CHEBI_subject_duckdb_query = "e.subject = 'CHEBI:17968'"
    # CHEBI_object_duckdb_query = "e.object = 'CHEBI:17968'"

    ## To query organismal traits by species and strain
    # Get all bugs that have organismal trait of produces/consumes metabolite
    ## NCBITaxon - produces/consumes - metabolite
    #! be specific about edge
    ## species_metabolite_edges
    if direction == "produces_consumes":
        direction_query = "AND (predicate LIKE '%produces' OR predicate LIKE '%consumes')"
    elif direction == "produces":
        direction_query = "AND (predicate LIKE '%produces')"
    elif direction == "consumes":
        direction_query = "AND (predicate LIKE '%consumes')"
    query = (
        f"""
        CREATE OR REPLACE TABLE species_metabolite_edges AS
        SELECT DISTINCT * FROM edges e
        WHERE ((subject LIKE 'NCBITaxon%' OR subject LIKE 'strain%') AND subject <> 'NCBITaxon:9606')
        AND {CHEBI_exact_object_duckdb_query}
        {direction_query};
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject, predicate, object FROM species_metabolite_edges)", output_dir + "/" + ORGANISMAL_TRAITS_ANNOTATIONS_FILE + ".tsv")


    # Instead, replace the strains with their NCBITaxon associated species
    ## strain - produces/consumes - metabolite
    #! be specific about edge
    ## strain_to_ncbi
    query = (
        f"""
        CREATE OR REPLACE TABLE strain_to_ncbi AS
        SELECT e.subject AS strain_subject,
            e.object AS ncbi_object
        FROM edges e
        WHERE split_part(e.subject, ':', 1) = 'strain'
        AND e.predicate = 'biolink:subclass_of'
        AND split_part(e.object, ':', 1) = 'NCBITaxon'
        AND object <> 'NCBITaxon:9606';
        """
    )

    print(query)
    conn.execute(query)

    query = (
        f"""
        CREATE OR REPLACE TABLE strain_metabolite_edges AS
        SELECT DISTINCT e.*, 
            split_part(subject, ':', 1) AS subject_prefix,
            COALESCE(s2n.ncbi_object, e.subject) AS updated_subject
        FROM edges e
        LEFT JOIN strain_to_ncbi s2n ON e.subject = s2n.strain_subject
        WHERE (split_part(subject, ':', 1) = 'NCBITaxon' OR split_part(subject, ':', 1) = 'strain')
        AND {CHEBI_exact_object_duckdb_query}
        {direction_query};
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "strain_metabolite_edges", output_dir + "/" + ORGANISMAL_TRAITS_STRAINS_ANNOTATIONS_FILE + ".tsv")

    # Get only NCBITaxon from organismal traits that have proteomes in genomic traits
    query = (
        f"""
        CREATE OR REPLACE TABLE strain_metabolite_edges_proteomes AS
        SELECT DISTINCT s.predicate, s.object, s.updated_subject
        FROM strain_metabolite_edges s
        JOIN edges e
        ON s.updated_subject = e.object
        WHERE e.subject LIKE '%UniprotKB:%'
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "strain_metabolite_edges_proteomes", output_dir + "/" + ORGANISMAL_TRAITS_STRAINS_PROTEOMES_ANNOTATIONS_FILE + ".tsv")

    #### To query genomic traits through GO
    ## UniprotKB - derives_from - Proteome - derives_from - NCBITaxon
    ## uniprot_to_ncbi
    query = (
        f"""
        CREATE OR REPLACE TABLE uniprot_to_ncbi AS
        SELECT subject AS uniprotkb, object AS ncbi_taxon
        FROM edges
        WHERE split_part(subject, ':', 1) = 'UniprotKB'
        AND split_part(object, ':', 1) = 'NCBITaxon'
        AND predicate = 'biolink:derives_from';
        """
    )

    #query = (
    #    f"""
    #    CREATE OR REPLACE TABLE uniprot_to_ncbi AS
    #    SELECT e1.subject AS uniprotkb, e2.object AS ncbi_taxon
    #    FROM edges e1
    #    JOIN edges e2
    #    ON e1.object = e2.subject
    #    WHERE split_part(e1.subject, ':', 1) = 'UniprotKB'
    #    AND e1.predicate = 'biolink:derives_from'
    #    AND split_part(e2.object, ':', 1) = 'NCBITaxon'
    #    AND e2.object <> 'NCBITaxon:9606'
    #    AND e2.predicate = 'biolink:derives_from';
    #    """
    #)

    print(query)
    conn.execute(query)

    result = conn.execute(
        f"""
        SELECT COUNT(*) FROM uniprot_to_ncbi;
        """
    ).fetchone()

    uniprot_to_ncbi = result[0]

    #! Including statement to get rid of proteins without proteome, which is bug in tar build only- NOT NULL
    ## NCBITaxon - participates_in - GO
    ## ncbi_to_go
    query = (
        f"""
        CREATE OR REPLACE TABLE ncbi_to_go AS
        SELECT DISTINCT u.ncbi_taxon AS subject,
                e.predicate,
                e.object
        FROM edges e
        LEFT JOIN uniprot_to_ncbi u ON e.subject = u.uniprotkb
        WHERE split_part(subject, ':', 1) = 'UniprotKB'
        AND {GO_duckdb_query}
        AND u.ncbi_taxon IS NOT NULL;
        """ 
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject, object FROM ncbi_to_go)", output_dir + "/NCBI_genomic_traits_GO.tsv")

    # Get all NCBITaxon IDs in organismal trait species and genomic results
    ## NCBITaxon - metabolite
    ## matching_species_organismal_go
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_species_organismal_go AS
        SELECT sme.updated_subject AS subject_id, sme.object AS metabolite_object, ncbi_to_go.subject
        FROM strain_metabolite_edges_proteomes sme
        JOIN ncbi_to_go ncbi_to_go ON sme.updated_subject = ncbi_to_go.subject;
        """
    )
    # species_metabolite_edges.subject

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject_id, metabolite_object FROM matching_species_organismal_go)", output_dir + "/NCBI_organismal_genomic_go_comparison_species.tsv")


    # Get all NCBITaxon IDs in organismal trait species and genomic results, updated_subject column is associated ncbitaxon
    ## NCBITaxon - metabolite
    ## matching_strain_organismal_go
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_strain_organismal_go AS
        SELECT sme.updated_subject AS subject_id, sme.object AS metabolite_object, ncbi_to_go.subject
        FROM strain_metabolite_edges_proteomes sme
        JOIN ncbi_to_go ncbi_to_go ON sme.updated_subject = ncbi_to_go.subject;
        """
    )
    # species_metabolite_edges.subject

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject_id, metabolite_object FROM matching_strain_organismal_go)", output_dir + "/NCBI_organismal_genomic_go_comparison_strain.tsv")

    #### To query genomic traits through RHEA
    # Create table of RHEA to CHEBI relevant terms
    ## rhea - has_input/output - chebi
    ## rhea_to_chebi
    query = (
        f"""
        CREATE OR REPLACE TABLE rhea_to_chebi AS
        SELECT *
        FROM edges e
        WHERE split_part(e.subject, ':', 1) = 'RHEA'
        AND {CHEBI_exact_object_duckdb_query};
        """
    )

    print(query)
    conn.execute(query)

    result = conn.execute(
        f"""
        SELECT COUNT(*) FROM rhea_to_chebi;
        """
    ).fetchone()

    rhea_to_chebi = result[0]

    # Create table of NCBITaxon nodes that have edges to those Chebi results
    ## ncbitaxon - rhea
    ## ncbitaxon_to_rhea
    #! Including statement to get rid of proteins without proteome, which is bug in this build only- NOT NULL
    query = (
        f"""
        CREATE OR REPLACE TABLE ncbitaxon_to_rhea AS
        SELECT 
        COALESCE(u.ncbi_taxon, e.subject) AS subject, 
        e.predicate, 
        e.object,
        u.ncbi_taxon, 
        u.uniprotkb
        FROM edges e
        LEFT JOIN uniprot_to_ncbi u
        ON e.subject = u.uniprotkb
        WHERE split_part(e.subject, ':', 1) = 'UniprotKB'
        AND split_part(e.object, ':', 1) = 'RHEA'
        AND u.ncbi_taxon IS NOT NULL; 
        """
    )

    print(query)
    conn.execute(query)

    result = conn.execute(
        f"""
        SELECT COUNT(*) FROM ncbitaxon_to_rhea;
        """
    ).fetchone()

    ncbitaxon_to_rhea = result[0]

    # Create table of ncbitaxon to Chebi
    ## ncbitaxon - chebi
    ## ncbitaxon_to_chebi
    query = (
        f"""
        CREATE OR REPLACE TABLE ncbitaxon_to_chebi AS
        SELECT DISTINCT ncbitaxon_to_rhea.subject AS ncbitaxon, rhea_to_chebi.subject AS rhea, rhea_to_chebi.predicate, rhea_to_chebi.object AS chebi, ncbitaxon_to_rhea.uniprotkb AS uniprotkb
        FROM ncbitaxon_to_rhea
        JOIN rhea_to_chebi ON ncbitaxon_to_rhea.object = rhea_to_chebi.subject;
        """
    )

    print(query)
    conn.execute(query)

    result = conn.execute(
        f"""
        SELECT COUNT(*) FROM ncbitaxon_to_chebi;
        """
    ).fetchone()

    ncbitaxon_to_chebi = result[0]

    output_table_to_file(conn, "(SELECT * FROM ncbitaxon_to_chebi)", output_dir + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv")

    # Get all NCBITaxon IDs in organismal trait species and RHEA genomic results
    ## NCBITaxon - metabolite
    ## matching_species_organismal_rhea
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_species_organismal_rhea AS
        SELECT DISTINCT sme.updated_subject, sme.object, ncbitaxon_to_chebi.ncbitaxon
        FROM strain_metabolite_edges_proteomes sme
        JOIN ncbitaxon_to_chebi ON sme.updated_subject = ncbitaxon_to_chebi.ncbitaxon;
        """
    )
    # species_metabolite_edges.subject

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT * FROM matching_species_organismal_rhea)", output_dir + "/" + ORGANISMAL_RHEA_CHEBI_OVERLAP_SPECIES_FILE + ".tsv")

    # Get all NCBITaxon IDs in organismal trait species and RHEA genomic results, updated_subject column is associated ncbitaxon
    ## NCBITaxon - metabolite
    ## matching_strain_organismal_rhea
    # Note the column was subject from sme for the final line
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_strain_organismal_rhea AS
        SELECT sme.updated_subject AS subject_id, sme.object AS metabolite_object, ncbitaxon_to_chebi.ncbitaxon
        FROM strain_metabolite_edges sme
        JOIN ncbitaxon_to_chebi ON sme.updated_subject = ncbitaxon_to_chebi.ncbitaxon;
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject_id, metabolite_object FROM matching_strain_organismal_rhea)", output_dir + "/" + ORGANISMAL_RHEA_CHEBI_OVERLAP_STRAINS_FILE + ".tsv")


    #### To query genomic traits to CHEBI
    # Create table of NCBITaxon to CHEBI relevant terms
    ## ncbitaxon - binds - chebi
    ## ncbitaxon_binds_chebi
    #! Including statement to get rid of proteins without proteome, which is bug in this build only- NOT NULL
    query = (
        f"""
        CREATE OR REPLACE TABLE ncbitaxon_binds_chebi AS
        SELECT DISTINCT
            e.subject AS subject, 
            e.predicate AS predicate, 
            u.ncbi_taxon AS object
        FROM edges e
        LEFT JOIN uniprot_to_ncbi u
        ON e.object = u.uniprotkb
        WHERE {CHEBI_exact_subject_duckdb_query}
        AND e.predicate LIKE '%binds%'
        AND u.ncbi_taxon IS NOT NULL;
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject, object FROM ncbitaxon_binds_chebi)", output_dir + "/NCBI_genomic_traits_CHEBI.tsv")

    # Get all NCBITaxon IDs in organismal trait species and CHEBI genomic results
    ## NCBITaxon - metabolite
    ## matching_species_organismal_chebi
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_species_organismal_chebi AS
        SELECT sme.updated_subject, sme.object, ncbitaxon_binds_chebi.object
        FROM strain_metabolite_edges_proteomes sme
        JOIN ncbitaxon_binds_chebi ON sme.updated_subject = ncbitaxon_binds_chebi.object;
        """
    )
    # species_metabolite_edges.subject

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT * FROM matching_species_organismal_chebi)", output_dir + "/NCBI_organismal_genomic_chebi_comparison_species.tsv")

    # Get all NCBITaxon IDs in organismal trait species and RHEA genomic results, updated_subject column is associated ncbitaxon
    ## NCBITaxon - metabolite
    ## matching_strain_organismal_rhea
    # Note updated subject in final line to updated_subject
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_strain_organismal_chebi AS
        SELECT sme.updated_subject AS subject_id, sme.object AS metabolite_object, ncbitaxon_binds_chebi.object
        FROM strain_metabolite_edges sme
        JOIN ncbitaxon_binds_chebi ON sme.updated_subject = ncbitaxon_binds_chebi.object;
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject_id, metabolite_object FROM matching_strain_organismal_chebi)", output_dir + "/NCBI_organismal_genomic_chebi_comparison_strain.tsv")


    #### To query genomic traits through RHEA to GO
    # Create table of RHEA to GO relevant terms
    ## rhea - enables - go
    ## rhea_to_go
    query = (
        f"""
        CREATE OR REPLACE TABLE rhea_to_go AS
        SELECT *
        FROM edges e
        WHERE split_part(e.subject, ':', 1) = 'RHEA'
        AND {GO_duckdb_query};
        """
    )

    print(query)
    conn.execute(query)

    # Create table of NCBITaxon to GO relevant terms
    ## ncbitaxon - rhea - go
    ## ncbitaxon_rhea_go
    #! Including statement to get rid of proteins without proteome, which is bug in this build only- NOT NULL
    query = (
        f"""
        CREATE OR REPLACE TABLE ncbitaxon_rhea_go AS
        SELECT DISTINCT ncbitaxon_to_rhea.subject AS ncbitaxon, rhea_to_go.predicate,
            rhea_to_go.object AS go
        FROM ncbitaxon_to_rhea
        JOIN rhea_to_go ON ncbitaxon_to_rhea.object = rhea_to_go.subject;
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT ncbitaxon, predicate, go FROM ncbitaxon_rhea_go)", output_dir + "/NCBI_genomic_traits_RHEA_GO.tsv")

    ####

    # Get all NCBITaxon IDs in organismal trait species and RHEA-GO genomic results
    ## matching_species_organismal_go
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_species_organismal_go AS
        SELECT sme.updated_subject, sme.object, ncbitaxon_rhea_go.ncbitaxon
        FROM strain_metabolite_edges_proteomes sme
        JOIN ncbitaxon_rhea_go ON sme.updated_subject = ncbitaxon_rhea_go.ncbitaxon;
        """
    )
    # species_metabolite_edges.subject

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT updated_subject FROM matching_species_organismal_go)", output_dir + "/NCBI_organismal_genomic_rhea_go_comparison_species.tsv")

    # Get all NCBITaxon IDs in organismal trait species and RHEA genomic results, updated_subject column is associated ncbitaxon
    ## NCBITaxon - metabolite
    ## matching_strain_organismal_rhea
    # Note changed subject to updated_subject for sme in final line
    query = (
        f"""
        CREATE OR REPLACE TABLE matching_strain_organismal_rhea AS
        SELECT sme.updated_subject AS subject_id, sme.object AS object, ncbitaxon_rhea_go.ncbitaxon
        FROM strain_metabolite_edges sme
        JOIN ncbitaxon_rhea_go ON sme.updated_subject = ncbitaxon_rhea_go.ncbitaxon;
        """
    )

    print(query)
    conn.execute(query)

    output_table_to_file(conn, "(SELECT subject_id FROM matching_strain_organismal_rhea)", output_dir + "/NCBI_organismal_genomic_rhea_go_comparison_strain.tsv")

    # Get uniprot microbes
    ncbitaxon_func_ids = get_ncbitaxon_with_uniprot(conn, "./Phylogeny_Search")

    return conn
    #genomic_ec_competency(conn, metabolite, direction, output_dir)

def process_congruency_competency_questions():

    directory = "./Intermediate_Files_Competencies/Congruency"

    # Open a text file in write mode
    with open(directory + '/Congruency_output.txt', 'w') as f:
        # Redirect standard output to the file
        sys.stdout = f

        # Read in 


def congruency_competencies():

    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    duckdb_load_table(conn, "./Input_Files/merged-kg/merged-kg_edges.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./Input_Files/merged-kg/merged-kg_nodes.tsv", "nodes", ["id", "name"])

    # Get all Uniprot Protein-RHEA-GO
    
    # Create table of uniprot to Rhea
    ## uniprot - participates_in - rhea
    ## uniprot_to_rhea
    query = (
        f"""
        CREATE OR REPLACE TABLE uniprot_to_rhea AS
        SELECT *
        FROM edges e
        WHERE split_part(e.subject, ':', 1) = 'UniprotKB'
        AND split_part(e.object, ':', 1) = 'RHEA';
        """
    )

    # print(query)
    conn.execute(query)

    query = (
        f"""
        CREATE OR REPLACE TABLE uniprot_rhea_go AS
        SELECT
        COALESCE(u.subject, e.subject) AS subject, 
        e.predicate, 
        e.object
        FROM edges e
        LEFT JOIN uniprot_to_rhea u
        ON e.subject = u.object
        WHERE split_part(e.subject, ':', 1) = 'RHEA'
        AND split_part(e.object, ':', 1) = 'GO'
        AND u.subject IS NOT NULL; 
        """
    )

    # print(query)
    conn.execute(query)

    result = conn.execute(
        f"""
        SELECT COUNT(*) FROM uniprot_rhea_go;
        """
    ).fetchone()

    uniprot_rhea_only = result[0]

    # result = conn.execute("SELECT COUNT(*) FROM uniprot_rhea_go").fetchone()
    # print(f"Number of rows: {result[0]}")

    # Get all Uniprot Protein-GO
    # Create table of Protein to GO
    ## uniprot - participates_in - GO
    ## uniprot_to_go
    query = (
        f"""
        CREATE OR REPLACE TABLE uniprot_to_go AS
        SELECT * FROM edges e
        WHERE split_part(e.subject, ':', 1) = 'UniprotKB'
        AND split_part(e.object, ':', 1) = 'GO';
        """
    )

    # print(query)
    conn.execute(query)

    result = conn.execute(
        f"""
        SELECT COUNT(*) FROM uniprot_to_go;
        """
    ).fetchone()

    uniprot_go_only = result[0]


    result = conn.execute(
        f"""
        SELECT COUNT(*) AS matching_count
        FROM uniprot_rhea_go
        INNER JOIN uniprot_to_go
        ON uniprot_rhea_go.subject = uniprot_to_go.subject
        AND uniprot_rhea_go.object = uniprot_to_go.object;
        """
    ).fetchone()

    rhea_go_overlap = result[0]

    print("Fraction of overlap out of all uniprot go only: ",rhea_go_overlap/uniprot_go_only)
    print("Fraction of overlap out of all uniprot rhea only: ",rhea_go_overlap/uniprot_rhea_only)

    # print(query)
    # conn.execute(query)

    # result = conn.execute("SELECT COUNT(*) FROM uniprot_to_go").fetchone()
    # print(f"Number of rows: {result[0]}")


    # Get all Uniprot Protein-RHEA-CHEBI

    # Get all Uniprot Protein-CHEBI

def get_chemical_direction_upa():

    upa_edges = pd.read_csv('./Input_Files/upa_edges.tsv', sep='\t')

    # Pivot the table
    pivot_df = upa_edges.pivot_table(
        index="object", 
        columns="predicate", 
        aggfunc="size", 
        fill_value=0
    ).reset_index()[["object","biolink:has_input","biolink:has_output"]]

    # Rename columns to match desired format
    pivot_df.columns.name = None  # Remove the pivot column name
    pivot_df.columns = ["chemical", "biolink:has_input", "biolink:has_output"]
    pivot_df = pivot_df.loc[pivot_df["chemical"].str.contains("CHEBI:")]


    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])
    chebi_labels = []
    # Convert CHEBI to labels
    for chebi in pivot_df["chemical"]:
        lab = get_node_label(conn, chebi)
        chebi_labels.append(lab)
    pivot_df["chemical_label"] = chebi_labels


    pivot_df = calc_significance(pivot_df)

    pivot_df.to_csv('./Intermediate_Files_Competencies/upa_chemical_directons.tsv',sep='\t',index=False)

    return pivot_df, conn

def calc_significance(pivot_df):

    # Set alpha level
    alpha = 0.05
    p_values = []

    # Collect p-values for each CHEBI ID
    for idx, row in pivot_df.iterrows():
        count_input = row['biolink:has_input']
        count_output = row['biolink:has_output']
        total = count_input + count_output
        # Binomial test for each CHEBI
        result = binomtest(count_input,n=total, p=0.5)
        p_values.append(result.pvalue)

    # Apply multiple testing correction
    # Bonferroni correction
    _, p_values_bonferroni, _, _ = multipletests(p_values, alpha=alpha, method='bonferroni')

    # Benjamini-Hochberg (FDR) correction
    _, p_values_bh, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')

    # Append results to DataFrame for clarity
    pivot_df['p_value'] = p_values
    pivot_df['p_value_bonferroni'] = p_values_bonferroni
    pivot_df['p_value_bh'] = p_values_bh

    # Output results
    return pivot_df

def get_reactions_with_sig_chemicals(conn,pivot_df, reactions, metabolite, predicted_direction_reactions_df):

    # Filter for significant results after correction
    significant_chebis = pivot_df[pivot_df['p_value_bh'] < 0.05]["chemical"].tolist()
    significant_chebis_duckdb_query =  '(' + ' OR '.join([f"e.object = '{term}'" for term in significant_chebis]) + ')'
    
    if reactions == "all":
        reactions_query = "split_part(e.subject, ':', 1) = 'RHEA'"
    else:
        reactions_query =  '(' + ' OR '.join([f"e.subject = '{term}'" for term in reactions]) + ')'
    #CREATE OR REPLACE TABLE rhea_to_sig_chebi AS
    query = (
        f"""
        SELECT *
        FROM edges e
        WHERE {reactions_query}
        AND predicate = 'biolink:has_participant'
        AND {significant_chebis_duckdb_query};
        """
    )

    # print(query)
    sig_rhea_df = conn.execute(query).fetchdf()

    chebi_labels = []
    # Convert CHEBI to labels
    for chebi in sig_rhea_df["object"]:
        lab = get_node_label(conn, chebi)
        chebi_labels.append(lab)
    sig_rhea_df["object_label"] = chebi_labels

    sig_rhea_df.to_csv('./Intermediate_Files_Competencies/reactions_with_sig_chemicals_' + metabolite + '.tsv', sep='\t')

    # Get total # of undirected reactions in graph for all
    # CREATE OR REPLACE TABLE all_rhea_undirected AS
    if reactions == 'all':
        query = (
        f"""
        SELECT *
        FROM edges e
        WHERE {reactions_query}
        AND predicate = 'biolink:has_participant';
        """
        )
        all_rhea_undirected_df = conn.execute(query).fetchdf()
        reactions = all_rhea_undirected_df['subject'].unique().tolist()

    # Fraction of reactions included in directionality result
    used_reactions = sig_rhea_df['subject'].unique().tolist()
    common_values = list(set(used_reactions).intersection(set(reactions)))
    unique_chemicals = sig_rhea_df['object_label'].unique().tolist()
    used_chemicals = ';'.join(str(value) for value in unique_chemicals)
    columns = ['Metabolite', 'Total_Reactions', 'Overlapping_Direction_Reactions', "Unique_Chemicals"]
    new_rows = pd.DataFrame([[metabolite, len(set(reactions)), len(common_values), used_chemicals]], columns = columns)
    predicted_direction_reactions_df = pd.concat([new_rows, predicted_direction_reactions_df], ignore_index=True)

    return predicted_direction_reactions_df

def get_all_reactions(metabolite):

    directory = "./Intermediate_Files_Competencies/" + metabolite
    df = pd.read_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", sep = "\t")

    rhea_ids = df["rhea"].unique().tolist()

    return rhea_ids

def equilibrator_reaction_direction(conn, metabolite,direction):

    directory = "./Intermediate_Files_Competencies/" + metabolite + "_" + direction

    metabolite_mappings_df = pd.read_csv(directory + "/term_mappings.csv")
    reaction_direction_dict_file = directory + "/reaction_direction_dict.json"
    if not os.path.exists(reaction_direction_dict_file):

        reaction_direction_dict = {}

        # Create a DuckDB connection
        #conn = duckdb.connect(":memory:")

        #print("Loading full nodes table.")

        #duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", "edges", ["subject", "predicate", "object"])
        #duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])

        df = pd.read_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", sep = "\t")
        df = df.drop(columns=['ncbitaxon', "uniprotkb"]).drop_duplicates()

        for i in tqdm.tqdm(range(len(df))):
            rhea = df.iloc[i].loc["rhea"]
            if rhea not in reaction_direction_dict.keys():
                query = (
                    f"""
                    SELECT *
                    FROM edges e
                    WHERE e.object = '{rhea}'
                    AND e.predicate = 'biolink:subclass_of'
                    AND split_part(e.subject, ':', 1) = 'RHEA';
                    """
                )

                result = conn.execute(query).fetchdf()
                # Get first rhea subclass
                directional_rhea = result.iloc[0].loc["subject"]
                query = (
                    f"""
                    SELECT *
                    FROM edges e
                    WHERE e.subject = '{directional_rhea}'
                    AND e.predicate = 'biolink:has_input'
                    AND split_part(e.object, ':', 1) = 'CHEBI';
                    """
                )

                result = conn.execute(query).fetchdf()
                inputs = result["object"].tolist()
                query = (
                    f"""
                    SELECT *
                    FROM edges e
                    WHERE e.subject = '{directional_rhea}'
                    AND e.predicate = 'biolink:has_output'
                    AND split_part(e.object, ':', 1) = 'CHEBI';
                    """
                )

                result = conn.execute(query).fetchdf()
                outputs = result["object"].tolist()
                # Check where metabolite of interest is
                metabs = metabolite_mappings_df["id"].tolist()
                if any(value in inputs for value in metabs):
                    metab_association = "input"
                elif any(value in outputs for value in metabs):
                    metab_association = "output"
                else:
                    metabolite_role = "unknown"
                    reaction_direction_dict[rhea] = metabolite_role
                    continue
                
                equation = f"{' + '.join(inputs)} = {' + '.join(outputs)}"
                try:
                    reaction = cc.parse_reaction_formula(
                        equation
                    )
                except Exception:
                    metabolite_role = "unknown"
                    reaction_direction_dict[rhea] = metabolite_role
                    continue
                gibbs = cc.standard_dg_prime(reaction).to('kilojoule/mole').magnitude
                gibbs = float(gibbs.nominal_value)
                
                metabolite_role = get_reaction_direction(gibbs, metab_association)
                reaction_direction_dict[rhea] = metabolite_role
                print(rhea,metab_association,metabolite_role,gibbs)

        with open(reaction_direction_dict_file, 'w') as json_file:
            json.dump(reaction_direction_dict, json_file, indent=4)

    else:
        with open(reaction_direction_dict_file, 'r') as json_file:
            reaction_direction_dict = json.load(json_file)
            reaction_direction_dict = defaultdict(list, reaction_direction_dict)

    conn.close()

    return reaction_direction_dict

def get_reaction_direction(gibbs, metab_association):

    # Positive gibbs means reaction procedes in reverse direction
    if gibbs > 0 and metab_association == "input":
        return "output"
    elif gibbs > 0 and metab_association == "output":
        return "input"
    # Negative gibbs means reaction procedes in forward direction
    elif gibbs < 0 and metab_association == "input":
        return "input"
    elif gibbs < 0 and metab_association == "output":
        return "output"
    # 0 gibbs means reaction is at equilibrium
    elif gibbs == 0:
        return "unknown"

def convert_ncbitaxon_label(ncbitaxon_df, label):

    try:
        id = ncbitaxon_df.loc[ncbitaxon_df["name"].str.lower() == label.lower(), "id"].values[0]
    except IndexError:
        id = "not found"

    return id

def exact_synonym_match_identification(ncbitaxon_df,label):

    id = "not found"
    ### Check for exact matches in synonym and return only those, don't assume nodes with special characters is a regex
    exact_synonym_matches = ncbitaxon_df[ncbitaxon_df["synonym"].str.contains(label, case=False, na=False)]
    if len(exact_synonym_matches) > 0:
        for i in range(len(exact_synonym_matches)):
            synonym_list = exact_synonym_matches.iloc[i].loc["synonym"].split("|")
            synonym_match = [i for i in synonym_list if i.lower() == label.lower()]
            if len(synonym_match) == 1:
                id = exact_synonym_matches.iloc[i][["id"]].values[0]

    return id

def get_ncbitaxon_df():

    #! TODO: Find a better way to get this path
    ncbitaxon_nodes_file = (
        "./Input_Files/ncbitaxon_nodes.tsv"
    )
    # Get NCBITaxon IDs from ontology nodes file
    if os.path.exists(ncbitaxon_nodes_file):
        ncbitaxon_df = pd.read_csv(ncbitaxon_nodes_file,sep='\t')
    else:
        ncbitaxon_df = pd.DataFrame()

    return ncbitaxon_df

def create_gs_file(metabolite, direction):

    gs_file = GOLD_STANDARD_FILES.get(metabolite + "_" + direction, None)
    if not gs_file:
        return None
    directory = "./Intermediate_Files_Competencies/" + metabolite + "_" + direction
    
    microbe_labels_ids_file = directory + "/gold_standard_ids.tsv"
    microbe_labels_ids_manual_file = directory + "/gold_standard_ids_manual.tsv"
    if not os.path.exists(microbe_labels_ids_file):

        # Get all NCBITaxon IDs
        ncbitaxon_label_dict = {}
        ncbitaxon_synonyms_dict = {}
        ncbitaxon_df = get_ncbitaxon_df()

        if len(ncbitaxon_df) < 0:
            raise RuntimeError("Ensure Input_Files/ncbitaxon_nodes.tsv exists")
        

        all_species_id = []
        all_family_id = []
        gs_df = pd.read_csv(gs_file,sep=",")
        # Convert taxa names to NCBITaxon IDs
        for i in range(len(gs_df)):
            family = gs_df.iloc[i].loc["Family"]
            family_id = convert_ncbitaxon_label(ncbitaxon_df, family)
            if family_id == "not found":
                # Get only exact synonym
                family_id = exact_synonym_match_identification(ncbitaxon_df, family)
            species = gs_df.iloc[i].loc["Name"]
            species_id = convert_ncbitaxon_label(ncbitaxon_df, species)
            if species_id == "not found":
            #     # Try with brackets around genus name
                species_brackets = re.sub(r"^(\w+)", r"[\1]", species)
                species_id = convert_ncbitaxon_label(ncbitaxon_df, species_brackets)
                if species_id == "not found":
                    species_no_sp = species.replace("sp. ","")
                    species_id = convert_ncbitaxon_label(ncbitaxon_df, species_no_sp)
                    if species_id == "not found":
                        species_id = exact_synonym_match_identification(ncbitaxon_df, species)
            all_family_id.append(family_id)
            all_species_id.append(species_id)

        new_df = gs_df.copy(deep=True)
        new_df["Family_ID"] = all_family_id
        new_df["Name_ID"] = all_species_id
        new_df.to_csv(microbe_labels_ids_file,sep = "\t", index=False)
    elif os.path.exists(microbe_labels_ids_manual_file):
        new_df = pd.read_csv(microbe_labels_ids_manual_file, sep = '\t')
    else:
        new_df = pd.read_csv(microbe_labels_ids_file, sep = '\t')

    return new_df

def create_gold_standard_venn_diagram(directory, venn_pairs, title, filename):

    fig, axes = plt.subplots(len(venn_pairs), 1, figsize=(12, 8))

    for i in range(len(venn_pairs)):
        # Create Venn diagram
        venn = venn2(
            subsets=(venn_pairs[i][1],venn_pairs[i][2],venn_pairs[i][3]),
            set_labels=(venn_pairs[i][0], "Gold Standard"),
            ax=axes[i]
        )
    
        # Make labels larger, shift the vertical position of overlap labels for better clarity
        for region_id in ['10', '01']:
            label = venn.get_label_by_id(region_id)
            if label:
                label.set_fontsize(14)
        label = venn.get_label_by_id('11')
        if label: label.set_fontsize(14); label.set_y(0.2)

    # Add title
    fig.suptitle(title,fontsize=18)
    # Adjust spacing manually if needed (optional)
    fig.subplots_adjust(top=0.85, bottom=0.02, left=0.04, right=0.98, wspace=0.1, hspace=0.1)
    fig.savefig(directory + "/" + filename)
    plt.close()

def gold_standard_comparison_species(metabolite, direction):
    
    conn = duckdb.connect(":memory:")
            
    print("Loading ncbitaxon relevant table.")
        
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])

    directory = "./Intermediate_Files_Competencies/" + metabolite + "_" + direction

    gs_df = create_gs_file(metabolite, direction)
    if gs_df is not None:
        gold_standard_list = gs_df["Name_ID"].unique().tolist()
        print("Len lit_comparison taxa")
        print(len(gold_standard_list))
    else:
        gs_df = pd.DataFrame()
    # Organismal
    organismal_strains_list = pd.read_csv(directory + "/" + ORGANISMAL_TRAITS_STRAINS_ANNOTATIONS_FILE + ".tsv",delimiter="\t").drop_duplicates(subset=["updated_subject"])["updated_subject"].tolist()
    # Functional
    rhea_chebi_strains_list = pd.read_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", delimiter="\t").drop_duplicates(subset=["ncbitaxon"])["ncbitaxon"].tolist()
    # Organismal-GS Overlap
    #organismal_gs_overlap = list(set(organismal_strains_list) & set(gold_standard_list))
    #gs_no_organismal_overlap = list(set(gold_standard_list) - set(organismal_gs_overlap))
    # Functional-GS Overlap
    #functional_gs_overlap = list(set(rhea_chebi_strains_list) & set(gold_standard_list))
    #kg_species_list = list(set(organismal_strains_list + rhea_chebi_strains_list))
    #kg_gs_overlap = (list(set(kg_species_list) & set(gold_standard_list)))

    
    # EC-GS Overlap
    if os.path.exists(directory + "/" + EC_ANNOTATIONS_FILE_SUBSTRING + "all.tsv"):
        ec_strains_list = pd.read_csv(directory + "/" + EC_ANNOTATIONS_FILE_SUBSTRING + "all.tsv", delimiter="\t").drop_duplicates(subset=["subject"])["subject"].tolist()
        #ec_gs_overlap = list(set(ec_strains_list) & set(gold_standard_list))
        # Get Pathway Number
        ec_df = pd.read_csv(directory + "/" + EC_ANNOTATIONS_FILE_SUBSTRING + "all.tsv", delimiter="\t")
        ec_dict = ec_df.groupby('subject')['Pathway'].apply(lambda x: [get_node_label(conn, v) for v in x]).to_dict()
    else:
        ec_strains_list = []
        ec_df = pd.DataFrame()
        ec_dict = {}

    #conn = load_graph()
    
    ncbitaxon_func_ids = get_ncbitaxon_with_uniprot(conn, "./Phylogeny_Search")
    print("Len Taxa with a functional annotation")
    print(len(ncbitaxon_func_ids))

    proteome_gs_overlap = (list(set(ncbitaxon_func_ids) & set(gold_standard_list)))
    print("proteome_gs_overlap")
    print(len(proteome_gs_overlap))
    
    print("orig len organismal, functional, ec")
    print(len(organismal_strains_list), len(rhea_chebi_strains_list), len(ec_strains_list))

    # all_values = set(gold_standard_list + rhea_chebi_strains_list + ec_strains_list + organismal_strains_list)

    # Organismal-GS Overlap
    organismal_gs_overlap = list(set(organismal_strains_list) & set(gold_standard_list))
    #gs_no_organismal_overlap = list(set(gold_standard_list) - set(organismal_gs_overlap))
    # Functional-GS Overlap
    functional_gs_overlap = list(set(rhea_chebi_strains_list) & set(gold_standard_list))
    kg_species_list = list(set(organismal_strains_list + rhea_chebi_strains_list))
    kg_gs_overlap = (list(set(kg_species_list) & set(gold_standard_list)))
    ec_gs_overlap = list(set(ec_strains_list) & set(gold_standard_list))

    print("unique len organismal, functional, ec")
    print(len(organismal_strains_list), len(rhea_chebi_strains_list), len(ec_strains_list))

    all_values = set(gold_standard_list + rhea_chebi_strains_list + ec_strains_list + organismal_strains_list)
    sorted_all_values = sorted(all_values)

    # Get protein name
    rhea_chebi_df = pd.read_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", delimiter="\t")

    gold_standard_overlap_file = directory + "/Gold_Standard_Species_Overlap_" + metabolite + "_" + direction + ".csv"
    if not os.path.exists(gold_standard_overlap_file):
        rhea_chebi_protein_dict = rhea_chebi_df.groupby('ncbitaxon')['uniprotkb'].apply(lambda x: [get_node_label(conn, v) for v in x]).to_dict()

        # Create a DataFrame with columns indicating membership
        df = pd.DataFrame({
            "Value": sorted_all_values,
            "Proteome": [1 if val in ncbitaxon_func_ids else 0 for val in sorted_all_values],
            "Gold_Standard": [1 if val in gold_standard_list else 0 for val in sorted_all_values],
            "Organismal": [1 if val in organismal_strains_list else 0 for val in sorted_all_values],
            "Functional": [1 if val in rhea_chebi_strains_list else 0 for val in sorted_all_values],
            "Functional_EC": [1 if val in ec_strains_list else 0 for val in sorted_all_values],
            "Functional_Protein_Name": ["|".join(rhea_chebi_protein_dict[val]) if val in rhea_chebi_strains_list else 0 for val in sorted_all_values],
            "Functional_EC_Pathway" :  ["|".join(ec_dict[val]) if val in ec_strains_list else 0 for val in sorted_all_values],
        })

        df.to_csv(directory + "/Gold_Standard_Species_Overlap_" + metabolite + "_" + direction + ".csv",index=False)
    else:
        df = pd.read_csv(gold_standard_overlap_file,sep = ',')

    venn_pairs = [("Organismal Trait", len(organismal_strains_list), len(gold_standard_list), len(organismal_gs_overlap)),("RHEA Trait", len(rhea_chebi_strains_list), len(gold_standard_list), len(functional_gs_overlap)),("Proteomes vs. RHEA Traits", len(proteome_gs_overlap), len(gold_standard_list), len(functional_gs_overlap))]

    venn_title = "Gold Standard Comparison Across Sources For Species, " + metabolite.capitalize() + ", " + direction.capitalize()
    create_gold_standard_venn_diagram(directory, venn_pairs, venn_title, "Gold_Standard_Species_Venn_Diagrams.png")

    # Show EC Venn
    rhea_chebi_and_ec_overlap = (list(set(rhea_chebi_strains_list) & set(ec_strains_list)))
    rhea_chebi_and_ec_strains_list = list(set(rhea_chebi_strains_list) | set(ec_strains_list))
    functional_and_ec_gs_overlap = (list(set(functional_gs_overlap) | set(ec_gs_overlap)))

    venn_pairs = [("EC Trait", len(ec_strains_list), len(gold_standard_list), len(ec_gs_overlap)),("EC, RHEA Trait", len(rhea_chebi_and_ec_strains_list), len(gold_standard_list), len(functional_and_ec_gs_overlap)),("Proteomes vs. EC, RHEA Traits", len(proteome_gs_overlap), len(gold_standard_list), len(functional_and_ec_gs_overlap))]

    venn_title = "Gold Standard Comparison Across Sources For Species, " + metabolite.capitalize() + ", " + direction.capitalize()
    create_gold_standard_venn_diagram(directory, venn_pairs, venn_title, "Gold_Standard_Species_Venn_Diagrams_EC.png")

    org_rhea_chebi_and_ec_strains_list = list(set(rhea_chebi_strains_list) | set(ec_strains_list) | set(organismal_strains_list))
    org_functional_and_ec_gs_overlap = (list(set(functional_gs_overlap) | set(ec_gs_overlap) | set(organismal_gs_overlap)))

    venn_pairs = [("Organismal, RHEA, EC Traits", len(org_rhea_chebi_and_ec_strains_list), len(gold_standard_list), len(org_functional_and_ec_gs_overlap)),("Proteomes vs. EC, RHEA Traits", len(proteome_gs_overlap), len(gold_standard_list), len(functional_and_ec_gs_overlap))]

    venn_title = "Gold Standard Comparison Across Sources For Species, " + metabolite.capitalize() + ", " + direction.capitalize()
    create_gold_standard_venn_diagram(directory, venn_pairs, venn_title, "Gold_Standard_Species_Venn_Diagrams_EC_RHEA_Org.png")

    # # Compare KG to GS
    # venn_pairs = [("All KG Traits", len(kg_species_list), len(gold_standard_list), len(kg_gs_overlap))]

    # venn_title = "Gold Standard Comparison Across Sources For Families, " + metabolite.capitalize() + ", " + direction.capitalize()
    # create_gold_standard_venn_diagram(directory, venn_pairs, venn_title, "Gold_Standard_KGSpecies_Venn_Diagram.png")

    return conn

def convert_to_species(conn, taxa_list, taxa_list_type):

    # Create a DuckDB connection
    '''
    conn = duckdb.connect(":memory:")
             
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_og.tsv", "edges", ["subject", "predicate", "object"])
    '''
    output_dir = "./Phylogeny_Search"
    os.makedirs(output_dir, exist_ok=True)

    ncbi_taxa_ranks_df = get_all_ranks(output_dir)

    microbes_species_dict = find_microbes_species(conn, ncbi_taxa_ranks_df, taxa_list, output_dir, taxa_list_type)

    # species_list = []
    # for i in taxa_list:
    #     species_list.append([key for key, values in microbes_species_dict.items() if i in values][0])

    return microbes_species_dict

def convert_to_family(conn, taxa_list, taxa_list_type):

    output_dir = "./Phylogeny_Search"
    os.makedirs(output_dir, exist_ok=True)

    ncbi_taxa_ranks_df = get_all_ranks(output_dir)

    # # For all phyla
    #conn = load_graph()
    # ncbitaxa = get_all_kg_taxa(conn)
    # print("Number of taxa in KG")
    # print(len(ncbitaxa))

    # ncbitaxon_func_ids, ncbitaxon_func_dict = get_ncbitaxon_with_functional_annotation(conn, output_dir)
    # print("Len Taxa with a functional annotation")
    # print(len(ncbitaxon_func_ids))

    # #Add functional annotated microbes to ncbitaxa
    # ncbitaxa = list(set(ncbitaxa) | set(ncbitaxon_func_ids))
    # print("Number of taxa in KG + func")
    # print(len(ncbitaxa))

    microbes_family_dict = find_microbes_family(conn, ncbi_taxa_ranks_df, taxa_list, output_dir, taxa_list_type)

    family_list = []
    taxa_without_family = []
    for i in taxa_list:
        # family_list.append(microbes_family_dict[i])
        matching_families = [key for key, values in microbes_family_dict.items() if i in values]
        if matching_families:
            family_list.append(matching_families[0])
        else:
            # Taxa without family mapping - use the taxon itself as placeholder
            family_list.append(i)
            taxa_without_family.append(i)

    if taxa_without_family:
        print(f"Warning: {len(taxa_without_family)} taxa without family mapping (using taxon ID as placeholder)")
        # Optionally log the first few for debugging
        if len(taxa_without_family) <= 5:
            print(f"  Taxa without family: {taxa_without_family}")
        else:
            print(f"  First 5 taxa without family: {taxa_without_family[:5]}")

    return family_list

def gold_standard_comparison_family(conn, metabolite, direction):

    directory = "./Intermediate_Files_Competencies/" + metabolite + "_" + direction

    gs_df = create_gs_file(metabolite, direction)
    if gs_df is not None:
        gold_standard_list = gs_df["Family_ID"].unique().tolist()
        # Organismal
        organismal_strains_list = pd.read_csv(directory + "/" + ORGANISMAL_TRAITS_STRAINS_ANNOTATIONS_FILE + ".tsv",delimiter="\t").drop_duplicates(subset=["updated_subject"])["updated_subject"].tolist()
        # Convert to family
        organismal_families_list = convert_to_family(conn, organismal_strains_list, "organismal_traits_" + metabolite + "_" + direction)
        # Functional
        rhea_chebi_strains_list = pd.read_csv(directory + "/" + RHEA_CHEBI_ANNOTATIONS_FILE + ".tsv", delimiter="\t").drop_duplicates(subset=["ncbitaxon"])["ncbitaxon"].tolist()
        functional_families_list = convert_to_family(conn, rhea_chebi_strains_list, "functional_traits_" + metabolite + "_" + direction)
        # Organismal-GS Overlap
        organismal_gs_overlap = list(set(organismal_families_list) & set(gold_standard_list))
        # Functional-GS Overlap
        functional_gs_overlap = list(set(functional_families_list) & set(gold_standard_list))
        kg_families_list = list(set(organismal_families_list + functional_families_list))
        kg_gs_overlap = (list(set(organismal_families_list + functional_families_list) & set(gold_standard_list)))

        # Combine all unique values from the three lists
        all_values = set(gold_standard_list + organismal_families_list + functional_families_list)
        sorted_all_values = sorted(all_values)

        # Create a DataFrame with columns indicating membership
        df = pd.DataFrame({
            "Value": sorted_all_values,
            "Gold_Standard": [1 if val in gold_standard_list else 0 for val in sorted_all_values],
            "Organismal": [1 if val in organismal_families_list else 0 for val in sorted_all_values],
            "Functional": [1 if val in functional_families_list else 0 for val in sorted_all_values],
        })

        df.to_csv(directory + "/Gold_Standard_Families_Overlap_" + metabolite + "_" + direction + ".csv",index=False)

        venn_pairs = [("Organismal Trait", len(organismal_strains_list), len(gold_standard_list), len(organismal_gs_overlap)),("Functional Trait", len(rhea_chebi_strains_list), len(gold_standard_list), len(functional_gs_overlap)),("All KG Traits", len(kg_families_list), len(gold_standard_list), len(kg_gs_overlap))]

        venn_title = "Gold Standard Comparison Across Sources For Families, " + metabolite.capitalize() + ", " + direction.capitalize()
        create_gold_standard_venn_diagram(directory, venn_pairs, venn_title, "Gold_Standard_Families_Venn_Diagrams.png")
