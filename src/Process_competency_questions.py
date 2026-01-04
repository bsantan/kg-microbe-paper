


import pandas as pd
from Competencies import combine_all_competency_dfs, congruency_competencies, create_metabolite_competency_df, equilibrator_reaction_direction, genomic_ec_competency, get_all_reactions, get_reactions_with_sig_chemicals, get_rhea_participants, gold_standard_comparison_family, gold_standard_comparison_species, organismal_genomic_competency, plot_competencies_barplot, plot_competencies_venn_diagrams, plot_competencies_venn_diagrams_with_proteomes, process_metabolite_competency_questions, get_chemical_direction_upa, visualize_all_competencies
from constants import ALL_DIRECTIONS, ALL_METABOLITES
from duckdb_utils import duckdb_load_table
import duckdb

def main():

    for direction in ALL_DIRECTIONS:
        all_metabolite_dfs = []
        combined_competency_kg_summary_dfs = []
        for metabolite in ALL_METABOLITES:
            c = organismal_genomic_competency(metabolite,direction)
            # # Create a DuckDB connection
            # conn = duckdb.connect(":memory:")
            
            # print("Loading og relevant table.")
             
            # duckdb_load_table(conn, "./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv", "edges", ["subject", "predicate", "object"])
            # duckdb_load_table(conn, "./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])
            #output_dir = "./data/Intermediate_Files_Competencies" + "/" + metabolite + "_" + direction
            reaction_direction_dict = equilibrator_reaction_direction(c, metabolite,direction)
            genomic_ec_competency(metabolite, direction)
            metabolite_df, combined_competency_kg_summary_df = create_metabolite_competency_df(metabolite, direction, reaction_direction_dict)
            all_metabolite_dfs.append(metabolite_df)
            combined_competency_kg_summary_dfs.append(combined_competency_kg_summary_df)
            conn = gold_standard_comparison_species(metabolite, direction)
            gold_standard_comparison_family(conn, metabolite, direction)
        final_df = combine_all_competency_dfs(all_metabolite_dfs,direction, "raw")
        final_summary_df = combine_all_competency_dfs(combined_competency_kg_summary_dfs,direction, "summary")
        visualize_all_competencies(final_summary_df, direction)
        # plot_competencies_barplot(final_df,direction)
    plot_competencies_venn_diagrams()
    plot_competencies_venn_diagrams_with_proteomes(conn)

    # # To look ats produces only - note need to update organismal_genomic_competency queries
    # butyrate_df = create_metabolite_competency_df("butyrate")
    # propionate_df = create_metabolite_competency_df("propionate")
    # acetate_df = create_metabolite_competency_df("acetate")
    # tryptophan_df = create_metabolite_competency_df("tryptophan")
    # indole_df = create_metabolite_competency_df("indole")


    # # tryptamine_df = create_metabolite_competency_df("tryptamine")
    # # ipa_df = create_metabolite_competency_df("3-(1H-indol-3-yl)propanoic acid")

    

    # get_rhea_participants("produces_consumes_butyrate")
    # # get_rhea_participants("produces_consumes_propionate")
    # # get_rhea_participants("produces_consumes_acetate")
    # # get_rhea_participants("produces_consumes_tryptophan")
    # # get_rhea_participants("produces_consumes_indole")

    # chemicals_sig_df, conn = get_chemical_direction_upa()

    # # To keep track of fraction of reactions with predicted direction
    # predicted_direction_reactions_df = pd.DataFrame()
    # # Get all significant reactions
    # predicted_direction_reactions_df = get_reactions_with_sig_chemicals(conn,chemicals_sig_df, 'all', 'all', predicted_direction_reactions_df)

    # # Get only competency reactions
    # # reactions = []
    # for metabolite in ['produces_consumes_butyrate', 
    #                    'produces_consumes_propionate',
    #                    'produces_consumes_acetate',
    #                    'produces_consumes_tryptophan',
    #                    'produces_consumes_indole']:
    #     reactions = get_all_reactions(metabolite)
    #     predicted_direction_reactions_df = get_reactions_with_sig_chemicals(conn,chemicals_sig_df, reactions, metabolite, predicted_direction_reactions_df)

    # predicted_direction_reactions_df.to_csv('./data/Intermediate_Files_Competencies/overlapping_direction_reactions_per_metabolite.tsv', sep='\t', index=False)


    # process_metabolite_competency_questions("butyrate")

    # organismal_genomic_competency("propionate")

    # process_metabolite_competency_questions("propionate")

    # organismal_genomic_competency("acetate")

    # process_metabolite_competency_questions("acetate")

    # organismal_genomic_competency("L-tryptophan zwitterion")

    # process_metabolite_competency_questions("L-tryptophan zwitterion")

    # congruency_competencies()

if __name__ == '__main__':
    main()
