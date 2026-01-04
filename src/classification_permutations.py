import os
import random
import pandas as pd
import numpy as np
from catboost import CatBoostClassifier, Pool, cv

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, balanced_accuracy_score
from sklearn.metrics import confusion_matrix

import matplotlib.pyplot as plt
from sklearn.utils import compute_class_weight

from classification_utils import add_closure, concatenate_features, create_feature_table, create_model_split, cross_validation, differentiate_edge_direction, get_feature_labels, include_isolation_sources, permute_outcomes, prepare_data, remove_conflicting_directionality, remove_duplicate_patterns, remove_sparse_features, select_top_important_features, shap_feature_dependence_plots, shap_feature_importance, subset_by_features, train_model
from constants import RANDOM_SEED

# Ensures that random numbers chosen for permuations are reproducible
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

def main():

    # For testing only
    # data_edges = pd.read_csv("/Users/brooksantangelo/Documents/LozuponeLab/FRMS_2024/duckdb/merged-kg_kg-microbe-host_subset/merged-kg_edges.tsv", header=0, sep="\t")
    # data_nodes = pd.read_csv("/Users/brooksantangelo/Documents/LozuponeLab/FRMS_2024/duckdb/merged-kg_kg-microbe-host_subset/merged-kg_nodes.tsv", header=0, sep="\t")
    ###

    print("RANDOM_SEED "+str(RANDOM_SEED))

    NUMBER_PERMUTATIONS = 100

    data_edges = pd.read_csv("./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", header=0, sep="\t")
    data_nodes = pd.read_csv("./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", header=0, sep="\t")
    # function_edges = pd.read_csv("./data/Intermediate_Files/NCBITaxon_to_GO.tsv", header=0, sep="\t")

    # Remove NaN rows
    data_edges = data_edges[~data_edges['predicate'].apply(lambda x: isinstance(x, float))]
    data_edges = data_edges[~data_edges['subject'].apply(lambda x: isinstance(x, float))]
    data_edges = data_edges[~data_edges['object'].apply(lambda x: isinstance(x, float))]

    # print("orig len edges before isolation sources")
    # print(len(data_edges))
    # # Add all isolation_source hierarachies
    # data_edges_test = include_isolation_sources(data_edges)
    # print("new len edges after isolation sources")
    # print(len(data_edges_test))

    # Optionally add closure
    # print("orig len edges before isolation sources")
    # print(len(data_edges))
    # data_edges = add_closure(['isolation_source:'], data_edges)
    # print("new len edges after isolation sources")
    # print(len(data_edges))
    # ['CHEBI:', 'GO:', 'ENVO:', 'UBERON:', 'EC:', 'PO:', 'PATO:', 'FOODON:', 'isolation_source:']

    # Remove all isolation source edges
    # print("Len data_edges with isolation source edges")
    # print(len(data_edges))
    # data_edges = data_edges[data_edges['predicate'] != "biolink:location_of"]
    # print("Len data_edges without isolation source edges")
    # print(len(data_edges))

    data = differentiate_edge_direction(data_edges)

    data_pairs = data[['subject','object']].drop_duplicates() 

    # function_edges.rename(columns={'NCBITaxon': 'subject', 'GO': 'object'}, inplace=True)
    # function_edges['object'] = function_edges.apply(
    #     lambda row: f"{'functional_annotation'}_{row['object']}",
    #     axis=1 
    # )

    generic_input_dir = "./data/Intermediate_Files"
    
    # Target Class
    data_subset = subset_by_features(generic_input_dir, data_pairs, "NCBITaxon:", True, "MONDO:0005011", "MONDO:0005265", "MONDO:0005101")

    # # Replace Crohns and UC with IBD label
    # data_pairs = data_pairs.replace("MONDO:0005011","MONDO:0005101", regex=True)
    # data_pairs = data_pairs.replace("MONDO:0005265","MONDO:0005101", regex=True)
    # data_subset = data_subset.replace("MONDO:0005011","MONDO:0005101", regex=True)
    # data_subset = data_subset.replace("MONDO:0005265","MONDO:0005101", regex=True)
    
    import pdb;pdb.set_trace()
    data_pairs_cleaned = remove_conflicting_directionality(data_subset)
    # 239 bugs leftover after this

    data_pairs_cleaned.to_csv(generic_input_dir + "/outcome_to_NCBITaxon_cleaned.tsv", sep="\t", header=True, index=False)


    for input_dir in ["./data/Intermediate_Files_traits"]:#,"./data/Intermediate_Files_traits_func", "./data/Intermediate_Files_func"]:

        os.makedirs(input_dir, exist_ok=True)
        if input_dir == "./data/Intermediate_Files_traits" or input_dir == "./data/Intermediate_Files_traits_func":
            # To include traits in feature table
            data_pairs_rest = concatenate_features(data_pairs, data_pairs_cleaned, "NCBITaxon:", "MONDO:0005101", input_dir, None, True)

            # To include func in feature table
            if input_dir == "./data/Intermediate_Files_traits_func":
                    data_pairs_rest = pd.concat([data_pairs_rest, function_edges], ignore_index=True)

            # print("orig len edges before isolation sources")
            # print(len(data_pairs_rest))
            # print("len of all data")
            # print(len(data))
            # # # Optionally add closure
            # # data_pairs_rest = add_closure(['isolation_source:'], data_pairs_rest, data)
            # # # ['CHEBI:', 'GO:', 'ENVO:', 'UBERON:', 'EC:', 'PO:', 'PATO:', 'FOODON:', 'isolation_source:']
            # print("new len edges after isolation sources")
            # print(len(data_pairs_rest))

            feature_table = create_feature_table(data_pairs_rest, data_pairs_cleaned, "MONDO:0005101", input_dir, True)

        elif input_dir == "./data/Intermediate_Files_func":

            # To include function only in feature table
            function_edges['Value'] = 1
            feature_table_temp = function_edges.pivot_table(index='subject', columns='object', values='Value', aggfunc='sum', fill_value=0)
            feature_table_temp = feature_table_temp.astype(int)
            data_pairs_cleaned.rename(columns={'subject': 'NCBITaxon'}, inplace=True)
            data_pairs_cleaned.rename(columns={'object': 'MONDO:0005101'}, inplace=True)
            data_pairs_cleaned.index = data_pairs_cleaned['NCBITaxon']
            data_pairs_cleaned = data_pairs_cleaned.drop(columns=['NCBITaxon'])
            feature_table = feature_table_temp.merge(data_pairs_cleaned, left_index=True, right_index=True, how='left')
            feature_table = feature_table[feature_table['MONDO:0005101'].notna()]
            # feature_table.to_csv(input_dir + '/feature_table_function.csv', index=True)

        # For testing only
        # feature_table = pd.read_csv(input_dir + "/feature_table.csv",index_col=0); feature_table.index.name = 'NCBITaxon'
        ###

        # Remove rows with only 1?
        feature_table = remove_sparse_features(feature_table)

        feature_table = remove_duplicate_patterns(feature_table, input_dir, True)

        feature_table.to_csv(input_dir + '/feature_table.csv', index=True)

        print("shape feature_table after removing dups and sparse features")
        print(feature_table.shape)

        # Test different threshold
        X, y = prepare_data(feature_table, "MONDO:0005101")

        X, y = get_feature_labels(X, y, data_edges, data_nodes)

        # Train original model with all features
        X_train_global, y_train_global, X_val_global, y_val_global, X_test_global, y_test_global, train_data_global, val_data_global, test_data_global = create_model_split(X, y, "associated_with_increased_likelihood_of_MONDO:0005101", RANDOM_SEED, test = True)

        classes = np.unique(y_train_global)
        class_weights_perm = compute_class_weight(
            class_weight='balanced',
            classes=classes,
            y=y_train_global
        )
        class_weights_perm_dict = dict(zip(classes, class_weights_perm))
        print("Class Weights perm:", class_weights_perm_dict)

        # Rename with dir name
        # Add val data to training because not doing calibration?
        model = train_model(train_data_global, val_data_global, test_data_global, y_test_global, input_dir, RANDOM_SEED, "all_features", True)
        print("y_test_global for model")
        print(len(y_test_global))

        explainer = shap_feature_importance(model, X_train_global, input_dir, "all_features")
        shap_feature_dependence_plots(model, X_train_global, y_train_global, "all_features", input_dir, explainer)

        cross_validation(X_train_global, y_train_global, "associated_with_increased_likelihood_of_MONDO:0005101", input_dir, "all_features")

        # Permutation knockoff
        # Recombine training and val data to get more data
        X_train_val_for_permutation = pd.concat([X_train_global, X_val_global])
        y_train_val_for_permutation = pd.concat([y_train_global, y_val_global])
        ranked_features_df = permute_outcomes(X_train_val_for_permutation, y_train_val_for_permutation, input_dir, NUMBER_PERMUTATIONS, class_weights_perm_dict)

        train_data_selected, val_data_selected, test_data_selected = select_top_important_features(ranked_features_df, X_train_global, X_val_global, X_test_global, y_train_global, y_val_global, y_test_global, input_dir)

        # Rerun model with selected features
        model = train_model(train_data_selected, val_data_selected, test_data_selected, y_test_global, input_dir, RANDOM_SEED, "selected_features", True)

        explainer = shap_feature_importance(model, X_train_val_for_permutation, input_dir, "selected_features")

        shap_feature_dependence_plots(model, X_train_val_for_permutation, y_train_val_for_permutation, "selected_features", input_dir, explainer)

        cross_validation(X_train_val_for_permutation, y_train_val_for_permutation, "associated_with_increased_likelihood_of_MONDO:0005101", input_dir, "selected_features")

if __name__ == '__main__':
    main()    

