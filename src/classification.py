import os
import pandas as pd
import numpy as np
from catboost import CatBoostClassifier, Pool, cv

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, balanced_accuracy_score
from sklearn.metrics import confusion_matrix

import matplotlib.pyplot as plt

from classification_utils import add_closure, concatenate_features, create_feature_table, cross_validation, differentiate_edge_direction, get_feature_labels, include_isolation_sources, prepare_data, remove_conflicting_directionality, remove_duplicate_patterns, shap_feature_importance, subset_by_features, train_model

def main():

    # For testing only
    # data_edges = pd.read_csv("/Users/brooksantangelo/Documents/LozuponeLab/FRMS_2024/duckdb/merged-kg_kg-microbe-host_subset/merged-kg_edges.tsv", header=0, sep="\t")
    # data_nodes = pd.read_csv("/Users/brooksantangelo/Documents/LozuponeLab/FRMS_2024/duckdb/merged-kg_kg-microbe-host_subset/merged-kg_nodes.tsv", header=0, sep="\t")
    ###

    data_edges = pd.read_csv("./data/kg-microbe-biomedical-function/merged-kg_edges.tsv", header=0, sep="\t")
    data_nodes = pd.read_csv("./data/kg-microbe-biomedical-function/merged-kg_nodes.tsv", header=0, sep="\t")
    function_edges = pd.read_csv("./data/Intermediate_Files/NCBITaxon_to_GO.tsv", header=0, sep="\t")

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
    print("orig len edges before isolation sources")
    print(len(data_edges))
    data_edges = add_closure(['isolation_source:'], data_edges)
    print("new len edges after isolation sources")
    print(len(data_edges))
    # ['CHEBI:', 'GO:', 'ENVO:', 'UBERON:', 'EC:', 'PO:', 'PATO:', 'FOODON:', 'isolation_source:']

    data = differentiate_edge_direction(data_edges)
    print("data columns")
    print(data.columns)

    data_pairs = data[['subject','object']].drop_duplicates() 

    function_edges.rename(columns={'NCBITaxon': 'subject', 'GO': 'object'}, inplace=True)
    function_edges['object'] = function_edges.apply(
    lambda row: f"{'functional_annotation'}_{row['object']}",
    axis=1 
    )
    
    # Target Class
    data_subset = subset_by_features(data_pairs, "NCBITaxon:", True, "MONDO:0005011", "MONDO:0005265", "MONDO:0005101")

    # Replace Crohns and UC with IBD label
    data_pairs = data_pairs.replace("MONDO:0005011","MONDO:0005101", regex=True)
    data_pairs = data_pairs.replace("MONDO:0005265","MONDO:0005101", regex=True)
    data_subset = data_subset.replace("MONDO:0005011","MONDO:0005101", regex=True)
    data_subset = data_subset.replace("MONDO:0005265","MONDO:0005101", regex=True)
    
    data_pairs_cleaned = remove_conflicting_directionality(data_subset)
    print("data_pairs_cleaned")
    print(len(data_pairs_cleaned))


    for input_dir in ["./data/Intermediate_Files_traits","./data/Intermediate_Files_traits_func", "./data/Intermediate_Files_func"]:

        os.makedirs(input_dir, exist_ok=True)
        if input_dir == "./data/Intermediate_Files_traits" or input_dir == "./data/Intermediate_Files_traits_func":
            # To include traits in feature table
            data_pairs_rest = concatenate_features(data_pairs, data_pairs_cleaned, "NCBITaxon:", "MONDO:0005101", input_dir, function_edges, True)

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

            print(feature_table.columns)
            print(feature_table.shape)

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
            feature_table.to_csv(input_dir + '/feature_table_function.csv', index=True)

        # For testing only
        # feature_table = pd.read_csv(input_dir + "/feature_table.csv",index_col=0); feature_table.index.name = 'NCBITaxon'
        ###

        # Remove rows with only 1?

        feature_table = remove_duplicate_patterns(feature_table, True)

        X, y = prepare_data(feature_table, "MONDO:0005101")

        X, y = get_feature_labels(X, y, data_edges, data_nodes)

        model, X_train, y_train = train_model(X, y, input_dir, True)

        shap_feature_importance(model, X_train, input_dir)

        cross_validation(X_train, y_train, "associated_with_increased_likelihood_of_MONDO:0005101", input_dir)

if __name__ == '__main__':
    main()    

