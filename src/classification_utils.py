import re
import sys
import pandas as pd
import numpy as np
from catboost import CatBoostClassifier, EFstrType, Pool

import shap
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import classification_report, f1_score, precision_score, recall_score

import matplotlib.pyplot as plt
import seaborn

import random
import os
import datetime

from constants import MODEL_PARAMETERS, RANDOM_SEED

current_datetime = datetime.datetime.now()
formatted_date = current_datetime.strftime("%Y-%m-%d_%H_%M_%S")

def differentiate_edge_direction(edges_data):
    """
    Replace subject and object columns with appended edge
    
    :param edges_data: A dataframe with at least the subject, predicate, object columns.
    :type input_dir: DataFrame

    :return: A df with the same columns and new values in the subject, object columns.
    """
    edges_data['predicate_object'] = edges_data.apply(
    lambda row: f"{row['predicate'].replace('biolink:', '')}_{row['object']}",
    axis=1
    )


    edges_data['subject_predicate'] = edges_data.apply(
        lambda row: f"{row['subject']}_{row[ 'predicate'].replace('biolink:', '')}",
        axis=1 
    )

    data_taxa_subject = edges_data[(edges_data["subject"].str.contains("NCBITaxon|strain:"))]
    data_taxa_subject = data_taxa_subject[~data_taxa_subject["object"].str.contains("NCBITaxon|strain:")]
    data_taxa_subject = data_taxa_subject[["subject","predicate_object"]]
    data_taxa_subject = data_taxa_subject.rename(columns = {"predicate_object" : "object"})

    data_taxa_object = edges_data[(edges_data["object"].str.contains("NCBITaxon|strain:"))]
    data_taxa_object = data_taxa_object[~data_taxa_object["subject"].str.contains("NCBITaxon|strain:")]
    data_taxa_object = data_taxa_object[["object","subject_predicate"]]
    data_taxa_object = data_taxa_object.rename(columns = {"object" : "subject", "subject_predicate": "object"})

    data_not_taxa = edges_data[~(edges_data["subject"].str.contains("NCBITaxon|strain:")) & ~(edges_data["object"].str.contains("NCBITaxon|strain:"))]
    data_not_taxa = data_not_taxa[["subject","object"]]

    data_only_taxa = edges_data[(edges_data["subject"].str.contains("NCBITaxon|strain:")) & (edges_data["object"].str.contains("NCBITaxon|strain:"))]
    data_only_taxa = data_only_taxa[["subject","object"]]

    new_data = pd.concat([data_taxa_subject, data_taxa_object], ignore_index=True)
    new_data = pd.concat([new_data, data_not_taxa], ignore_index=True)
    new_data = pd.concat([new_data, data_only_taxa], ignore_index=True)

    # # To append the new relation subject and objects to original df
    # edges_data = edges_data.drop(columns=['predicate_object', 'subject_predicate'])
    # new_data = pd.concat([new_data, edges_data], ignore_index=True).drop_duplicates(ignore_index=True)

    return new_data

def subset_by_features(input_dir, data_pairs, subject_prefix, features, intermediate_files=False):
    """
    Subset edges data to only those that include feature(s) of interest.
    :param input_dir: Directory to which intermediate files will be written.
    :type input_dir: Str
    :param data_pairs: A dataframe with only the subject, object columns.
    :type data_pairs: DataFrame
    :param subject_prefix: Prefix to subset the df by using the subject column (including :).
    :type subject_prefix: Str 
    :param features: One or many features to subset the df by using the object column.
    :type features: string(s)
    """
    if len(features) > 1:
        escaped_features = [re.escape(feature) for feature in features]
        # multiple_features = "|".join(f"({feature})" for feature in escaped_features)
        multiple_features = "|".join(escaped_features)
    elif len(features) == 1:
        multiple_features = features[0]

    # Subset the DataFrame based on the substring in subject
    data_pairs_clean = data_pairs[data_pairs['subject'].str.contains(subject_prefix)]
    # Subset the DataFrame based on the substring in object
    data_pairs_clean = data_pairs_clean[data_pairs_clean['object'].str.contains(multiple_features)] #feature_of_interest+":"
    # Replace Crohns and UC with IBD label
    data_pairs_clean = data_pairs_clean.replace("MONDO:0005011","MONDO:0005101", regex=True)
    data_pairs_clean = data_pairs_clean.replace("MONDO:0005265","MONDO:0005101", regex=True)
    data_pairs_clean = data_pairs_clean.replace("MONDO:0005011","MONDO:0005101", regex=True)
    data_pairs_clean = data_pairs_clean.replace("MONDO:0005265","MONDO:0005101", regex=True)

    data_pairs_clean = data_pairs_clean.drop_duplicates()
    # if intermediate_files:
    #     os.makedirs(input_dir, exist_ok=True)
    #     data_pairs_clean.to_csv(input_dir + "/outcome_to_NCBITaxon.tsv", sep="\t", header=True, index=False)

    return data_pairs_clean

def remove_conflicting_directionality(data_pairs):

    # Group by the 'subject' column and count unique values in the 'object' column for each group
    multiple_links = data_pairs.groupby('subject')['object'].nunique()

    # Filter the groups where the count is greater than 1
    multiple_links = multiple_links[multiple_links > 1]

    # Get the subjects with multiple links
    subjects_with_multiple_links = multiple_links.index.tolist()

    num_unique_values = data_pairs['subject'].nunique()
    r = len(subjects_with_multiple_links)/num_unique_values

    os.makedirs("./data/Intermediate_Files", exist_ok=True)
    with open('./data/Intermediate_Files/Used_data_results.txt', 'w') as file:
        file.write(f'Proportion of microbes both increased and decreased in IBD: {r}\n')

        #Remove bugs that are associated with both increased and decreased Crohns
        file.write(f'Length before removing conflicting pairs: {len(data_pairs)}\n')
        mask = ~data_pairs["subject"].isin(subjects_with_multiple_links)
        data_pairs = data_pairs[mask]
        file.write(f'Length after removing conflicting pairs: {len(data_pairs)}\n')

    return data_pairs

def concatenate_features(data_pairs, data_pairs_cleaned, subject_prefix, feature_of_interest, input_dir, function_edges=None, intermediate_files=False):

    #TODO add closure

    data_pairs_chem = data_pairs[data_pairs['subject'].str.contains(subject_prefix)]
    data_pairs_chem = data_pairs_chem[data_pairs_chem['object'].str.contains('CHEBI:')]

    #TODO add closure
    ###
    ### ESPECIALLY for Taxonomy subClassOf >> one hot
    ###

    data_pairs_go = data_pairs[data_pairs['subject'].str.contains(subject_prefix)]
    data_pairs_go = data_pairs_go[data_pairs_go['object'].str.contains('GO:')]
    #Remove triples with disbiome association
    data_pairs_go = data_pairs_go[~data_pairs_go['object'].str.contains('decreased_likelihood')] #_decreased
    data_pairs_go = data_pairs_go[~data_pairs_go['object'].str.contains('increased_likelihood')] #_increased

    data_pairs_rest_all = data_pairs[data_pairs['subject'].str.contains(subject_prefix)]
    data_pairs_rest = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('carbon_substrates:')]
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('pathways:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('production:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('CAS-RN:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('CHEBI:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('EC:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('GO:')]
    data_pairs_rest2 = data_pairs_rest2[~data_pairs_rest2['object'].str.contains('decreased_likelihood')]
    data_pairs_rest2 = data_pairs_rest2[~data_pairs_rest2['object'].str.contains('increased_likelihood')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)

    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('cell_shape:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('cell_length:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('cell_width:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('gc:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('gram_stain:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('pH_.*:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('temp_.*:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('pigment:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('sporulation:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('trophic_type:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('motility:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('medium:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('BSL:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)

    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('temperature:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('salinity:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('NaCl_.*:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('oxygen:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)


    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('pathogen:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('isolation_source:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('ENVO:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('UBERON:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all[data_pairs_rest_all['object'].str.contains('PO:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)

    data_pairs_rest_all2 = data_pairs[data_pairs['object'].str.contains(subject_prefix)]
    # Swap 'subject' and 'object' for the filtered DataFrame
    data_pairs_rest_all2_swapped = data_pairs_rest_all2.copy()
    data_pairs_rest_all2_swapped['subject'], data_pairs_rest_all2_swapped['object'] = data_pairs_rest_all2['object'], data_pairs_rest_all2['subject']
    data_pairs_rest2 = data_pairs_rest_all2_swapped[data_pairs_rest_all2_swapped['object'].str.contains('UBERON:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all2_swapped[data_pairs_rest_all2_swapped['object'].str.contains('FOODON:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all2_swapped[data_pairs_rest_all2_swapped['object'].str.contains('CHEBI:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all2_swapped[data_pairs_rest_all2_swapped['object'].str.contains('ENVO:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest2[~(data_pairs_rest2['object'].str.contains('related_to') & 
                                     data_pairs_rest2['object'].str.contains('NCBITaxon'))]
    data_pairs_rest2 = data_pairs_rest_all2_swapped[data_pairs_rest_all2_swapped['object'].str.contains('PATO:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest2 = data_pairs_rest_all2_swapped[data_pairs_rest_all2_swapped['object'].str.contains('assay:')]
    data_pairs_rest = pd.concat([data_pairs_rest, data_pairs_rest2], ignore_index=True)
    data_pairs_rest.shape

    data_df_pairs = pd.concat([data_pairs_chem, data_pairs_go], ignore_index=True)
    data_pairs_rest = pd.concat([data_pairs_rest, data_df_pairs], ignore_index=True)

    return data_pairs_rest

def add_closure(subject_prefixes, data):

    # Extract relevant object terms from data_pairs_rest
    relevant_objects = set()
    for obj in data['object']: # data_pairs_rest
        if any(obj.startswith(prefix) for prefix in subject_prefixes):
            relevant_objects.add(obj)

    # Filter data for biolink:subclass_of predicate and relevant objects
    filtered_data = data[(data['predicate'] == 'biolink:subclass_of') & (data['subject'].isin(relevant_objects))]

    # Build a dictionary of child to parent relationships
    subclass_dict = {}
    for _, row in filtered_data.iterrows():
        child = row['subject']
        parent = row['object']
        if child not in subclass_dict:
            subclass_dict[child] = []
        
        #print(f"child {child}\\tparent {parent}")
        subclass_dict[child].append(parent)


    print(f"subclass_dict {len(subclass_dict)}")

    # Function to get all parent terms following subclass relationships
    def get_parents(term, subclass_dict):
        parents = []
        current_term = term
        while current_term in subclass_dict:
            parent_terms = subclass_dict[current_term]
            if not parent_terms:
                break
            # Assume there is only one parent per term for simplicity
            parent = parent_terms[0]
            #print(parent)
            parents.append(parent)
            current_term = parent
        return parents

    # Create a new DataFrame for the closure
    data_pairs_rest_closure = data.copy(deep=True) # data_pairs_rest

    # Extend data_pairs_rest_closure with parent subclass edges
    new_edges = []
    for _, row in data.iterrows(): # data_pairs_rest
        subject = row['subject']
        obj = row['object']
        pred = row['predicate']
        #print(f"obj {obj}")
        if any(obj.startswith(prefix) for prefix in subject_prefixes):
            #print(f"obj2 {obj}")
            if obj in subclass_dict:
                parents = get_parents(obj, subclass_dict)
                #print(parents)
                for parent in parents:#[:-1]:  # Exclude the last parent term
                    new_edges.append({'subject': subject, 'predicate': pred, 'object': parent})
                #if parents:
                #    print(f"Last parent term for {obj}: {parents[-1]}")

    # Convert new_edges to DataFrame
    new_edges_df = pd.DataFrame(new_edges)
    # Concatenate the original and new edges DataFrames
    data_pairs_rest_closure = pd.concat([data_pairs_rest_closure, new_edges_df], ignore_index=True)
    data_pairs_rest_closure = data_pairs_rest_closure.drop_duplicates()

    return data_pairs_rest_closure

def create_feature_table(data_pairs_rest, data_pairs_cleaned, feature_of_interest, input_dir, intermediate_files=False):

    # Drop duplicate subject-object pairs to ensure only unique pairs
    data_pairs_rest = data_pairs_rest.drop_duplicates(subset=['subject', 'object'])

    data_pairs_rest['Value'] = 1

    # Step 2: Pivot the old DataFrame to form the new DataFrame structure
    data_df = data_pairs_rest.pivot_table(index='subject', columns='object', values='Value', aggfunc='sum', fill_value=0)
    # Step 3: Fill NaN values with 0 to indicate no relationship
    #data_df = data_df.fillna(0)

    # Optionally, convert the filled NaN values to integers if they were floats after pivot
    data_df = data_df.astype(int)

    data_pairs_cleaned.rename(columns={'subject': 'NCBITaxon'}, inplace=True)
    data_pairs_cleaned.rename(columns={'object': feature_of_interest}, inplace=True)
    data_pairs_cleaned.index = data_pairs_cleaned['NCBITaxon']
    data_pairs_cleaned = data_pairs_cleaned.drop(columns=['NCBITaxon'])

    data_df = data_df.merge(data_pairs_cleaned, left_index=True, right_index=True, how='left')

    # For testing only
    # data_df[feature_of_interest] = np.random.choice(['associated_with_increased_likelihood_of_MONDO:0005101', 'associated_with_decreased_likelihood_of_MONDO:0005101'], size=len(data_df))
    ###

    # Remove rows where feature of interest is na
    data_df = data_df[data_df[feature_of_interest].notna()]

    print("original data_df len")
    print(len(data_df))
    # Drop duplicate NCBITaxon:feature of interest pairs
    data_df = data_df.reset_index().drop_duplicates(subset=['index', feature_of_interest]).set_index('index')
    print("data_df len after drop dups")
    print(len(data_df))

    # if intermediate_files:
    #     os.makedirs(input_dir, exist_ok=True)
    #     # Write to CSV file
    #     data_df.to_csv(input_dir + '/feature_table.csv', index=True)

    return data_df

def remove_sparse_features(data_df_clean):

    # Select only numeric columns for the operation
    numeric_cols = data_df_clean.select_dtypes(include=['number']).columns

    # Remove columns with sum < 1, excluding any other non-numeric columns
    dimnow = data_df_clean.shape
    sum_less2 = data_df_clean[numeric_cols].sum(axis=0)
    cols_to_drop = data_df_clean[numeric_cols].columns[sum_less2 <= 1]
    data_df_clean = data_df_clean.drop(columns=cols_to_drop)

    #Remove cols which are all 1's
    numeric_cols = data_df_clean.select_dtypes(include=['number']).columns
    dimnow = data_df_clean.shape
    sum_less2 = data_df_clean[numeric_cols].sum(axis=0)
    cols_to_drop = data_df_clean[numeric_cols].columns[sum_less2 == dimnow[1]]
    data_df_clean = data_df_clean.drop(columns=cols_to_drop)


    # remove rows with sum < 1, ensuring to only sum over the updated numeric columns
    numeric_cols = data_df_clean.select_dtypes(include=['number']).columns
    dimnow = data_df_clean.shape
    sum_less2_row = data_df_clean[numeric_cols].sum(axis=1)
    rows_to_drop = data_df_clean.index[sum_less2_row <= 1]
    data_df_clean = data_df_clean.drop(index=rows_to_drop)

    #Remove rows which are all 1's
    numeric_cols = data_df_clean.select_dtypes(include=['number']).columns
    dimnow = data_df_clean.shape
    sum_less2_row = data_df_clean[numeric_cols].sum(axis=1)
    #print(sum_less2_row[sum_less2_row == dimnow[0]])
    rows_to_drop = data_df_clean.index[sum_less2_row == dimnow[0]]
    #print(rows_to_drop)
    data_df_clean = data_df_clean.drop(index=rows_to_drop)

    return data_df_clean

def include_isolation_sources(data_edges):

    location_of = "biolink:location_of"
    subclass_of = "biolink:subclass_of"

    # Keep track of the number of rows before and after each iteration
    previous_size = 0
    current_size = data_edges.shape[0]
    
    # Keep iterating until no more new rows are added
    while previous_size < current_size:
        previous_size = current_size
        
        # Find all location_of relationships
        location_of_rows = data_edges[data_edges['predicate'] == location_of]
        
        # Merge with itself to find related subclass_of relationships
        expanded_rows = location_of_rows.merge(
            data_edges[data_edges['predicate'] == subclass_of][['subject', 'object']],
            left_on='object', right_on='subject',
            suffixes=('', '_related')
        )

        
        # Create new rows from the expanded relationships
        new_rows = expanded_rows[['subject', 'predicate', 'object_related']]
        new_rows.rename(columns={'object_related': 'object'}, inplace=True)
        # Remove duplicate rows (to avoid adding the same row multiple times)
        data_edges = pd.concat([data_edges, new_rows]).drop_duplicates(ignore_index=True)
        # Update the size of the DataFrame after adding new rows
        current_size = data_edges.shape[0]

    return data_edges

def remove_duplicate_patterns(feature_table, input_dir, intermediate_files=False):

    patterns = {} # Same 0 and 1s in a column across taxa
    columns_to_drop = []
    retention_map = {}

    # Iterate over columns
    for col in feature_table.columns:
        pattern = tuple(feature_table[col])
        if pattern not in patterns:
            patterns[pattern] = col
            retention_map[col] = []  # Initialize the list of dropped columns for this pattern
        else:
            # Add the current column to the drop list and map it to the retained column
            columns_to_drop.append(col)
            retention_map[patterns[pattern]].append(col)

    # Drop duplicate columns
    feature_table = feature_table.drop(columns=columns_to_drop)

    # Prepare to write the mapping to a file
    retention_df = pd.DataFrame(
        [(retained, ','.join(duplicates)) for retained, duplicates in retention_map.items() if duplicates],
        columns=['Retained Column', 'Deleted Columns']
    )

    if intermediate_files:
        os.makedirs("data/Intermediate_Files", exist_ok=True)
        # Write to CSV file
        retention_df.to_csv(input_dir + '/feature_table_duplicate_column_mapping.csv', index=False)

    return feature_table

def prepare_data(feature_table, feature_of_interest):

    # Splitting the data into features and target labels
    # X is NCBITaxon x features (binary)
    # y is NCBITaxon x MONDO disease (categorical)

    # Make sure you have enough examples of annotated taxa - not used for disease classification
    # data_df_clean_filtered = feature_table.groupby(feature_of_interest).filter(lambda x : len(x)>10)
    X = feature_table.drop(feature_of_interest, axis=1)#data_pairs_clean[['subject']]
    y = feature_table[feature_of_interest]

    return X, y

def get_feature_labels(X, y, edges, nodes):
# def get_feature_labels(features, edges, nodes):

    unique_predicates = list(set(edges['predicate'].dropna()))
    unique_predicates = [s.replace("biolink:", '') for s in unique_predicates]
    unique_predicates.append("functional_annotation")

    # Sort predicates by length in descending order to ensure longest match first
    sorted_predicates = sorted(unique_predicates, key=len, reverse=True)

    # Remove uniprot to make lookup faster
    filtered_nodes = nodes[~nodes['id'].str.contains('UniprotKB:', na=False)]

    # Convert the input labels column to a categorical type if it isn't
    #X['subject'] = X['subject'].astype('category')

    # Convert categorical columns to integers
    #X['subject'] = X['subject'].cat.codes

    # Iterate over the current column names to get the new column names
    new_column_names = []
    # new_column_names = {}
    print("Getting feature labels...")
    for column in X.columns:
    # for column in features:
        #entity = [i for i in column.split("_") if ":" in i][0]
        predicate = [i for i in sorted_predicates if i in column][0]
        entity = column.replace(predicate,'').strip('_')
        name = filtered_nodes.loc[filtered_nodes["id"].str.lower() == entity.lower(),"name"].values[0]
        # print(entity,name)
        if any(name in s for s in new_column_names):
            new_column_names.append(column.replace(entity,name+"_"+column))
        else:
            new_column_names.append(column.replace(entity,name))
        # new_column_names[column] = column.replace(entity,name)

    # return new_column_names

    # Assign the new column names to the DataFrame
    X.columns = new_column_names
    y.columns = new_column_names

    return X, y

def create_model_split(feature_table, feature_of_interest_table, positive_class, seed, test = True):

    print("len cols all then unique")
    print(len(feature_table.columns))
    print(len(set(feature_table.columns)))

    X_train, X_temp, y_train, y_temp = train_test_split(feature_table, feature_of_interest_table, test_size=0.3, stratify=feature_of_interest_table, random_state=seed)
    if test:
        X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.33, stratify=y_temp, random_state=seed)

        # # Convert labels to binary formats
        # y_train_full_binary = y_train.apply(lambda x: 1 if x == positive_class else 0)
        # y_test_binary = y_test.apply(lambda x: 1 if x == positive_class else 0)

        train_data = Pool(data=X_train, label=y_train)
        val_data = Pool(data=X_val, label=y_val)
        test_data = Pool(data=X_test, label=y_test)

        return  X_train, y_train, X_val, y_val, X_test, y_test, train_data, val_data, test_data

    elif not test:
        train_data = Pool(data=X_train, label=y_train, cat_features=[0])
        val_data = Pool(data=X_temp, label=y_temp, cat_features=[0])

        return  X_train, y_train, X_temp, y_temp, train_data, val_data

def train_model(train_data, val_data, test_data, y_test, input_dir, seed, model_name, intermediate_files=False):

    model = CatBoostClassifier(**MODEL_PARAMETERS)
    
    os.makedirs(input_dir, exist_ok=True)
    # # Write the metrics to a text file
    with open(input_dir + '/' + model_name + '_Model_training_results.txt', 'w') as file:
        file.write(f'X_train shape: {train_data.shape}\n')
        file.write(f'X_val shape: {val_data.shape}\n')
        file.write(f'X_test shape: {test_data.shape}\n')
    
        # Redirect stdout to the file
        original_stdout = sys.stdout
        sys.stdout = file
    
        model.fit(train_data, 
            eval_set=val_data,
            early_stopping_rounds=10
            )#, plot=True)
        
    # Restore the original stdout
    sys.stdout = original_stdout

    # Predict on test data
    y_pred = model.predict(test_data)
    y_pred_proba = model.predict_proba(test_data)[:,1]  # Probabilities for the positive class
    # Generate classification report
    report = classification_report(y_test, y_pred, output_dict=True)
    # Convert the classification report to a DataFrame
    report_df = pd.DataFrame(report).transpose()
    # Filter rows where both precision and recall are greater than 0.9
    # Note: Precision and recall are not defined for the 'accuracy' row, so we exclude it from the filter
    #test_data_report = report_df[(report_df['precision'] > 0.2) & (report_df['recall'] > 0.2) & (report_df.index != 'accuracy')]
    # Print out the filtered rows
    report_df.to_csv(input_dir + '/' + model_name + '_test_data_report.csv', index=True)

    return model

def shap_feature_importance(model, X_train, input_dir, model_name):  #edges, nodes):
    explainer = shap.Explainer(model)

    shap_values = explainer(X_train)

    # # Order features by importance
    # # Compute mean absolute SHAP values for each feature
    # mean_abs_shap_values = np.abs(shap_values.values).mean(axis=0)

    # # Get the feature names in order of importance
    # sorted_feature_names = np.array(shap_values.feature_names)[np.argsort(mean_abs_shap_values)[::-1]]

    # # Get label of top 20 features
    # shap_values_labels = get_feature_labels(sorted_feature_names[0:20], edges, nodes)

    # # Replace values in X_train
    # X_train = X_train.rename(columns=shap_values_labels)

    # Pass the axis to the SHAP summary plot
    # shap.summary_plot(shap_values, X_train, show=False, plot_type="bar", plot_size=(12,5))
    shap.summary_plot(shap_values, X_train, show=False)

    # Save the plot to a file
    plt.savefig(input_dir + "/" + model_name + '_shap_summary_plot.png', bbox_inches='tight')
    plt.close()

    # Save the plot to a file
    plt.savefig(input_dir + '/shap_summary_plot.png', bbox_inches='tight')
    plt.close()

    shap.decision_plot(explainer.expected_value, shap_values.values, X_train.columns, show=False)

    # Save the plot to a file
    plt.savefig(input_dir + "/" + model_name + '_shap_decision_plot.png', bbox_inches='tight')
    plt.close()

    return explainer

def shap_feature_dependence_plots(model, X_train, y_train, model_name, input_dir, explainer):

    # Calculate feature interaction strengths
    interactions = model.get_feature_importance(type='Interaction', data=Pool(X_train, y_train))
    interaction_df = pd.DataFrame(interactions, columns=["Feature1", "Feature2", "InteractionStrength"])
    feature_names = X_train.columns.tolist()
    # Map indices to feature names
    interaction_df["Feature1"] = interaction_df["Feature1"].apply(lambda x: feature_names[int(x)])
    interaction_df["Feature2"] = interaction_df["Feature2"].apply(lambda x: feature_names[int(x)])

    interaction_df.to_csv(input_dir + '/' + model_name + '_feature_interaction.csv', index=True)

    # Feature dependence plots
    shap_interaction_values = explainer.shap_interaction_values(X_train)
    print(f"SHAP interaction values shape: {shap_interaction_values.shape}")

    shap.summary_plot(shap_interaction_values, X_train.iloc[:1000,:], show=False)
    # Save the plot to a file
    plt.savefig(input_dir + "/" + model_name + '_shap_interaction_values.png', bbox_inches='tight')
    plt.close()

    # Plot first n top interactions
    n = 2
    for i in range(2):
        feature1, feature2 = interaction_df['Feature1'].iloc[i], interaction_df['Feature2'].iloc[i]
        shap.dependence_plot(
            (feature1, feature2),
            shap_interaction_values, X_train.iloc[:1000,:],
            display_features=X_train.iloc[:1000,:],
            show=False
        )

        # Save the plot to a file
        plt.savefig(input_dir + "/" + model_name + '_shap_dependence_plot_top_features' + str(i) + '.png', bbox_inches='tight')
        plt.close()

def cross_validation(X_train, y_train, positive_class, input_dir, model_name):

    # Initialize StratifiedKFold
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=RANDOM_SEED)

    precision_scores = []
    recall_scores = []
    f1_scores = []

    # Perform cross-validation manually
    for train_index, val_index in kf.split(X_train, y_train):
        X_train_fold, X_val_fold = X_train.iloc[train_index], X_train.iloc[val_index]
        y_train_fold, y_val_fold = y_train.iloc[train_index], y_train.iloc[val_index]
        
        train_pool = Pool(X_train_fold, y_train_fold)
        val_pool = Pool(X_val_fold, y_val_fold)
        
        model = CatBoostClassifier(**MODEL_PARAMETERS)
        model.fit(train_pool, eval_set=val_pool, early_stopping_rounds=50, verbose=0)
    
        y_val_pred = model.predict(X_val_fold)
        
        precision_scores.append(precision_score(y_val_fold, y_val_pred, pos_label=positive_class))
        recall_scores.append(recall_score(y_val_fold, y_val_pred, pos_label=positive_class))
        f1_scores.append(f1_score(y_val_fold, y_val_pred, pos_label=positive_class))

    # Write the metrics to a text file
    with open(input_dir + '/' + model_name + '_Cross_validation_results.txt', 'w') as file:
        # Print the mean and standard deviation of each metric
        file.write(f'Cross-Validation Precision: {np.mean(precision_scores):.4f} ± {np.std(precision_scores):.4f}\n')
        file.write(f'Cross-Validation Recall: {np.mean(recall_scores):.4f} ± {np.std(recall_scores):.4f}\n')
        file.write(f'Cross-Validation F1-Score: {np.mean(f1_scores):.4f} ± {np.std(f1_scores):.4f}\n')

def permute_outcomes(feature_table_train, feature_of_interest_table_train, input_dir, num_permutations, class_weights_perm_dict):
    """_summary_

    Args:
        feature_table_train (_type_): Feature table without testing data
        feature_of_interest_table_train (_type_): Outcome table without testing data
        input_dir (_type_): _description_
        num_permutations (_type_): _description_

    Returns:
        _type_: _description_
    """
    RANDOM_SEEDS = [random.randint(0, 10000) for _ in range(num_permutations)]

    metrics_list = []
    shap_df_list = []
    test_data_report_list = []
    curn = 0

    # First document all seeds
    with open(input_dir + '/random_seeds.tsv', 'w') as file:
        for seed in RANDOM_SEEDS:
            file.write(f"{seed}\n")

    # Generate permuted features before splitting
    for seed in RANDOM_SEEDS:
        # List of series
        permuted_columns = []
        # Copy the full training data
        X_train_permuted_full = feature_table_train.copy()
        # print("doing "+str(curn))
        for col in feature_table_train.columns:
            permuted_col_name = f"{col}_perm"#_{seed}"
            permuted_col = feature_table_train[col].sample(frac=1, random_state=seed).reset_index(drop=True)
            permuted_col.name = permuted_col_name
            permuted_columns.append(permuted_col)

        # Concatenate permuted columns into a DataFrame
        permuted_df = pd.concat(permuted_columns, axis=1)
        # Concatenate permuted_df to X_train_permuted_full
        #X_train_permuted_full = pd.concat([X_train_permuted_full, permuted_df], axis=1)
        X_train_permuted_full = pd.concat([
            X_train_permuted_full.reset_index(drop=True),
            permuted_df.reset_index(drop=True)
        ], axis=1)

        print(f"After permutation - X_train_permuted_full shape: {X_train_permuted_full.shape}")
        print(f"After permutation - y_train_permuted_full shape: {feature_of_interest_table_train.shape}")

        # Just create train/val split
        # X_train, y_train, X_val, y_val, X_test, y_test, train_data, val_data, test_data = create_model_split(X_train_permuted_full, feature_of_interest_table_train, "associated_with_increased_likelihood_of_MONDO:0005101", RANDOM_SEED, test = True)
        X_train, y_train, X_val, y_val, train_data, val_data = create_model_split(X_train_permuted_full, feature_of_interest_table_train, "associated_with_increased_likelihood_of_MONDO:0005101", RANDOM_SEED, test = False)

        # Use this permutation val_data
        # shap_df_list, test_data_report_list = train_model_permutations(train_data, val_data, test_data, y_test, input_dir, RANDOM_SEED, metrics_list, shap_df_list, class_weights_perm_dict, test_data_report_list, False)
        shap_df_list = train_model_permutations(train_data, val_data, input_dir, RANDOM_SEED, metrics_list, shap_df_list, class_weights_perm_dict, False)

    # # Save report values DataFrame to CSV
    # test_data_report_df = pd.concat(test_data_report_list, ignore_index=True)
    # test_data_report_df.to_csv(input_dir + '/permuted_models_test_data_report.csv', index=True)
    ranked_features_df = rank_features_using_shap(shap_df_list, input_dir)

    return ranked_features_df

def train_model_permutations(train_data, val_data, input_dir, seed, metrics_list, shap_df_list, class_weights_perm_dict, intermediate_files=False): # train_data, val_data, test_data, y_test, input_dir, seed, metrics_list, shap_df_list, class_weights_perm_dict, test_data_report_list, intermediate_files=False):

    permuted_model_parameters = MODEL_PARAMETERS.copy() 
    permuted_model_parameters['random_seed'] = seed
    model = CatBoostClassifier(**permuted_model_parameters)
    #class_weights=class_weights_perm_dict)

    # Write the metrics to a metrics list
    model.fit(train_data, 
        eval_set=val_data,
        early_stopping_rounds=10
        )#, plot=True)

    # Predict on test data
    # y_pred = model.predict(test_data)
    # Generate classification report
    # report = classification_report(y_test, y_pred, output_dict=True)
    # Convert the classification report to a DataFrame
    # report_df = pd.DataFrame(report).transpose()
    # test_data_report_list.append(report_df)

    # # Save metrics to list- currently not output
    # metrics_list.append({
    #     'random_seed': seed,
    #     'accuracy': accuracy,
    #     'balanced_accuracy': balanced_accuracy,
    #     'auc_roc': auc_roc,
    #     'classification_report': report
    # })

    # Obtain SHAP values
    shap_values = model.get_feature_importance(data=train_data, type=EFstrType.ShapValues)
    shap_values = shap_values[:, :-1]  # Remove base value column

    # Convert SHAP values to DataFrame
    shap_df = pd.DataFrame(shap_values, columns=model.feature_names_)
    shap_df_list.append(shap_df)

    return shap_df_list #, test_data_report_list #metrics_list

def rank_features_using_shap(shap_df_list, input_dir):

    # Save SHAP values DataFrame to CSV
    shap_df = pd.concat(shap_df_list, ignore_index=True)
    shap_filename = "/all_features_shap_values.csv"
    shap_df.to_csv(input_dir + shap_filename, index=False)
    print(f"SHAP values saved to '{input_dir + shap_filename}'.")

    # Initialize a list to store mean SHAP values for each seed
    # List to store the mean SHAP values for each DataFrame
    # List to store the mean SHAP values for each DataFrame
    rank_shap_values_list = []
    shap_values_list = []

    # Process each DataFrame in shap_df_list
    for shap_df in shap_df_list:
        # Identify permuted feature names
        permuted_feature_indicator = ['_perm' in col for col in shap_df.columns]
        permuted_shap_values = shap_df.loc[:, permuted_feature_indicator]

        # Find the highest and lowest mean importance for permuted features
        max_permuted_importance = permuted_shap_values.mean(axis=0).max()
        min_permuted_importance = permuted_shap_values.mean(axis=0).min()
        print(f"min/max permuted importance {min_permuted_importance} {max_permuted_importance}")

        # Identify the original (non-permuted) features with mean importance above this threshold
        non_permuted_feature_indicator = ['_perm' not in col for col in shap_df.columns]
        original_shap_values = shap_df.loc[:, non_permuted_feature_indicator]

        # Calculate the mean SHAP values for non-permuted features
        mean_original_shap_values = original_shap_values.mean(axis=0)
        mean_original_shap_values.replace(0, np.nan, inplace=True)

        # Filter non-permuted features
        selected_features_positive = mean_original_shap_values[mean_original_shap_values > max_permuted_importance]
        selected_features_negative = mean_original_shap_values[mean_original_shap_values < min_permuted_importance]

        # Combine the selected features based on positive and negative criteria
        selected_features = pd.concat([selected_features_positive, selected_features_negative])

        # Rank the features
        ranked_features = selected_features.rank(method='average', ascending=True)

        # Append the ranked features to the list
        shap_values_list.append(mean_original_shap_values)
        rank_shap_values_list.append(ranked_features)

    # Combine the ranked features from all seeds into a single DataFrame
    shap_features_df = pd.concat(shap_values_list, axis=1)
    ranked_features_df = pd.concat(rank_shap_values_list, axis=1)

    # Define the file path for the TSV file
    shap_features_df_filename = "/shap_mean_across_seeds.tsv"
    ranked_features_df_filename = "/shap_rank_across_seeds.tsv"

    # Save the sorted DataFrame to a TSV file
    shap_features_df.to_csv(input_dir + shap_features_df_filename, sep='\t', header=True, index_label='Feature')
    ranked_features_df.to_csv(input_dir + ranked_features_df_filename, sep='\t', header=True, index_label='Feature')

    return ranked_features_df

def select_top_important_features(ranked_features_df, X_train, X_val, X_test, y_train, y_val, y_test, input_dir):
    """_summary_

    Args:
        ranked_features_df (DataFrame): ranked features from permutation test
        X_train (DataFrame): X_train from original model
        X_val (_type_): X_val from original model
        X_test (_type_): X_test from original model
        y_train (_type_): y_train from original model
        y_val (_type_): y_val from original model
        y_test (_type_): y_test from original model
    """

    # Drop features that are absent in more than a specified fraction of the columns
    min_non_na_count = 0#int(len(ranked_features_df.columns) * 0.75)  # Adjust the fraction as needed
    print("min_non_na_count")
    print(min_non_na_count)
    print("total cols in ranked_features_df")
    print(len(ranked_features_df.columns))
    ranked_features_df = ranked_features_df.dropna(thresh=min_non_na_count)

    # Calculate the mean rank for each feature
    mean_rank_across_seeds = ranked_features_df.mean(axis=1)

    # Sort the mean_rank_across_seeds DataFrame by the rank values
    sorted_mean_rank_across_seeds = mean_rank_across_seeds.sort_values()

    # Define the file path for the sorted mean rank TSV file
    sorted_mean_rank_filename = "/sorted_mean_rank_across_seeds.tsv"

    # Save the sorted DataFrame to a TSV file
    sorted_mean_rank_across_seeds.to_csv(input_dir + sorted_mean_rank_filename, sep='\t', header=True, index_label='Feature')

    print(f"Sorted mean rank values saved to {input_dir + sorted_mean_rank_filename}")

    # Convert the selected features to a list
    selected_feature_names = sorted_mean_rank_across_seeds.index.tolist()

    print("Selected Features with mean feature importance criteria:")
    print(len(selected_feature_names))
    print(selected_feature_names)

    X_train_selected = X_train[selected_feature_names]
    X_val_selected = X_val[selected_feature_names]
    X_test_selected = X_test[selected_feature_names]

    # Create Pool objects with the selected features only
    train_data_selected = Pool(data=X_train_selected, label=y_train)
    val_data_selected = Pool(data=X_val_selected, label=y_val)
    test_data_selected = Pool(data=X_test_selected, label=y_test)

    # Concatenate the selected features and labels for each set
    train_df_selected = pd.concat([X_train_selected, y_train.rename('label')], axis=1)
    val_df_selected = pd.concat([X_val_selected, y_val.rename('label')], axis=1)
    test_df_selected = pd.concat([X_test_selected, y_test.rename('label')], axis=1)

    train_file_path = 'taxa_to_media__data_df_clean__train_selected.tsv'
    val_file_path = 'taxa_to_media__data_df_clean__val_selected.tsv'
    test_file_path = 'taxa_to_media__data_df_clean__test_selected.tsv'

    # Save the DataFrames to TSV files
    train_df_selected.to_csv(input_dir + "/" + train_file_path, sep='\t', index=True, header=True)#, compression='gzip')
    val_df_selected.to_csv(input_dir + "/" + val_file_path, sep='\t', index=True, header=True)#, compression='gzip')
    test_df_selected.to_csv(input_dir + "/" + test_file_path, sep='\t', index=True, header=True)#, compression='gzip')

    return train_data_selected, val_data_selected, test_data_selected

