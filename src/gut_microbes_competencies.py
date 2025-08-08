## Ingests the download from HMP to list gut relevant microbes from stool, then queries kg-microbe for those microbes that 1. exist, 2. have at least 1 organismal trait, 3. have at least 1 functional annotation


import requests
import zipfile
import io
import os
import re
import pandas as pd
import tqdm
import duckdb

from Competencies import get_ncbitaxon_df, convert_ncbitaxon_label
from duckdb_utils import duckdb_load_table, query_with_edge_conditions, get_unique_values_with_substring, get_table_count
from ncbi_phylogeny_search import get_ncbitaxon_with_traits, get_ncbitaxon_with_uniprot

from constants import HMP_URL, HMP_ASSOCIATIONS_FILE, HMP_SHEET_NAME, HMP_STOOL_COLUMN_NAME, HMP_SITE_COLUMN_NAME, HMP_FEATURE_COLUMN_NAME, REPLACED_TAXA_NAMES, ORGANISMAL_TRAITS_EDGES

base_path = "Input_Files/hmp_supplementary/"

def load_data(input_dir, url):

    if not os.path.exists(input_dir) or not os.listdir(input_dir):
        print("hmp_supplementary directory is empty — downloading and extracting...")
        r = requests.get(url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(input_dir)

    print("extraction done")

def extract_hmp_data(input_dir, filename, sheet_name):

    supplementary_tables_df = pd.read_excel(input_dir + filename, skiprows=1, sheet_name=sheet_name)
    # Get only stool microbes
    subset_df = supplementary_tables_df[
        (supplementary_tables_df[HMP_SITE_COLUMN_NAME].str.lower() == HMP_STOOL_COLUMN_NAME) &
        (supplementary_tables_df[HMP_FEATURE_COLUMN_NAME].str.startswith("k__Bacteria"))
    ]

    unique_subset_df = subset_df[[HMP_FEATURE_COLUMN_NAME]].drop_duplicates() 

    return unique_subset_df

def create_ncbitaxon_dict():

    print("Creating NCBITaxon dictionary")
    # Get all NCBITaxon IDs
    ncbitaxon_label_dict = {}
    ncbitaxon_df = get_ncbitaxon_df()
    # Get NCBITaxon IDs from ontology nodes file
    for i, row in ncbitaxon_df.iterrows():
        ncbitaxon_label_dict[row["name"]] = row["id"]

    return ncbitaxon_label_dict

def get_ncbitaxon_id(df, ncbitaxon_label_dict, input_id):

    parts = input_id.strip("|").split("|")
    # Handle unclassified taxa
    if parts[-1].endswith("unclassified"):
        part2 = parts[-1].split("__")[1].replace("_", " ")
        part1 = parts[-2].split("__")[1].replace("_", " ")
        taxa = f"{part2} {part1}"
    else:
        taxa = input_id.split("|")[-1].replace("k__", "").replace("p__", "").replace("c__", "").replace("o__", "").replace("f__", "").replace("g__", "").replace("s__", "").replace("_", " ")
    # Replace name when necessary
    taxa = REPLACED_TAXA_NAMES.get(taxa, taxa)
    taxa_id = ncbitaxon_label_dict.get(taxa)
    if not taxa_id:
        # Try with brackets around genus name
        taxa_brackets = re.sub(r"^(\w+)", r"[\1]", taxa)
        taxa_id = ncbitaxon_label_dict.get(taxa_brackets)
        if not taxa_id:
            taxa_id = "not_found"
    return taxa_id, taxa

def get_taxa_rank(df, ncbitaxon_label_dict, input_id):

    rank_dict = {
        "k": "kingdom",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species"
    }

    parts = input_id.strip("|").split("|")
    # Handle unclassified taxa
    rank_symbol = parts[-1].split("__")[1]
    rank = rank_dict[rank_symbol]

    return rank

def main():

    output_dir = "./Intermediate_Files"
    
    load_data(base_path, HMP_URL)
    stool_microbes_df = extract_hmp_data(base_path, HMP_ASSOCIATIONS_FILE, HMP_SHEET_NAME)

    total_hmp_taxa = [len(stool_microbes_df)]

    ncbitaxon_label_dict = create_ncbitaxon_dict()
    
    # Get total taxa in graph
    hmp_microbes_names = []
    hmp_microbes_mapped_name = []
    hmp_microbes_mapped_id = []
    hmp_microbes_rank = []
    # Convert taxa names to NCBITaxon IDs
    for i in tqdm.tqdm(range(len(stool_microbes_df))):
        orig_taxa_id = stool_microbes_df.iloc[i].loc[HMP_FEATURE_COLUMN_NAME]
        hmp_microbes_names.append(orig_taxa_id)
        new_taxa_id, taxa_name = get_ncbitaxon_id(stool_microbes_df, ncbitaxon_label_dict, orig_taxa_id)
        taxa_rank = get_taxa_rank(stool_microbes_df, ncbitaxon_label_dict, orig_taxa_id)
        hmp_microbes_mapped_name.append(taxa_name)
        hmp_microbes_mapped_id.append(new_taxa_id)
        hmp_microbes_rank.append(taxa_rank)

    microbes_mapped_df = pd.DataFrame()
    microbes_mapped_df["HMP_Name"] = hmp_microbes_names
    microbes_mapped_df["NCBITaxon_Name"] = hmp_microbes_mapped_name
    microbes_mapped_df["NCBITaxon_ID"] = hmp_microbes_mapped_id
    microbes_mapped_df["NCBITaxon_Rank"] = hmp_microbes_rank

    total_hmp_mapped_taxa = [len(microbes_mapped_df[microbes_mapped_df["NCBITaxon_ID"] != "not_found"])]

    # Get total taxa with organismal traits
    conn = duckdb.connect(":memory:")
    duckdb_load_table(conn, "./Input_Files/kg-microbe-core/merged-kg_edges.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./Input_Files/kg-microbe-core/merged-kg_nodes.tsv", "nodes", ["id", "name"])

    ### To run using new traits search
    query_with_edge_conditions(conn, "edges", "organismal_traits_taxa", ORGANISMAL_TRAITS_EDGES)
    get_table_count(conn, "organismal_traits_taxa")
    total_hmp_organismal_taxa, ncbitaxon_traits_ids = get_unique_values_with_substring(conn, "organismal_traits_taxa", hmp_microbes_mapped_id)
    ###)
    print(len(ncbitaxon_traits_ids))

    ### To run using prior function; slower
    # ncbitaxon_traits_ids, ncbitaxon_traits_dict = get_ncbitaxon_with_traits(conn, output_dir)
    # total_hmp_organismal_taxa = len(set(ncbitaxon_traits_ids) & set(microbes_mapped_df["NCBITaxon_ID"].tolist()))
    ###

    hmp_organismal_taxa_query = microbes_mapped_df["NCBITaxon_ID"].isin(ncbitaxon_traits_ids).astype(str).str.upper()
    microbes_mapped_df["Has_Organismal_Trait"] = hmp_organismal_taxa_query

    # Get total taxa with functional annotations
    ncbitaxon_func_ids = get_ncbitaxon_with_uniprot(conn, "./Phylogeny_Search")
    total_hmp_func_microbes = len(set(ncbitaxon_func_ids) & set(microbes_mapped_df["NCBITaxon_ID"].tolist()))

    microbes_mapped_summary_df = pd.DataFrame(columns=["Total_HMP_Taxa", "Total_Taxa_Found", "Total_Taxa_Organismal_Traits", "Total_Taxa_Functional_Annotations"])
    microbes_mapped_summary_df["Total_HMP_Taxa"] = total_hmp_taxa
    microbes_mapped_summary_df["Total_Taxa_Found"] = total_hmp_mapped_taxa
    microbes_mapped_summary_df["Total_Taxa_Organismal_Traits"] = total_hmp_organismal_taxa
    microbes_mapped_summary_df["Total_Taxa_Functional_Annotations"] = total_hmp_func_microbes

    microbes_mapped_file = output_dir + '/HMP_Microbes_Mapped.csv'
    microbes_mapped_summary_file = output_dir + '/HMP_Microbes_Mapped_Summary.csv'
    microbes_mapped_df.to_csv(microbes_mapped_file,sep=",",index=False)
    microbes_mapped_summary_df.to_csv(microbes_mapped_summary_file,sep=",",index=False)

if __name__ == '__main__':
    main() 
