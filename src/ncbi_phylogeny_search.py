from collections import Counter
import json
import duckdb
from matplotlib import pyplot as plt
import numpy as np
from owlready2 import *
import pandas as pd
import tqdm

from duckdb_utils import duckdb_load_table, get_node_label

def get_all_ranks(output_dir):

    rank_file = './' + output_dir + '/ncbitaxon_rank.tsv'

    if not os.path.exists(rank_file):

        onto = get_ontology("http://purl.obolibrary.org/obo/ncbitaxon.owl")
        onto.load()
        data = []

        # Get the 'has_rank' property from the ontology
        has_rank = onto.search_one(iri="*has_rank")

        # Iterate through all classes
        for class_with_rank in onto.classes():
            rank = class_with_rank.has_rank  # Get the rank
            for rank in class_with_rank.has_rank:
                taxon_id = class_with_rank.iri.split("/")[-1].replace("_",":")  # Get the NCBITaxon ID
                rank_id = rank.iri.split("/")[-1].replace("NCBITaxon_","")  # Get the rank ID
                data.append([taxon_id, rank_id])  # Append to the data list

        df = pd.DataFrame(data, columns=['NCBITaxon_ID', 'Rank'])
        df.to_csv(rank_file,index=False)
    
    else:
        df = pd.read_csv(rank_file, index_col=False)

    return df

def load_graph():

    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    duckdb_load_table(conn, "./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", "edges", ["subject", "predicate", "object"])
    duckdb_load_table(conn, "./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_nodes.tsv", "nodes", ["id", "name"])

    return conn

def get_all_kg_taxa(conn):

    # Get NCBITaxon IDs from KG
    query = (
        f"""
        CREATE OR REPLACE TABLE all_taxa AS
        SELECT *
        FROM edges e
        WHERE e.subject LIKE '%NCBITaxon:%' OR e.object LIKE '%NCBITaxon:%'
        OR e.subject LIKE '%strain:%' OR e.object LIKE '%strain:%';
        SELECT * FROM all_taxa;
        """
    )
    # WHERE e.subject LIKE '%Proteomes:%' AND e.object LIKE '%NCBITaxon:%'

    result = conn.execute(query).fetchall()
    ncbitaxon_ids = [row[i] for row in result for i in [0, 2] if row[i].startswith(('NCBITaxon:','strain:'))]

    return list(set(ncbitaxon_ids))

def get_taxa_per_rank(ncbitaxa,rank,rank_df):

    # Get NCBITaxon IDs from owl file of given rank
    ranked_taxa = list(set(rank_df.loc[rank_df["Rank"] == rank, "NCBITaxon_ID"].tolist()))
    ranked_taxa = [i for i in ranked_taxa if i in ncbitaxa]

    return ranked_taxa


def search_lower_subclass_phylogeny(conn, microbe):

    query = (
        f"""
        CREATE TEMPORARY TABLE subclass_table AS
        SELECT *
        FROM edges
        WHERE subject LIKE '%NCBITaxon:%' AND predicate = 'biolink:subclass_of' AND object = '{microbe}';
        SELECT * FROM subclass_table;
        """
    )
    result = conn.execute(query).fetchall()
    conn.execute("DROP TABLE IF EXISTS subclass_table")
    # Get object from triple as parent taxa
    try:
        child_taxa = [row[0] for row in result]
    except IndexError:
        child_taxa = 'not found'

    return child_taxa

def search_subclass_phylogeny_parent(conn, microbe):

    query = (
        f"""
        CREATE TEMPORARY TABLE subclass_table AS
        SELECT *
        FROM edges
        WHERE subject = '{microbe}' AND predicate = 'biolink:subclass_of' AND object LIKE '%NCBITaxon:%';
        SELECT * FROM subclass_table;
        """
    )
    result = conn.execute(query).fetchall()
    conn.execute("DROP TABLE IF EXISTS subclass_table")
    # Get object from triple as parent taxa
    try:
        parent_taxa = [row[2] for row in result][0]
    except IndexError:
        parent_taxa = 'not found'


    return parent_taxa

# Get taxa from kg that is classified lower than phylum only
def find_relevant_taxa(phyla, ncbitaxa, ncbi_taxa_ranks_df):




    # Get only class, order, family, genus, species, strain
    kingdom = get_taxa_per_rank(ncbitaxa, "kingdom", ncbi_taxa_ranks_df)
    superkingdom = get_taxa_per_rank(ncbitaxa, "superkingdom", ncbi_taxa_ranks_df)
    relevant_ncbitaxa = list(set([taxon for taxon in ncbitaxa 
                    if taxon not in kingdom and 
                        taxon not in superkingdom and 
                        taxon not in phyla and 
                        taxon != "NCBITaxon:1"]))
    
    return relevant_ncbitaxa

def search_strains(conn, ncbi_taxa_ranks_df, microbe, strains_found, species_found):
    """Recursive function to search for strains."""
    # Get the subclasses of the current microbe
    child_taxa = search_lower_subclass_phylogeny(conn, microbe)

    for child in child_taxa:
        # Check if the child is a strain
        microbe_rank = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].get(child, None)
        if microbe_rank is None:
            continue
        if microbe_rank == "species":
            if child not in strains_found:
                species_found.append(child)
        if microbe_rank in ["subspecies","strain"]:
            if child not in strains_found:
                strains_found.append(child)
        else:
            # Continue searching subclasses
            search_strains(conn, ncbi_taxa_ranks_df, child, strains_found, species_found)

def find_all_strains(conn, ncbi_taxa_ranks_df, microbes,microbes_traits_strain, microbes_traits_species):
    """
    Finds all strains for a list of microbes by recursively searching subclasses.

    Parameters:
        conn: Database connection object.
        ncbi_taxa_ranks_df: DataFrame containing taxonomy ranks.
        microbes: List of original microbes to search.
        strain_rank: The rank to identify as a strain (default is "strain").

    Returns:
        dict: Dictionary where keys are original microbes and values are lists of strains.
    """

    # Perform search for each microbe
    for microbe in tqdm.tqdm(microbes, desc="Processing microbes"):
        print(f"Starting search for microbe: {microbe}")
        if microbe not in microbes_traits_strain.keys():
            microbe_rank = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].get(microbe, None)
            print(microbe_rank)
            # if microbe_rank not in ["phylum","class"] and microbe != "NCBITaxon:1280" and microbe != "NCBITaxon:1763":
            # if (microbe_rank in ["phylum","class"] or microbe == "NCBITaxon:1763" or microbe == "NCBITaxon:1280") and microbe != "NCBITaxon:1224" and microbe != "NCBITaxon:1236" and microbe != "NCBITaxon:1239":
            search_strains(conn, ncbi_taxa_ranks_df, microbe, microbes_traits_strain[microbe], microbes_traits_species[microbe])
            # else:
            #     microbes_traits_strain[microbe].extend([])            

    return microbes_traits_strain, microbes_traits_species


def get_microbe_strain_children(conn, ncbi_taxa_ranks_df, original_microbe, microbe, strain, microbes_strain):

    while strains_found == False:
        child_taxa = search_lower_subclass_phylogeny(conn, microbe)
        if any(child in strain for child in child_taxa):
            found_strains = [child for child in child_taxa if child in strain]
            microbes_strain[original_microbe].extend(found_strains)
            strains_found = True
        else:
            for i,c in enumerate(child_taxa): 
                result = get_microbe_strain_children(conn, ncbi_taxa_ranks_df, original_microbe, c, strain, microbes_strain)
        
    return microbes_strain



    # # microbe_list = [microbe]
    # print("top of function")
    # print("original_microbe: ",original_microbe,"microbe: ",microbe)
    # strain_found = False
    # while not strain_found:
    #     microbe_rank = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].get(microbe, None)
    #     if not microbe_rank in (["strain", "subspecies"]):
    #         child_taxa = search_lower_subclass_phylogeny(conn, microbe)
    #         if child_taxa == ["not found"] or len(child_taxa) == 0:
    #             continue
    #         else:
    #             for i,c in enumerate(child_taxa):
    #                 microbe_rank = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].get(c, None)
    #                 if microbe_rank in (["strain", "subspecies"]):
    #                     microbes_strain[original_microbe].append(c)
    #                 else:
    #                     child_taxa = search_lower_subclass_phylogeny(conn, microbe)

    #                 is_last = (i + 1 == len(child_taxa))
    #                 child_taxa = search_lower_subclass_phylogeny(conn, c)


        
    #     print(microbe_rank)
    #     child_taxa = search_lower_subclass_phylogeny(conn, microbe)
    #     print(microbe)
    #     print(child_taxa)


    #     #print("child taxa returned from subclass search")
    #     #print(child_taxa)
    #     #import pdb;pdb.set_trace()
    #     if child_taxa == ["not found"] or len(child_taxa) == 0 or not microbe_rank:
    #         strain_found = True
    #         microbes_strain[original_microbe].extend([])
    #         #print("strain not found, microbe: ",microbe)
    #         return microbes_strain

    #     #if child_taxa[0] in strain:
    #     if any(child in strain for child in child_taxa):
    #         found_strains = [child for child in child_taxa if child in strain]
    #         print("adding strains: ",len(found_strains))
    #         microbes_strain[original_microbe].extend(found_strains)
    #         #print(microbes_strain)
    #         if not more_to_search:
    #             strain_found = True
    #             #print(microbes_strain)
    #             #return microbes_strain
    #             #print("added strains for all microbe: ",original_microbe)
    #             #return microbes_strain
    #         return microbes_strain
    #         #pass

    #     strains_found_in_children = False   
    #     for i,c in enumerate(child_taxa):
    #         is_last = (i + 1 == len(child_taxa))
    #         print(i)
    #         print(is_last)
    #         #if i+1 < len(child_taxa):
    #         #print("iterating: ",i)
    #         #print("original_microbe: ",original_microbe,"c: ",c)
    #         result = get_microbe_strain_children(conn, ncbi_taxa_ranks_df, original_microbe, c, strain, microbes_strain,more_to_search = not is_last)
        
    #         unique_taxa = [taxon for taxon in result[original_microbe] if taxon not in microbes_strain[original_microbe]]
    #         if unique_taxa:
    #             strains_found_in_children = True
    #             microbes_strain[original_microbe].extend(unique_taxa)

    #         #if result and original_microbe in result:
    #         #    strains_found_in_children = True
    #         #    # Add the strains from the result to the current strains
    #         #    microbes_strain[original_microbe].extend(result[original_microbe])

    #     if strains_found_in_children:
    #         strain_found = True

    #     #import pdb;pdb.set_trace()
    # print("end of function")
    # print(microbes_strain)
    # return microbes_strain

def get_microbe_species(conn, microbe, species, microbes_species):

    microbe_list = [microbe]
    species_found = False
    while not species_found:
        parent_taxa = search_subclass_phylogeny_parent(conn, microbe)

        if parent_taxa in species or parent_taxa == 'not found':
            species_found = True
            microbes_species[parent_taxa].extend(microbe_list)
        # Keep track of each bug in the phylogeny
        else:
            microbe_list.append(parent_taxa)
            microbe = parent_taxa

    return microbes_species

def get_microbe_parent_rank(conn, microbe, all_of_rank, microbes_rank):

    microbe_list = [microbe]
    rank_found = False
    while not rank_found:
        parent_taxa = search_subclass_phylogeny_parent(conn, microbe)

        if parent_taxa in all_of_rank or parent_taxa == 'not found':
            rank_found = True
            microbes_rank[parent_taxa].extend(microbe_list)
        # Keep track of each bug in the phylogeny
        else:
            microbe_list.append(parent_taxa)
            microbe = parent_taxa

    return microbes_rank

def get_microbe_family(conn, microbe, family, microbes_family):

    microbe_list = [microbe]
    family_found = False
    while not family_found:
        parent_taxa = search_subclass_phylogeny_parent(conn, microbe)

        if parent_taxa in family or parent_taxa == 'not found':
            family_found = True
            microbes_family[parent_taxa].extend(microbe_list)
        # Keep track of each bug in the phylogeny
        else:
            microbe_list.append(parent_taxa)
            microbe = parent_taxa

    return microbes_family

def get_microbe_phylum(conn, microbe, phyla, microbes_phyla):

    microbe_list = [microbe]
    phyla_found = False
    while not phyla_found:
        parent_taxa = search_subclass_phylogeny_parent(conn, microbe)

        if parent_taxa in phyla or parent_taxa == 'not found':
            phyla_found = True
            microbes_phyla[parent_taxa].extend(microbe_list)
        # Keep track of each bug in the phylogeny
        else:
            microbe_list.append(parent_taxa)
            microbe = parent_taxa

    return microbes_phyla

def create_species_strains_dictionary(output_dir):

    microbes_traits_strain = defaultdict(list)
    # To also keep track of species if strains are not found
    microbes_traits_species = defaultdict(list)

    dictionary_map = {
        "_microbes_strain.json" : microbes_traits_strain,
        "_microbes_species.json" : microbes_traits_species
    }

    for file_substring, dic in dictionary_map.items():
        for root, _, files in os.walk(output_dir):
            for file in files:
                if file_substring in file:
                    filepath = os.path.join(root, file)
                    with open(filepath, 'r') as f:
                        try:
                            data = json.load(f)  # Load JSON data
                            for key, value in data.items():
                                dic[key].extend(value if isinstance(value, list) else [value])
                        except json.JSONDecodeError as e:
                            print(f"Error reading {filepath}: {e}")

    return microbes_traits_strain, microbes_traits_species

def find_microbes_strain(conn, ncbi_taxa_ranks_df, all_taxa, output_dir, feature_type):
    '''Takes in list of all relevant taxa or just 1 microbe'''
    # microbes_traits_strain_file = './' + output_dir + '/' + feature_type + '_microbes_strain_phyla_class_1280_1763.json'
    # microbes_traits_species_file = './' + output_dir + '/' + feature_type + '_microbes_species_phyla_class_1280_1763.json'

    microbes_traits_strain_file = './' + output_dir + '/' + feature_type + '_microbes_strain.json'
    microbes_traits_species_file = './' + output_dir + '/' + feature_type + '_microbes_species.json'

    print("Getting strain for all taxa: " + feature_type)
    if not os.path.exists(microbes_traits_strain_file):
        # microbes_traits_strain = defaultdict(list)
        # # To also keep track of species if strains are not found
        # microbes_traits_species = defaultdict(list)
        microbes_traits_strain, microbes_traits_species = create_species_strains_dictionary(output_dir)
        microbes_traits_strain, microbes_traits_species = find_all_strains(conn, ncbi_taxa_ranks_df, all_taxa, microbes_traits_strain, microbes_traits_species)

        # for microbe in tqdm.tqdm(all_taxa): #["NCBITaxon:853"]:#tqdm.tqdm(relevant_ncbitaxa):
        #     print("microbe from all taxa: ",microbe)
        #     microbes_traits_strain = get_microbe_strain_children(conn, ncbi_taxa_ranks_df, microbe, microbe, strain, microbes_traits_strain)
        #     print("after first in all taxa")
        #     print(microbes_traits_strain)
        #     import pdb;pdb.set_trace()
        with open(microbes_traits_strain_file, 'w') as json_file:
            json.dump(microbes_traits_strain, json_file, indent=4)
        with open(microbes_traits_species_file, 'w') as json_file:
            json.dump(microbes_traits_species, json_file, indent=4)
    else:
        with open(microbes_traits_strain_file, 'r') as json_file:
            data = json.load(json_file)
            microbes_traits_strain = defaultdict(list, data)
        with open(microbes_traits_species_file, 'r') as json_file:
            data = json.load(json_file)
            microbes_traits_species = defaultdict(list, data)

    return microbes_traits_strain, microbes_traits_species

def find_microbes_species(conn, ncbi_taxa_ranks_df, all_taxa, output_dir, feature_type):
    '''Takes in list of all relevant taxa or just 1 microbe'''
    microbes_traits_species_file = './' + output_dir + '/' + feature_type + '_microbes_species.json'

    print("Getting species for all taxa: " + feature_type)
    if not os.path.exists(microbes_traits_species_file):
        # Get only species
        # species = get_taxa_per_rank(all_taxa, "species", ncbi_taxa_ranks_df)
        species = list(set(ncbi_taxa_ranks_df.loc[ncbi_taxa_ranks_df["Rank"] == "species", "NCBITaxon_ID"].tolist()))
        microbes_traits_species = defaultdict(list)
        for microbe in tqdm.tqdm(all_taxa): #["NCBITaxon:853"]:#tqdm.tqdm(relevant_ncbitaxa):
            microbe_rank = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].get(microbe, None)
            if microbe_rank in ("subspecies", "strain"):
                microbes_traits_species = get_microbe_species(conn, microbe, species, microbes_traits_species,first_call=True)
            elif microbe_rank is None:
                print("Rank not found")
                print(microbe)
            else:
                continue
        with open(microbes_traits_species_file, 'w') as json_file:
            json.dump(microbes_traits_species, json_file, indent=4)
    else:
        with open(microbes_traits_species_file, 'r') as json_file:
            data = json.load(json_file)
            microbes_traits_species = defaultdict(list, data)

    return microbes_traits_species

def find_microbes_rank(conn, ncbi_taxa_ranks_df, all_taxa, output_dir, feature_type, rank):
    '''Takes in list of all relevant taxa or just 1 microbe'''
    microbes_traits_rank_file = './' + output_dir + '/' + feature_type + '_microbes_' + rank + '.json'

    print("Getting rank for all taxa: " + feature_type + ", " + rank)
    if not os.path.exists(microbes_traits_rank_file):
        # Get only family
        # rank = get_taxa_per_rank(all_taxa, "genus", ncbi_taxa_ranks_df)
        all_of_rank = list(set(ncbi_taxa_ranks_df.loc[ncbi_taxa_ranks_df["Rank"] == rank, "NCBITaxon_ID"].tolist()))
        microbes_traits_rank = defaultdict(list)
        for microbe in tqdm.tqdm(all_taxa): #["NCBITaxon:853"]:#tqdm.tqdm(relevant_ncbitaxa):
            microbes_traits_rank = get_microbe_parent_rank(conn, microbe, all_of_rank, microbes_traits_rank)
        with open(microbes_traits_rank_file, 'w') as json_file:
            json.dump(microbes_traits_rank, json_file, indent=4)
    else:
        with open(microbes_traits_rank_file, 'r') as json_file:
            data = json.load(json_file)
            microbes_traits_rank = defaultdict(list, data)

    return microbes_traits_rank


def find_microbes_family(conn, ncbi_taxa_ranks_df, all_taxa, output_dir, feature_type):
    '''Takes in list of all relevant taxa or just 1 microbe'''
    microbes_traits_family_file = './' + output_dir + '/' + feature_type + '_microbes_family.json'

    print("Getting family for all taxa: " + feature_type)
    if not os.path.exists(microbes_traits_family_file):
        # Get only family
        # family = get_taxa_per_rank(all_taxa, "family", ncbi_taxa_ranks_df)
        family = list(set(ncbi_taxa_ranks_df.loc[ncbi_taxa_ranks_df["Rank"] == "family", "NCBITaxon_ID"].tolist()))
        microbes_traits_family = defaultdict(list)
        for microbe in tqdm.tqdm(all_taxa): #["NCBITaxon:853"]:#tqdm.tqdm(relevant_ncbitaxa):
            microbes_traits_family = get_microbe_family(conn, microbe, family, microbes_traits_family)
        with open(microbes_traits_family_file, 'w') as json_file:
            json.dump(microbes_traits_family, json_file, indent=4)
    else:
        with open(microbes_traits_family_file, 'r') as json_file:
            data = json.load(json_file)
            microbes_traits_family = defaultdict(list, data)

    return microbes_traits_family

def find_microbes_phylum(conn, ncbi_taxa_ranks_df, all_taxa, output_dir, feature_type):
    '''Takes in list of all relevant taxa or just 1 microbe'''
    microbes_traits_phyla_file = './' + output_dir + '/' + feature_type + '_microbes_phyla.json'

    print("Getting phylum for all taxa")
    if not os.path.exists(microbes_traits_phyla_file):
        phyla = list(set(ncbi_taxa_ranks_df.loc[ncbi_taxa_ranks_df["Rank"] == "phylum", "NCBITaxon_ID"].tolist()))
        microbes_traits_phyla = defaultdict(list)
        for microbe in tqdm.tqdm(all_taxa): #["NCBITaxon:853"]:#tqdm.tqdm(relevant_ncbitaxa):
            microbes_traits_phyla = get_microbe_phylum(conn, microbe, phyla, microbes_traits_phyla)
        with open(microbes_traits_phyla_file, 'w') as json_file:
            json.dump(microbes_traits_phyla, json_file, indent=4)
    else:
        with open(microbes_traits_phyla_file, 'r') as json_file:
            data = json.load(json_file)
            microbes_traits_phyla = defaultdict(list, data)

    return microbes_traits_phyla

def get_ncbitaxon_with_traits(conn, output_dir):

    traits_prefixes = ['carbon_substrates', 'pathways', 'production', 'CAS-RN', 'CHEBI',
    'EC', 'GO', 'cell_shape', 'gc', 'gram_stain', 'pH_.%', 'temp_.%', 
    'pigment', 'sporulation', 'trophic_type', 'motility', 'temperature',
    'salinity', 'NaCl_.%', 'oxygen', 'pathogen', 'isolation_source', 
    'ENVO', 'UBERON', 'PO', 'PATO', 'cell_length', 'cell_width', 'FOODON', 'assay', 'medium'] #'strain',

    trait_prefixes_subject_query =  '(' + ' OR '.join([f"split_part(e.subject, \':\', 1) = '{term}'" for term in traits_prefixes]) + ')'
    trait_prefixes_object_query =  '(' + ' OR '.join([f"split_part(e.object, \':\', 1) = '{term}'" for term in traits_prefixes]) + ')'

    query = (
        f"""
        CREATE OR REPLACE TABLE ncbitaxon_to_trait AS
        SELECT *
        FROM edges e
        WHERE (split_part(e.subject, ':', 1) = 'NCBITaxon'
        AND {trait_prefixes_object_query})
        OR
        (split_part(e.subject, ':', 1) = 'strain'
        AND {trait_prefixes_object_query})
        OR
        (split_part(e.object, ':', 1) = 'NCBITaxon'
        AND {trait_prefixes_subject_query})
        OR
        (split_part(e.object, ':', 1) = 'strain'
        AND {trait_prefixes_subject_query});
        SELECT * FROM ncbitaxon_to_trait;
        """
    )

    result = conn.execute(query).fetchall()
    ncbitaxon_traits_ids = [row[i] for row in result for i in [0, 2] if row[i].startswith(('NCBITaxon:', 'strain:'))]
    ncbitaxon_traits_ids = list(set(ncbitaxon_traits_ids))

    ncbitaxon_traits_dict_file = output_dir + "/ncbitaxon_traits_dict.json"
    if not os.path.exists(ncbitaxon_traits_dict_file):
        # To create a dictionary of trait and each taxa id associated with it
        ncbitaxon_traits_dict = defaultdict(list)
        for subj, _, obj in tqdm.tqdm(result):
            if "NCBITaxon:" in subj or "strain" in subj:
                trait, ncbi_taxon = obj, subj
            elif "NCBITaxon:" in obj or "strain" in obj:
                trait, ncbi_taxon = subj, obj
            trait_label = get_node_label(conn, trait)
            ncbitaxon_traits_dict[trait_label].append(ncbi_taxon)
        with open(ncbitaxon_traits_dict_file, 'w') as json_file:
            json.dump(ncbitaxon_traits_dict, json_file, indent=4)
    else:
        with open(ncbitaxon_traits_dict_file, 'r') as json_file:
            data = json.load(json_file)
            ncbitaxon_traits_dict = defaultdict(list, data)
            # Get all unique taxa with a trait
            ncbitaxon_traits_ids = list(set(value for values in ncbitaxon_traits_dict.values() for value in values))

    return ncbitaxon_traits_ids, ncbitaxon_traits_dict

def get_ncbitaxon_with_uniprot(conn, output_dir):

    ncbitaxon_uniprot_file = output_dir + "/unique_ncbitaxon_uniprot_ids.txt"

    if not os.path.exists(ncbitaxon_uniprot_file):
        query = (
            f"""
            SELECT DISTINCT split_part(object, ':', 2) AS ncbi_taxon
            FROM edges
            WHERE split_part(subject, ':', 1) = 'UniprotKB'
            AND split_part(object, ':', 1) = 'NCBITaxon'
            AND predicate = 'biolink:derives_from';
            """
        )

        unique_ncbitaxon = conn.execute(query).fetchall()

        # Convert the results to a list
        unique_ncbitaxon_list = ["NCBITaxon:" + row[0] for row in unique_ncbitaxon]

        with open(ncbitaxon_uniprot_file, "w") as f:
            for taxon in unique_ncbitaxon_list:
                f.write(f"{taxon}\n")

    else:
        with open(ncbitaxon_uniprot_file, 'r') as f:
            unique_ncbitaxon_list = [line.strip() for line in f]

    return unique_ncbitaxon_list


def get_ncbitaxon_with_functional_annotation(conn, output_dir):

    ncbitaxon_func_dict_file = output_dir + "/ncbitaxon_func_dict.json"
    functional_mappings_df = pd.read_csv("./data/Intermediate_Files/NCBITaxon_to_GO.tsv", sep='\t',index_col=False)

    if not os.path.exists(ncbitaxon_func_dict_file):
        ncbitaxon_func_dict = functional_mappings_df.groupby('GO')['NCBITaxon'].apply(lambda x: list(set(x))).to_dict()
        ncbitaxon_func_dict = {get_node_label(conn, k): v for k, v in ncbitaxon_func_dict.items()}
        with open(ncbitaxon_func_dict_file, 'w') as json_file:
            json.dump(ncbitaxon_func_dict, json_file, indent=4)
    else:
        with open(ncbitaxon_func_dict_file, 'r') as json_file:
            data = json.load(json_file)
            ncbitaxon_func_dict = defaultdict(list, data)

    ncbitaxon_func_ids = functional_mappings_df["NCBITaxon"].unique().tolist()
    return ncbitaxon_func_ids, ncbitaxon_func_dict

def differentiate_taxa_by_phylum(conn, ncbitaxon_traits_ids, phyla_dict, output_dir, ncbitaxon_traits_dict, feature_type):

    print("Differentiating taxa by phylum")

    phyla_counts_file = './' + output_dir + '/' + feature_type + '_microbes_phyla_counts.json'
    phyla_traits_dict_file = './' + output_dir + '/' + feature_type + '_phyla_dict.json'

    if not os.path.exists(phyla_traits_dict_file):
        phyla_traits_dict = defaultdict(list)

    if not (os.path.exists(phyla_counts_file) and os.path.exists(phyla_traits_dict_file)):
        phyla_counts = defaultdict(int)
        for ncbitaxon in tqdm.tqdm(ncbitaxon_traits_ids):
            phylum = find_key_by_value(phyla_dict, ncbitaxon)
            phylum_label = get_node_label(conn, phylum)
            phyla_counts[phylum_label] += 1
            # Create counts of traits per phylum
            if not os.path.exists(phyla_traits_dict_file):
                traits = [key for key, values in ncbitaxon_traits_dict.items() if ncbitaxon in values]
                phyla_traits_dict[phylum_label].extend(traits)

        # Sort
        phyla_counts = dict(sorted(phyla_counts.items(), key=lambda item: item[1], reverse=True))

        with open(phyla_counts_file, 'w') as json_file:
            json.dump(phyla_counts, json_file, indent=4)
    
    else:
        with open(phyla_counts_file, 'r') as json_file:
            data = json.load(json_file)
            phyla_counts = dict(data)

    if os.path.exists(phyla_traits_dict_file):
        with open(phyla_traits_dict_file, 'r') as json_file:
            data = json.load(json_file)
            phyla_traits_dict = dict(data)
    else:
        with open(phyla_traits_dict_file, 'w') as json_file:
            json.dump(phyla_traits_dict, json_file, indent=4)

    return phyla_counts, phyla_traits_dict

def find_key_by_value(dict, value):
    key = next((k for k, v in dict.items() if value in v), "None")
    
    return key

def plot_taxa_by_phylum(count_dict, output_dir, feature_type):

    gut_phyla = ["Bacteroidota", "Bacillota", "Actinomycetota", "Pseudomonadota", "Fusobacteriota"]

    # Plot histogram
    plt.bar(count_dict.keys(), count_dict.values(), color='skyblue')
    plt.xlabel("Organism")
    plt.ylabel("Count (log)")
    plt.title("Microbial " + feature_type.replace("_"," ") + " by Phylum")
    plt.yscale('log')
    plt.xticks(rotation=90, fontsize = 6)
    plt.subplots_adjust(bottom=0.3)


    plt.savefig('./' + output_dir + '/' + feature_type + '_organismal_phyla_histogram.png', format='png', dpi=300)

def plot_taxa_by_phylum_and_feature(phyla_traits_dict, output_dir, feature_type):

    # Count occurrences of each phylum for each trait
    count_data = {key: Counter(values) for key, values in phyla_traits_dict.items()}

    # Convert counts to DataFrame, filling missing values with 0
    df = pd.DataFrame(count_data).fillna(0)
    # Categorize traits with maximum counts across all phyla < min_number_phylum_trait
    if feature_type == "Traits":
        min_number_phylum_trait = 10
    elif feature_type == "Functional_Annotations":
        min_number_phylum_trait = 2490
    to_group = df[df.max(axis=1) < min_number_phylum_trait]
    df.loc['other'] = to_group.sum()
    df = df[df.max(axis=1) >= min_number_phylum_trait]

    # Sort columns by sum in descending order
    column_totals = df.sum(axis=0)
    sorted_columns = column_totals.sort_values(ascending=False).index
    df = df[sorted_columns]

    # Combine colors from different tab colormaps
    colors = np.concatenate([
        plt.cm.tab20.colors,  # First 20 colors
        plt.cm.tab20b.colors,  # Next 20 colors
        plt.cm.tab20c.colors   # Additional 20 colors
    ])

    # Plot with the combined colors
    ax = df.T.plot(kind='bar', stacked=True, figsize=(18, 9), color=colors[:len(df.columns)])  # Limit to the number of columns
    # ax = df.T.plot(kind='bar', stacked=True, figsize=(18, 9), color=plt.cm.tab60.colors)

    # Adding titles and labels
    ax.set_title('Microbial ' + feature_type.replace("_"," ") + ' by Phylum')
    ax.set_xlabel('Phylum')
    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Counts (log)')
    plt.yscale('log')
    handles, labels = ax.get_legend_handles_labels()
    labels = [label.replace('_', ' ') for label in labels]  # Replace underscores with spaces
    ax.legend(handles, labels, title=feature_type, bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2, frameon=False)

    # plt.legend(title=feature_type, bbox_to_anchor=(1.05, 1), loc='upper left', ncol = 2, frameon=False)  # Place legend outside the plot
    # plt.subplots_adjust(left=0.05, right=0.6, top=0.9, bottom=0.3)
    plt.subplots_adjust(left=0.05, right=0.85, top=0.9, bottom=0.3)

    plt.tight_layout()
    plt.savefig('./' + output_dir + '/' + feature_type + '_organismal_phyla_stacked_barplot.png', format='png', dpi=300)

def main():

    output_dir = "./data/Phylogeny_Search"
    os.makedirs(output_dir, exist_ok=True)

    ncbi_taxa_ranks_df = get_all_ranks(output_dir)

    # # For all phyla
    conn = load_graph()
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

    # Get only phyla
    # phyla = get_taxa_per_rank([], "phylum", ncbi_taxa_ranks_df) #ncbitaxa

    ncbitaxon_traits_ids, ncbitaxon_traits_dict = get_ncbitaxon_with_traits(conn, output_dir)
    print("Len Taxa with a trait")
    print(len(ncbitaxon_traits_ids))

    # If desired to only get phylum assignments for specific subset of microbes
    microbes_traits_phyla = find_microbes_phylum(conn, ncbi_taxa_ranks_df, ncbitaxon_traits_ids, output_dir, "Traits")
    microbes_func_phyla = find_microbes_phylum(conn, ncbi_taxa_ranks_df, ncbitaxon_func_ids, output_dir, "Functional_Annotation")

    # If desired to produce microbes_phyla.json with all microbe phylum assignments
    # relevant_ncbitaxa = find_relevant_taxa(phyla, ncbitaxa, ncbi_taxa_ranks_df)
    # microbes_phyla = find_microbes_phylum(conn, phyla, relevant_ncbitaxa)

    phyla_traits_counts,phyla_traits_dict = differentiate_taxa_by_phylum(conn, ncbitaxon_traits_ids, microbes_traits_phyla, output_dir, ncbitaxon_traits_dict, "Traits")
    phyla_func_counts,phyla_func_dict = differentiate_taxa_by_phylum(conn, ncbitaxon_func_ids, microbes_func_phyla, output_dir, ncbitaxon_func_dict, "Functional_Annotations")

    plot_taxa_by_phylum(phyla_traits_counts, output_dir, "Traits")
    plot_taxa_by_phylum(phyla_func_counts, output_dir, "Functional_Annotations")

    plot_taxa_by_phylum_and_feature(phyla_traits_dict, output_dir, "Traits")
    plot_taxa_by_phylum_and_feature(phyla_func_dict, output_dir, "Functional_Annotations")

if __name__ == '__main__':
    main()
