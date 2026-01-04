

import os
import duckdb

from duckdb_utils import create_subject_object_pair_table, duckdb_load_table, get_table_count, get_total_unique_column, get_total_unique_pair, join_tables_subject_object, join_tables_unique_subject_object, output_table_to_file, query_with_multiple_conditions


def main():

    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    #duckdb_load_table(conn, "/Users/brooksantangelo/Documents/LozuponeLab/FRMS_2024/uniprot_transform_20240725/edges.tsv", "edges", ["subject", "object"])
    # duckdb_load_table(conn, "./data/merged-kg/merged-kg_edges.tsv", "edges", ["subject", "object"])
    duckdb_load_table(conn, "./data/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv", "edges", ["subject", "object"])
    conditions = [
        ('UniprotKB:%', 'NCBITaxon:%'),
        ('UniprotKB:%', 'GO:%')
    ]

    if not os.path.exists("./data/Intermediate_Files"):
        os.makedirs("./data/Intermediate_Files")
    
    print("Relevant edges loaded.")

    query_with_multiple_conditions(
        conn,
        table_name = "edges",
        conditions = conditions)

    print("Subset to UniprotKB-NCBITaxon pairs.")

    create_subject_object_pair_table(
        conn,
        table_name = "uniprot_ncbitaxon",
        base_table_name = "edges",
        subject = "UniprotKB",
        object = "NCBITaxon",
        subject_prefix = "UniprotKB:%",
        object_prefix = "NCBITaxon:%"
    )

    print("Subset to Uniprot-GO edges.")

    create_subject_object_pair_table(
        conn,
        table_name = "uniprot_go",
        base_table_name = "edges",
        subject = "UniprotKB",
        object = "GO",
        subject_prefix = "UniprotKB:%",
        object_prefix = "GO:%"
    )

    print("Subset to NCBITaxon-GO edges.")

    join_tables_unique_subject_object(
         conn, 
        base_table_name = "uniprot_go", 
        compared_table_name = "uniprot_ncbitaxon", 
        output_table_name = "ncbitaxon_go", 
        output_subject = "NCBITaxon", 
        output_object = "GO", 
        comparison = "UniprotKB"
    )

    get_table_count(conn, "ncbitaxon_go")

    output_table_to_file(conn, "ncbitaxon_go", "./data/Intermediate_Files/NCBITaxon_to_GO.tsv")

    num_unique_ncbitaxon_go_pairs = get_total_unique_pair(
        conn,
        table_name = "ncbitaxon_go",
        subject = "NCBITaxon",
        object = "GO",
        subject_prefix = "NCBITaxon:%",
        object_prefix = "GO:%"
    )
    
    print("Total NCBITaxon-GO Pairs: ", num_unique_ncbitaxon_go_pairs[0])
    
    num_unique_ncbitaxon = get_total_unique_column(
        conn,
        table_name = "ncbitaxon_go",
        column_name = "NCBITaxon",
    )

    print("Total unique NCBITaxon:", num_unique_ncbitaxon[0])

if __name__ == '__main__':
    main()
