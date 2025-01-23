

import os
import duckdb

from duckdb_utils import create_subject_object_pair_table, duckdb_load_table, get_table_count, get_total_unique_column, get_total_unique_pair, join_tables_subject_object, join_tables_unique_subject_object, output_table_to_file, query_with_multiple_conditions


def main():

    # Create a DuckDB connection
    conn = duckdb.connect(":memory:")

    print("Loading full table.")

    #duckdb_load_table(conn, "/Users/brooksantangelo/Documents/LozuponeLab/FRMS_2024/uniprot_transform_20240725/edges.tsv", "edges", ["subject", "object"])
    # duckdb_load_table(conn, "./Input_Files/merged-kg/merged-kg_edges.tsv", "edges", ["subject", "object"])
    duckdb_load_table(conn, "./Input_Files/kg-microbe-biomedical-function/merged-kg_edges.tsv", "edges", ["subject", "object"])
    conditions = [
        ('UniprotKB:%', 'Proteomes:%'),
        ('Proteomes:%', 'NCBITaxon:%'),
        ('UniprotKB:%', 'EC:%')
    ]

    if not os.path.exists("./Intermediate_Files"):
        os.makedirs("./Intermediate_Files")
    
    print("Relevant edges loaded.")

    query_with_multiple_conditions(
        conn,
        table_name = "edges",
        conditions = conditions)

    print("Subset to Proteomes-NCBITaxon pairs.")

    create_subject_object_pair_table(
        conn,
        table_name = "proteomes_ncbitaxon",
        base_table_name = "edges",
        subject = "Proteomes",
        object = "NCBITaxon",
        subject_prefix = "Proteomes:%",
        object_prefix = "NCBITaxon:%"
    )

    print("Subset to UniprotKB-Proteomes pairs.")

    create_subject_object_pair_table(
        conn,
        table_name = "uniprot_proteomes",
        base_table_name = "edges",
        subject = "UniprotKB",
        object = "Proteomes",
        subject_prefix = "UniprotKB:%",
        object_prefix = "Proteomes:%"
    )

    print("Subset to Uniprot-EC edges.")

    create_subject_object_pair_table(
        conn,
        table_name = "uniprot_ec",
        base_table_name = "edges",
        subject = "UniprotKB",
        object = "EC",
        subject_prefix = "UniprotKB:%",
        object_prefix = "EC:%"
    )

    print("Subset to UniprotKB-NCBITaxon edges.")

    join_tables_unique_subject_object(
         conn, 
        base_table_name = "proteomes_ncbitaxon", 
        compared_table_name = "uniprot_proteomes", 
        output_table_name = "uniprot_ncbitaxon", 
        output_subject = "UniprotKB", 
        output_object = "NCBITaxon", 
        comparison = "Proteomes"
    )

    print("Subset to UniprotKB-NCBITaxon edges.")

    join_tables_unique_subject_object(
         conn, 
        base_table_name = "uniprot_ec", 
        compared_table_name = "uniprot_ncbitaxon", 
        output_table_name = "ncbitaxon_ec", 
        output_subject = "NCBITaxon", 
        output_object = "EC", 
        comparison = "UniprotKB"
    )

    get_table_count(conn, "ncbitaxon_ec")

    output_table_to_file(conn, "ncbitaxon_ec", "./Intermediate_Files/NCBITaxon_to_EC.tsv")

    num_unique_ncbitaxon_ec_pairs = get_total_unique_pair(
        conn,
        table_name = "ncbitaxon_ec",
        subject = "NCBITaxon",
        object = "EC",
        subject_prefix = "NCBITaxon:%",
        object_prefix = "EC:%"
    )

    print("Total NCBITaxon-EC Pairs: ", num_unique_ncbitaxon_ec_pairs[0])

    num_unique_ncbitaxon = get_total_unique_column(
        conn,
        table_name = "ncbitaxon_ec",
        column_name = "NCBITaxon",
    )

    print("Total unique NCBITaxon:", num_unique_ncbitaxon[0])

if __name__ == '__main__':
    main()
