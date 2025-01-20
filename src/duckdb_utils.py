



def get_table_count(con, table):
    """Get the number of rows of a given duckdb table name."""
    # Execute the SQL query to count the rows
    result = con.execute(
        f"""
        SELECT COUNT(*) FROM {table};
    """
    ).fetchone()

    # Print the number of rows
    print(f"Number of rows in '{table}': {result[0]}")

def duckdb_load_table(con, file, table_name, columns):
    """Create a duckDB tables for any given graph."""
    columns_str = ", ".join(columns)

    # Read the subset file into a DuckDB table
    query = (
        f"""
    CREATE OR REPLACE TABLE {table_name} AS
    SELECT {columns_str}
    FROM read_csv_auto('{file}', delim='\t');
    """
    )

    print(query)
    con.execute(query)

def query_with_multiple_conditions(con, table_name, conditions):
    """
    Execute a query with multiple conditions connected by OR statements.

    Parameters:
    con (duckdb.DuckDBPyConnection): The DuckDB connection.
    table_name (str): The name of the table to query.
    conditions (list of tuples): Each tuple contains two elements representing the subject and object conditions.

    Returns:
    list: The query results.
    """

    # Build the WHERE clause
    where_clauses = []
    for s, o in conditions:
        where_clauses.append(f"(subject LIKE '{s}' AND object LIKE '{o}')")
    
    where_clause = " OR ".join(where_clauses)
    
    # Execute the query
    query = f"SELECT * FROM {table_name} WHERE {where_clause};"

    print(query)
    con.execute(query)

def create_subject_object_pair_table(con, table_name, base_table_name, subject, object, subject_prefix, object_prefix):

    query = (
        f"""
        CREATE TEMPORARY TABLE {table_name} AS
        SELECT subject AS {subject}, object AS {object}
        FROM {base_table_name}
        WHERE subject LIKE '{subject_prefix}' AND object LIKE '{object_prefix}';
        """
    )

    print(query)
    con.execute(query)

def create_subject_table(con, table_name, base_table_name, subject, subject_prefix):

    query = (
        f"""
        CREATE TEMPORARY TABLE {table_name} AS
        SELECT subject AS {subject}
        FROM {base_table_name}
        WHERE subject LIKE '{subject_prefix}';
        """
    )

    print(query)
    con.execute(query)

def join_tables_subject_object(con, base_table_name, compared_table_name, output_table_name, output_subject, output_object, comparison):

    query = (
        f"""
        CREATE TEMPORARY TABLE {output_table_name} AS
        SELECT {compared_table_name}.{output_subject}, {compared_table_name}.{comparison}, {base_table_name}.{output_object}
        FROM {base_table_name}
        JOIN {compared_table_name} ON {base_table_name}.{comparison} = {compared_table_name}.{comparison};

        DROP TABLE {compared_table_name};
        """
    )

    print(query)
    con.execute(query)

def join_tables_subject(con, base_table_name, compared_table_name, output_table_name, comparison):

    query = (
        f"""
        CREATE TEMPORARY TABLE {output_table_name} AS
        SELECT subject, object
        FROM {base_table_name}
        JOIN {compared_table_name} ON {base_table_name}.subject = {compared_table_name}.{comparison};

        DROP TABLE {compared_table_name};
        """
    )

    print(query)
    con.execute(query)

def join_tables_unique_subject_object(con, base_table_name, compared_table_name, output_table_name, output_subject, output_object, comparison):

    query = (
        f"""
        CREATE TEMPORARY TABLE {output_table_name} AS
        SELECT DISTINCT {compared_table_name}.{output_subject}, {base_table_name}.{output_object}
        FROM {base_table_name}
        JOIN {compared_table_name} ON {base_table_name}.{comparison} = {compared_table_name}.{comparison};

        DROP TABLE {compared_table_name};
        """
    )

    print(query)
    con.execute(query)

def get_total_unique_pair(con, table_name, subject, object, subject_prefix, object_prefix):

    query = (
        f"""
        SELECT COUNT(*) AS unique_pairs
        FROM (
            SELECT DISTINCT {subject}, {object}
            FROM {table_name}
            WHERE {subject} LIKE '{subject_prefix}' AND {object} LIKE '{object_prefix}'
        ) AS distinct_pairs;
        """
    )

    print(query)
    result = con.execute(query).fetchone()

    return result

def get_total_unique_column(con, table_name, column_name):

    query = (
        f"""
        SELECT COUNT(DISTINCT {column_name}) AS unique_count
        FROM {table_name};
        """
    )

    print(query)
    result = con.execute(query).fetchone()

    return result

def output_table_to_file(con, table_name, filename):

    query = (
        f"""
        COPY {table_name} TO '{filename}' (FORMAT CSV, DELIMITER '\t', HEADER);
        """
    ) 

    print(query)
    con.execute(query)

def get_node_label(con, node_id):

    query = (
        f"""
        SELECT name 
        FROM nodes 
        WHERE id = '{node_id}';
        """
    ) 

    # print(query)
    try:
        result = con.execute(query).fetchone()[0]
    except TypeError:
        return node_id

    return result
