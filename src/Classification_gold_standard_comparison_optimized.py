"""
Memory-optimized version of Classification_gold_standard_comparison.py

Key optimizations:
1. Use DuckDB to filter edges BEFORE loading into pandas
2. Only load disease-relevant edges (not entire 30GB KG)
3. Process one disease at a time, clearing memory between iterations
4. Use DuckDB queries instead of pandas string operations where possible
"""

from collections import Counter, defaultdict
import duckdb
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gc

from scipy.stats import chi2_contingency
from duckdb_utils import duckdb_load_table
from classification_utils import remove_conflicting_directionality
from constants import COMPETENCY_DISEASE_MAP
from ncbi_phylogeny_search import find_microbes_strain, get_all_ranks

def get_disease_pairs_optimized(conn, output_dir, disease_id, disease_name):
    """
    Get disease-microbe associations using DuckDB queries (memory efficient).

    Instead of loading entire KG into pandas, we:
    1. Query DuckDB for disease-specific edges only
    2. Load small filtered result into pandas
    """
    print(f"Querying edges for disease: {disease_name} ({disease_id})")

    # Query DuckDB to get only disease-relevant edges with NCBITaxon subjects
    # This replaces the massive pandas filtering operations
    # Note: disease_id comes from trusted COMPETENCY_DISEASE_MAP constants
    query = f"""
    SELECT DISTINCT
        CASE
            WHEN subject = 'NCBITaxon:Actinomyces sp. ICM47' THEN 'NCBITaxon:936548'
            ELSE subject
        END as subject,
        predicate,
        object,
        CONCAT(predicate, '_', object) as predicate_object
    FROM edges
    WHERE subject LIKE 'NCBITaxon:%'
      AND (object = '{disease_id}'
           OR object = 'MONDO:0005011'
           OR object = 'MONDO:0005265'
           OR object = 'MONDO:0005101')
      AND predicate IS NOT NULL
      AND subject IS NOT NULL
      AND object IS NOT NULL
    """

    data_edges = conn.execute(query).df()
    print(f"  Found {len(data_edges)} disease-relevant edges")

    if len(data_edges) == 0:
        print(f"  WARNING: No edges found for {disease_name}")
        return pd.DataFrame(columns=['subject', 'object'])

    # Replace Crohn's and UC with IBD label
    data_edges['object'] = data_edges['object'].replace({
        'MONDO:0005011': 'MONDO:0005101',
        'MONDO:0005265': 'MONDO:0005101'
    })

    # Keep only relevant disease
    data_edges = data_edges[data_edges['object'] == disease_id]

    # Create subject-object pairs
    data_pairs = data_edges[['subject', 'predicate_object']].rename(
        columns={'predicate_object': 'object'}
    ).drop_duplicates()

    print(f"  Found {len(data_pairs)} unique subject-object pairs")

    # Remove conflicting directionality
    data_pairs_cleaned = remove_conflicting_directionality(data_pairs)
    print(f"  After removing conflicts: {len(data_pairs_cleaned)}")

    # Save results
    output_file = f"{output_dir}/outcome_to_NCBITaxon_cleaned_{disease_name}.tsv"
    data_pairs_cleaned.to_csv(output_file, sep="\t", header=True, index=False)

    return data_pairs_cleaned


def combine_labels(taxonomic_labels_counts, keep_label, new_label):
    try:
        taxonomic_labels_counts[keep_label] = taxonomic_labels_counts[keep_label] + taxonomic_labels_counts[new_label]
        del taxonomic_labels_counts[new_label]
    except KeyError:
        return taxonomic_labels_counts
    return taxonomic_labels_counts


def main():
    output_dir = "./data/Intermediate_Files"

    # Initialize DuckDB connection and load edges ONCE
    print("Initializing DuckDB and loading edges (this may take a few minutes)...")
    conn = duckdb.connect(":memory:")

    # Load only the full edges file into DuckDB (not pandas!)
    # DuckDB handles large files efficiently without loading entirely into RAM
    duckdb_load_table(
        conn,
        "./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges.tsv",
        "edges",
        ["subject", "predicate", "object"]
    )
    print("  DuckDB edges table loaded")

    # Load taxonomy edges for phylogeny searches
    conn_taxonomy = duckdb.connect(":memory:")
    duckdb_load_table(
        conn_taxonomy,
        "./data/Input_Files/kg-microbe-biomedical-function-cat/merged-kg_edges_ncbitaxon.tsv",
        "ncbitaxon_edges",
        ["subject", "predicate", "object"]
    )
    print("  DuckDB taxonomy edges loaded")

    # Load phylogeny ranks (small file)
    phylogeny_output_dir = "./data/Phylogeny_Search"
    ncbi_taxa_ranks_df = get_all_ranks(phylogeny_output_dir)
    print(f"  Loaded {len(ncbi_taxa_ranks_df)} taxonomy ranks")

    # Load butyrate producers (small file)
    butyrate_production_output_dir = "./data/Intermediate_Files_Competencies/butyrate_produces"
    gs_analysis_microbes_file = f"{butyrate_production_output_dir}/Gold_Standard_Species_Overlap_butyrate_produces.csv"
    gs_analysis_microbes_df = pd.read_csv(gs_analysis_microbes_file)
    kg_analysis_microbes = gs_analysis_microbes_df[
        (gs_analysis_microbes_df['Organismal'] == 1) |
        (gs_analysis_microbes_df['Functional'] == 1) |
        (gs_analysis_microbes_df['Functional_EC'] == 1)
    ]["Value"].tolist()
    print(f"  Loaded {len(kg_analysis_microbes)} butyrate producers")

    # Load taxa with proteomes (protein annotations in KG)
    # From Methods: "From the collection of protein annotations to these taxa, we then characterized them as being butyrate producers"
    proteomes_file = f"{phylogeny_output_dir}/unique_ncbitaxon_uniprot_ids.txt"
    proteomes_df = pd.read_csv(proteomes_file, header=None, names=['taxon'])
    ncbitaxa_with_proteomes = set(proteomes_df['taxon'].tolist())
    print(f"  Loaded {len(ncbitaxa_with_proteomes)} taxa with proteomes (UniprotKB annotations)")

    # Storage for plotting
    original_taxon_labels_counts = []
    original_taxon_labels = []
    representative_taxon_labels_counts = []
    representative_taxon_labels = []

    # Process each disease ONE AT A TIME (memory efficient)
    for disease_name, disease_ids in COMPETENCY_DISEASE_MAP.items():
        # Use the canonical disease ID (last in list) - this is what remains after replacement
        # e.g., IBD: Crohn's/UC are replaced with MONDO:0005101 (IBD), which is the last element
        disease_id = disease_ids[-1]
        print(f"\n{'='*60}")
        print(f"Processing disease: {disease_name}")
        print(f"{'='*60}")

        # Get disease-microbe pairs using DuckDB (memory efficient!)
        disease_microbes_df = get_disease_pairs_optimized(conn, output_dir, disease_id, disease_name)

        if len(disease_microbes_df) == 0:
            print(f"Skipping {disease_name} - no data")
            continue

        disease_microbes = disease_microbes_df["subject"].unique().tolist()
        print(f"Total disease-associated microbes: {len(disease_microbes)}")

        # Find overlap with butyrate producers
        gs_analysis_disease_overlap = set(kg_analysis_microbes) & set(disease_microbes)
        print(f"Overlap with butyrate producers: {len(gs_analysis_disease_overlap)}")

        # Build phylogeny relationships
        print("Building phylogeny relationships...")
        microbes_strain_dict, microbes_species_dict = find_microbes_strain(
            conn_taxonomy,
            ncbi_taxa_ranks_df,
            disease_microbes,
            output_dir,
            f"classification_butyrate_produces_{disease_name}"
        )

        # Analyze each microbe
        disease_microbes_ranks = []
        disease_microbes_disease_relationship = []
        disease_microbes_species_numbers = []
        disease_microbes_species_butyrate_producers = []
        disease_microbes_strains_numbers = []
        disease_microbes_strains_butyrate_producers = []

        # Create lookup dict for faster access
        disease_relationship_dict = dict(zip(
            disease_microbes_df["subject"],
            disease_microbes_df["object"]
        ))
        rank_dict = ncbi_taxa_ranks_df.set_index("NCBITaxon_ID")["Rank"].to_dict()

        print(f"Analyzing {len(disease_microbes)} microbes...")
        for microbe in disease_microbes:
            # Get relationship and rank
            ibd_relationship = disease_relationship_dict[microbe]
            disease_microbes_disease_relationship.append(ibd_relationship)

            microbe_rank = rank_dict.get(microbe, "not_found")
            disease_microbes_ranks.append(microbe_rank)

            # Get species
            if microbe_rank != "species":
                ibd_species = microbes_species_dict.get(microbe, [])
                if ibd_species == ['not found']:
                    ibd_species = []
            else:
                ibd_species = [microbe]

            # Filter species to only those with proteomes (protein annotations in KG)
            ibd_species_with_proteomes = [s for s in ibd_species if s in ncbitaxa_with_proteomes]

            species_with_trait = set(kg_analysis_microbes) & set(ibd_species_with_proteomes)
            disease_microbes_species_numbers.append(len(ibd_species_with_proteomes))
            disease_microbes_species_butyrate_producers.append(len(species_with_trait))

            # Get strains
            if microbe_rank not in ["strain", "subspecies"]:
                ibd_strains = microbes_strain_dict.get(microbe, [])
                if ibd_strains == ['not found']:
                    ibd_strains = []
            else:
                ibd_strains = [microbe]

            # Filter strains to only those with proteomes (protein annotations in KG)
            ibd_strains_with_proteomes = [s for s in ibd_strains if s in ncbitaxa_with_proteomes]

            strains_with_trait = set(kg_analysis_microbes) & set(ibd_strains_with_proteomes)
            disease_microbes_strains_numbers.append(len(ibd_strains_with_proteomes))
            disease_microbes_strains_butyrate_producers.append(len(strains_with_trait))

        # Create results dataframe
        ranks_df = pd.DataFrame({
            "Name": disease_microbes,
            "Disease_Relationship": disease_microbes_disease_relationship,
            "Rank": disease_microbes_ranks,
            "Num_Species": disease_microbes_species_numbers,
            "Num_Species_Butyrate_Producers": disease_microbes_species_butyrate_producers,
            "Num_Strains": disease_microbes_strains_numbers,
            "Num_Strains_Butyrate_Producers": disease_microbes_strains_butyrate_producers
        })

        ranks_df.to_csv(
            f"{output_dir}/outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_{disease_name}.csv",
            sep=",",
            index=False
        )

        # Calculate chi-square statistics
        num_increased_total = 0
        num_decreased_total = 0
        num_increased_disease_butyrate_producers = 0
        num_decreased_disease_butyrate_producers = 0
        total_strains = 0
        total_species = 0

        for _, row in ranks_df.iterrows():
            if row["Num_Strains"] > 0:
                num_total = row["Num_Strains"]
                num_producers = row["Num_Strains_Butyrate_Producers"]
                total_strains += num_total
            else:
                num_total = row["Num_Species"]
                num_producers = row["Num_Species_Butyrate_Producers"]
                total_species += num_total

            if "increased" in row["Disease_Relationship"]:
                num_increased_total += num_total
                num_increased_disease_butyrate_producers += num_producers
            elif "decreased" in row["Disease_Relationship"]:
                num_decreased_total += num_total
                num_decreased_disease_butyrate_producers += num_producers

        # Chi-square test
        contingency_table = [
            [num_increased_disease_butyrate_producers, num_increased_total - num_increased_disease_butyrate_producers],
            [num_decreased_disease_butyrate_producers, num_decreased_total - num_decreased_disease_butyrate_producers]
        ]
        chi2, pval, dof, expected = chi2_contingency(contingency_table)

        print(f"\nChi-square results for {disease_name}:")
        print(f"  Chi2: {chi2:.4f}")
        print(f"  P-value: {pval:.4e}")
        print(f"  Increased total: {num_increased_total}, producers: {num_increased_disease_butyrate_producers}")
        print(f"  Decreased total: {num_decreased_total}, producers: {num_decreased_disease_butyrate_producers}")

        # Save summary
        classification_summary_df = pd.DataFrame({
            "Num_Decreased_disease_Total": [num_decreased_total],
            "Num_Increased_disease_Total": [num_increased_total],
            "Num_Decreased_disease_Butyrate_Producers": [num_decreased_disease_butyrate_producers],
            "Num_Increased_disease_Butyrate_Producers": [num_increased_disease_butyrate_producers],
            "chi2": [chi2],
            "P_Val": [pval]
        })

        classification_summary_file = f"{output_dir}/{disease_name}_Classification_butyrate_producers_summary.csv"
        classification_summary_df.to_csv(classification_summary_file, sep=",", index=False)

        # Track taxonomy labels for plotting
        taxonomic_labels_counts = Counter(disease_microbes_ranks)
        taxonomic_labels_counts = combine_labels(taxonomic_labels_counts, "strain", "subspecies")
        taxonomic_labels_counts = combine_labels(taxonomic_labels_counts, "genus", "NCBITaxon#_species_group")
        taxonomic_labels = list(taxonomic_labels_counts.keys())
        original_taxon_labels_counts.append(taxonomic_labels_counts)
        original_taxon_labels.append(taxonomic_labels)

        new_labels_counts = {'species': total_species, 'strain': total_strains}
        representative_taxon_labels_counts.append(new_labels_counts)
        representative_taxon_labels.append(list(new_labels_counts.keys()))

        # Clean up memory after each disease
        del disease_microbes_df, disease_microbes, ranks_df
        del microbes_strain_dict, microbes_species_dict
        gc.collect()
        print(f"Completed {disease_name}")

    # Create pie charts
    print("\nCreating visualization...")
    combined_original_taxon_labels = list(set([item for sublist in original_taxon_labels for item in sublist]))
    order_list = ['phylum', 'class', 'order', 'family', 'genus', 'species', 'strain', 'not_found']
    combined_original_taxon_labels = sorted(
        combined_original_taxon_labels,
        key=lambda x: (x not in order_list, order_list.index(x) if x in order_list else float('inf'))
    )

    colors = plt.cm.Paired(range(len(combined_original_taxon_labels)))
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))

    for i, disease_name in enumerate(COMPETENCY_DISEASE_MAP.keys()):
        # Original taxa
        wedges, texts = ax[0, i].pie(
            [original_taxon_labels_counts[i].get(label, 0) for label in combined_original_taxon_labels],
            labels=None,
            autopct=None,
            startangle=90,
            colors=colors,
            wedgeprops={'edgecolor': 'black'}
        )
        ax[0, i].axis('equal')
        ax[0, i].set_title(
            f"{disease_name}: {sum(original_taxon_labels_counts[i].values())} total",
            fontsize=18
        )

        # Representative taxa
        wedges, texts = ax[1, i].pie(
            [representative_taxon_labels_counts[i].get(label, 0) for label in combined_original_taxon_labels],
            labels=None,
            autopct=None,
            startangle=90,
            colors=colors,
            wedgeprops={'edgecolor': 'black'}
        )
        ax[1, i].axis('equal')
        ax[1, i].set_title(
            f"{disease_name}: {sum(representative_taxon_labels_counts[i].values())} total",
            fontsize=18
        )

    fig.legend(wedges, combined_original_taxon_labels, title="Rank", loc="upper left", fontsize=12, bbox_to_anchor=(0.8, 0.8))
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    fig.savefig(f"{output_dir}/All_Disease_Comparison_Ranks.png")
    plt.close()

    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

    # Close DuckDB connections
    conn.close()
    conn_taxonomy.close()


if __name__ == '__main__':
    main()
