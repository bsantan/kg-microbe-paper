#!/usr/bin/env python3
"""
Figure 6C Competency Queries
=============================

Helper functions to answer competency questions for the comprehensive
validation report.
"""

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from pathlib import Path
from collections import Counter

# Base directories
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
INTERMEDIATE_DIR = DATA_DIR / "Intermediate_Files"
COMPETENCIES_DIR = DATA_DIR / "Intermediate_Files_Competencies" / "butyrate_produces"


def query_butyrate_producers():
    """
    Competency Q1: Total butyrate producers identified

    Returns dict with counts by evidence source.
    """
    gs_file = COMPETENCIES_DIR / "Gold_Standard_Species_Overlap_butyrate_produces.csv"
    df = pd.read_csv(gs_file)

    # Count by evidence type
    organismal_count = (df['Organismal'] == 1).sum()
    functional_count = (df['Functional'] == 1).sum()
    functional_ec_count = (df['Functional_EC'] == 1).sum()

    # Count unique producers (any evidence)
    any_evidence = ((df['Organismal'] == 1) |
                   (df['Functional'] == 1) |
                   (df['Functional_EC'] == 1)).sum()

    # Count overlap
    org_and_func = ((df['Organismal'] == 1) & (df['Functional'] == 1)).sum()
    func_and_ec = ((df['Functional'] == 1) & (df['Functional_EC'] == 1)).sum()
    org_and_ec = ((df['Organismal'] == 1) & (df['Functional_EC'] == 1)).sum()
    all_three = ((df['Organismal'] == 1) & (df['Functional'] == 1) &
                 (df['Functional_EC'] == 1)).sum()

    return {
        'total_entries': len(df),
        'organismal': organismal_count,
        'functional': functional_count,
        'functional_ec': functional_ec_count,
        'any_evidence': any_evidence,
        'org_and_func': org_and_func,
        'func_and_ec': func_and_ec,
        'org_and_ec': org_and_ec,
        'all_three': all_three
    }


def query_disease_taxa(disease='IBD'):
    """
    Competency Q2: Disease-associated taxa counts

    Args:
        disease: 'IBD' or 'PD'

    Returns dict with taxa counts by direction.
    """
    ranks_file = INTERMEDIATE_DIR / f"outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_{disease}.csv"
    df = pd.read_csv(ranks_file)

    # Count by direction
    increased = (df['Disease_Relationship'].str.contains('increased')).sum()
    decreased = (df['Disease_Relationship'].str.contains('decreased')).sum()

    return {
        'disease': disease,
        'total_taxa': len(df),
        'increased': increased,
        'decreased': decreased
    }


def query_rank_distribution(disease='IBD'):
    """
    Competency Q3: Taxonomic rank distribution

    Args:
        disease: 'IBD' or 'PD'

    Returns dict with rank counts and percentages.
    """
    ranks_file = INTERMEDIATE_DIR / f"outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_{disease}.csv"
    df = pd.read_csv(ranks_file)

    rank_counts = df['Rank'].value_counts().to_dict()
    total = len(df)

    rank_percentages = {rank: (count/total)*100 for rank, count in rank_counts.items()}

    return {
        'counts': rank_counts,
        'percentages': rank_percentages,
        'total': total
    }


def trace_taxa_expansion(disease='IBD', n=5):
    """
    Competency Q4: Worked expansion examples

    Args:
        disease: 'IBD' or 'PD'
        n: Number of examples to return

    Returns list of dicts with expansion details.
    """
    ranks_file = INTERMEDIATE_DIR / f"outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_{disease}.csv"
    df = pd.read_csv(ranks_file)

    # Select diverse examples by rank
    examples = []

    # Get one of each major rank type if available
    for rank in ['species', 'genus', 'family', 'order', 'class']:
        rank_df = df[df['Rank'] == rank]
        if len(rank_df) > 0:
            # Pick an interesting example (has producers, moderate expansion)
            producer_examples = rank_df[
                (rank_df['Num_Species_Butyrate_Producers'] > 0) |
                (rank_df['Num_Strains_Butyrate_Producers'] > 0)
            ]

            if len(producer_examples) > 0:
                example = producer_examples.iloc[0]
            else:
                example = rank_df.iloc[0]

            # Determine which count is used (strain priority)
            if example['Num_Strains'] > 0:
                used_total = example['Num_Strains']
                used_producers = example['Num_Strains_Butyrate_Producers']
                count_type = 'strains'
            else:
                used_total = example['Num_Species']
                used_producers = example['Num_Species_Butyrate_Producers']
                count_type = 'species'

            examples.append({
                'taxon_id': example['Name'],
                'rank': example['Rank'],
                'relationship': example['Disease_Relationship'],
                'num_species': example['Num_Species'],
                'species_producers': example['Num_Species_Butyrate_Producers'],
                'num_strains': example['Num_Strains'],
                'strain_producers': example['Num_Strains_Butyrate_Producers'],
                'used_count_type': count_type,
                'used_total': used_total,
                'used_producers': used_producers,
                'direction': 'decreased' if 'decreased' in example['Disease_Relationship'] else 'increased'
            })

            if len(examples) >= n:
                break

    return examples


def explain_count_differences():
    """
    Competency Q5: Taxa count discrepancy explanation

    Returns dict comparing current vs published counts.
    """
    # Current counts from summary files
    ibd_summary = pd.read_csv(INTERMEDIATE_DIR / "IBD_Classification_butyrate_producers_summary.csv")
    pd_summary = pd.read_csv(INTERMEDIATE_DIR / "PD_Classification_butyrate_producers_summary.csv")

    current = {
        'IBD': {
            'decreased_total': int(ibd_summary['Num_Decreased_disease_Total'].iloc[0]),
            'decreased_producers': int(ibd_summary['Num_Decreased_disease_Butyrate_Producers'].iloc[0]),
            'increased_total': int(ibd_summary['Num_Increased_disease_Total'].iloc[0]),
            'increased_producers': int(ibd_summary['Num_Increased_disease_Butyrate_Producers'].iloc[0])
        },
        'PD': {
            'decreased_total': int(pd_summary['Num_Decreased_disease_Total'].iloc[0]),
            'decreased_producers': int(pd_summary['Num_Decreased_disease_Butyrate_Producers'].iloc[0]),
            'increased_total': int(pd_summary['Num_Increased_disease_Total'].iloc[0]),
            'increased_producers': int(pd_summary['Num_Increased_disease_Butyrate_Producers'].iloc[0])
        }
    }

    # Published counts.
    # Reference: revisions2/KG-Microbe_Responses2_mpj2.docx ("Summary of the
    # second computational reproducibility review") confirms paper chi²/p-values
    # for PD as 1317 / 2e-288 and notes IBD was "consistent" between paper and
    # reviewer reproduction. The per-cell totals/butyrate counts below are not
    # restated in that responses doc; they are from Figure 6 panels (transcribed
    # in data/NERSC_vs_new_comparison_reports/figure6_publication_results.txt).
    published = {
        'IBD': {
            'increased_total': 15_398,
            'increased_producers': 514,
            'decreased_total': 39_554,
            'decreased_producers': 221
        },
        'PD': {
            'increased_total': 2_990,
            'increased_producers': 299,
            'decreased_total': 22_531,
            'decreased_producers': 152
        }
    }

    # Calculate differences
    differences = {}
    for disease in ['IBD', 'PD']:
        curr = current[disease]
        pub = published[disease]

        curr_total = curr['decreased_total'] + curr['increased_total']
        curr_prod = curr['decreased_producers'] + curr['increased_producers']
        pub_total = pub['decreased_total'] + pub['increased_total']
        pub_prod = pub['decreased_producers'] + pub['increased_producers']

        differences[disease] = {
            'current_total': curr_total,
            'published_total': pub_total,
            'total_diff_pct': ((curr_total - pub_total) / pub_total) * 100,
            'current_producers': curr_prod,
            'published_producers': pub_prod,
            'producers_diff_pct': ((curr_prod - pub_prod) / pub_prod) * 100
        }

    return {
        'current': current,
        'published': published,
        'differences': differences
    }


def recalculate_chi_square():
    """
    Competency Q6: Recalculate chi-square independently

    Returns dict with calculated vs reported values.
    """
    results = {}

    for disease in ['IBD', 'PD']:
        summary_file = INTERMEDIATE_DIR / f"{disease}_Classification_butyrate_producers_summary.csv"
        df = pd.read_csv(summary_file)

        # Extract values
        dec_total = int(df['Num_Decreased_disease_Total'].iloc[0])
        dec_prod = int(df['Num_Decreased_disease_Butyrate_Producers'].iloc[0])
        inc_total = int(df['Num_Increased_disease_Total'].iloc[0])
        inc_prod = int(df['Num_Increased_disease_Butyrate_Producers'].iloc[0])

        reported_chi2 = float(df['chi2'].iloc[0])
        reported_pval = float(df['P_Val'].iloc[0])

        # Calculate contingency table
        contingency = [
            [inc_prod, inc_total - inc_prod],
            [dec_prod, dec_total - dec_prod]
        ]

        # Calculate chi-square
        chi2, pval, dof, expected = chi2_contingency(contingency)

        # Calculate proportions and effect sizes
        inc_prop = inc_prod / inc_total
        dec_prop = dec_prod / dec_total
        proportion_ratio = dec_prop / inc_prop

        # Calculate odds ratio
        odds_ratio = (dec_prod * (inc_total - inc_prod)) / (inc_prod * (dec_total - dec_prod))

        results[disease] = {
            'contingency': contingency,
            'reported_chi2': reported_chi2,
            'calculated_chi2': chi2,
            'chi2_match': abs(chi2 - reported_chi2) < 1.0,
            'reported_pval': reported_pval,
            'calculated_pval': pval,
            'pval_match': abs(np.log10(pval) - np.log10(reported_pval)) < 1.0,
            'increased_proportion': inc_prop,
            'decreased_proportion': dec_prop,
            'proportion_ratio': proportion_ratio,
            'odds_ratio': odds_ratio
        }

    return results


def validate_directionality(disease='IBD'):
    """
    Competency Q8-Q9: Manual taxa tracing and literature cross-validation

    Returns dict with validation results.
    """
    ranks_file = INTERMEDIATE_DIR / f"outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_{disease}.csv"
    df = pd.read_csv(ranks_file)

    # Known protective butyrate producers (should be in decreased group)
    known_protective = {
        'NCBITaxon:853': 'Faecalibacterium prausnitzii',
        'NCBITaxon:172901': 'Roseburia hominis',
        'NCBITaxon:39491': 'Eubacterium rectale',
        'NCBITaxon:815': 'Bacteroidaceae family'
    }

    # Known pathogenic taxa (should be in increased group)
    known_pathogenic = {
        'NCBITaxon:562': 'Escherichia coli',
        'NCBITaxon:573': 'Klebsiella pneumoniae'
    }

    protective_validation = []
    for taxon_id, name in known_protective.items():
        matches = df[df['Name'] == taxon_id]
        if len(matches) > 0:
            relationship = matches.iloc[0]['Disease_Relationship']
            is_decreased = 'decreased' in relationship
            protective_validation.append({
                'taxon_id': taxon_id,
                'name': name,
                'relationship': relationship,
                'in_decreased_group': is_decreased,
                'expected': True,
                'validated': is_decreased
            })

    pathogenic_validation = []
    for taxon_id, name in known_pathogenic.items():
        matches = df[df['Name'] == taxon_id]
        if len(matches) > 0:
            relationship = matches.iloc[0]['Disease_Relationship']
            is_increased = 'increased' in relationship
            pathogenic_validation.append({
                'taxon_id': taxon_id,
                'name': name,
                'relationship': relationship,
                'in_increased_group': is_increased,
                'expected': True,
                'validated': is_increased
            })

    # Check for ambiguous relationships
    ambiguous = []
    for idx, row in df.iterrows():
        rel = row['Disease_Relationship'].lower()
        if 'increased' in rel and 'decreased' in rel:
            ambiguous.append(row['Name'])

    return {
        'protective_taxa': protective_validation,
        'pathogenic_taxa': pathogenic_validation,
        'ambiguous_relationships': ambiguous,
        'all_protective_validated': all(v['validated'] for v in protective_validation),
        'all_pathogenic_validated': all(v['validated'] for v in pathogenic_validation)
    }


def run_counterfactual():
    """
    Competency Q10: Label swap analysis

    Returns dict comparing current vs counterfactual interpretation.
    """
    results = {}

    for disease in ['IBD', 'PD']:
        summary_file = INTERMEDIATE_DIR / f"{disease}_Classification_butyrate_producers_summary.csv"
        df = pd.read_csv(summary_file)

        # Current values
        dec_total = int(df['Num_Decreased_disease_Total'].iloc[0])
        dec_prod = int(df['Num_Decreased_disease_Butyrate_Producers'].iloc[0])
        inc_total = int(df['Num_Increased_disease_Total'].iloc[0])
        inc_prod = int(df['Num_Increased_disease_Butyrate_Producers'].iloc[0])
        chi2 = float(df['chi2'].iloc[0])
        pval = float(df['P_Val'].iloc[0])

        # Current interpretation
        current_inc_prop = inc_prod / inc_total
        current_dec_prop = dec_prod / dec_total
        current_ratio = current_dec_prop / current_inc_prop
        current_interpretation = "Protective" if current_ratio > 1 else "Risk Factor"

        # Counterfactual (swapped labels)
        counter_inc_prop = dec_prod / dec_total  # Swapped
        counter_dec_prop = inc_prod / inc_total  # Swapped
        counter_ratio = counter_dec_prop / counter_inc_prop
        counter_interpretation = "Protective" if counter_ratio > 1 else "Risk Factor"

        results[disease] = {
            'current': {
                'increased_total': inc_total,
                'increased_producers': inc_prod,
                'increased_proportion': current_inc_prop,
                'decreased_total': dec_total,
                'decreased_producers': dec_prod,
                'decreased_proportion': current_dec_prop,
                'ratio': current_ratio,
                'interpretation': current_interpretation,
                'chi2': chi2,
                'pval': pval
            },
            'counterfactual': {
                'increased_total': dec_total,  # Swapped
                'increased_producers': dec_prod,  # Swapped
                'increased_proportion': counter_inc_prop,
                'decreased_total': inc_total,  # Swapped
                'decreased_producers': inc_prod,  # Swapped
                'decreased_proportion': counter_dec_prop,
                'ratio': counter_ratio,
                'interpretation': counter_interpretation,
                'chi2': chi2,  # Same!
                'pval': pval  # Same!
            }
        }

    return results


if __name__ == '__main__':
    # Test all queries
    print("Testing competency queries...")

    print("\nQ1: Butyrate Producers")
    print(query_butyrate_producers())

    print("\nQ2: Disease Taxa (IBD)")
    print(query_disease_taxa('IBD'))

    print("\nQ3: Rank Distribution (IBD)")
    print(query_rank_distribution('IBD'))

    print("\nQ4: Taxa Expansion Examples (IBD)")
    examples = trace_taxa_expansion('IBD', n=3)
    for ex in examples:
        print(f"  {ex['taxon_id']} ({ex['rank']}): {ex['used_total']} {ex['used_count_type']}, {ex['used_producers']} producers")

    print("\nQ5: Count Differences")
    diffs = explain_count_differences()
    print(f"  IBD Current Total: {diffs['differences']['IBD']['current_total']}")
    print(f"  IBD Published Total: {diffs['differences']['IBD']['published_total']}")

    print("\nQ6: Chi-Square Recalculation")
    chi2_results = recalculate_chi_square()
    for disease in ['IBD', 'PD']:
        print(f"  {disease}: Calculated={chi2_results[disease]['calculated_chi2']:.2f}, Reported={chi2_results[disease]['reported_chi2']:.2f}, Match={chi2_results[disease]['chi2_match']}")

    print("\nQ8-Q9: Directionality Validation (IBD)")
    direction_results = validate_directionality('IBD')
    print(f"  Protective taxa validated: {direction_results['all_protective_validated']}")
    print(f"  Pathogenic taxa validated: {direction_results['all_pathogenic_validated']}")
    print(f"  Ambiguous relationships: {len(direction_results['ambiguous_relationships'])}")

    print("\nQ10: Counterfactual Analysis")
    counter_results = run_counterfactual()
    for disease in ['IBD', 'PD']:
        print(f"  {disease} Current: {counter_results[disease]['current']['interpretation']} ({counter_results[disease]['current']['ratio']:.2f}x)")
        print(f"  {disease} Counterfactual: {counter_results[disease]['counterfactual']['interpretation']} ({counter_results[disease]['counterfactual']['ratio']:.2f}x)")

    print("\nAll competency queries completed successfully!")
