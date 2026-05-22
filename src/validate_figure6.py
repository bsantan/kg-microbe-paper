#!/usr/bin/env python3
"""
Figure 6 Validation Script
===========================

This script validates the statistical results for Figure 6 by:
1. Loading the classification data from JSON files
2. Loading the gold standard butyrate producer annotations
3. Recalculating aggregations independently
4. Reconstructing chi-square tests
5. Comparing with published and current results
6. Auditing directionality preservation

This addresses the user's concern about counterintuitive patterns where
butyrate producers appear in BOTH increased and decreased disease groups.
"""

import json
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from pathlib import Path
from collections import defaultdict, Counter

# Base directories
BASE_DIR = Path(__file__).parent
INTERMEDIATE_DIR = BASE_DIR / "Intermediate_Files"
COMPETENCIES_DIR = BASE_DIR / "Intermediate_Files_Competencies" / "butyrate_produces"
DATA_DIR = BASE_DIR.parent / "data"

# Published results (from bioRxiv paper Figure 6)
PUBLISHED_RESULTS = {
    'IBD': {
        'increased_total': 15_398,
        'increased_producers': 514,
        'decreased_total': 39_554,
        'decreased_producers': 221,
        'chi2': 647,
        'pval': 1e-142
    },
    'PD': {
        'increased_total': 2_990,
        'increased_producers': 299,
        'decreased_total': 22_531,
        'decreased_producers': 152,
        'chi2': 1317,
        'pval': 2e-288
    }
}

# Current repository results (from comparison report)
CURRENT_RESULTS = {
    'IBD': {
        'decreased_total': 3_036,
        'increased_total': 7_764,
        'decreased_producers': 478,
        'increased_producers': 186,
        'chi2': 671.68,
        'pval': 4.30e-148
    },
    'PD': {
        'decreased_total': 1_240,
        'increased_total': 7_401,
        'decreased_producers': 302,
        'increased_producers': 150,
        'chi2': 1063.60,
        'pval': 2.70e-233
    }
}

def load_gold_standard():
    """Load gold standard butyrate producer annotations."""
    print("="*80)
    print("LOADING GOLD STANDARD BUTYRATE PRODUCERS")
    print("="*80)

    gs_file = DATA_DIR / "Gold_Standard_Species_Overlap_butyrate_produces.csv"
    if not gs_file.exists():
        # Try alternative location
        gs_file = COMPETENCIES_DIR / "Gold_Standard_Species_Overlap_butyrate_produces.csv"

    if not gs_file.exists():
        print(f"❌ ERROR: Gold standard file not found at {gs_file}")
        return set()

    gs_df = pd.read_csv(gs_file)
    print(f"✓ Loaded gold standard file: {len(gs_df)} total entries")
    print(f"  Columns: {list(gs_df.columns)}")

    # Define producers as having ANY of the three evidence types
    # Organismal = 1, Functional = 1, or Functional_EC = 1
    if all(col in gs_df.columns for col in ['Organismal', 'Functional', 'Functional_EC', 'Value']):
        producers = gs_df[
            (gs_df['Organismal'] == 1) |
            (gs_df['Functional'] == 1) |
            (gs_df['Functional_EC'] == 1)
        ]['Value'].tolist()
    elif 'Gold_Standard' in gs_df.columns and 'Value' in gs_df.columns:
        producers = gs_df[gs_df['Gold_Standard'] == 1]['Value'].tolist()
    else:
        print(f"❌ ERROR: Cannot identify producer columns in gold standard file")
        print(f"  Available columns: {list(gs_df.columns)}")
        return set()

    producer_set = set(producers)
    print(f"✓ Identified {len(producer_set)} unique butyrate producers")
    print(f"  Example producers: {list(producer_set)[:5]}")

    return producer_set

def load_classification_json(disease):
    """Load classification JSON files for a disease."""
    print(f"\n{'='*80}")
    print(f"LOADING CLASSIFICATION DATA FOR {disease}")
    print(f"{'='*80}")

    species_file = INTERMEDIATE_DIR / f"classification_butyrate_produces_{disease}_microbes_species.json"
    strain_file = INTERMEDIATE_DIR / f"classification_butyrate_produces_{disease}_microbes_strain.json"

    data = {'species': None, 'strain': None}

    if species_file.exists():
        with open(species_file, 'r') as f:
            data['species'] = json.load(f)
        print(f"✓ Loaded species data: {len(data['species'])} entries")
    else:
        print(f"⚠ Species file not found: {species_file}")

    if strain_file.exists():
        with open(strain_file, 'r') as f:
            data['strain'] = json.load(f)
        print(f"✓ Loaded strain data: {len(data['strain'])} entries")
    else:
        print(f"⚠ Strain file not found: {strain_file}")

    return data

def analyze_directionality(classification_data, disease):
    """Analyze the directionality patterns in classification data."""
    print(f"\n{'='*80}")
    print(f"ANALYZING DIRECTIONALITY FOR {disease}")
    print(f"{'='*80}")

    # Combine species and strain data
    all_entries = []

    if classification_data['species']:
        for microbe, info in classification_data['species'].items():
            all_entries.append({
                'microbe': microbe,
                'type': 'species',
                'relationship': info.get('Disease_Relationship', ''),
                'data': info
            })

    if classification_data['strain']:
        for microbe, info in classification_data['strain'].items():
            all_entries.append({
                'microbe': microbe,
                'type': 'strain',
                'relationship': info.get('Disease_Relationship', ''),
                'data': info
            })

    # Count direction patterns
    direction_counts = Counter()
    for entry in all_entries:
        rel = entry['relationship'].lower()
        if 'increased' in rel:
            direction_counts['increased'] += 1
        elif 'decreased' in rel:
            direction_counts['decreased'] += 1
        else:
            direction_counts['other'] += 1

    print(f"Direction distribution:")
    print(f"  Increased likelihood: {direction_counts['increased']}")
    print(f"  Decreased likelihood: {direction_counts['decreased']}")
    print(f"  Other/Unknown: {direction_counts['other']}")

    # Sample some entries for manual inspection
    print(f"\nSample entries with 'increased' relationship:")
    increased_samples = [e for e in all_entries if 'increased' in e['relationship'].lower()][:5]
    for i, entry in enumerate(increased_samples, 1):
        print(f"  {i}. {entry['microbe']} ({entry['type']}): {entry['relationship']}")

    print(f"\nSample entries with 'decreased' relationship:")
    decreased_samples = [e for e in all_entries if 'decreased' in e['relationship'].lower()][:5]
    for i, entry in enumerate(decreased_samples, 1):
        print(f"  {i}. {entry['microbe']} ({entry['type']}): {entry['relationship']}")

    return all_entries

def recalculate_aggregations(classification_data, producer_set, disease):
    """
    Independently recalculate aggregations from classification data.

    This mimics the logic in Classification_gold_standard_comparison.py lines 196-210.
    """
    print(f"\n{'='*80}")
    print(f"RECALCULATING AGGREGATIONS FOR {disease}")
    print(f"{'='*80}")

    # Counters for aggregation
    increased_total = 0
    increased_producers = 0
    decreased_total = 0
    decreased_producers = 0

    # Track what we're counting
    counted_microbes = []

    # Process species and strain data
    for level in ['species', 'strain']:
        if not classification_data[level]:
            continue

        for microbe, info in classification_data[level].items():
            relationship = info.get('Disease_Relationship', '')

            # For this microbe, how many species/strains does it represent?
            # This matches the Classification script logic
            num_species = info.get('Num_Species', 0)
            num_strains = info.get('Num_Strains', 0)

            # Use strains if available, otherwise species
            if num_strains > 0:
                num_total = num_strains
                num_prod = info.get('Num_Strains_Butyrate_Producers', 0)
            else:
                num_total = num_species
                num_prod = info.get('Num_Species_Butyrate_Producers', 0)

            # Skip if no count
            if num_total == 0:
                continue

            # Classify by direction using substring matching (as in original code)
            is_producer = microbe in producer_set

            if "increased" in relationship:
                increased_total += num_total
                increased_producers += num_prod
                counted_microbes.append({
                    'microbe': microbe,
                    'direction': 'increased',
                    'total': num_total,
                    'producers': num_prod,
                    'is_producer': is_producer,
                    'relationship': relationship
                })
            elif "decreased" in relationship:
                decreased_total += num_total
                decreased_producers += num_prod
                counted_microbes.append({
                    'microbe': microbe,
                    'direction': 'decreased',
                    'total': num_total,
                    'producers': num_prod,
                    'is_producer': is_producer,
                    'relationship': relationship
                })

    results = {
        'increased_total': increased_total,
        'increased_producers': increased_producers,
        'decreased_total': decreased_total,
        'decreased_producers': decreased_producers,
        'counted_microbes': counted_microbes
    }

    print(f"\nRecalculated counts:")
    print(f"  Increased disease - Total: {increased_total}, Producers: {increased_producers}")
    print(f"  Decreased disease - Total: {decreased_total}, Producers: {decreased_producers}")
    print(f"  Grand total: {increased_total + decreased_total}")
    print(f"  Unique microbes counted: {len(counted_microbes)}")

    # Calculate proportions
    if increased_total > 0:
        increased_prop = increased_producers / increased_total
        print(f"  Proportion producers in increased: {increased_prop:.4f} ({increased_prop*100:.2f}%)")

    if decreased_total > 0:
        decreased_prop = decreased_producers / decreased_total
        print(f"  Proportion producers in decreased: {decreased_prop:.4f} ({decreased_prop*100:.2f}%)")

    return results

def compute_chi_square(counts):
    """Compute chi-square test from counts."""
    print(f"\n{'='*80}")
    print(f"COMPUTING CHI-SQUARE TEST")
    print(f"{'='*80}")

    # Construct contingency table (matching Classification script lines 213-216)
    contingency = [
        [counts['increased_producers'], counts['increased_total'] - counts['increased_producers']],
        [counts['decreased_producers'], counts['decreased_total'] - counts['decreased_producers']]
    ]

    print(f"Contingency table:")
    print(f"                      | Producers | Non-Producers |")
    print(f"  Increased in Disease| {contingency[0][0]:9d} | {contingency[0][1]:13d} |")
    print(f"  Decreased in Disease| {contingency[1][0]:9d} | {contingency[1][1]:13d} |")

    chi2, pval, dof, expected = chi2_contingency(contingency)

    print(f"\nChi-square results:")
    print(f"  χ² = {chi2:.2f}")
    print(f"  p-value = {pval:.2e}")
    print(f"  degrees of freedom = {dof}")

    print(f"\nExpected frequencies:")
    print(f"                      | Producers | Non-Producers |")
    print(f"  Increased in Disease| {expected[0][0]:9.1f} | {expected[0][1]:13.1f} |")
    print(f"  Decreased in Disease| {expected[1][0]:9.1f} | {expected[1][1]:13.1f} |")

    return chi2, pval, contingency

def compare_with_published(recalc, disease):
    """Compare recalculated values with published and current results."""
    print(f"\n{'='*80}")
    print(f"COMPARISON WITH PUBLISHED AND CURRENT RESULTS - {disease}")
    print(f"{'='*80}")

    pub = PUBLISHED_RESULTS[disease]
    curr = CURRENT_RESULTS[disease]

    print(f"\n{'Metric':<30} | {'Published':<12} | {'Current Repo':<14} | {'Recalculated':<14}")
    print(f"{'-'*30} | {'-'*12} | {'-'*14} | {'-'*14}")

    metrics = [
        ('Increased Total', 'increased_total'),
        ('Increased Producers', 'increased_producers'),
        ('Decreased Total', 'decreased_total'),
        ('Decreased Producers', 'decreased_producers'),
    ]

    for label, key in metrics:
        pub_val = pub.get(key, 'N/A')
        curr_val = curr.get(key, 'N/A')
        recalc_val = recalc.get(key, 'N/A')

        pub_str = f"{pub_val:,}" if isinstance(pub_val, int) else str(pub_val)
        curr_str = f"{curr_val:,}" if isinstance(curr_val, (int, float)) else str(curr_val)
        recalc_str = f"{recalc_val:,}" if isinstance(recalc_val, (int, float)) else str(recalc_val)

        print(f"{label:<30} | {pub_str:<12} | {curr_str:<14} | {recalc_str:<14}")

    # Compare chi-square
    print(f"\n{'Chi-square':<30} | {pub['chi2']:<12.2f} | {curr['chi2']:<14.2f} | {'[computed]':<14}")
    print(f"{'P-value':<30} | {pub['pval']:<12.2e} | {curr['pval']:<14.2e} | {'[computed]':<14}")

def trace_sample_microbes(counted_microbes, producer_set, n=10):
    """Trace a sample of microbes to verify correct classification."""
    print(f"\n{'='*80}")
    print(f"TRACING SAMPLE MICROBES")
    print(f"{'='*80}")

    # Sample from both directions
    increased = [m for m in counted_microbes if m['direction'] == 'increased']
    decreased = [m for m in counted_microbes if m['direction'] == 'decreased']

    print(f"\nSample from INCREASED disease group:")
    for i, microbe in enumerate(increased[:n//2], 1):
        is_prod = "✓ PRODUCER" if microbe['is_producer'] else "✗ Non-producer"
        print(f"  {i}. {microbe['microbe']}")
        print(f"     {is_prod} | Total: {microbe['total']}, Producers: {microbe['producers']}")
        print(f"     Relationship: {microbe['relationship']}")

    print(f"\nSample from DECREASED disease group:")
    for i, microbe in enumerate(decreased[:n//2], 1):
        is_prod = "✓ PRODUCER" if microbe['is_producer'] else "✗ Non-producer"
        print(f"  {i}. {microbe['microbe']}")
        print(f"     {is_prod} | Total: {microbe['total']}, Producers: {microbe['producers']}")
        print(f"     Relationship: {microbe['relationship']}")

def check_string_matching_robustness(counted_microbes):
    """Check if substring matching could cause misclassification."""
    print(f"\n{'='*80}")
    print(f"STRING MATCHING ROBUSTNESS CHECK")
    print(f"{'='*80}")

    # Check for any relationships that contain both "increased" and "decreased"
    ambiguous = []
    for microbe in counted_microbes:
        rel = microbe['relationship'].lower()
        if 'increased' in rel and 'decreased' in rel:
            ambiguous.append(microbe)

    if ambiguous:
        print(f"⚠ WARNING: Found {len(ambiguous)} microbes with ambiguous relationships:")
        for m in ambiguous[:5]:
            print(f"  - {m['microbe']}: {m['relationship']}")
    else:
        print(f"✓ No ambiguous relationships found (no relationship contains both 'increased' and 'decreased')")

    # Check for unexpected relationship strings
    unique_rels = set(m['relationship'] for m in counted_microbes)
    print(f"\nUnique relationship strings found ({len(unique_rels)} total):")
    for rel in sorted(unique_rels)[:10]:
        count = sum(1 for m in counted_microbes if m['relationship'] == rel)
        print(f"  - '{rel}' (n={count})")

def investigate_counterintuitive_pattern(recalc, disease):
    """
    Investigate the counterintuitive pattern where butyrate producers
    appear in BOTH increased and decreased groups.
    """
    print(f"\n{'='*80}")
    print(f"INVESTIGATING COUNTERINTUITIVE PATTERN - {disease}")
    print(f"{'='*80}")

    increased_prod = recalc['increased_producers']
    increased_total = recalc['increased_total']
    decreased_prod = recalc['decreased_producers']
    decreased_total = recalc['decreased_total']

    print(f"\nKey Question: Why are there butyrate producers in BOTH groups?")
    print(f"  Increased disease: {increased_prod} producers out of {increased_total} total")
    print(f"  Decreased disease: {decreased_prod} producers out of {decreased_total} total")

    print(f"\nInterpretation:")
    print(f"  1. The analysis is NOT counting individual producer species")
    print(f"  2. It's counting ASSOCIATIONS between microbes and disease")
    print(f"  3. A single producer can have multiple associations (increased OR decreased)")

    if increased_total > 0 and decreased_total > 0:
        inc_prop = increased_prod / increased_total
        dec_prop = decreased_prod / decreased_total

        print(f"\n  Proportion of producers:")
        print(f"    In increased-disease: {inc_prop:.4f} ({inc_prop*100:.2f}%)")
        print(f"    In decreased-disease: {dec_prop:.4f} ({dec_prop*100:.2f}%)")

        if dec_prop > inc_prop:
            print(f"\n  ✓ EXPECTED PATTERN: Producer proportion is HIGHER in decreased-disease")
            print(f"    This means producers are DEPLETED in disease (protective effect)")
            print(f"    Ratio: {dec_prop/inc_prop:.2f}x more producers in protective group")
        else:
            print(f"\n  ⚠ UNEXPECTED PATTERN: Producer proportion is HIGHER in increased-disease")
            print(f"    This would suggest producers are ENRICHED in disease (risk factor)")
            print(f"    Ratio: {inc_prop/dec_prop:.2f}x more producers in disease group")

def main():
    """Main validation workflow."""
    print("""
    ╔════════════════════════════════════════════════════════════════════════════╗
    ║                    FIGURE 6 VALIDATION AND AUDIT                           ║
    ║                                                                            ║
    ║  Validating statistical results for butyrate producer disease             ║
    ║  associations in Inflammatory Bowel Disease (IBD) and Parkinson's          ║
    ║  Disease (PD).                                                             ║
    ║                                                                            ║
    ║  Goal: Verify counts, directionality, and statistical tests               ║
    ╚════════════════════════════════════════════════════════════════════════════╝
    """)

    # Load gold standard producers
    producer_set = load_gold_standard()

    if not producer_set:
        print("\n❌ FATAL ERROR: Cannot proceed without gold standard data")
        return

    # Analyze each disease
    for disease in ['IBD', 'PD']:
        print(f"\n\n")
        print(f"{'#'*80}")
        print(f"{'#'*80}")
        print(f"###{' '*74}###")
        print(f"###   ANALYZING {disease:<64}###")
        print(f"###{' '*74}###")
        print(f"{'#'*80}")
        print(f"{'#'*80}")

        # Load classification data
        classification_data = load_classification_json(disease)

        if not classification_data['species'] and not classification_data['strain']:
            print(f"\n⚠ WARNING: No classification data found for {disease}, skipping...")
            continue

        # Analyze directionality patterns
        all_entries = analyze_directionality(classification_data, disease)

        # Recalculate aggregations
        recalc = recalculate_aggregations(classification_data, producer_set, disease)

        # Compute chi-square
        chi2, pval, contingency = compute_chi_square(recalc)
        recalc['chi2'] = chi2
        recalc['pval'] = pval

        # Compare with published
        compare_with_published(recalc, disease)

        # Trace sample microbes
        trace_sample_microbes(recalc['counted_microbes'], producer_set)

        # Check string matching robustness
        check_string_matching_robustness(recalc['counted_microbes'])

        # Investigate counterintuitive pattern
        investigate_counterintuitive_pattern(recalc, disease)

    print(f"\n\n")
    print(f"{'='*80}")
    print(f"VALIDATION COMPLETE")
    print(f"{'='*80}")
    print(f"\nNext steps:")
    print(f"  1. Review the recalculated counts vs published/current values")
    print(f"  2. Check if directionality is preserved correctly")
    print(f"  3. Verify chi-square test construction")
    print(f"  4. Investigate any discrepancies found")
    print(f"\nOutput: This analysis should help identify if there are any")
    print(f"        logic inversions or data processing errors causing the")
    print(f"        counterintuitive pattern observed in Figure 6.")

if __name__ == '__main__':
    main()
