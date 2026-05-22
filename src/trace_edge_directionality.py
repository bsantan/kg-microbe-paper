#!/usr/bin/env python3
"""
Edge Directionality Tracing Script
===================================

This script traces individual edges through the differentiate_edge_direction()
transformation to verify that directionality is preserved correctly.

Use this to manually verify that "increased_likelihood_of" and
"decreased_likelihood_of" maintain their semantic meaning after transformation.

Usage:
    python trace_edge_directionality.py --disease IBD
    python trace_edge_directionality.py --disease PD --sample-size 20
"""

import argparse
import pandas as pd
import duckdb
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from classification_utils import differentiate_edge_direction

BASE_DIR = Path(__file__).parent
INPUT_DIR = BASE_DIR / "Input_Files"

# Try alternative input locations
if not INPUT_DIR.exists():
    INPUT_DIR = BASE_DIR.parent / "data" / "Input_Files"

KG_EDGES_FILE = INPUT_DIR / "kg-microbe-biomedical-function-cat" / "merged-kg_edges.tsv"

DISEASE_MAP = {
    'IBD': 'MONDO:0005265',  # Inflammatory Bowel Disease
    'PD': 'MONDO:0005180'    # Parkinson's Disease
}

def load_disease_edges(disease_id, sample_size=None):
    """
    Load disease-microbe edges from KG using DuckDB for efficient filtering.

    Args:
        disease_id: MONDO ID for the disease
        sample_size: Number of edges to sample (None = all)

    Returns:
        DataFrame with subject, predicate, object columns
    """
    print(f"\nLoading edges for {disease_id}...")

    if not KG_EDGES_FILE.exists():
        print(f"❌ ERROR: KG edges file not found at {KG_EDGES_FILE}")
        print(f"   Please ensure Input_Files directory is set up correctly")
        return None

    conn = duckdb.connect(":memory:")

    # Filter to disease-relevant edges only
    query = f"""
    SELECT subject, predicate, object
    FROM read_csv_auto('{KG_EDGES_FILE}', delim='\\t', null_padding=true)
    WHERE ((subject LIKE 'NCBITaxon:%' AND object = '{disease_id}')
           OR (object LIKE 'NCBITaxon:%' AND subject = '{disease_id}'))
      AND predicate IS NOT NULL
      AND subject IS NOT NULL
      AND object IS NOT NULL
    """

    if sample_size:
        query += f" LIMIT {sample_size}"

    try:
        edges_df = conn.execute(query).df()
        print(f"✓ Loaded {len(edges_df)} edges")
    except Exception as e:
        print(f"❌ ERROR loading edges: {e}")
        return None
    finally:
        conn.close()

    return edges_df

def trace_single_edge(row, edge_num):
    """
    Trace a single edge through the transformation.

    Args:
        row: DataFrame row with subject, predicate, object
        edge_num: Edge number for display

    Returns:
        Dict with before/after information
    """
    subject = row['subject']
    predicate = row['predicate']
    obj = row['object']

    # Determine original direction
    if subject.startswith('NCBITaxon:'):
        direction = "Microbe → Disease"
        microbe = subject
        disease = obj
    else:
        direction = "Disease → Microbe"
        microbe = obj
        disease = subject

    # Interpret predicate meaning
    if 'increased' in predicate.lower():
        if direction == "Microbe → Disease":
            meaning = "Microbe INCREASES risk of disease"
        else:
            meaning = "Disease is MORE common when microbe present"
    elif 'decreased' in predicate.lower():
        if direction == "Microbe → Disease":
            meaning = "Microbe DECREASES risk of disease (protective)"
        else:
            meaning = "Disease is LESS common when microbe present (protective)"
    else:
        meaning = "Unknown predicate type"

    return {
        'edge_num': edge_num,
        'original_subject': subject,
        'original_predicate': predicate,
        'original_object': obj,
        'direction': direction,
        'microbe': microbe,
        'disease': disease,
        'meaning': meaning
    }

def apply_transformation(edges_df):
    """
    Apply differentiate_edge_direction() transformation.

    Args:
        edges_df: DataFrame with edges

    Returns:
        Transformed DataFrame
    """
    print(f"\nApplying differentiate_edge_direction() transformation...")

    try:
        transformed_df = differentiate_edge_direction(edges_df)
        print(f"✓ Transformation complete")
        return transformed_df
    except Exception as e:
        print(f"❌ ERROR during transformation: {e}")
        return None

def trace_after_transformation(original_info, transformed_row):
    """
    Trace the edge after transformation.

    Args:
        original_info: Dict with original edge info
        transformed_row: DataFrame row after transformation

    Returns:
        Dict with after-transformation info
    """
    new_subject = transformed_row['subject']
    new_object = transformed_row['object']

    # Check if subject/object were swapped
    swapped = new_subject != original_info['original_subject']

    # Extract predicate from object string
    if '_' in new_object:
        # Object is now "predicate_disease"
        parts = new_object.split('_', 1)
        if len(parts) >= 2:
            new_predicate = parts[0]
        else:
            new_predicate = new_object
    else:
        new_predicate = "unknown"

    # Determine new direction classification
    if 'increased' in new_object.lower():
        new_direction = "increased_likelihood_of"
        classification = "INCREASED-DISEASE (risk)"
    elif 'decreased' in new_object.lower():
        new_direction = "decreased_likelihood_of"
        classification = "DECREASED-DISEASE (protective)"
    else:
        new_direction = "unknown"
        classification = "UNKNOWN"

    return {
        'new_subject': new_subject,
        'new_object': new_object,
        'swapped': swapped,
        'new_direction': new_direction,
        'classification': classification
    }

def check_directionality_preserved(original_info, after_info):
    """
    Check if directionality semantic meaning is preserved.

    Args:
        original_info: Dict with original edge info
        after_info: Dict with transformed edge info

    Returns:
        Bool indicating if directionality is preserved
    """
    original_is_protective = 'decreased' in original_info['original_predicate'].lower()
    after_is_protective = after_info['classification'] == "DECREASED-DISEASE (protective)"

    return original_is_protective == after_is_protective

def print_trace(original_info, after_info, directionality_preserved):
    """Print detailed trace of edge transformation."""

    print(f"\n{'='*80}")
    print(f"Edge #{original_info['edge_num']}")
    print(f"{'='*80}")

    print(f"\n📋 BEFORE TRANSFORMATION:")
    print(f"  Subject:   {original_info['original_subject']}")
    print(f"  Predicate: {original_info['original_predicate']}")
    print(f"  Object:    {original_info['original_object']}")
    print(f"\n  Direction: {original_info['direction']}")
    print(f"  Microbe:   {original_info['microbe']}")
    print(f"  Disease:   {original_info['disease']}")
    print(f"  Meaning:   {original_info['meaning']}")

    print(f"\n🔄 AFTER TRANSFORMATION:")
    print(f"  Subject:   {after_info['new_subject']}")
    print(f"  Object:    {after_info['new_object']}")
    print(f"  Swapped:   {after_info['swapped']}")
    print(f"  Direction: {after_info['new_direction']}")
    print(f"  Classification: {after_info['classification']}")

    print(f"\n✓ VALIDATION:")
    if directionality_preserved:
        print(f"  ✅ Directionality PRESERVED")
        print(f"     Original meaning matches transformed classification")
    else:
        print(f"  ❌ Directionality INVERTED")
        print(f"     Original meaning DOES NOT match transformed classification")
        print(f"     ⚠️ THIS IS A BUG!")

def main():
    """Main tracing workflow."""
    parser = argparse.ArgumentParser(
        description='Trace edge directionality through transformation pipeline'
    )
    parser.add_argument(
        '--disease',
        choices=['IBD', 'PD'],
        required=True,
        help='Disease to analyze (IBD or PD)'
    )
    parser.add_argument(
        '--sample-size',
        type=int,
        default=10,
        help='Number of edges to trace (default: 10)'
    )
    parser.add_argument(
        '--show-all',
        action='store_true',
        help='Show all edges (not just sample)'
    )

    args = parser.parse_args()

    disease_name = args.disease
    disease_id = DISEASE_MAP[disease_name]
    sample_size = None if args.show_all else args.sample_size

    print("="*80)
    print(f"EDGE DIRECTIONALITY TRACING - {disease_name}".center(80))
    print("="*80)

    # Load edges
    edges_df = load_disease_edges(disease_id, sample_size)

    if edges_df is None or len(edges_df) == 0:
        print("❌ No edges to trace. Exiting.")
        return

    # Apply transformation
    transformed_df = apply_transformation(edges_df)

    if transformed_df is None:
        print("❌ Transformation failed. Exiting.")
        return

    # Trace each edge
    print(f"\n{'='*80}")
    print(f"TRACING {len(edges_df)} EDGES")
    print(f"{'='*80}")

    preserved_count = 0
    inverted_count = 0

    for idx, (_, original_row) in enumerate(edges_df.iterrows(), 1):
        # Trace original
        original_info = trace_single_edge(original_row, idx)

        # Find corresponding transformed row
        # Match by microbe ID (should be subject after transformation)
        transformed_row = transformed_df[
            transformed_df['subject'] == original_info['microbe']
        ]

        if len(transformed_row) == 0:
            print(f"\n⚠️ Edge #{idx}: Could not find transformed version")
            continue

        transformed_row = transformed_row.iloc[0]

        # Trace after transformation
        after_info = trace_after_transformation(original_info, transformed_row)

        # Check if directionality preserved
        directionality_preserved = check_directionality_preserved(
            original_info, after_info
        )

        if directionality_preserved:
            preserved_count += 1
        else:
            inverted_count += 1

        # Print trace
        print_trace(original_info, after_info, directionality_preserved)

    # Summary
    print(f"\n\n{'='*80}")
    print(f"SUMMARY")
    print(f"{'='*80}")
    print(f"Total edges traced: {len(edges_df)}")
    print(f"Directionality preserved: {preserved_count} ✅")
    print(f"Directionality inverted: {inverted_count} ❌")

    if inverted_count == 0:
        print(f"\n✅ ALL EDGES: Directionality correctly preserved!")
    else:
        print(f"\n❌ WARNING: {inverted_count} edges have inverted directionality!")
        print(f"   This indicates a bug in the edge transformation logic.")

    print(f"\n{'='*80}")

if __name__ == '__main__':
    main()
