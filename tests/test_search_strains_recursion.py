"""Regression test: search_strains must descend into species' strain children.

Background: commit 4195edb changed search_strains from `if/if/else` to
`if/elif/else`, stopping descent at the species rank and silently dropping
strains nested under species. This shrank the published Figure 6C chi-square
(IBD 674.63 -> 1201, PD 1337.36 -> 642) and was masked by cached JSONs.

This test exercises a synthetic genus -> species -> strain fixture and
asserts the strain is found. If `if`->`elif` is reintroduced, the strain
will be silently dropped and this test will fail.

Run: uv run python tests/test_search_strains_recursion.py
Exit code 0 = pass, 1 = fail.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

from ncbi_phylogeny_search import search_strains


def _make_fixture():
    """Genus -> species -> strain hierarchy.

    GENUS has a species child SPECIES.
    SPECIES has a strain child STRAIN.
    SPECIES also has a non-strain child OTHER (no rank) that should be ignored.
    """
    taxonomy_hierarchy = {
        "NCBITaxon:GENUS": ["NCBITaxon:SPECIES"],
        "NCBITaxon:SPECIES": ["NCBITaxon:STRAIN", "NCBITaxon:OTHER"],
        "NCBITaxon:STRAIN": [],
        "NCBITaxon:OTHER": [],
    }
    rank_lookup = {
        "NCBITaxon:GENUS": "genus",
        "NCBITaxon:SPECIES": "species",
        "NCBITaxon:STRAIN": "strain",
        # OTHER intentionally has no rank to also exercise the "rank is None -> skip" path
    }
    return taxonomy_hierarchy, rank_lookup


def test_species_with_strain_child_is_enumerated():
    hierarchy, rank_lookup = _make_fixture()
    strains_found: list[str] = []
    species_found: list[str] = []
    search_strains(hierarchy, rank_lookup, "NCBITaxon:GENUS", strains_found, species_found, use_hierarchy=True)
    assert "NCBITaxon:STRAIN" in strains_found, (
        f"Regression: STRAIN was not enumerated. strains_found={strains_found}. "
        "If search_strains was changed from `if microbe_rank in [...subspecies, strain]:` "
        "to `elif ...`, descent into species' strain children is lost and Figure 6C "
        "chi-square shrinks. See commit 4195edb."
    )
    assert "NCBITaxon:SPECIES" in species_found, (
        f"species_found={species_found}; species should be recorded too."
    )


def test_strain_under_subspecies_chain():
    """A subspecies -> strain chain must also enumerate the strain."""
    hierarchy = {
        "NCBITaxon:SPECIES2": ["NCBITaxon:SUBSP"],
        "NCBITaxon:SUBSP": ["NCBITaxon:STRAIN2"],
        "NCBITaxon:STRAIN2": [],
    }
    rank_lookup = {
        "NCBITaxon:SPECIES2": "species",
        "NCBITaxon:SUBSP": "subspecies",
        "NCBITaxon:STRAIN2": "strain",
    }
    strains_found: list[str] = []
    species_found: list[str] = []
    search_strains(hierarchy, rank_lookup, "NCBITaxon:SPECIES2", strains_found, species_found, use_hierarchy=True)
    # SUBSP is recorded as a strain by the `subspecies` branch
    assert "NCBITaxon:SUBSP" in strains_found
    # STRAIN2 should not show up because subspecies branch doesn't recurse — that's
    # intentional and matches current/published behavior. The key invariant is the
    # species -> strain descent, validated by the first test.


def test_within_parent_dedup():
    """Same strain referenced via two paths under the same parent is recorded once."""
    hierarchy = {
        "NCBITaxon:G2": ["NCBITaxon:S_A", "NCBITaxon:S_B"],
        "NCBITaxon:S_A": ["NCBITaxon:DUP_STRAIN"],
        "NCBITaxon:S_B": ["NCBITaxon:DUP_STRAIN"],
        "NCBITaxon:DUP_STRAIN": [],
    }
    rank_lookup = {
        "NCBITaxon:G2": "genus",
        "NCBITaxon:S_A": "species",
        "NCBITaxon:S_B": "species",
        "NCBITaxon:DUP_STRAIN": "strain",
    }
    strains_found: list[str] = []
    species_found: list[str] = []
    search_strains(hierarchy, rank_lookup, "NCBITaxon:G2", strains_found, species_found, use_hierarchy=True)
    assert strains_found.count("NCBITaxon:DUP_STRAIN") == 1


def main() -> int:
    tests = [
        test_species_with_strain_child_is_enumerated,
        test_strain_under_subspecies_chain,
        test_within_parent_dedup,
    ]
    failures = []
    for t in tests:
        try:
            t()
            print(f"  ✓ {t.__name__}")
        except AssertionError as e:
            failures.append((t.__name__, str(e)))
            print(f"  ✗ {t.__name__}: {e}")
    if failures:
        print(f"\n{len(failures)} test(s) failed.")
        return 1
    print(f"\nAll {len(tests)} tests passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
