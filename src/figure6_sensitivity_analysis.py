"""Sensitivity analysis for Figure 6C chi-square (Codex finding 1).

The published Figure 6C contingency builds a 2x2 by summing per-parent strain
(or species fallback) occurrences within each disease direction
(increased/decreased likelihood), then runs chi2_contingency. The same
descendant strain can appear under multiple disease parents AND in both
likelihood directions, so the summed counts are not independent observations
under the chi-square model.

This script does NOT alter the published outputs. It re-runs the test under
five scenarios so the response can disclose the multiplicity transparently:

  A0 - Published baseline: replays the per-parent sum exactly (self-check).
  A1 - Globally-unique strains per disease: each strain counted once;
       cross-bin overlaps assigned to the first-seen direction (with the
       drop-count reported).
  A2 - Unique strains, conflicting overlaps removed entirely. NOTE: dropping
       conflicts is not statistically neutral — the conflict bucket has a
       producer fraction much closer to the increased-only bin than the
       decreased-only bin, so excluding it inflates the apparent
       decreased-likelihood enrichment. A2 also dedups exact NCBITaxon IDs
       only; sibling/descendant phylogenetic clustering still violates
       chi-square independence. See A2_keep_overlaps for the conservative
       inclusion and A3 for the cluster-aware test.
  A2_keep_overlaps - Same unique-descendant unit as A2, but retains cross-
       direction conflict descendants in BOTH bins. Discloses the bias
       introduced by A2's drop-conflicts choice.
  A3 - Cluster permutation test: permutes parent -> direction labels while
       keeping each parent's strain bundle intact, then reports the empirical
       p-value vs the observed A0 chi-square. This is cluster-aware in the
       sense that it respects the (parent -> strain set) hierarchy, which
       is what A0/A1/A2/A2_keep_overlaps' chi-square cannot address.

Diagnostics CSV reports the multiplicity numbers (total occurrences vs
unique strains, cross-bin overlap counts, cross-parent fan-out histogram).

Run: uv run python src/figure6_sensitivity_analysis.py
Outputs are written to ./data/Intermediate_Files/figure6c_sensitivity_<DIS>.csv
and ./data/Intermediate_Files/figure6c_multiplicity_diagnostics_<DIS>.csv
"""

from __future__ import annotations

import json
import os
import random
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

INTERMEDIATE_DIR = "./data/Intermediate_Files"
COMPETENCIES_DIR = "./data/Intermediate_Files_Competencies/butyrate_produces"
GOLD_STANDARD_CSV = COMPETENCIES_DIR + "/Gold_Standard_Species_Overlap_butyrate_produces.csv"

# Historical A0 self-check anchors. The PUBLISHED Figure 6C chi-square is the
# deduplicated A2 contingency (IBD 497.26, PD 918.80; see
# data/Intermediate_Files/{IBD,PD}_Classification_butyrate_producers_summary.csv).
# A0 below is the per-parent baseline used for reproducibility self-check only:
# A2 derives from the same strain/species enumeration that produces A0, so if
# A0 here does not match these values to 10 decimals, the sensitivity script is
# reading enumeration outputs different from what A2 depends on.
A0_SELF_CHECK_BASELINE = {
    "IBD": {"chi2": 674.6307631130468},
    "PD": {"chi2": 1337.3569687912984},
}

DISEASES = ["IBD", "PD"]
PERMUTATION_N = 1_000_000
PERMUTATION_SEED = 12  # match RANDOM_SEED in constants.py


def load_producers() -> set:
    df = pd.read_csv(GOLD_STANDARD_CSV)
    return set(
        df[(df["Organismal"] == 1) | (df["Functional"] == 1) | (df["Functional_EC"] == 1)]["Value"]
    )


def load_disease_inputs(disease: str):
    ranks_csv = f"{INTERMEDIATE_DIR}/outcome_to_NCBITaxon_cleaned_ranks_butyrate_production_{disease}.csv"
    strain_json = f"{INTERMEDIATE_DIR}/classification_butyrate_produces_{disease}_microbes_strain.json"
    species_json = f"{INTERMEDIATE_DIR}/classification_butyrate_produces_{disease}_microbes_species.json"
    ranks_df = pd.read_csv(ranks_csv)
    with open(strain_json) as f:
        strain_dict = json.load(f)
    with open(species_json) as f:
        species_dict = json.load(f)
    return ranks_df, strain_dict, species_dict


def parent_direction(rel: str) -> str:
    if "increased" in rel:
        return "increased"
    if "decreased" in rel:
        return "decreased"
    return "unknown"


def per_parent_counts_from_ranks(ranks_df: pd.DataFrame):
    """Per-parent (n_total, n_prod, is_increased) using the published
    Num_Strains>0-else-Num_Species fallback. Summing these reproduces the A0
    contingency exactly, so the A3 permutation test acts on the same statistic
    as the figure (not on descendant-list lengths, which can differ slightly)."""
    n_total, n_prod, is_inc = [], [], []
    for _, row in ranks_df.iterrows():
        if row["Num_Strains"] > 0:
            t = int(row["Num_Strains"]); p = int(row["Num_Strains_Butyrate_Producers"])
        else:
            t = int(row["Num_Species"]); p = int(row["Num_Species_Butyrate_Producers"])
        d = parent_direction(row["Disease_Relationship"])
        if d not in ("increased", "decreased"):
            continue
        n_total.append(t); n_prod.append(p); is_inc.append(d == "increased")
    return np.array(n_total), np.array(n_prod), np.array(is_inc, dtype=bool)


def effect_sizes(inc_total, inc_prod, dec_total, dec_prod, chi2):
    """Sample-size-independent effect sizes for the 2x2 producer x direction table.
    Enrichment (risk) ratio and odds ratio express how concentrated producers are
    among decreased-likelihood taxa; Cramér's V is the chi-square-based effect size
    (small even when chi2 is huge, because chi2 scales with N)."""
    inc_np = inc_total - inc_prod
    dec_np = dec_total - dec_prod
    N = inc_total + dec_total
    p_inc = inc_prod / inc_total if inc_total else float("nan")
    p_dec = dec_prod / dec_total if dec_total else float("nan")
    enrichment_ratio = (p_dec / p_inc) if p_inc else float("nan")
    odds_ratio = ((dec_prod * inc_np) / (inc_prod * dec_np)
                  if inc_prod and dec_np else float("nan"))
    cramers_v = (chi2 / N) ** 0.5 if N else float("nan")
    return {
        "producer_frac_increased": round(p_inc, 6),
        "producer_frac_decreased": round(p_dec, 6),
        "enrichment_ratio_dec_over_inc": round(enrichment_ratio, 4),
        "odds_ratio_dec_vs_inc": round(odds_ratio, 4),
        "cramers_v": round(cramers_v, 4),
    }


def compute_a0(ranks_df: pd.DataFrame):
    """Replay the per-parent sum exactly as Classification_gold_standard_comparison.py."""
    inc_total = dec_total = inc_prod = dec_prod = 0
    for _, row in ranks_df.iterrows():
        if row["Num_Strains"] > 0:
            n_total = int(row["Num_Strains"])
            n_prod = int(row["Num_Strains_Butyrate_Producers"])
        else:
            n_total = int(row["Num_Species"])
            n_prod = int(row["Num_Species_Butyrate_Producers"])
        d = parent_direction(row["Disease_Relationship"])
        if d == "increased":
            inc_total += n_total
            inc_prod += n_prod
        elif d == "decreased":
            dec_total += n_total
            dec_prod += n_prod
    table = [[inc_prod, inc_total - inc_prod], [dec_prod, dec_total - dec_prod]]
    chi2, pval, dof, _ = chi2_contingency(table)
    out = {
        "scenario": "A0_published_baseline",
        "unit": "per-parent occurrence (not deduped)",
        "increased_total": inc_total,
        "increased_producers": inc_prod,
        "decreased_total": dec_total,
        "decreased_producers": dec_prod,
        "chi2": chi2,
        "pval": pval,
        "dof": dof,
        "contingency": str(table),
    }
    out.update(effect_sizes(inc_total, inc_prod, dec_total, dec_prod, chi2))
    return out


def per_parent_descendants(ranks_df, strain_dict, species_dict):
    """For each parent return (direction, descendants_used). Strains take priority
    over species, matching the published Num_Strains>0 fallback rule."""
    out = []
    for _, row in ranks_df.iterrows():
        parent = row["Name"]
        d = parent_direction(row["Disease_Relationship"])
        if row["Rank"] in ("strain", "subspecies"):
            # A parent that is itself a strain/subspecies represents itself
            # (consistent with how ranks_df / the A0 contingency count it).
            descs = [parent]
            unit = "strain"
        elif row["Num_Strains"] > 0:
            descs = list(strain_dict.get(parent, []))
            unit = "strain"
        else:
            descs = [parent] if row["Rank"] == "species" else list(species_dict.get(parent, []))
            unit = "species"
        out.append({"parent": parent, "direction": d, "descendants": descs, "unit": unit})
    return out


def compute_a1(per_parent, producers, disease: str):
    """Globally-unique descendants per disease, cross-bin overlap assigned to the
    first-seen direction in iteration order (deterministic given ranks_df order)."""
    seen = {}
    dropped_from_other = Counter()
    for entry in per_parent:
        d = entry["direction"]
        for desc in entry["descendants"]:
            if desc not in seen:
                seen[desc] = d
            elif seen[desc] != d:
                dropped_from_other[d] += 1
    inc = [k for k, v in seen.items() if v == "increased"]
    dec = [k for k, v in seen.items() if v == "decreased"]
    inc_prod = sum(1 for k in inc if k in producers)
    dec_prod = sum(1 for k in dec if k in producers)
    table = [[inc_prod, len(inc) - inc_prod], [dec_prod, len(dec) - dec_prod]]
    chi2, pval, dof, _ = chi2_contingency(table)
    return {
        "scenario": "A1_unique_strains_first_seen_wins",
        "unit": "unique descendant per disease",
        "increased_total": len(inc),
        "increased_producers": inc_prod,
        "decreased_total": len(dec),
        "decreased_producers": dec_prod,
        "chi2": chi2,
        "pval": pval,
        "dof": dof,
        "contingency": str(table),
        "dropped_from_increased_due_to_overlap": dropped_from_other.get("increased", 0),
        "dropped_from_decreased_due_to_overlap": dropped_from_other.get("decreased", 0),
    }


def compute_a2(per_parent, producers, disease: str):
    """Unique descendants per disease, dropping any descendant that appeared in
    both directions.

    Note: dropping conflicts is *not* statistically neutral here. The conflict
    bucket has a producer fraction (~1.4 %% IBD, ~0.9 %% PD) much closer to the
    increased-only bin than the decreased-only bin, so excluding it inflates
    the apparent decreased-likelihood enrichment. The A2_keep_overlaps scenario
    below retains conflicts in both bins so reviewers can see the chi-square
    under the more conservative inclusion. A2 also dedups exact NCBITaxon IDs
    only — sibling/descendant phylogenetic clustering still violates
    chi-square independence; see A3 for the cluster-aware permutation test.
    """
    direction_sets = defaultdict(set)
    for entry in per_parent:
        for desc in entry["descendants"]:
            direction_sets[entry["direction"]].add(desc)
    inc_set = direction_sets["increased"]
    dec_set = direction_sets["decreased"]
    conflicting = inc_set & dec_set
    inc = inc_set - conflicting
    dec = dec_set - conflicting
    inc_prod = sum(1 for k in inc if k in producers)
    dec_prod = sum(1 for k in dec if k in producers)
    table = [[inc_prod, len(inc) - inc_prod], [dec_prod, len(dec) - dec_prod]]
    chi2, pval, dof, _ = chi2_contingency(table)
    return {
        "scenario": "A2_unique_strains_overlaps_excluded",
        "unit": "unique descendant, conflicts excluded",
        "increased_total": len(inc),
        "increased_producers": inc_prod,
        "decreased_total": len(dec),
        "decreased_producers": dec_prod,
        "chi2": chi2,
        "pval": pval,
        "dof": dof,
        "contingency": str(table),
        "conflicting_descendants_excluded": len(conflicting),
    }


def compute_a2_keep_overlaps(per_parent, producers, disease: str):
    """Unique descendants per disease, retaining cross-direction conflicts in
    BOTH bins (no dedup of overlaps). Discloses how non-neutral A2's
    drop-conflicts choice is: conflict descendants have a producer fraction
    between the two pure bins, so dropping them strengthens the apparent
    enrichment in the decreased-likelihood direction. Keeping them in both bins
    is the conservative inclusion of the same unique-descendant unit.
    """
    direction_sets = defaultdict(set)
    for entry in per_parent:
        for desc in entry["descendants"]:
            direction_sets[entry["direction"]].add(desc)
    inc = direction_sets["increased"]
    dec = direction_sets["decreased"]
    overlap = inc & dec
    inc_prod = sum(1 for k in inc if k in producers)
    dec_prod = sum(1 for k in dec if k in producers)
    overlap_prod = sum(1 for k in overlap if k in producers)
    table = [[inc_prod, len(inc) - inc_prod], [dec_prod, len(dec) - dec_prod]]
    chi2, pval, dof, _ = chi2_contingency(table)
    return {
        "scenario": "A2_keep_overlaps_unique_strains_both_bins",
        "unit": "unique descendant, conflicts retained in both bins",
        "increased_total": len(inc),
        "increased_producers": inc_prod,
        "decreased_total": len(dec),
        "decreased_producers": dec_prod,
        "chi2": chi2,
        "pval": pval,
        "dof": dof,
        "contingency": str(table),
        "conflict_descendants_retained": len(overlap),
        "conflict_producers_retained": overlap_prod,
    }


def _yates_chi2_vec(inc_prod, inc_total, total_prod, N):
    """Vectorized 2x2 Pearson chi-square with Yates' continuity correction.
    total_prod (column total) and N are fixed across permutations; only the
    increased/decreased split varies. Matches scipy.stats.chi2_contingency default."""
    a = inc_prod.astype(np.int64)          # increased producers
    r1 = inc_total.astype(np.int64)        # increased total
    total_np = N - total_prod              # non-producers (fixed)
    b = r1 - a                             # increased non-producers
    c = total_prod - a                     # decreased producers
    d = total_np - b                       # decreased non-producers
    num = np.abs(a * d - b * c) - N / 2.0
    num = np.where(num < 0, 0.0, num)      # Yates floor
    denom = r1.astype(float) * (N - r1) * total_prod * total_np
    return np.where(denom > 0, N * num ** 2 / denom, 0.0)


def compute_a3(ranks_df, n_perm=PERMUTATION_N, seed=PERMUTATION_SEED, batch=50_000):
    """Cluster-level permutation test. The unit of permutation is the disease-
    associated PARENT taxon: each parent carries one direction and contributes its
    entire (n_total, n_prod) bundle to that direction. We permute the parent->
    direction labels (preserving the number of increased/decreased parents),
    recompute the same Yates chi-square the figure uses, and report the fraction
    of permutations whose chi-square is >= the observed value. This is valid under
    the actual clustered (and overlapping) assignment structure; reading the
    observed chi-square against the chi-square(1) distribution is not.

    Operates on the CSV-derived per-parent counts so the observed statistic equals
    the published A0 value exactly."""
    n_total, n_prod, is_inc = per_parent_counts_from_ranks(ranks_df)
    P = len(n_total)
    k = int(is_inc.sum())                  # number of increased parents (fixed)
    total_prod = int(n_prod.sum())
    N = int(n_total.sum())

    obs_inc_prod = np.array([int(n_prod[is_inc].sum())])
    obs_inc_total = np.array([int(n_total[is_inc].sum())])
    observed_chi2 = float(_yates_chi2_vec(obs_inc_prod, obs_inc_total, total_prod, N)[0])

    rng = np.random.default_rng(seed)
    ge = 0
    done = 0
    while done < n_perm:
        bsz = min(batch, n_perm - done)
        # each row selects k of P parents as 'increased' (relabeling, counts fixed)
        rnd = rng.random((bsz, P))
        order = np.argpartition(rnd, k, axis=1)
        mask = np.zeros((bsz, P), dtype=bool)
        rows = np.arange(bsz)[:, None]
        mask[rows, order[:, :k]] = True
        inc_prod = mask @ n_prod
        inc_total = mask @ n_total
        chi2v = _yates_chi2_vec(inc_prod, inc_total, total_prod, N)
        ge += int((chi2v >= observed_chi2 - 1e-9).sum())
        done += bsz

    emp_p = (ge + 1) / (n_perm + 1)
    return {
        "scenario": "A3_cluster_permutation",
        "unit": "per-parent bundle (parent->direction permuted, bundle intact)",
        "increased_total": int(obs_inc_total[0]),
        "increased_producers": int(obs_inc_prod[0]),
        "decreased_total": N - int(obs_inc_total[0]),
        "decreased_producers": total_prod - int(obs_inc_prod[0]),
        "chi2": observed_chi2,
        "pval": emp_p,
        "dof": "",
        "contingency": "",
        "n_permutations": n_perm,
        "permutation_seed": seed,
        "n_permutations_ge_observed": ge,
        "n_parents": P,
    }


def diagnostics(per_parent, disease: str):
    occurrences = 0
    unique = set()
    direction_sets = defaultdict(set)
    fanout = Counter()
    for entry in per_parent:
        occurrences += len(entry["descendants"])
        for desc in entry["descendants"]:
            unique.add(desc)
            direction_sets[entry["direction"]].add(desc)
            fanout[desc] += 1
    overlap = direction_sets["increased"] & direction_sets["decreased"]
    fan_hist = Counter(fanout.values())
    rows = [
        {"diagnostic": "n_parents", "value": len(per_parent)},
        {"diagnostic": "total_descendant_occurrences", "value": occurrences},
        {"diagnostic": "unique_descendants", "value": len(unique)},
        {"diagnostic": "ratio_occurrences_to_unique", "value": round(occurrences / len(unique), 4) if unique else 0},
        {"diagnostic": "descendants_in_both_directions", "value": len(overlap)},
        {"diagnostic": "increased_unique", "value": len(direction_sets['increased'])},
        {"diagnostic": "decreased_unique", "value": len(direction_sets['decreased'])},
    ]
    for k in sorted(fan_hist.keys()):
        rows.append({"diagnostic": f"fanout_eq_{k}_parents", "value": fan_hist[k]})
    return pd.DataFrame(rows)


def run_disease(disease: str, producers: set) -> None:
    ranks_df, strain_dict, species_dict = load_disease_inputs(disease)
    per_parent = per_parent_descendants(ranks_df, strain_dict, species_dict)

    a0 = compute_a0(ranks_df)

    # Self-check vs historical A0 baseline (10-decimal match). Note: A0 is the
    # reproducibility anchor here, not the published Figure 6C chi-square (that
    # is A2, derived from the same enumeration A0 sums per-parent).
    expected = A0_SELF_CHECK_BASELINE[disease]["chi2"]
    if abs(a0["chi2"] - expected) > 1e-9:
        print(f"  ⚠ A0 chi2 self-check FAILED for {disease}: "
              f"computed {a0['chi2']:.10f} vs A0 baseline {expected:.10f}")
    else:
        print(f"  ✓ A0 chi2 self-check OK for {disease}: {a0['chi2']:.10f}")

    a1 = compute_a1(per_parent, producers, disease)
    a2 = compute_a2(per_parent, producers, disease)
    a2_keep = compute_a2_keep_overlaps(per_parent, producers, disease)
    a3 = compute_a3(ranks_df)
    print(f"  A3 cluster permutation ({disease}): chi2={a3['chi2']:.4f}, "
          f"{a3['n_permutations_ge_observed']}/{a3['n_permutations']} >= observed, "
          f"empirical p={a3['pval']:.3e}")

    rows = [a0, a1, a2, a2_keep, a3]
    cols = [
        "scenario", "unit", "increased_total", "increased_producers",
        "decreased_total", "decreased_producers", "chi2", "pval", "dof",
        "contingency",
        "producer_frac_increased", "producer_frac_decreased",
        "enrichment_ratio_dec_over_inc", "odds_ratio_dec_vs_inc", "cramers_v",
        "dropped_from_increased_due_to_overlap",
        "dropped_from_decreased_due_to_overlap",
        "conflicting_descendants_excluded",
        "conflict_descendants_retained", "conflict_producers_retained",
        "n_permutations", "permutation_seed", "n_permutations_ge_observed",
    ]
    sens_df = pd.DataFrame(rows)
    for c in cols:
        if c not in sens_df.columns:
            sens_df[c] = ""
    sens_df = sens_df[cols]
    sens_path = f"{INTERMEDIATE_DIR}/figure6c_sensitivity_{disease}.csv"
    sens_df.to_csv(sens_path, index=False)
    print(f"  wrote {sens_path}")

    diag_df = diagnostics(per_parent, disease)
    diag_path = f"{INTERMEDIATE_DIR}/figure6c_multiplicity_diagnostics_{disease}.csv"
    diag_df.to_csv(diag_path, index=False)
    print(f"  wrote {diag_path}")


def main():
    print("Figure 6C sensitivity analysis")
    producers = load_producers()
    print(f"  butyrate producers in Gold Standard: {len(producers)}")
    for disease in DISEASES:
        print(f"--- {disease} ---")
        run_disease(disease, producers)
    print("Done.")


if __name__ == "__main__":
    main()
