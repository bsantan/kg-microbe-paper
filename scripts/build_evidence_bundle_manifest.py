"""Build a manifest.json and diff_report.txt for a 3-way evidence bundle.

Expects directory layout:
  data/evidence_bundles/<bundle_id>/
    shipped/<artifact files>
    fix_from_scratch/<artifact files>
    bug_4195edb/<artifact files>

Writes:
  data/evidence_bundles/<bundle_id>/manifest.json
  data/evidence_bundles/<bundle_id>/diff_report.txt

Usage: uv run python scripts/build_evidence_bundle_manifest.py <bundle_dir>
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import pandas as pd


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_dir(path: str) -> dict:
    out = {}
    if not os.path.isdir(path):
        return out
    for name in sorted(os.listdir(path)):
        p = os.path.join(path, name)
        if os.path.isfile(p):
            out[name] = {
                "sha256": sha256_file(p),
                "size": os.path.getsize(p),
            }
    return out


def load_summary(path: str) -> dict | None:
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path)
        if len(df) == 0:
            return None
        return df.iloc[0].to_dict()
    except Exception:
        return None


def disease_metrics(run_dir: str, disease: str) -> dict:
    """Pull the key numbers for one disease from the run directory."""
    summary = load_summary(os.path.join(run_dir, f"{disease}_Classification_butyrate_producers_summary.csv"))
    strain_json = os.path.join(run_dir, f"classification_butyrate_produces_{disease}_microbes_strain.json")
    species_json = os.path.join(run_dir, f"classification_butyrate_produces_{disease}_microbes_species.json")
    strain_occurrences = 0
    strain_keys = 0
    species_occurrences = 0
    species_keys = 0
    if os.path.exists(strain_json):
        with open(strain_json) as f:
            d = json.load(f)
        strain_keys = len(d)
        strain_occurrences = sum(len(v) for v in d.values())
    if os.path.exists(species_json):
        with open(species_json) as f:
            d = json.load(f)
        species_keys = len(d)
        species_occurrences = sum(len(v) for v in d.values())
    return {
        "summary": summary,
        "strain_keys": strain_keys,
        "strain_occurrences": strain_occurrences,
        "species_keys": species_keys,
        "species_occurrences": species_occurrences,
    }


def main(bundle_dir: str) -> int:
    if not os.path.isdir(bundle_dir):
        print(f"error: bundle dir not found: {bundle_dir}")
        return 1

    runs = ["shipped", "fix_from_scratch", "bug_4195edb"]
    manifest = {
        "bundle_dir": bundle_dir,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runs": {},
    }
    for r in runs:
        sub = os.path.join(bundle_dir, r)
        manifest["runs"][r] = {
            "files": hash_dir(sub),
            "metrics": {d: disease_metrics(sub, d) for d in ("IBD", "PD")},
        }

    manifest_path = os.path.join(bundle_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print(f"wrote {manifest_path}")

    # ---- diff_report.txt ----
    lines: list[str] = []
    lines.append(f"3-way evidence bundle: {bundle_dir}")
    lines.append(f"Built: {manifest['created_utc']}")
    lines.append("")
    lines.append("== Metrics by run/disease ==")
    for disease in ("IBD", "PD"):
        lines.append(f"\n--- {disease} ---")
        cols = [
            "run", "chi2", "P_Val",
            "Num_Increased_disease_Total", "Num_Increased_disease_Butyrate_Producers",
            "Num_Decreased_disease_Total", "Num_Decreased_disease_Butyrate_Producers",
            "strain_keys", "strain_occurrences", "species_keys", "species_occurrences",
        ]
        lines.append("  " + "  ".join(c.ljust(18) for c in cols))
        for r in runs:
            m = manifest["runs"][r]["metrics"][disease]
            s = m["summary"] or {}
            row = [
                r,
                f"{s.get('chi2', ''):.6f}" if s.get("chi2") not in (None, "") else "—",
                f"{s.get('P_Val', '')}" if s.get("P_Val") not in (None, "") else "—",
                str(s.get("Num_Increased_disease_Total", "—")),
                str(s.get("Num_Increased_disease_Butyrate_Producers", "—")),
                str(s.get("Num_Decreased_disease_Total", "—")),
                str(s.get("Num_Decreased_disease_Butyrate_Producers", "—")),
                str(m["strain_keys"]),
                str(m["strain_occurrences"]),
                str(m["species_keys"]),
                str(m["species_occurrences"]),
            ]
            lines.append("  " + "  ".join(c.ljust(18) for c in row))

    lines.append("\n== Byte-level file comparison ==")
    # Compare each run's files against shipped baseline
    shipped_hashes = {k: v["sha256"] for k, v in manifest["runs"]["shipped"]["files"].items()}
    for r in ["fix_from_scratch", "bug_4195edb"]:
        lines.append(f"\n--- {r} vs shipped ---")
        run_hashes = {k: v["sha256"] for k, v in manifest["runs"][r]["files"].items()}
        all_keys = sorted(set(shipped_hashes) | set(run_hashes))
        for name in all_keys:
            sh = shipped_hashes.get(name)
            rh = run_hashes.get(name)
            if sh is None:
                lines.append(f"  {name}: only in {r} (sha {rh[:12] if rh else '—'}...)")
            elif rh is None:
                lines.append(f"  {name}: missing in {r}")
            elif sh == rh:
                lines.append(f"  {name}: IDENTICAL")
            else:
                lines.append(f"  {name}: DIFFERS  shipped={sh[:12]}...  {r}={rh[:12]}...")

    report_path = os.path.join(bundle_dir, "diff_report.txt")
    with open(report_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {report_path}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: build_evidence_bundle_manifest.py <bundle_dir>")
        sys.exit(2)
    sys.exit(main(sys.argv[1]))
