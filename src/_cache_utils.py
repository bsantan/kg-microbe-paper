"""Cache validation helpers for find_microbes_strain / find_microbes_species.

Why this exists: a regression in search_strains (commit 4195edb) was hidden
because the cached *_microbes_strain.json was loaded by default and there was
no way to detect that the cached output came from a different code version.
The helpers here add a sibling .manifest.json that records the input
fingerprint, output sha256s, and a hand-managed code_version_marker; loading
the cache is gated on all three matching.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from typing import Iterable

MANIFEST_SCHEMA_VERSION = 2


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


_SOURCE_FILE_SHA_CACHE: dict[str, tuple[float, int, str]] = {}


def source_file_sha256(path: str) -> str:
    """sha256 a large source file, cached per process by (path, mtime, size).

    Used to fingerprint inputs that drive cache validity but are too big to
    re-hash on every call (e.g. ``merged-kg_edges_ncbitaxon.tsv``, ~45 MB).
    Recomputes if mtime or size changes. Returns "missing" if the file
    is absent so the fingerprint changes meaningfully rather than crashing
    callers that catch and fall back.
    """
    try:
        st = os.stat(path)
    except OSError:
        return "missing"
    cached = _SOURCE_FILE_SHA_CACHE.get(path)
    if cached and cached[0] == st.st_mtime and cached[1] == st.st_size:
        return cached[2]
    sha = sha256_file(path)
    _SOURCE_FILE_SHA_CACHE[path] = (st.st_mtime, st.st_size, sha)
    return sha


def get_git_head(repo_dir: str | None = None) -> tuple[str, bool]:
    """Return (short_sha, dirty). On any failure returns ('unknown', False)."""
    cwd = repo_dir or os.getcwd()
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=cwd, stderr=subprocess.DEVNULL
        ).decode().strip()
        status = subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=cwd, stderr=subprocess.DEVNULL
        ).decode()
        return sha, bool(status.strip())
    except Exception:
        return "unknown", False


def compute_input_fingerprint(
    all_taxa,
    feature_type: str,
    ranks_df,
    extra: dict | None = None,
    edge_table_path: str | None = None,
) -> str:
    """Hash the inputs that, if changed, should invalidate the cache.

    Includes: feature_type, taxa count + sorted-taxa sha, ranks_df shape + a
    sha of the sorted (taxon, rank) pairs that intersect all_taxa, and (when
    ``edge_table_path`` is supplied) a sha of the underlying NCBITaxon
    subclass edge file that drives descendant enumeration. Output sha256 is
    independent (validated separately on load).
    """
    taxa_list = sorted(str(t) for t in (all_taxa or []))
    taxa_sha = _sha256_bytes(("\n".join(taxa_list)).encode("utf-8"))

    try:
        sub = ranks_df[ranks_df["NCBITaxon_ID"].isin(taxa_list)]
        rank_pairs = sub[["NCBITaxon_ID", "Rank"]].astype(str).values.tolist()
        rank_pairs.sort()
        ranks_sha = _sha256_bytes(json.dumps(rank_pairs, separators=(",", ":")).encode("utf-8"))
        ranks_shape = list(getattr(ranks_df, "shape", (0, 0)))
    except Exception:
        ranks_sha = ""
        ranks_shape = [0, 0]

    payload = {
        "feature_type": feature_type,
        "n_taxa": len(taxa_list),
        "taxa_sha256": taxa_sha,
        "ranks_shape": ranks_shape,
        "ranks_sample_sha256": ranks_sha,
    }
    if edge_table_path is not None:
        payload["edge_table_sha256"] = source_file_sha256(edge_table_path)
    if extra:
        payload["extra"] = extra
    return _sha256_bytes(json.dumps(payload, sort_keys=True).encode("utf-8"))


def atomic_write_json(path: str, obj) -> None:
    """Write JSON atomically: tempfile in same dir, fsync, os.replace."""
    directory = os.path.dirname(os.path.abspath(path)) or "."
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp.", suffix=".json", dir=directory)
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, indent=4)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def write_manifest(
    manifest_path: str,
    kind: str,
    feature_type: str,
    output_paths: Iterable[str],
    input_fingerprint: str,
    code_version_marker: str,
    repo_dir: str | None = None,
    extra: dict | None = None,
) -> dict:
    """Build a manifest from the freshly-written outputs and write it last."""
    git_head, git_dirty = get_git_head(repo_dir)
    outputs = []
    for p in output_paths:
        outputs.append({
            "name": os.path.basename(p),
            "sha256": sha256_file(p),
            "size": os.path.getsize(p),
        })
    payload = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "kind": kind,
        "feature_type": feature_type,
        "outputs": outputs,
        "input_fingerprint": input_fingerprint,
        "code_version_marker": code_version_marker,
        "git_head": git_head,
        "git_dirty": git_dirty,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    if extra:
        payload["extra"] = extra
    atomic_write_json(manifest_path, payload)
    return payload


def read_manifest(manifest_path: str) -> dict | None:
    try:
        with open(manifest_path, "r") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return None


def cache_is_valid(
    output_paths: Iterable[str],
    manifest_path: str,
    expected_fingerprint: str,
    expected_code_version_marker: str,
    feature_type: str,
) -> tuple[bool, str]:
    """Return (valid, reason). Reason is a short token suitable for logging."""
    output_paths = list(output_paths)
    for p in output_paths:
        if not os.path.exists(p):
            return False, f"missing_output:{os.path.basename(p)}"
    if not os.path.exists(manifest_path):
        return False, "no_manifest"
    manifest = read_manifest(manifest_path)
    if manifest is None:
        return False, "manifest_unparseable"
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        return False, "schema_version_mismatch"
    if manifest.get("feature_type") != feature_type:
        return False, "feature_type_mismatch"
    if manifest.get("input_fingerprint") != expected_fingerprint:
        return False, "input_fingerprint_mismatch"
    if manifest.get("code_version_marker") != expected_code_version_marker:
        return False, "code_version_marker_mismatch"
    declared = {o["name"]: o for o in manifest.get("outputs", [])}
    for p in output_paths:
        name = os.path.basename(p)
        if name not in declared:
            return False, f"output_missing_from_manifest:{name}"
        if sha256_file(p) != declared[name].get("sha256"):
            return False, f"output_sha256_mismatch:{name}"
    return True, "valid"
