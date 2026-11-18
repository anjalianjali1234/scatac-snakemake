"""Shared input contracts. No dependencies outside the Python standard library."""
import csv
import hashlib
import json
import re
from pathlib import Path

IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path, required):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required) <= set(reader.fieldnames):
            raise ValueError(f"{path}: missing required columns {sorted(set(required) - set(reader.fieldnames or []))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: table is empty")
    if any(None in row or any(row[k] is None or row[k] == "" for k in required) for row in rows):
        raise ValueError(f"{path}: malformed or empty required fields")
    return rows


def load_samples(path, mode, genome):
    required = ["sample", "donor", "condition", "genome"]
    if mode == "atlas":
        required += ["counts_h5", "metadata_csv", "fragments", "fragments_index"]
    rows = read_tsv(path, required)
    ids = [r["sample"] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("Sample identifiers must be unique")
    for row in rows:
        for key in ("sample", "donor", "condition"):
            if not IDENTIFIER.fullmatch(row[key]):
                raise ValueError(f"Unsafe {key} identifier: {row[key]!r}; use letters, numbers, dot, underscore or hyphen")
        if row["genome"] != genome:
            raise ValueError("Mixed or mismatched genome assemblies; renaming chromosomes is not liftOver")
    if mode == "atlas" and len(rows) != 1:
        raise ValueError("Atlas mode currently handles one library. Recount a common peak universe before joint analysis")
    return rows


def validate_settings(config):
    mode = config["mode"]
    if mode == "atlas":
        a = config["atlas"]
        if not isinstance(a["max_cells"], int) or a["max_cells"] < 0:
            raise ValueError("atlas.max_cells must be a nonnegative integer; zero means all")
        if not 0 < a.get("expected_doublet_rate", 0) < 1:
            raise ValueError("Specify an expected doublet rate between zero and one for the original capture")
        q = a["qc"]
        for key in ("reviewed", "rationale", "min_fragments", "max_fragments", "min_tss", "max_nucleosome", "min_frip", "max_blacklist_fraction"):
            if key not in q:
                raise ValueError(f"Missing QC setting: {key}")
        if type(q["reviewed"]) is not bool or not isinstance(q["rationale"], str):
            raise ValueError("QC reviewed must be boolean and rationale must be text")
        if not 0 <= q["min_frip"] <= 1 or not 0 <= q["min_fragments"] < q["max_fragments"]:
            raise ValueError("Invalid FRiP or fragment bounds")
        if q["min_tss"] < 0 or q["max_nucleosome"] <= 0:
            raise ValueError("Invalid TSS/nucleosome bounds")
        if not 0 <= q["max_blacklist_fraction"] <= 1:
            raise ValueError("Blacklist fraction must be between zero and one")
        lsi = a["lsi"]
        for key in ("reviewed", "rationale", "dimensions", "n_components", "resolution"):
            if key not in lsi:
                raise ValueError(f"Missing LSI setting: {key}")
        dims = lsi["dimensions"]
        if type(lsi["reviewed"]) is not bool or not isinstance(lsi["rationale"], str):
            raise ValueError("LSI reviewed must be boolean and rationale must be text")
        if not dims or any(type(x) is not int or x < 1 or x > lsi["n_components"] for x in dims) or len(set(dims)) != len(dims):
            raise ValueError("LSI dimensions must be distinct, positive and within n_components")
        if lsi["resolution"] <= 0:
            raise ValueError("Clustering resolution must be positive")
    if mode == "pseudobulk":
        p = config["pseudobulk"]
        if p["control"] == p["treatment"] or p["min_pairs"] < 3 or p["min_cells"] < 1:
            raise ValueError("Require distinct conditions, >=3 independent donor pairs and >=1 cell per group")
        if not p["rationale"].strip():
            raise ValueError("Document the contrast, lineage selection and cell-support rationale")


def input_paths(config, samples):
    if config["mode"] == "demo":
        return [config["demo"]["counts"], config["demo"]["cells"]]
    if config["mode"] == "atlas":
        return [samples[0][key] for key in ("counts_h5", "metadata_csv", "fragments", "fragments_index")]
    p = config["pseudobulk"]
    return [p[key] for key in ("matrix", "features", "barcodes", "cells")]


def verify_checksums(paths, checksum_file):
    expected = json.loads(Path(checksum_file).read_text())
    records = []
    for name in paths:
        if name not in expected or not re.fullmatch(r"[0-9a-f]{64}", expected[name]):
            raise ValueError(f"Missing real SHA256 for {name}; use scripts/pin_inputs.py after inspecting provenance")
        actual = sha256(name)
        if actual != expected[name]:
            raise ValueError(f"Checksum mismatch: {name}; do not silently refresh the input lock")
        records.append({"path": name, "sha256": actual, "bytes": Path(name).stat().st_size})
    return records


def validate_paired_cells(cells, samples, settings):
    """Check sample attribution and all donor pairs BEFORE selecting a lineage."""
    sample_map = {row["sample"]: row for row in samples}
    seen = set()
    for row in cells:
        cell = row["cell_id"]
        if cell in seen:
            raise ValueError("Cell IDs must be globally unique; preserve a barcode-to-sample mapping")
        seen.add(cell)
        if row["sample"] not in sample_map:
            raise ValueError(f"Unknown sample for {cell}")
        sample = sample_map[row["sample"]]
        for key in ("donor", "condition"):
            if row[key] != sample[key]:
                raise ValueError(f"Cell/sample {key} mismatch for {cell}")
    conditions = {settings["control"], settings["treatment"]}
    if {r["condition"] for r in samples} != conditions:
        raise ValueError("Sample sheet must contain exactly the chosen two conditions")
    donors = {r["donor"] for r in samples}
    if len(donors) < settings["min_pairs"]:
        raise ValueError("Not enough independent donors; splitting cells cannot create replication")
    counts = {(d, c): 0 for d in donors for c in conditions}
    for row in cells:
        if row["celltype"] == settings["celltype"]:
            counts[(row["donor"], row["condition"])] += 1
    deficient = {f"{d}/{c}": n for (d, c), n in counts.items() if n < settings["min_cells"]}
    if deficient:
        raise ValueError(f"Missing or insufficient lineage support; no donors silently dropped: {deficient}")
    return counts
