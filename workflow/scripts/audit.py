import argparse
import importlib.metadata
import json
import platform
import subprocess
from pathlib import Path
from common import input_paths, load_samples, read_tsv, sha256, validate_paired_cells, validate_settings, verify_checksums


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--settings", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    config = json.loads(args.settings)
    validate_settings(config)
    samples = load_samples(config["samples"], config["mode"], config["genome"])
    records = verify_checksums(input_paths(config, samples) + [config["samples"]], config["checksums"])
    if config["mode"] == "pseudobulk":
        p = config["pseudobulk"]
        cells = read_tsv(p["cells"], ["cell_id", "sample", "donor", "condition", "celltype"])
        validate_paired_cells(cells, samples, p)
    source_files = sorted(p for folder in ("workflow", "schemas") for p in Path(folder).rglob("*")
                          if p.is_file() and "__pycache__" not in p.parts)
    source_files += [Path("requirements.lock")]
    result = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
    git_sha = result.stdout.strip() if result.returncode == 0 else None
    dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True)
    doc = {
        "mode": config["mode"], "synthetic": config["mode"] == "demo" or config.get("pseudobulk", {}).get("synthetic", False),
        "parameters": config, "inputs": records,
        "source": [{"path": str(p), "sha256": sha256(p)} for p in source_files],
        "python": platform.python_version(), "platform": platform.platform(),
        "snakemake": importlib.metadata.version("snakemake"),
        "git_commit": git_sha, "git_dirty": bool(dirty.stdout.strip()) if dirty.returncode == 0 else None,
        "checksum_origin": "Locally recorded content identity; not a claim of publisher-authenticated checksums"
    }
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n")
    print(f"Audited {len(records)} input files")


if __name__ == "__main__":
    main()
