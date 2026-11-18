"""Record actual file identities after the analyst reviews their provenance."""
import argparse
import json
import sys
from pathlib import Path
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
from common import input_paths, load_samples, sha256, validate_settings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    validate_settings(cfg)
    rows = load_samples(cfg["samples"], cfg["mode"], cfg["genome"])
    paths = sorted(set(input_paths(cfg, rows) + [cfg["samples"]]))
    output = Path(cfg["checksums"])
    if output.exists():
        raise SystemExit(f"Refusing to overwrite {output}. Review and version an input change explicitly.")
    checksums = {p: sha256(p) for p in paths}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        json.dump(checksums, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Recorded {len(paths)} real SHA256 values in {output}; this does not authenticate the publisher.")


if __name__ == "__main__":
    main()
