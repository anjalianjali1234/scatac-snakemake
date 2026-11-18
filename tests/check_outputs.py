"""Integration assertions after Snakemake has completed the demo."""
import csv
import json
from pathlib import Path

out = Path("results/demo")
summary = json.loads((out / "summary.json").read_text())
assert summary["synthetic"] is True
assert summary["panel_counts"] == [8, 4, 4, 7, 4, 5]
assert summary["total_counts"] == 32 and summary["zero_entries"] == 16
with (out / "region_summary.csv").open() as handle:
    regions = list(csv.DictReader(handle))
assert regions[-1]["region"] == "unobserved_region"
assert float(regions[-1]["detection_fraction"]) == 0
manifest = json.loads((out / "provenance.json").read_text())
assert manifest["snakemake"] == "9.27.2"
assert len(manifest["inputs"]) == 3
assert all(len(r["sha256"]) == 64 for r in manifest["inputs"])
assert "Synthetic accessibility panel" in (out / "report.html").read_text()
assert all((out / p).is_file() for p in summary["artifacts"])
print("Demo output and provenance checks passed")
