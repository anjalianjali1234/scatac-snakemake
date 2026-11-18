"""Optional real edgeR execution against a synthetic, paired experimental design."""
import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class PseudobulkIntegration(unittest.TestCase):
    def test_paired_workflow_with_known_synthetic_shift(self):
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not installed")
        check = subprocess.run(["Rscript", "--vanilla", "-e", 'quit(status=if(requireNamespace("edgeR",quietly=TRUE) && requireNamespace("jsonlite",quietly=TRUE) && requireNamespace("Matrix",quietly=TRUE)) 0 else 1)'], capture_output=True)
        if check.returncode:
            self.skipTest("Optional edgeR/Matrix/jsonlite environment not installed")
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for name in ("workflow", "schemas", "config"):
                shutil.copytree(ROOT / name, work / name, ignore=shutil.ignore_patterns("__pycache__"))
            shutil.copy(ROOT / "requirements.lock", work / "requirements.lock")
            subprocess.run(["Rscript", "--vanilla", str(ROOT / "tests/make_paired_fixture.R"), "data/synthetic"], cwd=work, check=True, capture_output=True)
            settings = {"mode": "pseudobulk", "genome": "hg19", "samples": "data/synthetic/samples.tsv", "checksums": "data/synthetic/checksums.json", "seed": 42,
                        "pseudobulk": {"matrix": "data/synthetic/counts.mtx", "features": "data/synthetic/features.tsv", "barcodes": "data/synthetic/barcodes.tsv", "cells": "data/synthetic/cells.tsv",
                        "celltype": "synthetic_lineage", "control": "control", "treatment": "treated", "min_cells": 6, "min_pairs": 4, "synthetic": True,
                        "rationale": "Fictional four-pair integration fixture; the first 20 features have a planted fourfold shift."}}
            # JSON is also a valid Snakemake config file and avoids another serialization dependency.
            (work / "config/test.json").write_text(json.dumps(settings))
            subprocess.run([sys.executable, str(ROOT / "scripts/pin_inputs.py"), "--config", "config/test.json"], cwd=work, check=True, capture_output=True)
            run = subprocess.run([sys.executable, "-m", "snakemake", "--cores", "1", "--configfile", "config/test.json"], cwd=work, capture_output=True, text=True)
            if run.returncode:
                logs = "\n".join(p.read_text() for p in (work / "logs").rglob("*.log"))
                self.fail(run.stdout + run.stderr + logs)
            out = work / "results/pseudobulk"
            summary = json.loads((out / "summary.json").read_text())
            self.assertTrue(summary["synthetic"])
            self.assertTrue(json.loads((out / "provenance.json").read_text())["synthetic"])
            self.assertEqual(summary["independent_donor_pairs"], 4)
            self.assertEqual(summary["pseudobulk_groups"], 8)
            self.assertEqual(summary["cells"], 48)
            with (out / "groups.tsv").open() as h:
                groups = list(csv.DictReader(h, delimiter="\t"))
            self.assertEqual([int(g["n_cells"]) for g in groups], [6] * 8)
            with (out / "differential_accessibility.tsv").open() as h:
                result = list(csv.DictReader(h, delimiter="\t"))
            planted = [r for r in result if int(r["feature_id"].removeprefix("peak")) <= 20]
            self.assertEqual(len(planted), 20)
            self.assertGreater(sum(float(r["logFC"]) > 1 for r in planted), 15)
            self.assertTrue(all(0 <= float(r["FDR"]) <= 1 for r in result))
            self.assertTrue((out / "report.html").is_file())


if __name__ == "__main__":
    unittest.main()
