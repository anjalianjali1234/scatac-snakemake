import copy
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from common import load_samples, sha256, validate_paired_cells, verify_checksums
from demo import load_panel, summarize


class TeachingDataTests(unittest.TestCase):
    def test_known_quantities_and_reproducibility(self):
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "a", Path(tmp) / "b"
            result = summarize(ROOT / "tests/data/counts.csv", ROOT / "tests/data/cells.tsv", first)
            self.assertEqual(result["panel_counts"], [8, 4, 4, 7, 4, 5])
            self.assertEqual(result["total_counts"], 32)
            self.assertEqual(result["zero_entries"], 16)
            self.assertAlmostEqual(result["zero_fraction"], 16 / 30)
            summarize(ROOT / "tests/data/counts.csv", ROOT / "tests/data/cells.tsv", second)
            self.assertEqual({p.name: sha256(p) for p in first.iterdir()},
                             {p.name: sha256(p) for p in second.iterdir()})

    def test_reordered_barcodes_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cells.tsv"
            lines = (ROOT / "tests/data/cells.tsv").read_text().splitlines()
            path.write_text("\n".join([lines[0]] + list(reversed(lines[1:]))) + "\n")
            with self.assertRaisesRegex(ValueError, "exactly match"):
                load_panel(ROOT / "tests/data/counts.csv", path)

    def test_changed_input_not_silently_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "input.txt"; path.write_text("original\n")
            lock = Path(tmp) / "lock.json"; lock.write_text(json.dumps({str(path): sha256(path)}))
            verify_checksums([str(path)], lock)
            path.write_text("changed\n")
            with self.assertRaisesRegex(ValueError, "Checksum mismatch"):
                verify_checksums([str(path)], lock)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is needed for cross-language parity")
    def test_r_python_agree(self):
        with tempfile.TemporaryDirectory() as tmp:
            py, r = Path(tmp) / "py", Path(tmp) / "r"
            summarize(ROOT / "tests/data/counts.csv", ROOT / "tests/data/cells.tsv", py)
            subprocess.run(["Rscript", "--vanilla", str(ROOT / "workflow/scripts/module01.R"),
                            str(ROOT / "tests/data/counts.csv"), str(ROOT / "tests/data/cells.tsv"), str(r)], check=True, capture_output=True)
            for name in ("cell_summary.csv", "region_summary.csv"):
                with (py / name).open() as a, (r / name).open() as b:
                    pa, rb = list(csv.DictReader(a)), list(csv.DictReader(b))
                self.assertEqual(len(pa), len(rb))
                for x, y in zip(pa, rb):
                    self.assertEqual(x.keys(), y.keys())
                    for key in x:
                        try:
                            self.assertAlmostEqual(float(x[key]), float(y[key]))
                        except ValueError:
                            self.assertEqual(x[key], y[key])


class ReplicationTests(unittest.TestCase):
    def setUp(self):
        self.settings = {"celltype": "monocyte", "control": "control", "treatment": "treated", "min_cells": 2, "min_pairs": 3}
        self.samples, self.cells = [], []
        for donor in ("d1", "d2", "d3"):
            for condition in ("control", "treated"):
                sample = f"{donor}_{condition}"
                self.samples.append(dict(sample=sample, donor=donor, condition=condition))
                for i in range(2):
                    self.cells.append(dict(cell_id=f"{sample}_{i}", sample=sample, donor=donor, condition=condition, celltype="monocyte"))

    def test_independent_pairs_accepted(self):
        self.assertEqual(set(validate_paired_cells(self.cells, self.samples, self.settings).values()), {2})

    def test_missing_pair_rejected(self):
        cells = [r for r in self.cells if r["sample"] != "d3_treated"]
        with self.assertRaisesRegex(ValueError, "no donors silently dropped"):
            validate_paired_cells(cells, self.samples, self.settings)

    def test_many_cells_do_not_create_donors(self):
        samples = [r for r in self.samples if r["donor"] == "d1"]
        cells = [r for r in self.cells if r["donor"] == "d1"]
        with self.assertRaisesRegex(ValueError, "splitting cells cannot create replication"):
            validate_paired_cells(cells, samples, self.settings)

    def test_inconsistent_sample_donor_rejected(self):
        cells = copy.deepcopy(self.cells); cells[0]["donor"] = "d2"
        with self.assertRaisesRegex(ValueError, "donor mismatch"):
            validate_paired_cells(cells, self.samples, self.settings)

    def test_duplicate_cell_id_rejected(self):
        with self.assertRaisesRegex(ValueError, "globally unique"):
            validate_paired_cells(self.cells + [self.cells[0]], self.samples, self.settings)

    def test_genome_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, "mismatched genome"):
            load_samples(ROOT / "config/samples.demo.tsv", "demo", "hg38")


if __name__ == "__main__":
    unittest.main()
