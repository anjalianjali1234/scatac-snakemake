"""Deterministic, synthetic Module 1 exercise. Not a biological ATAC pipeline."""
import argparse
import csv
import html
import json
from pathlib import Path
from common import read_tsv


def load_panel(counts_file, cells_file):
    cells = read_tsv(cells_file, ["cell_id", "group", "donor", "sample"])
    ids = [row["cell_id"] for row in cells]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate cell IDs")
    with open(counts_file, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header[0] != "region" or header[1:] != ids:
            raise ValueError("Counts columns must exactly match metadata cell IDs and order")
        names, matrix = [], []
        for row in reader:
            if len(row) != len(header) or not row[0]:
                raise ValueError("Malformed count matrix row")
            if any(not x.isascii() or not x.isdecimal() for x in row[1:]):
                raise ValueError("Raw counts must be nonnegative integers")
            names.append(row[0])
            matrix.append([int(x) for x in row[1:]])
    if not names or len(names) != len(set(names)):
        raise ValueError("Region IDs must be unique and nonempty")
    return cells, names, matrix


def write_csv(path, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(counts_file, cells_file, output):
    cells, names, matrix = load_panel(counts_file, cells_file)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    totals = [sum(row[j] for row in matrix) for j in range(len(cells))]
    for j, row in enumerate(cells):
        row["panel_counts"] = totals[j]
        row["detected_regions"] = sum(x[j] > 0 for x in matrix)
    regions = [{"region": name, "total_count": sum(row),
                "detection_fraction": sum(x > 0 for x in row) / len(cells)}
               for name, row in zip(names, matrix)]
    write_csv(output / "cell_summary.csv", cells)
    write_csv(output / "region_summary.csv", regions)
    nzero = sum(x == 0 for row in matrix for x in row)
    summary = {"synthetic": True, "n_cells": len(cells), "n_regions": len(names),
               "panel_counts": totals, "total_counts": sum(totals), "zero_entries": nzero,
               "zero_fraction": nzero / (len(cells) * len(names)),
               "interpretation": "Synthetic panel counts only. No measured RNA, TF binding, causal effects or independent real donors.",
               "artifacts": ["cell_summary.csv", "region_summary.csv", "panel_counts.svg"]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    # An exact SVG bar chart keeps the smoke test dependency-free.
    width, height, base = 800, 410, 320
    maximum = max(totals) or 1
    bar_width = min(70, 600 / len(cells) * 0.65)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
           '<title id="title">Synthetic accessibility panel counts</title>',
           '<desc id="desc">Counts across five teaching regions, not whole-library depth.</desc>',
           '<rect width="800" height="410" fill="white"/>',
           '<g font-family="sans-serif" fill="#153542"><text x="70" y="35" font-size="23">Synthetic accessibility panel</text>',
           '<text x="70" y="62" font-size="14">Observed counts in the teaching panel</text>']
    for i, (cell, total) in enumerate(zip(cells, totals)):
        x = 85 + i * 600 / len(cells)
        h = total / maximum * 210
        color = "#347d9c" if cell["group"] == "A" else "#b17f3c"
        svg += [f'<rect x="{x:.2f}" y="{base-h:.2f}" width="{bar_width:.2f}" height="{h:.2f}" fill="{color}"/>',
                f'<text x="{x+bar_width/2:.2f}" y="{base-h-9:.2f}" text-anchor="middle">{total}</text>',
                f'<text x="{x+bar_width/2:.2f}" y="345" text-anchor="middle">{html.escape(cell["cell_id"])}</text>']
    svg += ['<text x="70" y="385" font-size="14">Invented data for the course; no biological results inferred.</text></g></svg>']
    (output / "panel_counts.svg").write_text("\n".join(svg) + "\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.counts, args.cells, args.output), sort_keys=True))
