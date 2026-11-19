# Reproducible scATAC-seq analysis with Snakemake

A practical tutorial for exploring chromatin accessibility in human peripheral blood mononuclear cells (PBMCs), using **R, Signac, Seurat and Snakemake**.

**[Start the step-by-step tutorial →](docs/tutorial.md)**

Follow one public dataset from processed files to quality-control plots, doublet scores, a low-dimensional representation, accessibility clusters and marker-based interpretation. Each step explains what the command does, which files it creates and what to check before continuing.

## Public dataset

[10x Genomics: 10k Human PBMCs, ATAC v2, Chromium Controller](https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard). Human; one healthy donor; **10,246 recovered nuclei**; **hg38/GRCh38**; Cell Ranger ATAC **2.1.0**. The dataset is freely available under **CC BY 4.0**. The exact processed-file URLs are recorded in [the download manifest](config/downloads.atlas.tsv).

The tutorial selects **1,000 cells with seed 42** to reduce downstream computation. This is a subset of real measurements. All four source files still need to be downloaded, and the full sparse H5 matrix is read before subsetting. Use `max_cells: 0` for the full capture when resources permit.

## Before you start

- Linux x86_64 is the reference platform. Linux/WSL2 is recommended; other platforms need their own environment validation.
- Install Git, curl, Python **3.12.14**, and Conda using the [setup instructions](docs/tutorial.md#1-prepare-your-computer).
- Budget **16 GB RAM and 20 GB free disk** for the tutorial. These are planning allowances, not a guarantee for every machine. See [dataset and resource details](docs/datasets.md).
- Basic familiarity with opening a terminal is enough to start. Commands below run from the repository root.

This repository is private: clone/open it while signed into the GitHub account that owns it or has been granted access.

## Run the tutorial

```bash
git clone https://github.com/anjalianjali1234/scatac-snakemake.git
cd scatac-snakemake
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock

# Download the four public files. Interrupted downloads can resume.
python scripts/download_atlas.py

# Generate metrics first; inspect them before removing cells.
snakemake atlas_qc --configfile config/atlas.yaml --cores 2 --sdm conda
```

The [walkthrough](docs/tutorial.md) explains input checksum verification and the next two steps:

1. Review `results/atlas/qc.png` and `qc_metrics.csv`. Record justified filtering choices in `config/atlas.yaml`, then run `atlas_lsi`.
2. Review `lsi_depth.png` and `depth_correlations.csv`. Record the dimensions to use, then run the complete workflow to generate `results/atlas/report.html`.

**Do not mark a decision reviewed just to bypass an error.** The two checkpoints make analysis choices visible and reproducible. See [the validation record](docs/validation.md) for the exact execution status and limitations.

## What you get

| Output in `results/atlas/` | What it tells you |
|---|---|
| `qc.png`, `qc_metrics.csv` | Fragment depth, TSS enrichment, nucleosome signal, vendor FRiP, blacklist overlap and doublet scores |
| `cell_filtering.csv` | Which cells passed each individual filter |
| `lsi_depth.png`, `depth_correlations.csv` | How strongly each LSI component tracks sequencing depth |
| `umap.png`, `cells.csv` | Accessibility clusters and their cell-level coordinates/metadata |
| `marker_activity.png` | Accessibility-derived marker evidence for provisional cell labels |
| `cluster_annotation_template.csv` | A table for documenting proposed labels and conflicting evidence |
| `atlas.rds` | The analyzed Seurat/Signac object for further inspection |
| `report.html`, `provenance.json`, `session_*.txt` | Plots, parameters, actual input hashes and software versions |

## Check your installation without downloading data

```bash
snakemake --cores 2
python tests/check_outputs.py
```

The default `demo` mode is only a small **synthetic smoke test**. It verifies workflow mechanics and is not the PBMC analysis. Open `results/demo/report.html` to inspect it. An optional R/Python agreement check is available with `snakemake smoke_test_r --cores 1 --sdm conda`.

## Reproducibility

The workflow uses Snakemake **9.27.2**, R **4.5.1**, Signac **1.16.0**, Seurat **5.3.0**, SeuratObject **5.2.0**, and scDblFinder **1.24.0**. These are fixed compatibility targets, not a claim to use the newest releases. See [the environment recipes](workflow/envs), [reproducibility notes](docs/reproducibility.md), and [validation record](docs/validation.md).

Raw input files and large analysis objects stay out of Git. Input identities, parameter choices, scripts and environment specifications are versioned. GitHub Actions runs the small-data checks; a passing CI run does not establish a biological result.

## Interpretation and scope

Accessibility is not RNA expression. Gene activity provides supporting evidence for a label, not a direct expression measurement. A UMAP is a visualization, not a statistical test. This single-donor dataset supports exploratory analysis; it cannot establish population-level differences between conditions. No cell identities are assigned automatically.

A separate [paired-donor pseudobulk extension](docs/datasets.md#paired-donor-extension) is available for correctly prepared replicated data. Motif analysis, trajectories, paired multiome integration and raw FASTQ processing are outside the main tutorial.

## Documentation

- [Step-by-step tutorial](docs/tutorial.md)
- [Public data, file formats and computing requirements](docs/datasets.md)
- [Reproducibility and troubleshooting](docs/reproducibility.md)
- [Validation record](docs/validation.md)

Method references: [versioned Signac workflow](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette), [scDblFinder ATAC guidance](https://plger.github.io/scDblFinder/articles/scATAC.html), and [Snakemake deployment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). Please cite the original data provider and methods when using this work.
