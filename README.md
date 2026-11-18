# scatac-snakemake

A 20-week single-cell ATAC-seq course with a Snakemake workflow for repeatable teaching exercises and a documented path to real-data analysis.

Start with **[the complete course](docs/course.md)**: weekly objectives, lessons, exercises, quiz answer keys, assessments, four capstone briefs, and Module 1 taught in full. The primary analysis route is R, Signac and Seurat. Advanced ArchR, Python, motifs, trajectories and multiome material is in the course; those topics are not all automated in this first workflow.

## Choose a route

| Route | Input | Output | What it establishes |
|---|---|---|---|
| `demo` (default) | Six fictional cells, five regions; included in Git | CSV summaries, a chart, HTML report, provenance | Installation and workflow mechanics; no biological conclusions |
| `atlas` | A processed 10x PBMC ATAC library, hg38 | QC, doublet scores, LSI diagnostics, UMAP/clusters, gene activity, R objects | An exploratory atlas, with analyst-reviewed thresholds and dimensions |
| `pseudobulk` | Common-peak raw counts and reviewed cell/sample metadata | Donor-condition sums, paired edgeR model, effect sizes, FDR and diagnostics | A replicated two-condition contrast within one annotated cell type |

**Execution status:** see [validation notes](docs/validation.md). Environment recipes and code are distinct from evidence of running a biological analysis. No public PBMC dataset or capstone result is included as a claimed completed analysis.

## Run the small offline example

Use Linux, macOS, or WSL2, Git, and Python **3.12.14**. Run these commands from the repository root after cloning. Installing packages requires the internet; the default analysis itself uses only included files.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock
snakemake --cores 2 --dry-run
snakemake --cores 2 --printshellcmds
python tests/check_outputs.py
```

Open `results/demo/report.html`. Its table, exact chart, parameters and SHA256 input identities are generated locally. The expected synthetic totals are `8, 4, 4, 7, 4, 5`, with 16 zero entries out of 30. These are known teaching values, not measured PBMC results.

Run the same Module 1 exercise in **base R 4.5.1**:

```bash
# With Rscript already on PATH:
snakemake --cores 1 module1_r
# Or with a working Conda installation:
snakemake module1_r --cores 1 --sdm conda
```

The R exercise creates CSV summaries, a PDF chart and `sessionInfo.txt`. `python -m unittest discover -s tests -v` compares the R/Python summaries when R is available, and tests checksum failures, barcode mismatches and replication safeguards.

## Real PBMC accessibility analysis

Read [datasets and input contracts](docs/datasets.md) and [reproducibility](docs/reproducibility.md) first. The exact 10x 2.1.0 file links are in `config/downloads.atlas.tsv`. Downloads are a separate, explicit step and can occupy several GB.

```bash
python scripts/download_atlas.py
python scripts/pin_inputs.py --config config/atlas.yaml
snakemake atlas_qc --configfile config/atlas.yaml --cores 2 --sdm conda
```

Inspect `results/atlas/qc.pdf` and `qc_metrics.csv`. Edit the thresholds, the capture-level expected doublet rate, and the written justification in `config/atlas.yaml`; set `atlas.qc.reviewed: true`. The supplied numbers are starting values to inspect, not universal cutoffs.

```bash
snakemake atlas_lsi --configfile config/atlas.yaml --cores 2 --sdm conda
```

Inspect `depth_correlations.csv` and `lsi_depth.pdf`. Decide which components to retain using technical associations and biology. Set the dimensions and rationale, then `atlas.lsi.reviewed: true`.

```bash
snakemake --configfile config/atlas.yaml --cores 2 --sdm conda
```

The final report links QC, embedding/activity plots, cell metadata and `atlas.rds`. Clusters are deliberately not assigned biological identities automatically. Use multiple marker loci, gene activity and a reviewed reference as taught in the course; record uncertainty. Gene activity is named `ACTIVITY`, because it is not an RNA measurement.

This route handles one capture. It uses EnsDb.Hsapiens.v86 (hg38) for annotation, which need not match the exact annotation release used by the original alignment reference. TSS enrichment and activity depend on that choice. Whole-file input hashes are checked on each audit; changing settings currently reruns dependent stages, including QC.

## Paired donor comparison

`config/samples.paired.tsv` identifies the eight GSE178428 libraries for four matched donors, control versus LPS at 6 h, without Golgi inhibitor. The workflow starts **after** joint peak definition, raw recounting, QC/doublet handling and cell annotation. It does not pretend that separately called peak matrices can be concatenated.

Prepare the documented files under `data/paired/`, review the selected lineage and group support, and edit `config/pseudobulk.yaml`:

```bash
python scripts/pin_inputs.py --config config/pseudobulk.yaml
snakemake --configfile config/pseudobulk.yaml --cores 2 --sdm conda --dry-run
snakemake --configfile config/pseudobulk.yaml --cores 2 --sdm conda
```

The design is `~ donor + condition`; positive logFC means greater relative accessibility in the configured treatment. Technical libraries can contribute to the same donor-condition sum. Individual cells never become independent biological replicates. Missing donor pairs, duplicate barcodes, conflicting sample metadata and assembly mismatches fail instead of being silently repaired.

This is a deliberately bounded paired design. Additional covariates, unpaired studies, repeated time courses and interactions require a reviewed model change, rank/confounding checks and an updated report. A condition completely confounded with processing batch cannot be rescued by a model formula.

## Files and versions

| Path | Purpose |
|---|---|
| `docs/course.md` | The complete learning pathway and Module 1 |
| `workflow/Snakefile`, `workflow/rules/` | Dependency graph and execution rules |
| `workflow/scripts/` | Analysis, validation and reporting code |
| `config/` | Parameters, sample sheets, download sources and input locks |
| `workflow/envs/` | R/Conda environment recipes and any actual solved locks |
| `requirements.lock` | Exact Python package versions recorded from the installed environment |
| `tests/` | Tiny synthetic inputs and scientific-contract tests |
| `.github/workflows/ci.yaml` | Small-data checks on GitHub; no multi-GB downloads |

The engine is **Snakemake 9.27.2**. The atlas targets **R 4.5.1, Signac 1.16.0, Seurat 5.3.0 and SeuratObject 5.2.0**. The paired model targets **edgeR 4.8.2**. These are deliberate teaching targets, not a claim to use the newest APIs. Do not upgrade Signac independently and assume old QC or chromVAR examples still work.

Input data, large objects, results and credentials are excluded from Git. Keep a durable, access-controlled copy of the actual inputs as well as their hashes. A checksum records identity; it cannot recover a deleted file. Commit approved configs and locks with the code.

## Interpretation boundaries

- Accessibility is not expression; motif accessibility is not demonstrated TF binding.
- Co-accessibility, peak-to-gene links and pseudotime are associations, not proof of regulatory causation or observed developmental time.
- A one-donor atlas cannot establish population-level condition effects.
- Motif families, enzyme sequence bias, sparse coverage, peak selection and footprinting assumptions can change an interpretation.
- Review effect sizes, uncertainty, donor consistency and the testing universe, alongside adjusted P values.

The [course](docs/course.md) develops these distinctions and the four capstone projects in detail. The optional FASTQ lab remains a separately provisioned advanced exercise; the default workflow begins with processed matrices/fragments.

## Troubleshooting

- **Checksum mismatch:** stop and investigate file identity. Do not regenerate the lock just to make a run pass.
- **Missing data:** download/prepare the requested inputs; a dry-run does not create them.
- **QC or LSI review stop:** inspect the named diagnostics and record the scientific decision in the config.
- **Conda solve failure:** preserve the log and check the supported platform and pins. A YAML recipe is not a complete solved environment.
- **Out of memory:** use the synthetic route or reduce retained cells on a capable host. Reading the full H5 and scanning full fragments still consumes resources.
- **Too few donor pairs/cells:** revisit the study and lineage support. Do not create pseudo-replicates by partitioning cells.

Official references: [Snakemake deployment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html), [versioned Signac tutorial](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette), [scDblFinder](https://plger.github.io/scDblFinder/), [edgeR 3.22 release](https://bioconductor.org/packages/3.22/bioc/html/edgeR.html). Full readings and dataset citations are in the course.
