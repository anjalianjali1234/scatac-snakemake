# A practical scATAC-seq tutorial: public PBMC data to interpretable clusters

[Repository overview](../README.md) · [Dataset details](datasets.md) · [Reproducibility](reproducibility.md) · [Validation](validation.md)

**Question:** which groups of nuclei have similar chromatin-accessibility patterns, and which broad immune-cell identities are consistent with those patterns?

This tutorial follows one analysis from start to finish. You will download a real public dataset, inspect data quality, construct a low-dimensional representation, visualize clusters and document marker evidence. Snakemake connects the steps so that inputs, parameters and outputs are explicit.

The main analysis uses a fixed random subset of 1,000 nuclei from 10x Genomics' human PBMC ATAC v2 dataset. It starts with processed files, so you do not need to run an alignment pipeline. No biological result is prescribed: inspect the actual figures and retain uncertainty.

## How to use this page

Run **bash** blocks in a terminal from the repository root. Run **R** blocks in the pinned R environment when the text says to inspect an existing object. Those R blocks are optional inspection commands; Snakemake already executes the main analysis scripts. YAML blocks describe configuration fields to edit, not terminal commands.

Allow roughly 2–3 hours of active reading and interpretation. Installation, downloads and fragment-file scans add elapsed time and depend strongly on the machine. The [validation record](validation.md) separates actual checks from resource estimates.

## 1. Prepare your computer

Use Linux x86_64 or WSL2 with about 16 GB RAM and 20 GB free disk. No GPU is needed. On a smaller machine, the synthetic installation check works with much less memory; reducing the biological subset does not remove the full download or H5-import requirement.

Install [Git](https://git-scm.com/downloads), curl, and [Miniforge](https://github.com/conda-forge/miniforge#install). Follow their official instructions for your operating system, then open a terminal where `conda` is available. Do not paste Linux installation commands into Windows Command Prompt; use a WSL terminal.

One way to obtain the exact Python version is:

```bash
conda create --name scatac-tools -c conda-forge python=3.12.14 pip -y
conda activate scatac-tools
```

Clone the repository while authenticated to an account with access:

```bash
git clone https://github.com/anjalianjali1234/scatac-snakemake.git
cd scatac-snakemake
python -m pip install -r requirements.lock
python --version
snakemake --version
conda --version
```

The Python version should be `3.12.14`, and Snakemake should be `9.27.2`. Snakemake is the workflow manager: it examines the requested output and runs the commands needed to create it. Python runs the workflow and utility scripts; R performs the biological analysis.

The `--sdm conda` option used below asks Snakemake to create and use the declared R environment automatically. It does not install R into your ordinary Python environment. Keep target names, such as `atlas_qc`, before that option.

**Optional installation check:**

```bash
snakemake --cores 2 --dry-run
snakemake --cores 2
python tests/check_outputs.py
```

This uses six fictional cells and five regions included in the repository. It should create `results/demo/report.html`. The known counts are software-test expectations, not public-data results. Return to the `atlas` commands below for the real analysis.

## 2. Understand the data you will use

The [10x dataset page](https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard) describes one healthy human donor, 10,246 recovered PBMC nuclei and Cell Ranger ATAC 2.1.0 processing. Data are available under CC BY 4.0; attribute 10x Genomics when sharing derived figures. The exact files match the [Signac 1.16 example](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette).

ATAC-seq identifies DNA regions that were accessible to the transposase. A promoter may be accessible in cells that use its gene, but accessibility alone does not measure the amount of RNA produced. A barcode identifies a captured partition; a doublet can contain material from more than one nucleus.

| File saved under `data/pbmc10k/` | Meaning | Why the workflow needs it |
|---|---|---|
| `filtered_peak_bc_matrix.h5` | Sparse peak-by-barcode counts | Starts the count matrix and the cell list |
| `singlecell.csv` | Vendor barcode-level QC metadata | Provides usable-fragment and in-peak-fragment counts |
| `fragments.tsv.gz` | Sorted, bgzip-compressed genomic intervals with barcodes | Supports TSS/nucleosome metrics and gene activity |
| `fragments.tsv.gz.tbi` | Tabix index | Allows efficient interval queries without loading all fragments into RAM |

All files must come from this same processing release and assembly. The matrix has peaks as rows and cell barcodes as columns. Sparse storage keeps only nonzero entries: a zero can reflect incomplete sampling, not proof that a locus is closed. BED-style fragment coordinates are zero-based, half-open; do not manually add or subtract bases in the supplied files.

## 3. Download and verify the inputs

```bash
python scripts/download_atlas.py
```

This reads `config/downloads.atlas.tsv`, downloads from the official 10x HTTPS host, and saves temporary `.partial` files until each download completes. Run the same command again after an interrupted transfer; completed files are left unchanged. The download is several GB, even though the analysis uses 1,000 cells.

The repository's input lock is `config/checksums.atlas.json`. When present, it contains SHA256 identities calculated from the actual downloaded files, including the sample sheet. Snakemake verifies the inputs before analysis. If you are deliberately adapting this project to different files, review their source and create a separate lock with `scripts/pin_inputs.py`; do not overwrite expected hashes simply to silence an error.

```bash
# Prints the planned jobs without executing the analysis.
snakemake atlas_qc --configfile config/atlas.yaml --cores 2 --dry-run
```

Check that inputs point to `data/pbmc10k/` and outputs to `results/atlas/`. A missing file error means a required download or input lock is absent. It is not a reason to create an empty placeholder.

The sample sheet records `pbmc10k`, `donor1`, `healthy` and `hg38`. These identify one capture from one donor; 1,000 selected cells are not 1,000 independent biological replicates.

## 4. Calculate quality-control metrics

```bash
snakemake atlas_qc --configfile config/atlas.yaml --cores 2 --sdm conda
```

The first run installs the R dependencies. Later runs reuse the environment. `--cores 2` limits workflow concurrency; most of this analysis is deliberately run with one numerical-library thread.

The rule imports the sparse count matrix, matches metadata by barcode, selects 1,000 barcodes with seed 42, retains standard nuclear chromosomes and attaches hg38 gene annotations from EnsDb.Hsapiens.v86. It then calculates QC metrics and scores possible doublets using scDblFinder's ATAC feature aggregation.

Open `results/atlas/qc.png` and `qc_metrics.csv`. The table contains one row per selected barcode before the explicit QC filtering step.

| Metric | Plain-language interpretation | What to inspect |
|---|---|---|
| `passed_filters` | Vendor count of usable fragments | Both a low-count tail and unusually high-depth cells |
| `TSS.enrichment` | Enrichment near transcription-start sites relative to surrounding background | Weak signal, outliers and its relationship to fragment depth |
| `nucleosome_signal` | Relative abundance of nucleosome-sized versus short fragments | A high tail that may indicate poorer nuclear signal |
| `frip_vendor` | `peak_region_fragments / passed_filters` | Low in-peak signal; it depends on the vendor's peak set |
| `blacklist_fraction` | Fraction of peak-matrix counts overlapping the supplied blacklist | Unusually high signal in problematic genomic regions |
| `doublet_score`, `doublet_class` | Evidence for mixed profiles from scDblFinder | Scores near the decision boundary and sensitivity to the assumed rate |

`nCount_peaks` is a sum over matrix peaks; it is not the same denominator as all usable fragments. The blacklist metric is also matrix-based, so describe its definition when reporting it. Do not compare differently defined metrics as if they were interchangeable.

The initial doublet prior is 0.08 for the original capture. It is an explicit starting assumption, not a measured truth, and it should not be divided by ten because only 1,000 cells were retained. A subset can alter doublet detection and miss rare populations. Prefer the full capture and a sensitivity check for a substantive biological analysis.

**What the important R commands do:** `Read10X_h5()` imports counts; `CreateChromatinAssay()` connects counts, genomic ranges, fragments and annotations; `TSSEnrichment()` and `NucleosomeSignal()` compute the metrics supported by the pinned Signac version. Do not update Signac alone and assume these APIs remain unchanged.

## 5. Choose and document filters

Read the distributions before deciding which cells to retain. The starting values in `config/atlas.yaml` are:

| Parameter | Initial value | Definition |
|---|---:|---|
| `min_fragments` | 1,000 | Minimum usable fragments |
| `max_fragments` | 100,000 | Maximum usable fragments |
| `min_tss` | 2 | Minimum TSS enrichment |
| `max_nucleosome` | 4 | Maximum nucleosome signal |
| `min_frip` | 0.20 | Minimum vendor FRiP |
| `max_blacklist_fraction` | 0.01 | Maximum matrix-based blacklist fraction |

These are not universal acceptance criteria. Ask whether the tails represent poor quality, plausible cell-type differences, or uncertainty. Compare retained fractions and marker evidence after changing an influential cutoff. The workflow also retains only cells called singlets; that classification remains imperfect.

Edit the `atlas.qc` section in `config/atlas.yaml`. Replace the empty `rationale` with your actual observations and set `reviewed: true` after inspection. Keep YAML indentation, use a decimal for fractions, and put the explanation in quotes. Record which plots you checked and why the retained cells are suitable for the intended analysis. A blank justification intentionally stops the next step.

```bash
snakemake atlas_lsi --configfile config/atlas.yaml --cores 2 --sdm conda
```

Open `cell_filtering.csv`: each column records whether a particular check passed; `retained` requires all checks. Missing QC values are treated as failures. Summarize failures separately because one cell may fail several checks. Very few retained cells should trigger investigation, not automatic threshold relaxation.

## 6. Understand TF-IDF and LSI

The retained peak matrix is sparse and strongly affected by sequencing depth. TF-IDF adjusts the representation so that common, broadly accessible regions do not dominate every comparison. Singular-value decomposition compresses that representation into LSI components: combinations of peaks that describe major patterns across cells.

The workflow calls `RunTFIDF()`, `FindTopFeatures()` and `RunSVD()`. Raw counts remain stored in the assay. TF-IDF values and LSI coordinates are representations for comparing cells, not raw counts for a differential-accessibility model.

Open `lsi_depth.png` and `depth_correlations.csv`. Each bar is a Spearman correlation between one LSI component and log10 usable-fragment depth. A strong association can indicate a technical influence, but depth and biology can be related. Do not remove a component solely because of its number.

The initial config lists components 2–30. Inspect the actual correlations before accepting or changing that choice. Edit `atlas.lsi.dimensions`, write your observation in `atlas.lsi.rationale`, and set `atlas.lsi.reviewed: true`. Record ambiguous decisions and compare an alternative component set if it changes the biological interpretation.

## 7. Build a graph, cluster cells and create the report

```bash
snakemake --configfile config/atlas.yaml --cores 2 --sdm conda
```

`FindNeighbors()` connects cells with similar selected LSI coordinates. `FindClusters()` groups that graph; resolution controls the granularity of the grouping. `RunUMAP()` produces a two-dimensional visualization of the same selected representation. UMAP is not used as the input for clustering.

Open `results/atlas/report.html` in a browser, or inspect `umap.png` and `embedding_and_activity.pdf` directly. Keep the HTML in the same folder as its linked images and tables. A PDF gives a convenient portable view of the embedding and marker panels.

Check whether clusters correspond strongly to fragment depth, whether isolated small groups have unusual QC values, and whether the groupings remain interpretable under a reasonable change in resolution. Cluster numbers are arbitrary labels. Distance between disconnected UMAP islands is not a calibrated measure of biological difference; a smooth shape is not evidence of developmental time.

## 8. Interpret marker evidence and propose cell labels

The workflow counts fragments overlapping gene bodies and nearby promoters with `GeneActivity()`, stores the result in an assay called **ACTIVITY**, and normalizes it for visualization. This remains accessibility-derived evidence, not an RNA measurement.

`marker_activity.png` displays relative marker activity across clusters. In this dot plot, color summarizes scaled average activity and dot size summarizes the fraction with detected activity. The color scale is relative within a gene; it is not a cross-gene comparison of transcript abundance. Sparse detection, gene length and annotation choice influence these scores.

| Candidate broad identity | Markers to investigate together | A useful caution |
|---|---|---|
| T cells | CD3D, CD3E | Check supporting loci before attempting T-cell subtypes |
| B cells | MS4A1, CD79A | One weak marker is insufficient |
| Monocytes | LYZ, CD14, FCGR3A | Subtype boundaries need more than one gene |
| NK / cytotoxic lymphocytes | NKG7, GNLY | Cytotoxic T cells can share this program |
| Dendritic-like cells | FCER1A, CST3 | Rare groups and overlapping myeloid markers need caution |
| Platelet / megakaryocytic signal | PPBP | Inspect contamination and low-complexity profiles before labeling |

These are candidates to investigate, not the expected answer or a complete marker dictionary. Check multiple positive markers, contradictory evidence, sample quality and local accessibility. Keep a cluster `Unassigned` when evidence is insufficient. Reference mapping can provide additional evidence in an extended analysis, but it is not executed here.

Copy the generated template before editing so that rerunning the workflow does not overwrite your interpretation:

```bash
cp results/atlas/cluster_annotation_template.csv config/cluster_annotations.csv
```

Fill `proposed_celltype`, `supporting_markers`, `conflicting_evidence` and `confidence` for the actual cluster IDs in that file. Do not reuse numeric labels from a different run. Record the configuration/Git revision and retained cell count alongside your annotation table.

Optional **R inspection**, after the workflow has completed and with the pinned R environment active:

```r
library(Seurat)
library(Signac)
pbmc <- readRDS("results/atlas/atlas.rds")
DefaultAssay(pbmc) <- "ACTIVITY"
FeaturePlot(pbmc, features = c("CD3D", "MS4A1", "LYZ", "NKG7"), reduction = "umap")
table(Idents(pbmc))
head(pbmc[[]])
```

Optional **R application of your reviewed labels**:

```r
labels <- read.csv("config/cluster_annotations.csv", colClasses = "character")
stopifnot(!anyDuplicated(labels$cluster))
cluster_id <- as.character(Idents(pbmc))
stopifnot(setequal(unique(cluster_id), labels$cluster))
stopifnot(all(nzchar(labels$proposed_celltype)))
label_map <- setNames(labels$proposed_celltype, labels$cluster)
pbmc$provisional_celltype <- unname(label_map[cluster_id])
pdf("results/atlas/provisional_celltypes.pdf", width = 10, height = 7)
print(DimPlot(pbmc, group.by = "provisional_celltype", label = TRUE))
dev.off()
saveRDS(pbmc, "results/atlas/atlas_provisional_labels.rds")
```

These optional outputs are generated by the inspection code, not by the default Snakemake target. Commit the reviewed annotation table and a script containing your final labeling code if you use this extension. Manual labels are hypotheses, not independently validated cell identities.

## 9. Verify reproducibility

```bash
snakemake --configfile config/atlas.yaml --cores 2 --sdm conda --dry-run
```

With unchanged inputs, code, environment and configuration, the completed workflow should report no pending work. If a job is scheduled, inspect the reported reason before deleting files or changing timestamps.

- `provenance.json` records the effective config, input hashes/sizes, source hashes, platform and available Git revision.
- `session_qc.txt`, `session_lsi.txt` and `session_embedding.txt` record the R packages actually loaded.
- `logs/atlas/` stores command output and errors.
- `benchmarks/atlas/` records elapsed time and resource measurements for each major rule.
- `cell_filtering.csv` and the two config rationales make exclusions and component choices auditable.

Changing the QC justification or thresholds should rerun the filtering/LSI step and its downstream outputs, while an LSI dimension choice should rerun the embedding step. Full input identities are still audited. A fixed seed improves repeatability; it does not guarantee bitwise identical numerical results across every platform and library build.

When sharing, include your config, input lock, session information, annotation evidence and appropriate data/method attribution. Keep original data and large R objects outside ordinary Git history.

## 10. What this analysis can and cannot conclude

You can describe accessibility variation, QC decisions, clusters and marker evidence within this capture. You cannot use these cells as independent donors in a condition comparison. The separate pseudobulk extension requires independent biological replicates and common raw peak counts.

Accessibility is not expression; motif accessibility is not TF binding; proximity or correlation is not proof that a region regulates a gene. This tutorial does not perform motif inference, peak-to-gene linking, trajectory analysis or a condition test. Adding those analyses requires additional assumptions and validation.

A useful final account states the dataset/assembly, selected and retained cells, peak definition, QC criteria, doublet method, LSI dimensions, clustering settings, label evidence and limitations. Report observations from your run, including uncertainty, rather than matching the appearance of another person's UMAP.

## Troubleshooting

| Symptom | What to check |
|---|---|
| GitHub returns 404 | The repository is private; check that you are signed into an account with access |
| `snakemake` or `conda` not found | Activate the intended terminal environment; confirm the version commands in step 1 |
| Download interrupted | Run the downloader again; retain `.partial` files for resumption |
| Checksum mismatch | Check the URL, release, file completeness and sample sheet; do not replace the expected hash automatically |
| Missing fragment index | Ensure the `.tbi` belongs to this exact fragment file and is adjacent to it |
| QC/LSI review error | Inspect the named diagnostics, edit the relevant rationale, then mark that section reviewed |
| Missing/very weak TSS signal | Check assembly, chromosome naming, annotation and fragment quality before changing cutoffs |
| Too few retained cells | Inspect individual failures, doublet scores and choices; do not silently fit an invalid LSI rank |
| Out of memory | Close other jobs, use a capable machine, or reduce `max_cells`; full-file import/download requirements remain |
| Cell labels disagree across markers | Keep a broader or unassigned label and document the conflicting evidence |

For installation logs and environment-lock details, see [reproducibility](reproducibility.md). For what has actually been executed, see [validation](validation.md).
