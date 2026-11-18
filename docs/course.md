# Single cell ATAC sequencing analysis course

Navigation: [Course guide](#course-guide) · [Weekly curriculum](#complete-weekly-curriculum) · [Quiz keys](#quiz-answer-keys-with-reasoning) · [Four projects](#four-capstone-projects) · [Module 1 in full](#module-1-in-full)

## Course guide

This course takes a beginner from chromatin biology and basic programming to an independently documented single-cell ATAC-seq analysis. The main route is R, Signac and Seurat. ArchR is introduced for disk-backed projects and an alternative analysis workflow; SnapATAC2 provides a focused Python comparison. The goal is to explain and defend analytical decisions, not merely reproduce a UMAP.

**Duration:** 20 weeks, approximately 6–8 active study hours per week, about 140 hours overall. Four subsequent capstones require an additional 58–78 active hours if all are completed. Installation, downloads and unattended computation can add elapsed time. There are no paid course subscriptions required; computing services, if used, may cost money.

**Documentation checked:** 18 November 2026. Dataset identity, assay, assembly and available file types were checked against official dataset records and software tutorials. This is an authored curriculum, not a completed biological analysis. No real dataset was processed and no biological results are claimed. Module 1 uses explicitly synthetic data. Its R script is intended to run in base R but has not been executed in R in this authoring environment. Later weekly labs are guided workflow specifications with command examples, not complete scripts runnable from an empty session. Their full lessons and scripts are delivered one module at a time.

### Entry requirements and exit outcomes

You should recognize DNA, RNA, genes and cells, be able to create a folder, and understand fractions and a graph's axes. No previous genomics analysis or command-line experience is required. Calculus is unnecessary. Weeks 1–3 contain the preparation; they are part of the 20 weeks, not an extra prerequisite course.

By the end you should be able to:

1. Explain which molecules and events ATAC-seq measures and distinguish accessibility from expression and binding.
2. Audit samples, donors, genome assemblies, barcodes, fragment conventions and matrix dimensions before analysis.
3. Build a reproducible Signac workflow with justified QC, LSI components, neighborhoods, clusters and cell annotations.
4. Construct a common peak universe and perform donor-level differential accessibility with a defensible design.
5. Prioritize candidate motifs, regulators and regulatory regions while describing uncertainty and alternative explanations.
6. Compare separate-modality integration with experimentally paired multiome analysis.
7. Critique trajectories, computational resource use and apparent biological findings.
8. Produce a runnable report with provenance, versions, sensitivity analyses and appropriately limited conclusions.

### Weekly learning routine

Use four sessions: 90 minutes of lessons and reading; 150 minutes of guided lab; 90 minutes of exercises; and 90 minutes for the quiz, assignment and troubleshooting. That is seven hours. Six-hour weeks shorten the lab and reading; eight-hour weeks add a supervised sensitivity analysis. Each lab follows **predict → run → inspect → explain → change one decision**. Quiz retries are encouraged after writing a correction in your own words.

Keep one course project with `data_raw/`, `data_derived/`, `scripts/`, `figures/`, `reports/` and `metadata/`. Store input files unchanged. Keep a weekly decision log with the question, input version, parameter, reason, observed consequence and unresolved issue. Do not put large raw data into Git.

### Preparation map

| Skill | Where it is learned | Readiness task |
|---|---|---|
| R objects, vectors, indexing, tables, plots, missing values and functions | Weeks 1–2 | Read a CSV, check column types, calculate a summary, save a plot |
| Statistics: variability, sampling, independence, effect size, uncertainty, testing and FDR | Week 2, revisited in 11–13 | Identify the experimental unit and explain why a small p-value need not imply a large effect |
| Shell: paths, working directories, pipes, compressed text and safe scripts | Week 3 | Inspect the first records of a compressed file without loading the entire file |
| Python: variables, lists, dictionaries, loops, file IO and environments | Week 3, applied in 18 | Parse a small TSV using `csv` and count unique barcode strings |
| FASTQ, BAM, BED, GTF/GFF, fragments, MTX/HDF5 and sparse matrices | Weeks 3–5 | State which files provide sequence, intervals, annotation, counts and provenance |

An optional 4–6-hour warm-up before Week 1 can cover the R and shell setup sections in R1–R3. It is unnecessary if you are comfortable learning them during the course.

### Software and version policy

The teaching environment deliberately follows a versioned Signac workflow. Current documentation has changed: the Signac 1.17 release family deprecates `TSSEnrichment()` and `NucleosomeSignal()` in favor of `ATACqc()` and removes `RunChromVAR()`. Do not combine an old tutorial with newly installed defaults and assume compatibility. The live Seurat documentation now displays version 6.0.0; this course's primary API target remains the documented Seurat 5 workflow. See R7–R10 and R16.

| Component | Course target or recorded requirement | Purpose |
|---|---|---|
| R | 4.5.1 | Main teaching environment; matches the archived Signac tutorial's R version |
| Bioconductor | 3.22 on R 4.5 | Explicit package repository generation; record all resolved dependency versions |
| Signac | 1.16.0 | Versioned primary tutorials, including the chromVAR wrapper |
| Seurat and SeuratObject | 5.3.0 and 5.2.0 | Versions recorded in the archived Signac PBMC tutorial |
| chromVAR | 1.32.0, Bioconductor 3.22 | Bias-corrected deviations across motif-associated peaks |
| edgeR, scDblFinder, GenomicRanges, rtracklayer, Rsamtools | Bioconductor 3.22-compatible releases, exact resolved versions saved at installation | Pseudobulk inference, ATAC-aware doublets, intervals and indexed files |
| Genome and gene annotation resources | hg38 for datasets A/B; hg19 for C/D; exact FASTA/annotation release recorded per project | Avoid mixing assemblies or silently changing TSS definitions |
| ArchR | 1.0.3 teaching branch, separate project using its official release lockfile and matching R | Arrow files, doublets, iterative LSI and resource comparison; do not merge its dependencies into the Signac project |
| Python and SnapATAC2 | Python 3.12; SnapATAC2 2.10.0 | Backed fragment processing and a spectral-embedding alternative |
| MACS | MACS3 3.1.0 in its own Python environment for the Signac peak-calling exercise | Record the executable path and command; the Signac argument can still be named `macs2.path` |
| Raw-data optional lab | Cell Ranger ATAC 2.2 documentation; original dataset A was processed with 2.1.0 | Reprocessing is a new analysis, not an exact reproduction of the deposited output |
| bgzip/tabix | HTSlib 1.25 documentation | Position-sorted, block-compressed genomic interval files |
| Editor and reporting | RStudio Desktop or another R editor; Quarto optional; Git recommended | No particular IDE version is required |

These are **documentation targets, not a tested complete dependency lockfile**. The R environment must pass installation checks before Week 4. `renv` pins packages after successful installation but does not install R or operating-system libraries. On Windows the main R route is possible, but use Linux/WSL2 or a Linux server for ArchR, MACS and the optional raw-data pipeline. ArchR does not support native Windows. No GPU is required.

Illustrative setup procedure for a new, isolated R 4.5.1 project; run only after R is installed. Package installation can require system compilers and several GB of downloads. If a dependency fails, resolve that package first and preserve the error log instead of repeatedly upgrading everything.

```r
# Installation recipe, not executed here. Use a fresh course project.
stopifnot(getRversion() == numeric_version("4.5.1"))
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("renv")
renv::init(bare = TRUE)
install.packages(c("BiocManager", "remotes"))
BiocManager::install(version = "3.22", ask = FALSE)
renv::settings$bioconductor.version("3.22")

remotes::install_version("SeuratObject", "5.2.0", upgrade = "never")
remotes::install_version("Seurat", "5.3.0", upgrade = "never")
remotes::install_version("Signac", "1.16.0", upgrade = "never")
BiocManager::install(
  c("GenomicRanges", "rtracklayer", "Rsamtools", "edgeR",
    "scDblFinder", "chromVAR", "SingleCellExperiment"),
  ask = FALSE, update = FALSE
)
library(Signac)
library(Seurat)
stopifnot(as.character(packageVersion("Signac")) == "1.16.0")
stopifnot(as.character(packageVersion("Seurat")) == "5.3.0")
stopifnot(as.character(packageVersion("SeuratObject")) == "5.2.0")
BiocManager::valid()  # inspect and resolve mismatched Bioconductor packages
renv::snapshot(type = "all", prompt = FALSE)
writeLines(capture.output(sessionInfo()), "sessionInfo_install.txt")
```

Add genome packages and motif databases when needed, recording their versions. For the pinned PBMC tutorials, `EnsDb.Hsapiens.v86` is an hg38 annotation option and `EnsDb.Hsapiens.v75` is an hg19 option. These are not automatically identical to the annotation used by the original alignment reference. State that distinction when comparing TSS enrichment and gene activity. Verify the reference header rather than inferring an assembly from a filename.

### Datasets and computing plan

All biological datasets below are human. A/B are vendor demonstration datasets, not replicated condition-comparison cohorts. C is the primary condition-comparison dataset. D contains developmental structure useful for a cautious trajectory exercise.

**A — Introductory PBMC ATAC atlas.** [10k Human PBMCs ATAC v2 Chromium Controller](https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard), with exact downloads in the [versioned Signac PBMC tutorial](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette). Organism: *Homo sapiens*. Assay: single-nucleus ATAC, one healthy donor. Assembly: GRCh38/hg38, reference package GRCh38-2020-A in the tutorial. Approximately 10,246 recovered nuclei and about 165,000 retained peaks at the tutorial's initial object stage; your filtering can change these numbers. Required: filtered peak-by-barcode H5, `singlecell.csv`, `fragments.tsv.gz`, its `.tbi`, and a matching annotation resource. Its well-documented workflow and familiar immune lineages make it suitable for learning. It cannot establish a population-level disease or treatment effect.

The exact prefix is `10k_pbmc_ATACv2_nextgem_Chromium_Controller`; the tutorial links four files ending in `_filtered_peak_bc_matrix.h5`, `_singlecell.csv`, `_fragments.tsv.gz` and `_fragments.tsv.gz.tbi` under `samples/cell-atac/2.1.0/`. Start with the matrix, metadata and a reproducible 1,000–2,000-cell subset. Keep the original indexed fragments on disk; subsetting the Seurat object does not shrink that file.

**B — Paired PBMC RNA–ATAC multiome.** [10x dataset description](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0) and [Signac's exact ARC 1.0.0 downloads](https://stuartlab.org/signac/articles/pbmc_multiomic). Organism: *Homo sapiens*. Assay: RNA and ATAC measured in the same nuclei from one donor. Assembly: GRCh38/hg38. Use the archived `pbmc_granulocyte_sorted_10k` ARC **1.0.0** files used by Signac: approximately 11,909 input barcodes, 108,377 ATAC features and 36,601 RNA features in that example. Required: `_filtered_feature_bc_matrix.h5`, `_atac_fragments.tsv.gz` and `.tbi`; obtain QC metadata and reference details from the same processing version when available. The current 10x page describes ARC **2.0.0** reprocessing and reports 11,898 cells. Do not mix versions or present the two cell counts as conflicting measurements of the same processed object. This dataset supports paired WNN and peak–RNA association exercises, not a replicated treatment test.

**C — Donor-matched PBMC stimulation.** [GSE178428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178428). Organism: *Homo sapiens*. Assay: Bio-Rad SureCell scATAC-seq, **hg19/GRCh37**. The study contains 42 libraries and approximately 67,581 stimulated/control ATAC profiles described in the associated study. Use only the following eight libraries for the main comparison: LPS versus DMSO control at 6 h, without Golgi inhibitor, four matched donors. The exact retained cell count must be calculated from the metadata and your QC; it is not prescribed.

| Donor | Control at 6 h | LPS at 6 h |
|---|---|---|
| 1 | [GSM5390937](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390937) | [GSM5390957](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390957) |
| 2 | [GSM5390938](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390938) | [GSM5390958](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390958) |
| 3 | [GSM5390939](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390939) | [GSM5390959](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390959) |
| 4 | [GSM5390940](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390940) | [GSM5390960](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5390960) |

Required: each selected sample's fragments, the series files `GSE178428_PBMCs_stim_peaks.bed.gz` and `GSE178428_PBMCs_stim_scATACseq_cellMeta.txt.gz`, a manually verified sample sheet, and hg19 annotation. The series lists a 6.7 GB fragment archive for all samples; the two series-level small files are each about 1.6 MB. Individual selected files are much smaller: donor 1 control at 6 h is listed as 46.4 MB. Download the selected sample files, not the full archive, for the course. Recompress with bgzip and index if necessary. Inspect the actual barcode schema; capture IDs must be joined correctly to sample and donor. Never blindly strip barcode suffixes. The GEO cell table must not be assumed to contain trustworthy cell-type labels; annotate and validate your chosen lineage. The associated study's computational RNA–ATAC pairing is not experimentally paired multiome data.

The reported peak-construction method and a BED file's coordinate convention must be checked independently. Inspect widths and endpoints; retain any discrepancy in the provenance log rather than silently changing coordinates to match a reported width. The deposited peaks are useful for practice; a balanced common peak set rebuilt from the eight selected samples is preferable for the main project if computing allows.

**D — Hematopoietic differentiation and alternative workflows.** [GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785), [Signac CD34 trajectory tutorial](https://stuartlab.org/signac/1.16.0/articles/monocle), and the smaller [ArchR Hematopoiesis tutorial download](https://www.archrproject.com/articles/Articles/tutorial.html). Organism: *Homo sapiens*. Assay: scATAC-seq; assembly hg19. The Signac example selects 9,913 cells from `CD34_Progenitors_Rep1`. Required: the CD34 matrix, peak table, cell-barcode metadata, `GSM3722029_CD34_Progenitors_Rep1_fragments.tsv.gz` and a locally created index. The source tutorial warns that fragments are inside the large series archive. Prefer ArchR's approximately 0.5 GB teaching subset, obtained with `getTutorialData("Hematopoiesis")`, for the resource-limited route. Its three fragment files cover BMMC, CD34-enriched BMMC and PBMC samples; it is not a matched three-donor condition study. Record the downloaded subset's cell count rather than claiming it equals the full Signac CD34 example.

**Tiny practice route.** Signac's bundled `atac_small` object contains 100 cells and 323 features, hg19. Use it for object inspection and command mechanics. The package also supplies small example fragments via `system.file("extdata", "fragments.tsv.gz", package="Signac")`; do not assume they cover all cells/features in the object. This is not adequate for atlas discovery, realistic QC thresholds, footprints or donor-level inference. Module 1 needs neither Signac nor a biological download.

**Planning estimates, not measured benchmarks.** Assume a recent 4–8-core laptop or workstation with SSD storage. RAM is peak working memory, disk includes working copies, packages and temporary files; estimates depend strongly on retained peaks, concurrency and method. Times below exclude downloading and installation.

| Route | RAM target | Free disk allowance | Approximate elapsed computation |
|---|---:|---:|---|
| Module 1 synthetic lab | 2 GB | <100 MB | Seconds |
| Bundled tiny example | 4–8 GB | 5–10 GB including packages | Seconds to minutes |
| A, 1,000–2,000 cells after sparse import | 8–16 GB | 15–30 GB including the full fragments | 15–90 minutes for core processing |
| A, full dataset | 16–32 GB | 25–50 GB | 1–4 hours; motifs/links may add 1–4 hours |
| B, 2,000 paired cells | 16 GB | 20–40 GB | 30 minutes–3 hours |
| B, full paired dataset | 32–64 GB | 40–80 GB | 2–8 hours including selected motif/link calculations |
| C, eight libraries, one lineage at a time | 16–32 GB; 64 GB comfortable for joint preprocessing | 20–50 GB | 2–8 hours preprocessing/recounting; minutes for a small pseudobulk model |
| D, ArchR teaching subset | 8–16 GB | 10–20 GB | 30 minutes–3 hours |
| Optional full A FASTQ processing | 64–128 GB, 8–16 cores | 250–500 GB | Roughly 8–24+ hours |

On an 8 GB machine use the tiny object for mechanics, a 500–1,000-cell real-data subset for embeddings, one chromosome for regional recounting, and one lineage for downstream analysis. For C retain all four donor pairs: reduce cells or peaks per sample, never reduce replication and pretend the design is equivalent. A chromosome-only analysis has a chromosome-only testing universe. Large-H5 import may itself exceed laptop memory; use a larger machine once to save a sparse subset or use backed processing. The course does not supply an already downloaded biological subset.

### Free resource list

Read the specified sections, not every linked manual cover to cover. All links are sources or author-maintained documentation, not paid course requirements.

| ID | Resource and use |
|---|---|
| R1 | [R for Data Science, second edition](https://r4ds.hadley.nz/): workflow, import, transformation, visualization and joins |
| R2 | [R's introductory manual](https://cran.r-project.org/doc/manuals/r-release/R-intro.html): vectors, factors, matrices, indexing and functions; use your installed help for version-specific behavior |
| R3 | [Software Carpentry Unix Shell](https://swcarpentry.github.io/shell-novice/): paths, files, pipes and scripts |
| R4 | [Software Carpentry Python](https://swcarpentry.github.io/python-novice-inflammation/): variables, loops, functions and errors |
| R5 | [OpenIntro Statistics](https://www.openintro.org/book/os/): free PDF, sampling, study design, uncertainty and inference |
| R6 | [OpenStax gene regulation](https://openstax.org/books/biology-2e/pages/16-3-eukaryotic-epigenetic-gene-regulation) and [original ATAC-seq paper, free manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC3959825/): biology and assay rationale |
| R7 | [Signac 1.16 PBMC workflow](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette): primary analysis reference |
| R8 | [Signac installation](https://stuartlab.org/signac/articles/install) and [current changelog](https://stuartlab.org/signac/news/): version-aware setup |
| R9 | [Seurat installation](https://satijalab.org/seurat/articles/install) and [changelog](https://satijalab.org/seurat/news/): verify API changes |
| R10 | [Signac motif tutorial](https://stuartlab.org/signac/1.16.0/articles/motif_vignette), [chromVAR 3.22 package](https://bioconductor.org/packages/3.22/bioc/html/chromVAR.html) and [original chromVAR paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5623146/): motif statistics and deviations |
| R11 | [UCSC file formats](https://genome.ucsc.edu/FAQ/FAQformat.html), [10x fragments](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/outputs/fragments-file) and [tabix manual](https://www.htslib.org/doc/tabix.html): coordinates, fragment schema and indexing |
| R12 | [Cell Ranger ATAC count](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/running-pipelines/single-library-analysis-count): optional raw pipeline |
| R13 | [Signac peak calling](https://stuartlab.org/signac/1.16.0/articles/peak_calling), [FeatureMatrix](https://stuartlab.org/signac/reference/featurematrix), [MACS installation](https://macs3-project.github.io/MACS/docs/INSTALL.html): feature construction and recounting |
| R14 | [scDblFinder ATAC vignette](https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scATAC.html): use ATAC-specific feature aggregation, not RNA defaults |
| R15 | [edgeR package and user's guide](https://bioconductor.org/packages/3.22/bioc/html/edgeR.html) and [Squair et al. on biological replication](https://pmc.ncbi.nlm.nih.gov/articles/PMC8479118/): the latter benchmarks RNA, while its independence principle informs the ATAC design |
| R16 | [Signac paired multiome](https://stuartlab.org/signac/articles/pbmc_multiomic), [Seurat WNN](https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis) and [LinkPeaks reference](https://stuartlab.org/signac/reference/linkpeaks): paired modalities, separate preprocessing and candidate-link statistics |
| R17 | [Signac ATAC integration](https://stuartlab.org/signac/1.16.0/articles/integrate_atac): common features and integration diagnostics |
| R18 | [Signac Cicero](https://stuartlab.org/signac/1.16.0/articles/cicero) and [footprinting](https://stuartlab.org/signac/1.16.0/articles/footprint): hypothesis generation and limitations |
| R19 | [Signac trajectories](https://stuartlab.org/signac/1.16.0/articles/monocle): HSPC example and documented input files |
| R20 | [ArchR tutorial](https://www.archrproject.com/articles/Articles/tutorial.html), [installation](https://www.archrproject.com/) and [release environments](https://www.archrproject.com/bookdown/using-renv-to-manage-dependencies.html): isolate dependencies and record the release |
| R21 | [SnapATAC2 installation](https://snapatac2.scverse.org/install.html), [release notes](https://snapatac2.scverse.org/changelog.html) and [PBMC workflow](https://snapatac2.scverse.org/tutorials/pbmc.html): version 2.10 requires Python 3.12 or newer and uses `pp.import_fragments`, not the removed `pp.import_data` |
| R22 | [renv introduction](https://rstudio.github.io/renv/articles/renv.html): snapshot, restore and limitations |
| R23 | [GSE178428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178428) and [associated original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9534481/): experimental design and provenance |

### Assessment and progression

The 20-week course grade is 20% quizzes, 50% weekly assignments, 15% the Week 10 analysis audit, and 15% the Week 20 report/design defense. Every quiz has five questions: give 1 point for a correct answer and 1 point for its reasoning, 10 points total. The separate answer-key section below supports self-study; attempt questions before reading it.

Weekly assignment rubric, 100 points: computational correctness 30; interpretation 25; justified choices and sensitivity checks 20; reproducibility/provenance 15; clear communication 10. A passing assignment earns at least 70 and has no unresolved critical error. Critical errors include mixing genome assemblies, treating cells as independent donors, testing an integrated/TF-IDF matrix as raw counts, or claiming causation from a correlation. Correct these and resubmit; they are teaching checkpoints rather than permanent penalties.

Progress gates: after Week 3, independently read and validate a small data file; after Week 6, defend QC and barcode integrity; after Week 10, reproduce an annotated atlas checkpoint; after Week 12, explain the experimental unit and model; after Week 16, explain paired versus unpaired integration; after Week 20, defend a complete analysis proposal and report. A null result or an unassigned cell population can receive full marks when the analysis supports it.

The four capstones are graded separately using their project rubrics. Completing all four constitutes the full portfolio; completing only one should be described as completion of that project, not of the full portfolio.

## Complete weekly curriculum

The following are complete module specifications: objectives, prerequisite, lessons, guided lab, three exercises, five quiz prompts, pitfalls, assignment and study budget. Module 1 is taught in full later in this workbook. Subsequent full lessons are supplied on request. **Command examples here are workflow guides:** objects such as `pbmc`, `peaks`, `sce`, `design` and `dims_keep` must first be created and checked in the relevant lab; these snippets are not standalone scripts. Do not fabricate intermediate results to make a downstream command run.

### Week 1 — Chromatin biology and the first data table

**Prerequisite:** basic DNA/RNA knowledge; no programming prerequisite. **Objectives:** explain nucleosomes, promoters, enhancers and transcription factors; describe what Tn5-accessibility measurements mean; distinguish cells from donors; create a simple R table and plot.

**Short lessons:** DNA packaging influences which regions are accessible to proteins. ATAC-seq samples accessible DNA through transposition and sequencing; a peak is evidence of enriched accessibility, not proof of enhancer function. B and T lymphocytes illustrate how cells with broadly shared genomes can maintain different regulatory programs.

**Guided lab:** complete the fully specified synthetic lab below. Use `data.frame()` to store fictional cell metadata, `matrix()` to represent region counts, `colSums()` to summarize observed counts, `barplot()` to inspect them, and `write.csv()` to save a reproducible deliverable. Read dataset A's description without downloading its large files.

**Exercises:** (1) Guided: retrieve one named region and calculate its detection rate. (2) Partly guided: compare two cell groups after dividing counts by their observed panel totals. (3) Independent: write three supported and three unsupported conclusions from the toy table.

**Quiz:** Q1 What molecule does ATAC-seq sequence? Q2 Is an accessible promoter proof of high RNA expression? Q3 What does a zero in a sparse peak matrix mean? Q4 Are 1,000 cells from one donor 1,000 independent biological replicates? Q5 Does a distal accessible region identify its target gene?

**Mistakes and troubleshooting:** confusing a region with a gene; interpreting a barcode as a guaranteed single cell; assuming zero means permanently closed; confusing a named R object with a quoted string. Use `str()` and `dim()` before trying to interpret an object.

**Assignment and completion:** submit the R script, two CSV summaries, one labeled plot and a 300–500-word interpretation. Pass when the script runs from a clean session, row/column orientation is correct, and all four main inference limitations are stated. **Time:** 7 h: lessons 1.5, lab 2, exercises 1.5, quiz 0.5, assignment 1.5. **Read:** R2 introductory vectors/matrices; R6 assay overview; dataset A description.

### Week 2 — R and statistical thinking

**Prerequisite:** Week 1 table, vector and assay concepts. **Objectives:** manipulate metadata, handle missing values, use factors, distinguish variability from uncertainty, and identify sampling and experimental units.

**Short lessons:** a mean summarizes a distribution but hides differences among cells and donors. A confidence interval describes uncertainty under a model; it is not the range containing 95% of cells. Sampling more cells from one donor improves characterization of that donor but does not create more donors.

**Guided lab:** extend the toy metadata with fictional donor and condition columns, clearly labeling them as simulated. Inspect `str()` and `table(donor, condition)`. Use `aggregate()` for donor summaries, `factor()` for categorical variables, `is.na()` to detect missing values, and `set.seed()` for a repeatable sampling exercise. Compare resampling cells with resampling donor summaries; do not interpret the toy exercise as a real treatment test.

**Exercises:** (1) Guided: compute mean, median and standard deviation with one missing value handled explicitly. (2) Partly guided: join a sample table and detect an intentionally unmatched sample. (3) Independent: design a balanced six-donor paired stimulation experiment and explain the unit of inference.

**Quiz:** Q1 When can the median be more informative than the mean? Q2 Why can a missing QC value differ from zero? Q3 What is a biological replicate here? Q4 Does a small p-value guarantee a useful effect? Q5 Why set a random seed?

**Mistakes and troubleshooting:** R silently recycles incompatible-length vectors in some operations; sample order can change in joins; `na.rm=TRUE` can conceal systematic missingness. Check row counts, unique IDs and missingness before and after every join.

**Assignment and completion:** a metadata-cleaning script, missingness report, donor-level summary and 250-word study design. Pass when joins preserve intended sample identity and the distinction between cells and donors is explicit. **Time:** 7 h. **Read:** R1 data transformation/joins; R2 factors; R5 study design and sampling. **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 3 — Shell, Python and genomic coordinates

**Prerequisite:** Week 2 data tables. **Objectives:** navigate files safely, parse a small TSV in Python, explain FASTQ/BAM/BED/GTF/HDF5, and convert between interval conventions without off-by-one errors.

**Short lessons:** a reference assembly is the coordinate system for a genome; changing `1` to `chr1` changes naming, not assembly. BED uses zero-based, half-open intervals, whereas R `GRanges` uses one-based, closed intervals. A fragment record describes a sequenced DNA fragment and a barcode, while a peak matrix summarizes counts over selected regions.

**Guided lab:** create a tiny synthetic BED and fragment TSV. Inspect them with `pwd`, `ls`, `head` and `wc -l`; use `gzip -cd` for compressed text. Parse with Python's `csv.reader(..., delimiter="\t")`, keep barcodes as strings, and check widths as `end-start`. In R import a three-column BED with `makeGRangesFromDataFrame(..., starts.in.df.are.0based=TRUE)`; compare widths before/after conversion. Use `tabix -p bed` only on sorted bgzip data.

**Exercises:** (1) Guided: convert BED `[100,150)` into a 1-based interval. (2) Partly guided: identify a mismatched chromosome naming style and explain what it does not fix. (3) Independent: write a small Python validator rejecting negative starts, end≤start and empty barcode fields.

**Quiz:** Q1 What is the width of BED `[100,150)`? Q2 What does a FASTQ quality string describe? Q3 Does changing `chr1` to `1` perform liftOver? Q4 Why is an indexed bgzip file useful? Q5 Should a barcode such as `00123-1` be numeric?

**Mistakes and troubleshooting:** double-shifting coordinates, shell commands entered at the R prompt, treating gzip as bgzip, or querying a stale index. Confirm your prompt, file type and coordinate contract before analysis.

**Assignment and completion:** BED/GRanges conversion table, Python validation script, file-format glossary and shell transcript. Pass when lengths are conserved and malformed records are reported. **Time:** 8 h. **Read:** R3 navigation/pipes; R4 loops/errors; R11 BED and fragments. **Budget:** 2 h lessons, 3 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 4 — Import real data and establish provenance

**Prerequisite:** Weeks 1–3 and a working versioned R environment. **Objectives:** load sparse counts, link metadata by barcode, create a chromatin assay and validate genome annotation.

**Short lessons:** zero entries dominate single-cell ATAC matrices, so sparse storage is essential. Matrix rows are genomic regions; metadata rows describe cells. A Seurat object stores a fragment-file reference, so moving the file can break downstream calculations even if the object still loads.

**Guided lab:** download dataset A's processed files from R7. Use `Read10X_h5()` and `read.csv(..., row.names=1)`; confirm `all(colnames(counts) %in% rownames(metadata))` and reorder metadata with explicit indexing. Construct `CreateChromatinAssay()` and `CreateSeuratObject()`. Attach matching gene annotations after checking seqnames and assembly. Inspect `dim()`, `class()`, `object.size()` and `Fragments()`. Save the imported checkpoint and a seeded cell subset without converting the full matrix to dense form.

**Exercises:** (1) Guided: report cell/peak counts and sparse storage size. (2) Partly guided: fix shuffled metadata safely. (3) Independent: create a provenance manifest containing source URL, processing version, assembly, download date, checksum and local path.

**Quiz:** Q1 Why is `as.matrix()` dangerous here? Q2 What identifies a matrix column? Q3 Can a fragment path be discarded after import? Q4 What must agree between fragments and annotation? Q5 Why save input-file checksums?

**Mistakes and troubleshooting:** confusing raw and filtered H5 files; silently dropping barcode suffixes; attaching hg19 annotations to hg38 peaks. Use assertions and inspect mismatches; do not weaken validation until the object happens to load.

**Assignment and completion:** saved sparse object, provenance manifest and dimension/barcode audit. Pass with zero unexplained barcode mismatches and documented assembly. **Time:** 7 h. **Read:** R7 import section; R11 fragments; R22 project isolation. **Budget:** 1 h lessons, 3 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 5 — From reads to fragments and peaks

**Prerequisite:** Week 4 input audit. **Objectives:** describe alignment, barcode correction, deduplication, Tn5 adjustment, cell calling and peak calling; explain what information processed files have already lost.

**Short lessons:** paired reads locate the ends of a DNA fragment. PCR can create several sequenced copies of one molecule, so duplicate handling must respect the cell barcode and protocol. The same genomic coordinates in different cells are not automatically PCR duplicates. A peak caller detects enrichment in pooled signal; it does not reliably call a separate full peak set from every sparse cell.

**Guided lab:** audit the vendor web summary and fragment schema rather than rerun FASTQs. Query a locus with `tabix`; compare unique fragment records with their support-count field. Recount a small region list using `FeatureMatrix(fragments=..., features=..., cells=...)`. Trace which alignment and cell-calling decisions remain fixed in the provided fragments. Sketch the separate optional FASTQ workflow at the end of this curriculum.

**Exercises:** (1) Guided: classify duplicates within versus between barcodes. (2) Partly guided: compare counts from one fragment overlap and from endpoint counts, documenting the convention. (3) Independent: audit a fragment header and pipeline report for build, software and filtering provenance.

**Quiz:** Q1 Does the support-count field necessarily represent independent molecules? Q2 Why must deduplication be barcode-aware? Q3 Should deposited Tn5-corrected fragments be shifted again? Q4 Why pool signal for peak calling? Q5 Can a peak matrix recover off-peak fragments?

**Mistakes and troubleshooting:** multiplying unique-fragment counts by PCR support; applying +4/−5 a second time; treating all fragment formats as exactly five columns. Newer formats may include a strand column; inspect the documented version.

**Assignment and completion:** a read-to-matrix process diagram, provenance audit and a small recounted matrix with its counting rule. Pass when fragment versus insertion versus read-pair counting are distinguished. **Time:** 6 h, excluding optional FASTQ computation. **Read:** R11–R13. **Budget:** 1.5 h lessons, 2 lab, 1 exercises, 0.5 quiz, 1 report.

### Week 6 — Quality control and doublets

**Prerequisite:** Week 5 count definitions and Week 4 annotations/fragments. **Objectives:** compute and explain fragment depth, TSS enrichment, nucleosome signal and FRiP; detect likely doublets; defend dataset-specific filters.

**Short lessons:** depth measures evidence quantity, while TSS enrichment measures one aspect of signal concentration near annotated starts. Nucleosome signal compares fragment-length classes and depends on the implementation. FRiP depends on the peak universe and denominator; changing peaks can change FRiP without changing the cells. High depth can indicate a doublet but can also occur in a good singlet.

**Guided lab:** in the pinned environment run `TSSEnrichment()` and `NucleosomeSignal()`. Use vendor fragment/QC metadata where appropriate; distinguish `nCount_peaks` from total unique fragments. Plot QC distributions and relationships per library. Build a `SingleCellExperiment` from raw peak counts and apply scDblFinder's ATAC-specific `aggregateFeatures=TRUE` workflow per capture/library, saving score and class. Compare permissive and strict filters while retaining an unfiltered checkpoint. Specify whether FRiP counts reads, insertions or unique fragments, and count a fragment once if using a union-of-peaks fragment definition.

**Exercises:** (1) Guided: interpret four contrasting QC profiles. (2) Partly guided: compare retention under two justified rules. (3) Independent: inspect whether one donor or putative rare lineage is disproportionately removed.

**Quiz:** Q1 Why is there no universal TSS cutoff? Q2 Can low FRiP reflect an incomplete peak set? Q3 Is every high-depth barcode a doublet? Q4 Why run doublet detection per capture? Q5 Can a peak-only matrix supply a reliable whole-library FRiP denominator?

**Mistakes and troubleshooting:** copying tutorial thresholds; running RNA default doublet settings on ATAC; forcing the predicted doublet rate to equal an expectation; considering doublet calls certain. Check annotation/build if TSS enrichment is implausibly low across all cells and check for homotypic doublets missed by the classifier.

**Assignment and completion:** before/after QC dashboard, retention table by library and written rule justification with sensitivity analysis. Pass when filters are explainable and no population loss is hidden. **Time:** 8 h. **Read:** R7 QC; R8 API changes; R14. **Budget:** 1.5 h lessons, 3 lab, 1.5 exercises, 0.5 quiz, 1.5 report.

### Week 7 — TF-IDF and latent semantic indexing

**Prerequisite:** retained QC object and basic matrix vocabulary. **Objectives:** explain TF-IDF intuitively, compute LSI, detect depth-associated components and preserve raw counts.

**Short lessons:** TF reduces dominance by deeply sequenced cells; IDF changes the influence of broadly versus selectively accessible regions. LSI summarizes patterns across many sparse peaks through a singular-value decomposition. A component can reflect depth, biology or both; neither keeping nor removing component 1 is a universal rule.

**Guided lab:** set `DefaultAssay(pbmc) <- "peaks"`, then use `FindTopFeatures()`, `RunTFIDF()` and `RunSVD()`. Document the feature cutoff and TF-IDF method. Examine `DepthCor()` and correlations between each LSI coordinate, depth and QC. Create and justify `dims_keep`; compare neighbor structure with and without a suspect dimension, keeping counts untouched.

**Exercises:** (1) Guided: calculate a simple toy TF and explain the intuition behind IDF; do not assume every package implements the same formula. (2) Partly guided: compare two feature thresholds. (3) Independent: write a component-selection memo supported by depth and marker evidence.

**Quiz:** Q1 What does term frequency normalize? Q2 What is the purpose of IDF? Q3 Does LSI operate directly on untransformed raw counts in this workflow? Q4 Must LSI 1 always be removed? Q5 Should TF-IDF values enter an edgeR raw-count model?

**Mistakes and troubleshooting:** densifying the full matrix; treating singular-vector signs as biological directions; choosing dimensions solely by visual appeal. Start with fewer features/cells if memory fails and preserve the exact reduction settings.

**Assignment and completion:** LSI checkpoint, depth-correlation figure and two-setting comparison. Pass with an explicit component choice and retained raw counts. **Time:** 7 h. **Read:** R7 dimension reduction; [RunTFIDF documentation](https://stuartlab.org/signac/reference/runtfidf). **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 8 — Neighborhoods, UMAP and clustering

**Prerequisite:** Week 7 LSI and selected dimensions. **Objectives:** construct a graph, compare clustering resolutions and interpret UMAP without overclaiming geometry.

**Short lessons:** a neighborhood graph connects cells that are similar in the selected representation. Clustering partitions this graph; it does not discover an immutable number of biological types. UMAP is a useful view of local organization, but distances between separated islands and island areas should not be interpreted literally.

**Guided lab:** use `FindNeighbors(pbmc, reduction="lsi", dims=dims_keep)`, `FindClusters()` and `RunUMAP(..., reduction="lsi", dims=dims_keep, seed.use=...)`. Specify the clustering algorithm and resolution. Repeat with another seed, two resolutions and one alternative dimension choice. Compare cluster membership using a contingency table and an agreement measure; retain QC overlays.

**Exercises:** (1) Guided: identify a cluster dominated by low-quality cells. (2) Partly guided: compare membership across resolutions. (3) Independent: propose a stable, interpretable clustering with a reasoned exception for an unresolved population.

**Quiz:** Q1 On what representation is the graph built? Q2 Does higher resolution necessarily reveal more real cell types? Q3 Does a large UMAP gap prove a developmental barrier? Q4 Why retain random seeds? Q5 What evidence besides UMAP supports a cluster?

**Mistakes and troubleshooting:** choosing the prettiest map; clustering the two-dimensional UMAP instead of the chosen latent representation; assuming a fixed cluster number across versions. Inspect marker coherence and QC before assigning biological names.

**Assignment and completion:** paired embedding panels, cluster comparison table and a one-page clustering decision. Pass when the chosen representation and stability limitations are reported. **Time:** 6 h. **Read:** R7 clustering/visualization. **Budget:** 1 h lessons, 2 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 9 — Cell-type annotation with uncertainty

**Prerequisite:** Week 8 clusters; Weeks 1 and 6 biological/QC interpretation. **Objectives:** combine marker accessibility, gene-activity scores and reference labels; keep unresolved labels when evidence is weak.

**Short lessons:** a gene-activity score summarizes nearby DNA accessibility, not measured RNA. CD3D/CD3E, MS4A1/CD79A and LST1/LYZ are starting marker hypotheses for T, B and myeloid populations; their accessibility must be checked in the actual dataset. Reference mapping predicts resemblance to a reference and may force unfamiliar states into familiar labels.

**Guided lab:** use `GeneActivity()` and store its matrix in a separate assay named `ACTIVITY`, avoiding the misleading name `RNA`. Normalize that score assay for visualization. Inspect `CoveragePlot()` and marker panels, including negative evidence. Follow `FindTransferAnchors()`/`TransferData()` in R7 using its documented RNA reference, retain prediction scores, and compare predictions with manual annotations. Verify compatible gene identifiers and reference provenance.

**Exercises:** (1) Guided: support a broad annotation using three markers. (2) Partly guided: investigate a cluster with conflicting marker sets. (3) Independent: assign broad labels and confidence levels, allowing `unassigned`.

**Quiz:** Q1 Is gene activity RNA expression? Q2 Why use multiple markers? Q3 What can a low prediction score mean? Q4 Is a reference label ground truth? Q5 Why can a broad label be preferable to a precise subtype label?

**Mistakes and troubleshooting:** treating predicted labels as observations; circularly validating a reference with the same genes used for transfer; interpreting marker dropout as exclusion. Check ID mapping, fragment paths and negative markers before relabeling cells.

**Assignment and completion:** annotated atlas checkpoint, marker/coverage panels and an annotation evidence table. Pass when every named population has independent supporting evidence or is marked uncertain. **Time:** 7 h. **Read:** R7 annotation; R16 for future paired-data comparison. **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 10 — Peak universes and reproducible atlas checkpoint

**Prerequisite:** Week 9 broad cell types; Week 5 peak calling; Week 3 coordinate contracts. **Objectives:** construct common features, compare peak definitions and prepare multi-donor data without erasing provenance.

**Short lessons:** counting each sample against a different set of peaks makes comparisons ambiguous. A simple union can create very broad merged intervals; a strict intersection can discard real condition-specific regions. A defensible consensus policy balances representation across samples, cell types and conditions and records how competing peaks are resolved.

**Guided lab:** on A, call a small cell-type-specific peak set with `CallPeaks(..., group.by=..., macs2.path=...)` and compare it with vendor peaks. Use `rtracklayer::import()` and `GenomicRanges` to inspect/reconcile overlaps, then `FeatureMatrix()` to recount a shared universe. Begin C in a separate hg19 project: verify its eight-library sample sheet, identify a sufficiently represented lineage and define a balanced peak-discovery plan. Peak calling itself may run unattended or use provided peaks for the low-resource route.

**Exercises:** (1) Guided: show how `reduce()` can merge adjacent/overlapping windows. (2) Partly guided: compare union and intersection feature counts and retention. (3) Independent: propose a sample-balanced consensus rule that does not require every peak to appear in both conditions.

**Quiz:** Q1 Why need common matrix rows? Q2 What does strict intersection risk losing? Q3 Why can union/reduction create undesirable features? Q4 Does condition-balanced peak discovery guarantee an unbiased test? Q5 Can hg38 A peaks be used directly with hg19 C fragments?

**Mistakes and troubleshooting:** selecting peaks only from the condition expected to change; allowing high-depth samples to dominate discovery; hidden cell-ID collisions. Retain original barcodes plus globally unique sample-prefixed IDs and a mapping between them.

**Assignment and completion:** reproducible A atlas audit and C design/feature manifest. The checkpoint earns its separate course component: 30% input/QC integrity, 25% annotation evidence, 20% sensitivity analysis, 15% reproducibility and 10% communication. Pass with common coordinates, explicit feature construction and a replayable script. **Time:** 8 h. **Read:** R13 and R23. **Budget:** 1.5 h lessons, 3 lab, 1.5 exercises, 0.5 quiz, 1.5 audit.

### Week 11 — Differential accessibility and valid questions

**Prerequisite:** Week 10 common features and Week 2 independence. **Objectives:** distinguish descriptive cell-type markers from condition inference, define a contrast and interpret effect size/FDR separately.

**Short lessons:** a peak can differ between cell types without changing within a type after treatment. A condition comparison on all PBMCs can reflect changed cell composition. Testing many peaks produces chance discoveries, so the testing family, filtering rule and FDR procedure must be explicit.

**Guided lab:** use A to obtain exploratory cell-type markers with `FindMarkers(test.use="LR", latent.vars="nCount_peaks")`, checking the active assay and implementation. This is a descriptive ranking, not a donor-replicated treatment test. Report detection rates and effect estimates alongside adjusted p-values. In C, plot donor-by-condition cell counts and QC for the chosen lineage; write the contrast “LPS minus control at 6 h within this lineage” before fitting any model.

**Exercises:** (1) Guided: interpret a large versus small log2 fold change. (2) Partly guided: adjust a small list of p-values with `p.adjust(..., method="BH")`. (3) Independent: explain why adding depth as a cell-level covariate does not solve pseudoreplication.

**Quiz:** Q1 Is a cell-type marker necessarily a treatment-responsive peak? Q2 What is pseudoreplication? Q3 Does FDR 0.05 make each selected peak 95% likely to be real? Q4 Why report effect sizes? Q5 Does an integrated embedding supply raw counts for testing?

**Mistakes and troubleshooting:** treating `FindMarkers()` adjusted p-values as automatically BH-adjusted—check its documented method; interpreting nonsignificance as equivalence; testing only peaks selected by the same extreme contrast. Preserve full testing tables and define filters before looking at outcomes.

**Assignment and completion:** exploratory marker table, a contrast specification and a critique of a deliberately cell-pseudoreplicated analysis. Pass when descriptive and population-level claims are separated. **Time:** 7 h. **Read:** R7 marker examples; R15 inference and replication. **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 12 — Donor-level pseudobulk inference

**Prerequisite:** C sample/lineage audit, Week 11 contrast and raw common-peak counts. **Objectives:** aggregate within donor-condition-lineage groups, fit a paired design, check confounding and report uncertainty.

**Short lessons:** pseudobulk sums raw counts within an experimental sample and cell type; it does not average normalized expression. For four donors measured twice, there are eight sample columns but four independent donor blocks. A paired donor term accounts for baseline differences; it cannot repair perfect condition/batch confounding.

**Guided lab:** sum sparse counts into sample-by-lineage groups, retaining cell numbers and exclusion reasons. Align pseudobulk matrix columns exactly to the sample table. Set control as the reference level and use `model.matrix(~ donor + condition, data=sample_table)`. Check full rank with `qr(design)$rank`. In edgeR use `DGEList()`, `filterByExpr(..., design=design)`, `calcNormFactors()`, `estimateDisp()`, `glmQLFit()` and `glmQLFTest()` for the explicitly named LPS coefficient. Report BH FDR over the stated tested peak family. Draw paired donor plots and repeat after leaving each donor out. Inspect compositional-normalization assumptions, particularly if widespread accessibility changes are plausible.

**Exercises:** (1) Guided: aggregate a small integer matrix by donor-condition. (2) Partly guided: diagnose a rank-deficient batch/condition design. (3) Independent: justify a low-count filter, minimum cell-support rule and cautious interpretation for n=4.

**Quiz:** Q1 What values are summed? Q2 Why include donor in a paired model? Q3 Can a perfectly confounded batch coefficient be estimated separately? Q4 How many independent donor blocks are present? Q5 What does a leave-one-donor-out reversal suggest?

**Mistakes and troubleshooting:** summing TF-IDF values; entering zero columns for missing cell types as if they were true zero-accessibility samples; adding too many covariates for eight libraries; treating downsampled pseudo-replicates as new donors. Age/sex fixed within donor are usually absorbed by the donor fixed effects, not independently estimable alongside them.

**Assignment and completion:** pseudobulk counts/sample table, model design audit, full results, paired donor plot and sensitivity report. Pass with correct raw counts, full-rank design, stated FDR family and no invented replication. **Time:** 8 h. **Read:** R15 edgeR guide and R23. **Budget:** 1.5 h lessons, 3 lab, 1.5 exercises, 0.5 quiz, 1.5 report.

### Week 13 — Motifs, chromVAR and transcription-factor hypotheses

**Prerequisite:** Week 10 peak universe; Week 11/12 interpretation; matching genome sequence. **Objectives:** compare motif overrepresentation with chromVAR deviations, control sequence/depth biases and avoid claiming binding from motifs.

**Short lessons:** a motif is a sequence pattern compatible with binding by one or several TFs. Enrichment asks whether a selected peak set contains a motif unusually often relative to an appropriate background. chromVAR asks whether motif-associated peaks show accessibility deviations after bias correction; it does not measure TF protein abundance or occupancy.

**Guided lab:** select and record a motif-database release and motif IDs, then use `AddMotifs()` with the correct `BSgenome`. Inspect foreground/background GC content, width and accessibility before `FindMotifs()`. In the pinned Signac 1.16 environment use `RunChromVAR()` with chromVAR 1.32.0. Compare enrichment rankings and per-cell deviation summaries. For condition claims summarize deviations per donor-condition; do not treat thousands of cell scores as independent replicates or submit non-integer deviations to edgeR.

**Exercises:** (1) Guided: calculate a small motif foreground/background table. (2) Partly guided: show how unmatched GC backgrounds alter interpretation. (3) Independent: nominate a TF family with competing explanations and a required validation assay.

**Quiz:** Q1 Does motif presence demonstrate binding? Q2 Why can different TFs share a motif? Q3 Why match motif backgrounds? Q4 Is a chromVAR z-score a binding probability? Q5 Does RNA evidence for a TF prove its protein is active?

**Mistakes and troubleshooting:** using hg38 sequence for hg19 regions; confusing motif-family labels with unique TF identity; failing after upgrading to a Signac version without `RunChromVAR()`. Restore the pinned environment or explicitly redesign and document the replacement workflow.

**Assignment and completion:** motif enrichment table, background diagnostic, deviation figure and one-page TF hypothesis with uncertainty. Pass with motif identifiers, genome/database versions, appropriate replication and restrained wording. **Time:** 7 h. **Read:** R10 and R8. **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 14 — Gene scores, co-accessibility and footprinting limits

**Prerequisite:** Weeks 9 and 13; D subset for optional co-accessibility. **Objectives:** distinguish promoter/gene-body scores, co-accessibility and peak–RNA association; critique regulatory links and footprints.

**Short lessons:** nearby regions can covary because of a shared cell state, depth or technical effects. Co-accessibility suggests coordinated accessibility but does not establish physical contact or a regulatory direction. Footprints depend on sequencing depth, Tn5 sequence preference, nucleosomes, TF residence behavior and cell mixtures; a dip alone does not prove occupancy.

**Guided lab:** compare gene-activity definitions and coverage at one locus. Follow the one-chromosome Cicero workflow in R18: create its cell-data object, supply coordinates/genome sizes and compute a small network. Compare distance and co-accessibility evidence; leave genome-wide computation optional. Use `Footprint()` and `PlotFootprint()` on well-supported motif instances as an optional pooled diagnostic, checking expected bias profiles. Prepare a table of what extra evidence would be needed for each candidate link. Peak–RNA `LinkPeaks()` is deferred until genuinely paired data in Week 16.

**Exercises:** (1) Guided: distinguish a nearest-gene assignment from a tested association. (2) Partly guided: compare a link after changing aggregation/window settings. (3) Independent: design a perturbation experiment that could distinguish a candidate enhancer from a correlated passenger region.

**Quiz:** Q1 Is co-accessibility physical chromatin contact? Q2 Must the closest gene be the target? Q3 Can a pooled footprint be read as single-cell occupancy? Q4 How can Tn5 bias mimic a footprint? Q5 What would provide stronger causal evidence than correlation?

**Mistakes and troubleshooting:** building dense all-by-all peak correlations; mixing assemblies in genome-size tables; linking computationally paired RNA as if it came from the same cells. Restrict region and chromosome scope explicitly; treat unstable links as findings about uncertainty.

**Assignment and completion:** annotated locus view, a small link evidence table and a footprint critique. Pass when association, contact and causation are distinguished. **Time:** 7 h. **Read:** R18; R16 link interpretation. **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 15 — Batch effects and integration across experiments

**Prerequisite:** common features, LSI, annotation and donor-aware design. **Objectives:** diagnose batch effects, distinguish technical mixing from biological conservation, and compare separate RNA–ATAC mapping with paired assays.

**Short lessons:** integration aims to make comparable biological states comparable across datasets. It can erase true condition effects if those effects resemble the variation the algorithm removes. An embedding that mixes donors perfectly may also overcorrect; a fully confounded design cannot be rescued by an attractive UMAP.

**Guided lab:** on a small A/B ATAC subset, verify hg38 compatibility and recount shared features. Follow R17's LSI-based integration using `FindIntegrationAnchors(..., reduction="rlsi")` and `IntegrateEmbeddings()` as documented. Compare uncorrected and corrected graphs with donor/dataset mixing within each broad cell type and preservation of marker separation. Revisit RNA label transfer from Week 9; distinguish transferred labels from joint measurement. Keep original counts and sample IDs for all differential analyses.

**Exercises:** (1) Guided: recognize a batch-associated component. (2) Partly guided: compare mixing within cell type rather than globally. (3) Independent: reject or accept an integration based on biology-preservation evidence.

**Quiz:** Q1 Does global mixing always mean success? Q2 What is overcorrection? Q3 Can integration repair perfect donor-condition confounding? Q4 Are unpaired RNA and ATAC observations the same cells? Q5 Which counts should enter pseudobulk testing after integration?

**Mistakes and troubleshooting:** integrating incomparable lineages; optimizing a mixing score with no biology check; merging hg19 C directly with hg38 A. Keep this lab restricted to compatible A/B data and check features/barcode namespaces before integrating.

**Assignment and completion:** before/after embeddings, within-type batch diagnostic and marker-preservation report. Pass with a defensible decision to integrate or retain the uncorrected representation. **Time:** 7 h. **Read:** R17 and R7 reference mapping. **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 16 — Paired RNA–ATAC multiome

**Prerequisite:** Weeks 7–9, 14–15 and dataset B. **Objectives:** maintain identical paired barcodes, perform modality-specific QC/reductions, build a WNN graph and evaluate peak–gene associations.

**Short lessons:** paired assays measure two modalities in the same nucleus, but neither modality is complete. WNN learns a per-cell combination of modality information; its weights are not causal fractions or proof of one assay's biological importance. Peak–RNA correlations can still reflect cell identity or other shared factors.

**Guided lab:** use `Read10X_h5()` to load B's `Gene Expression` and `Peaks` matrices, align paired barcodes and create RNA/ATAC assays. Perform RNA normalization/PCA and ATAC TF-IDF/LSI separately, selecting dimensions and QC for each. Use `FindMultiModalNeighbors(reduction.list=list("pca","lsi"), dims.list=list(rna_dims,atac_dims))`, then WNN UMAP/clustering. Use `RegionStats()` and `LinkPeaks(peak.assay="ATAC", expression.assay="RNA", genes.use=...)` on a small gene list with matching annotations. Audit its candidate-distance/background assumptions and p-value/multiplicity behavior; document any additional FDR family explicitly.

**Exercises:** (1) Guided: verify barcode intersection and RNA/ATAC dimensions. (2) Partly guided: compare RNA-only, ATAC-only and WNN annotation. (3) Independent: evaluate one candidate peak–gene relationship within a cell type and across the full dataset, explaining confounding.

**Quiz:** Q1 What makes this dataset experimentally paired? Q2 Does RNA PCA replace ATAC LSI? Q3 Does high RNA WNN weight imply RNA causes accessibility? Q4 Can cell-type composition create a peak–RNA association? Q5 Does this one-donor dataset support general population treatment inference?

**Mistakes and troubleshooting:** retaining unmatched cells; mixing ARC 1.0.0 matrices with ARC 2.0.0 metadata; assuming every `LinkPeaks()` p-value is already adjusted; inferring causality. Check barcode identity, feature builds and the exact testing output before interpretation. Do not apply BH only to links that have already passed a p-value/effect cutoff: retain all tested candidate p-values and define their family, or label the filtered links exploratory when that information is unavailable.

**Assignment and completion:** paired QC table, three embeddings, modality-weight panel, selected locus plot and link-evidence table. Pass with genuine pairing, separate preprocessing and bounded conclusions. **Time:** 8 h. **Read:** R16. **Budget:** 1.5 h lessons, 3 lab, 1.5 exercises, 0.5 quiz, 1.5 report.

### Week 17 — Trajectories and regulatory programs

**Prerequisite:** Week 8 graph interpretation and D HSPC data; motif concepts. **Objectives:** infer a plausible ordering, justify its root, assess branch stability and distinguish pseudotime from measured time.

**Short lessons:** pseudotime is a data-driven ordering along a chosen model, not a clock. Distinct mature PBMC lineages should not automatically be forced into a continuous differentiation trajectory. Root choice, sampling gaps and the embedding can create apparent directional stories.

**Guided lab:** reproduce a restricted D lineage from R19 using `as.cell_data_set()`, Monocle3 clustering/graph learning and `order_cells()` with an explicitly chosen root. If using the smaller ArchR teaching subset, keep its data provenance separate. Examine marker accessibility and motif deviations across pseudotime bins; vary root and neighborhood settings. Use donor-aware models if adequate donors exist; a one-replicate practice trajectory is descriptive. Graph-association tests are exploratory when ordering is learned from the same data.

**Exercises:** (1) Guided: identify a biologically defensible root. (2) Partly guided: compare trajectories with another plausible root. (3) Independent: formulate a regulatory-program hypothesis and name longitudinal or lineage-tracing evidence needed to test it.

**Quiz:** Q1 Is pseudotime elapsed time? Q2 Why does root selection matter? Q3 Do branches prove lineage ancestry? Q4 Can a mature PBMC mixture imply a valid differentiation path? Q5 What is the danger of testing features used to construct the same trajectory?

**Mistakes and troubleshooting:** treating disconnected populations as one continuum; interpreting smoothed curves as precise timing; allowing root labels derived from the outcome of interest. Limit inference to the sampled process and report sensitivity.

**Assignment and completion:** trajectory/root justification, feature-program plot, alternative-root comparison and limitations paragraph. Pass even if the evidence rejects a unique trajectory. **Time:** 7 h. **Read:** R19 and R18. **Budget:** 1.5 h lessons, 2.5 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 18 — ArchR and a focused Python alternative

**Prerequisite:** established R workflow; Week 3 Python basics. **Objectives:** explain disk-backed computation, run a small alternative workflow and compare scientific conclusions across implementations.

**Short lessons:** on-disk storage reduces memory pressure but introduces IO and file-lifecycle requirements. Different feature sets, counting definitions and reductions can produce different neighborhoods without either workflow being intrinsically wrong. Compare corresponding cells and biological evidence, not numeric cluster IDs.

**Guided lab:** in a separate release-pinned ArchR project use `getTutorialData("Hematopoiesis")`, `addArchRGenome("hg19")`, `createArrowFiles()`, `addDoubletScores()` and `addIterativeLSI()` on a small subset; record the Arrow paths and sample mapping. In a separate Python 3.12 environment install SnapATAC2 2.10.0 and apply `snap.pp.import_fragments(..., file="subset.h5ad", sorted_by_barcode=False)` to a documented small fragment subset, using its correct genome sizes. Follow its tile/feature selection, `snap.tl.spectral`, nearest-neighbor and clustering steps. Do not call this spectral embedding “the same LSI.” Full installation/download waiting is outside active lab time; use official cached example outputs for a labeled interpretation exercise if a machine cannot support the run.

**Exercises:** (1) Guided: inspect a backed file and its metadata. (2) Partly guided: compare resource use at two thread settings on the same small input. (3) Independent: compare label agreement and marker evidence across two workflows, separating feature changes from algorithm changes.

**Quiz:** Q1 Does backed storage eliminate all RAM requirements? Q2 Are cluster 0 in two tools necessarily the same population? Q3 Is SnapATAC2 spectral embedding identical to Signac LSI? Q4 Why isolate ArchR's environment? Q5 Which SnapATAC2 import API replaces `pp.import_data`?

**Mistakes and troubleshooting:** opening old/new H5AD formats interchangeably; oversubscribing cores; deleting Arrow files before saving/exporting their dependent project; running native Windows where unsupported. Record storage backend, versions, working directories and peak RAM.

**Assignment and completion:** two small workflow records, a method-comparison table and measured or explicitly estimated resource costs. Pass when claims identify which workflows were actually executed. **Time:** 8 h. **Read:** R20–R21. **Budget:** 1.5 h lessons, 3 lab, 1.5 exercises, 0.5 quiz, 1.5 report.

### Week 19 — Reproducibility, efficiency and debugging

**Prerequisite:** at least one saved end-to-end analysis checkpoint. **Objectives:** rebuild an analysis, track dependencies/inputs and diagnose memory, path, assembly and API failures systematically.

**Short lessons:** a saved R object is helpful but is not a complete reproducible analysis. Reproduction also needs inputs, sample definitions, reference files, software, parameters and an executable sequence. More parallel workers can increase memory use because data may be copied to each worker.

**Guided lab:** write an entry script with early file/version/barcode assertions. Use `renv::snapshot()` and `sessionInfo()`, and check file checksums. Measure representative steps using `system.time()`; report `object.size()` only as object storage, not total process peak RAM. Restore in a fresh project where feasible. Inject a moved fragment path, shuffled metadata and incompatible assembly into a small copy, then diagnose them. Keep reusable intermediate checkpoints and avoid dense conversion, full-genome footprints and all-by-all peak correlations when unnecessary.

**Exercises:** (1) Guided: repair one broken fragment path while preserving the source file identity. (2) Partly guided: locate the first failing assertion in a scripted pipeline. (3) Independent: produce a minimal reproducible example for an API error without publishing sensitive data.

**Quiz:** Q1 Does an RDS file alone guarantee reproducibility? Q2 What does a checksum verify? Q3 Does `object.size()` equal peak process RAM? Q4 Can more parallel workers increase memory use? Q5 Why retain an unfiltered checkpoint?

**Mistakes and troubleshooting:** relying on a hidden workspace; treating a seed as complete reproducibility; removing provenance to save space; claiming bitwise identical stochastic outputs across platforms. Define tolerances for numerical results and stable biological conclusions.

**Assignment and completion:** clean-session rerun log, dependency record, manifest, timing table and troubleshooting notes. Pass when another learner can locate inputs and follow the workflow, with unresolved failures stated. **Time:** 6 h. **Read:** R22; software troubleshooting sections in R8/R20. **Budget:** 1 h lessons, 2 lab, 1.5 exercises, 0.5 quiz, 1 report.

### Week 20 — Scientific reporting and independent design

**Prerequisite:** Weeks 1–19 deliverables. **Objectives:** write a concise methods/results narrative, audit claims against evidence and design the four capstones without prespecifying their findings.

**Short lessons:** a useful report says what was observed, how uncertainty was assessed and what remains a hypothesis. Negative or inconclusive results are legitimate outcomes. Publication-quality figures expose relevant sample-level variation, not just attractive cell-level patterns.

**Guided lab:** assemble a report from existing scripts/outputs, optionally with Quarto. Generate an evidence table linking each biological claim to a figure, analysis, sample count and limitation. Re-run one key figure from its raw or versioned intermediate input. Compare the report against the pseudobulk design and data provenance before proposing mechanisms.

**Exercises:** (1) Guided: rewrite three overclaimed statements. (2) Partly guided: audit a figure with missing donor information. (3) Independent: write a capstone proposal with a primary contrast, alternative outcomes, resource plan and stopping criteria.

**Quiz:** Q1 Can a null result receive full credit? Q2 What should accompany a UMAP in a condition-comparison report? Q3 What separates a mechanistic hypothesis from a demonstrated mechanism? Q4 Why report exclusions and failed sensitivity checks? Q5 What is the independent replication unit in C's paired design?

**Mistakes and troubleshooting:** presenting expected markers as findings before running the data; hiding failed QC; implying all samples share the same genome assembly; selecting only supporting loci. Report the prespecified analysis family and contradictory evidence.

**Assignment and completion:** a 6–10-page report or equivalent notebook, a 5-minute explanation and four capstone briefs. The Week 20 assessment is 30% valid design, 25% evidence/interpretation, 20% reproducibility, 15% limitations and 10% clarity. Pass when every major claim is traceable and none relies on pseudoreplication or causal overstatement. **Time:** 7 h. **Read:** R1 communication; revisit R15/R23. **Budget:** 1 h lessons, 2 lab, 1.5 exercises, 0.5 quiz, 2 report/defense.

### Separate advanced lab — Optional FASTQ processing

**When:** after Week 5, preferably after Week 19's resource planning. **Extra time:** 4–6 active hours plus download and approximately 8–24+ hours of computation for a full 10k library. **Input:** dataset A FASTQs and the exact compatible human reference from its official page; retain sample/read-layout information. Use the 2.1.0 software/reference combination for historical reproduction if available, or record a new 2.2 run as reprocessing. Do not run 10x Cell Ranger ATAC on Bio-Rad SureCell dataset C.

The following is **pseudocode with placeholders**, not an executable command:

```bash
cellranger-atac count \
  --id=YOUR_RUN_ID \
  --reference=/absolute/path/to/validated_reference \
  --fastqs=/absolute/path/to/matching_FASTQs \
  --sample=EXACT_FASTQ_SAMPLE_PREFIX \
  --localcores=8 \
  --localmem=64
```

`--id` names a new output run; `--reference` selects the alignment reference; `--fastqs` and `--sample` identify the correct reads; the local resource flags cap requested cores and memory. Verify the CLI help for the installed release before execution.

Milestones: (1) verify FASTQ identity, assay read layout and checksums; (2) run alignment/barcode processing/cell calling; (3) audit mapping, duplication, valid-barcode and fragment/QC summaries; (4) compare output matrices and fragments with deposited outputs using matched definitions. Keep intermediate processing decisions explicit. The pipeline handles protocol-specific alignment, barcode correction and duplicate logic; do not add a second unconditional duplicate removal or Tn5 shift. Downsampling must retain synchronized read/barcode records from the same library. A tiny downsample is useful for pipeline mechanics, not for validating final biological sensitivity or cell calling.

Deliver a command/configuration record, version/reference manifest, pipeline report and comparison memo. A different retained-cell count is not automatically an error; explain software, reference, depth and threshold changes. This lab is optional and outside the core 6–8-hour weekly budget.

## Quiz answer keys with reasoning

These answers are intentionally separated from the weekly questions. Award one point for the answer and one for its explanation. Alternative wording is welcome when the reasoning is sound.

### Week 1 key

1. DNA fragments. ATAC assays transposase-accessible chromatin; it does not sequence RNA transcripts directly.
2. No. Accessibility can be permissive or poised, and expression also depends on transcription, RNA processing and decay.
3. No count was observed for that region/cell under this assay and feature definition. Sparse sampling prevents interpreting every zero as biological closure.
4. No: one donor. The cells are observations nested within a shared biological source.
5. No. Distance and accessibility identify candidates, not a validated target or regulatory effect.

### Week 2 key

1. With skewed distributions or strong outliers. The median is less influenced by extremes; show the distribution as well.
2. Missing means unknown/uncomputed/failed; zero is a recorded numerical value. Substitution changes the meaning.
3. An independently sampled biological source, usually a donor in these PBMC studies. Paired conditions from that donor remain related.
4. No. Precision and sample size can make a very small effect statistically detectable.
5. To make stochastic choices repeatable in a fixed environment. It does not replace dependency and input records.

### Week 3 key

1. 50 bases; the 1-based closed interval is 101–150. The BED end is excluded.
2. Base-call uncertainty for each sequenced base, encoded as quality scores. It is not a cell-quality score.
3. No. It only changes chromosome naming; assembly conversion requires coordinate mapping and validation.
4. It supports rapid region queries without reading the entire file. Index and data must match.
5. No. Numeric conversion can destroy leading zeros and punctuation required for matching.

### Week 4 key

1. A dense peak-by-cell matrix can consume many GB even when almost all entries are zero.
2. Its barcode/cell identifier in the correct sample namespace. The same barcode string may occur in another library.
3. Usually not if later QC, coverage or recounting needs it. The object may store only the path, not fragment bytes.
4. Organism, assembly, chromosome conventions and compatible genomic coordinates. Annotation release also affects definitions.
5. To detect file changes and confirm input identity. A checksum does not establish biological correctness.

### Week 5 key

1. No. It can include the unique pair and PCR duplicate pairs supporting one fragment.
2. Independent cells can contain genuine fragments at the same locus. Collapsing across barcodes destroys real signal.
3. No. A second shift changes already adjusted coordinates and can bias regional analyses.
4. Individual cells are sparse; pooled biological groups provide enough signal for enrichment detection.
5. No. Off-peak information was not retained in that summary, so fragments or another source are needed.

### Week 6 key

1. Chemistry, annotation, depth, cell type and implementation change the distribution and score scale.
2. Yes. A peak universe missing a lineage's regulatory regions can reduce its apparent FRiP.
3. No. High-quality singlets and different cell states can have high counts; combine evidence.
4. Doublets arise within capture/library events, and expected rates and data distributions differ among captures.
5. No. Peak counts alone omit off-peak fragments, and peak-count conventions may differ from fragment counts.

### Week 7 key

1. Within-cell counts relative to depth under the selected TF-IDF implementation.
2. To change the weight of features according to how broadly represented they are, helping expose informative variation.
3. No. This workflow decomposes a TF-IDF-transformed representation.
4. No. Examine depth and biological associations; the first component is not universally technical.
5. No. Those are transformed scores, not the raw count observations modeled by edgeR.

### Week 8 key

1. The selected LSI coordinates in this route, not automatically the plotted UMAP coordinates.
2. No. It makes a finer graph partition that may split continuous states or technical variation.
3. No. UMAP's two-dimensional geometry can distort global distances and disconnected layouts.
4. To record stochastic choices and support reproducibility of a given configuration.
5. Coherent marker accessibility, stability across reasonable settings, QC and biological plausibility.

### Week 9 key

1. No. It is derived from DNA accessibility near a gene under a specified scoring rule.
2. Individual loci can be sparse or nonspecific; coherent positive and negative evidence is stronger.
3. Low-quality data, a missing reference state, ambiguous identity or imperfect modality alignment; inspect rather than force.
4. No. It is a model-dependent prediction constrained by the reference's coverage and labels.
5. It can match the resolution supported by the data and avoid false specificity.

### Week 10 key

1. Each row must represent the same genomic feature for between-sample counts to be comparable.
2. Genuine cell-type- or condition-specific regions that are not detected in every sample.
3. Overlapping/adjacent intervals can form broad merged regions, obscuring individual regulatory elements and changing counts.
4. No. It reduces one selection imbalance but does not remove depth, feature-selection or other biases automatically.
5. No. Their coordinates refer to different assemblies; common names do not make them compatible.

### Week 11 key

1. No. Cell identity differences and within-lineage treatment responses are different contrasts.
2. Treating correlated observations, such as cells from one donor, as independent experimental replicates.
3. No. FDR describes expected false-discovery proportion among a selected set under assumptions, not individual posterior truth.
4. They express magnitude and direction, which statistical significance alone does not convey.
5. No. An embedding is a representation; inference here uses original sample-indexed raw feature counts.

### Week 12 key

1. Raw integer peak counts within the same sample/condition and cell type, with a documented counting convention.
2. It accounts for stable baseline differences between paired donors, comparing conditions within donor blocks.
3. No. The design cannot identify two perfectly aligned effects separately; software cannot invent missing comparisons.
4. Four. Eight condition-specific samples come from four biological donor blocks.
5. The conclusion may be driven by one donor or limited information; report instability and inspect confounding/QC.

### Week 13 key

1. No. A compatible sequence can be unbound, inaccessible or bound by a related factor.
2. Related DNA-binding domains recognize similar sequence patterns, limiting TF-specific assignment.
3. Sequence composition, width and baseline accessibility can bias motif frequencies independently of the hypothesis.
4. No. It is a standardized bias-corrected deviation under a specified model/background, not an occupancy probability.
5. No. Translation, localization, modification and cofactors can alter TF protein activity.

### Week 14 key

1. No. It is statistical covariation, potentially induced by shared state or technical effects.
2. No. Enhancers can act across intervening genes and nearest-gene assignment can be wrong.
3. No. It aggregates heterogeneous cells and motif instances; it does not identify which cell had occupancy.
4. Unequal sequence preferences can create insertion patterns around motifs independent of protein protection.
5. A well-controlled region/TF perturbation with an appropriate response readout, replicates and controls; orthogonal binding/contact assays add complementary evidence.

### Week 15 key

1. No. Mixing unrelated cell types or removing true condition differences can also improve a global mixing metric.
2. Removal or distortion of meaningful biological variation while correcting technical differences.
3. No. Confounding is a design limitation; mixing does not identify the missing biological contrast.
4. No. Computational matching is not experimental pairing.
5. Original raw counts aggregated by biological sample and cell type, not integrated assay values.

### Week 16 key

1. RNA and ATAC measurements share the same experimentally assayed nuclei/barcodes.
2. No. Each modality needs appropriate normalization, features and reduction before combination.
3. No. It is an algorithmic information weight for neighborhood construction, not a causal effect.
4. Yes. A lineage can have both high accessibility and high RNA without a direct regulatory relationship between the selected region and gene.
5. No. Many cells from one donor do not provide replicated treatment or population evidence.

### Week 17 key

1. No. It is a relative ordering conditional on the model and sampled cells.
2. It defines the starting direction of that ordering and can change downstream narratives.
3. No. Similarity-based branches require independent evidence for lineage ancestry.
4. Not by itself. A mixture of mature lineages need not sample a continuous process.
5. Data reuse can make associations circular or optimistic; use sensitivity checks and independent evidence/validation.

### Week 18 key

1. No. Embeddings, graphs, metadata and temporary arrays still require RAM, and disk speed matters.
2. No. Cluster IDs are arbitrary labels, so match cells and evidence rather than numeric labels.
3. No. They use different mathematical representations and can also differ in feature construction.
4. Its dependency and R-version requirements differ; isolation helps avoid breaking the primary teaching workflow.
5. `pp.import_fragments`. Check the installed version, schema and current tutorial rather than copying obsolete code.

### Week 19 key

1. No. It does not necessarily include raw fragments, reference files, package versions or the analysis sequence.
2. Byte-level identity/change relative to a recorded digest; it does not validate sample labels or biology.
3. No. It describes the object, excluding many process allocations, temporary objects and worker copies.
4. Yes. Parallel processes can duplicate large objects and have their own intermediate allocations.
5. It enables new QC decisions and auditing without losing removed cells or rereading all original inputs.

### Week 20 key

1. Yes. Valid design, transparent uncertainty and accurate reporting matter more than obtaining a preferred effect.
2. Donor/sample-level QC, cell counts, effect-size summaries, model details and appropriate sample-level evidence.
3. Demonstration requires evidence that rules out plausible alternatives, commonly intervention plus suitable controls and replication.
4. They affect generalizability and reveal whether conclusions depend on selectively retained observations or settings.
5. The donor: four independent paired blocks, each supplying a control and an LPS sample.

## Four capstone projects

These begin after the curriculum, although earlier checkpoints can be reused. The estimates are additional active work for analysis, interpretation and reporting; unattended computing is extra. No project has a prescribed number of clusters, significant peaks, enriched motifs or successful links. Full credit is possible for well-supported null or inconclusive findings.

### Project 1 — Build and annotate a PBMC accessibility atlas

**Biological question:** Which broad PBMC populations and accessibility programs can be supported in this donor, and which annotations remain uncertain?

**Dataset:** A; 1,000–2,000-cell practice subset followed by the full dataset if resources allow. A is a one-donor atlas, so the scope is this sample. **Prerequisites:** Weeks 1–10 and 19. **Extra time:** 12–16 h; typically 1–4 h core unattended computation.

**Workflow and milestones:**

1. **Input contract, 2 h:** manifest, genome/annotation audit, barcode consistency, sparse import and resource budget.
2. **QC, 3–4 h:** multiple QC metrics and doublet evidence; compare at least two justified filter rules; record retention.
3. **Atlas construction, 3–4 h:** TF-IDF/LSI, justified dimension selection, graph, two resolutions and sensitivity to one seed/feature setting.
4. **Annotation/report, 4–6 h:** marker accessibility, gene activity, reference predictions and uncertainty table; replay from saved inputs.

**Required figures:** QC distributions with chosen cutoffs; depth versus TSS/LSI diagnostic; UMAP colored separately by cluster, QC and broad label; marker/activity panel; coverage at two informative loci; a cluster/annotation stability comparison. Figures must state the retained cell count, reference build and whether a signal is activity or measured expression.

**Interpretation questions:** Which labels rely on several independent observations? Does a suspicious cluster survive QC and parameter changes? What rare populations could subsetting miss? Which conclusions are confined to one donor? Does reference mapping add evidence beyond the markers already used?

**Deliverables:** notebook/script, manifest and dependency file, saved atlas checkpoint, annotation evidence table and a 4–6-page report. **Rubric, 100 points:** input integrity 15; QC justification 20; reduction/clustering robustness 20; annotation evidence 25; reproducibility 10; interpretation/clarity 10. Genome mismatch or calling gene activity measured RNA requires revision.

**Optional extension:** call cell-type-specific peaks and compare annotation stability with the vendor peak universe, without prespecifying an improvement.

### Project 2 — Compare accessibility with biological replication

**Biological question:** Within a sufficiently represented PBMC lineage, is accessibility associated with 6-hour LPS exposure relative to matched DMSO control, and how consistent are estimated effects across donors?

**Dataset:** C, exactly the eight libraries listed above, four matched donors, hg19. Choose a lineage after auditing representation; do not promise that a rare type will have sufficient data. **Prerequisites:** Weeks 1–12 and 19–20. **Extra time:** 18–24 h; typically 2–8 h unattended preprocessing/recounting.

**Workflow and milestones:**

1. **Design and provenance, 3 h:** verify donor, time, condition and absence of Golgi inhibitor; define a primary contrast, lineage inclusion policy and covariate rationale before inspecting differential peaks.
2. **Comparable cells/features, 5–7 h:** annotate broad types with uncertainty; audit QC and lineage coverage in each library; build or justify a shared peak universe. Prefer balanced peak discovery across donor-condition samples; document reuse of author peaks if chosen.
3. **Inference, 5–6 h:** sum raw counts within sample/lineage; fit full-rank paired donor-plus-condition edgeR design; filter low-information peaks; report effect direction/magnitude and BH FDR for the full tested family.
4. **Robustness/report, 5–8 h:** donor-specific plots, leave-one-donor-out analysis, alternate QC/peak-set comparison and a limitation-aware narrative.

**Required figures:** donor/condition design table; QC and retained lineage-cell counts by sample; pseudobulk MDS/PCA colored by donor and condition; MA or volcano plot with effect sizes and FDR; paired donor-level accessibility plots for selected loci chosen by a declared rule; coverage split by donor/condition; sensitivity concordance plot.

**Interpretation questions:** Are changes consistent across donors or driven by one? Could cell-state proportions within the chosen lineage explain an effect? What do n=4 and low cell counts limit? Are normalization assumptions plausible if LPS affects accessibility broadly? Which covariates are estimable? How does author-peak selection affect interpretation? Absence of significant peaks is not evidence of equivalence without a suitable precision/equivalence analysis.

**Deliverables:** sample sheet, raw pseudobulk count matrix, design matrix, complete tested-peak table, scripts and a 5–8-page report. **Rubric:** design/replication 25; feature/QC comparability 20; statistical implementation 25; effect interpretation/sensitivity 20; reproducibility/clarity 10. Cell-level pseudoreplication, dropping problematic donors without disclosure, or testing transformed counts requires revision.

**Optional extension:** repeat the prespecified contrast in a second sufficiently sampled lineage and explicitly expand or separate the multiplicity family; compare direction and uncertainty without declaring cell-type specificity solely because one p-value crosses 0.05 and another does not.

### Project 3 — Prioritize candidate TFs and regulatory regions

**Biological question:** Which sequence motifs, candidate TF families and accessible regions are associated with a chosen cell type or response, and what evidence would be needed to establish their roles?

**Dataset:** A for a cell-identity study, or C plus Project 2 for a replicated condition-response study. State which question is being asked; do not quietly convert a one-donor atlas into a condition experiment. **Prerequisites:** Weeks 9–14 and 19; Week 12 for condition claims. **Extra time:** 12–16 h, plus roughly 1–4 h motif calculations depending on scope.

**Workflow and milestones:**

1. **Question/foreground, 2 h:** choose the cell type and define foreground peaks from a documented contrast; retain the full eligible background universe.
2. **Motif evidence, 3–4 h:** add versioned motifs, compare matched background composition and calculate enrichment with multiple-testing accounting.
3. **Complementary evidence, 3–4 h:** compute chromVAR deviations; compare donor summaries where relevant; inspect coverage, gene-activity or real RNA evidence with correct labels; examine a few candidate links.
4. **Candidate evaluation, 4–6 h:** compare two plausible TF-family interpretations, critique footprints if used, and propose a perturbation plus controls and measurable outcomes.

**Required figures:** foreground/background GC/accessibility diagnostic; motif-enrichment effect/significance plot; chromVAR distribution by type or donor-condition; coverage and motif-position annotation at 2–3 candidates; evidence matrix distinguishing sequence, accessibility, expression, contact and perturbation. A bias-aware pooled footprint is optional, not a requirement for inferring TF binding.

**Interpretation questions:** Can multiple TFs explain this motif? Could background selection explain the ranking? Is motif-associated accessibility consistent with expression of the candidate TF, and what does consistency fail to prove? Could a link arise from shared cell identity? Which experiment could disprove the favored explanation?

**Deliverables:** complete motif/background tables, candidate-region table with coordinates/build, analysis notebook and a 4–6-page hypothesis report. **Rubric:** foreground/background design 25; motif/deviation implementation 20; integration of evidence 20; limitations and validation design 25; reproducibility/clarity 10. Claiming TF occupancy or enhancer causality from ATAC alone requires revision.

**Optional extension:** use an independently generated, context-matched TF-binding or perturbation dataset after verifying its metadata and availability; compare against a negative/control factor and describe mismatched tissue or condition limitations.

### Project 4 — Integrate paired RNA and ATAC

**Biological question:** How do RNA and chromatin accessibility jointly support cell identity, where do they disagree, and which peak–gene relationships merit further testing?

**Dataset:** B, using the consistent archived ARC 1.0.0 file set. Start with 2,000 paired cells and extend if resources permit. **Prerequisites:** Weeks 7–9 and 14–16, plus 19. **Extra time:** 16–22 h; 2–8 h unattended full-data processing, less on a subset.

**Workflow and milestones:**

1. **Pairing/QC, 3–4 h:** verify identical barcodes and sample identity, modality-specific QC and retained intersections; distinguish cells missing a modality from true zero measurements.
2. **Modality processing, 4–5 h:** RNA normalization/PCA, ATAC TF-IDF/LSI and independent graphs; choose dimensions from evidence.
3. **Joint analysis, 4–6 h:** WNN graph/UMAP, labels/confidence and modality-weight diagnostics; compare with single-modality solutions using the same cells.
4. **Regulatory/reporting work, 5–7 h:** a limited, declared peak–gene test universe; compare across-type and within-type associations; inspect effect sizes, backgrounds, multiple testing and negative examples; record sensitivity to QC/dimension settings.

**Required figures:** two-modality QC scatter/distributions; RNA-only, ATAC-only and WNN embeddings; modality weights against QC; annotation-agreement table; coverage plus paired RNA at selected genes; peak–gene association table/plot with uncertainty and testing scope. Include discordant or inconclusive examples chosen using a stated rule.

**Interpretation questions:** Does pairing resolve the ambiguity present in unpaired mapping? Which cells depend strongly on one modality, and could quality explain that? Do candidate links persist within a coherent cell type? Why is same-cell correlation still not causal? Which conclusions cannot generalize beyond one donor?

**Deliverables:** reproducible paired object, preprocessing/feature manifest, complete candidate-link table and a 5–8-page report. **Rubric:** pairing/provenance 20; modality-specific processing 20; integration validation 20; regulatory interpretation 25; reproducibility/communication 15. False pairing, using activity as RNA, or equating WNN weights with causal effects requires revision.

**Optional extension:** mask one modality in a held-out subset and evaluate how well a reference-based prediction preserves broad labels, keeping the held-out data out of fitting and separating prediction from experimental measurement.

## Module 1 in full

### What you will learn

By the end of this module, you will be able to explain why cells with similar DNA can have different identities; describe chromatin, nucleosomes, promoters, enhancers and transcription factors; state what ATAC-seq measures; read a region-by-cell matrix; write and rerun a small R script; and distinguish a defensible observation from an unsupported regulatory claim.

You need basic DNA/RNA vocabulary, an R installation and a text editor or RStudio. No Signac, Seurat, Python or large dataset download is needed. The core lab uses base R functions available in R 4.x; the later analysis environment has the specific versions in the course guide. The data below are **invented for teaching**, not measured PBMC data. Expected values are arithmetic checks, not output from an executed R analysis here.

**Seven-hour plan:** 90 minutes of lessons/reading; 120 minutes of guided coding; 90 minutes of exercises; 30 minutes for the quiz; 90 minutes for the assignment, checking and reflection. Take a short break between sessions.

### Lesson 1 — DNA sequence and cell identity

Most nucleated cells in an individual share broadly the same genome, yet a lymphocyte and a neuron do very different jobs. Their differences involve which genes are used, when, and in what amounts, together with many other regulatory processes. DNA sequence alone does not tell you the full current state of a cell.

**Chromatin** is DNA together with its associated proteins. A **nucleosome** contains DNA wrapped around a histone protein complex. Nucleosomes package DNA and influence which regions are physically available to regulatory proteins. Chromatin is dynamic: accessibility depends on nucleosome positions, bound proteins, sequence context, cell state and the experimental preparation.

Imagine a book with some pages easy to reach and others partly covered. Access to a page makes reading possible, but does not show that someone is reading it now. Likewise, accessible DNA can support gene regulation without guaranteeing active transcription. The analogy has limits: chromatin is a molecular system, and “open” and “closed” are useful shorthand for probabilistic measurements, not two permanent physical states.

**Pause and explain:** If two immune cells share a genome, what could differ in their ATAC-seq profiles? A good answer mentions regional accessibility and regulatory state, rather than assuming every difference must be a DNA mutation.

### Lesson 2 — Promoters, enhancers and transcription factors

A **promoter** lies around a transcription start site and participates in transcription initiation. A **transcription start site**, or TSS, is the genomic position where an RNA transcript starts; genes can have multiple transcripts and start sites.

An **enhancer** is a DNA region that can regulate transcription in a particular biological context. A candidate enhancer may lie far from the gene it affects, within another gene or on either side of its target. Calling a distal accessible peak an enhancer is a hypothesis unless supported by additional functional evidence.

A **transcription factor**, or TF, is a protein involved in regulating transcription. Many TFs recognize particular sequence patterns, called **motifs**, often with help from other proteins. Related TFs can recognize very similar motifs. Therefore, finding a motif says that a compatible sequence is present; it does not prove that a particular TF was expressed, active or bound at the time of the experiment.

For a later PBMC example, you might inspect accessibility around several B-cell-associated genes and several T-cell-associated genes. A coherent pattern can support a cell-identity hypothesis. One accessible promoter alone is weak evidence, and an RNA marker's usefulness does not guarantee its promoter is a specific ATAC marker.

**Pause and explain:** “A distal peak contains a TF motif, so that TF activates the nearest gene.” Identify the three unsupported steps: motif to binding, nearest gene to target, and association to activation.

### Lesson 3 — What ATAC-seq measures

ATAC-seq stands for **Assay for Transposase-Accessible Chromatin using sequencing**. The enzyme Tn5 inserts sequencing adapters into DNA in a way influenced by accessibility. The resulting fragments are amplified and sequenced. The locations of those fragments or insertion sites provide evidence about accessible chromatin, subject to technical biases and sparse sampling. See the original assay paper in R6 for the experimental basis.

The simplified workflow is: prepare cells or nuclei; expose DNA to Tn5; create and sequence a library; assign reads to genomic locations and barcodes; remove or account for duplicates; summarize accessible regions. Later modules inspect each computational step and its failure modes.

In **bulk ATAC-seq**, DNA from many cells is summarized together. A signal change can arise because accessibility changes within a cell type, because the mixture of cell types changes, or both. In **single-cell/single-nucleus ATAC-seq**, barcodes allow measurements to be assigned to individual capture events. That reveals heterogeneity but yields sparse evidence for each cell. An empty droplet or a doublet can also receive a barcode, so a barcode is not a guarantee of one good cell.

| Assay | Direct material measured | Typical summary | Important limitation |
|---|---|---|---|
| Bulk ATAC-seq | DNA associated with accessible-chromatin sampling | Region-by-sample counts | Cell mixtures can conceal or mimic within-type changes |
| Single-cell ATAC-seq | Barcode-associated accessible DNA fragments/insertion events | Peak-by-cell counts | Sparse sampling and capture artifacts |
| Single-cell RNA-seq | RNA converted to sequenceable material | Gene-by-cell counts | RNA abundance is not TF binding or protein activity |
| Paired RNA–ATAC multiome | RNA and accessible DNA from the same nucleus | Two linked matrices | Pairing improves association analysis but does not establish causation |

### Lesson 4 — A matrix is a summary of evidence

A **peak** is a genomic region identified as enriched for accessible-chromatin signal under a particular peak-calling procedure. Peak definitions depend on the available data and parameters. A peak is not the same thing as a gene; a gene may have many nearby peaks, and many peaks are distant from annotated genes.

In our R pathway, rows are regions and columns are cells. Matrix entry `counts[region, cell]` is the count assigned to that combination under a declared counting convention. Real software may count fragment overlaps or insertion events differently; later labs check this before comparing outputs.

A zero means **no count was observed**. It might reflect low coverage, a region with lower accessibility, a feature-definition issue, or other technical/biological causes. It is insufficient evidence that the region can never be accessible in that cell type.

**Depth** describes how much evidence was recovered, but it also depends on how you define the total. A sum over five teaching regions is just a five-region panel total. It is not the total fragments in the whole library, and it cannot provide a real FRiP denominator. We will keep that distinction visible in our variable names.

### Lesson 5 — Design and replication come before significance

Suppose you measure 5,000 cells from one donor before stimulation and 5,000 afterward. You have many cellular observations but only one donor's response. Cells share genetics, handling and sample-level conditions, so treating all cells as independent donors would overstate evidence about a population.

If four independent donors each supply control and stimulated samples, you have four donor pairs. Pairing helps distinguish a treatment-associated change from baseline donor differences. Technical replicates and repeated sequencing of the same library do not add biological donors.

Good experimental design records donors, conditions, time points, batches and sample handling. Balance conditions across batches when possible, randomize processing order where appropriate and retain sample identity throughout analysis. If every control was processed in one batch and every stimulated sample in another, data alone may be unable to separate batch from condition.

The introductory vendor dataset is useful for learning an atlas workflow. Its one-donor design is why we will use a different dataset for the replicated condition-comparison project.

### Lesson 6 — Four rules for responsible interpretation

1. **Accessibility is not expression.** A region can be accessible without its nearby gene having high RNA abundance.
2. **Motif-associated accessibility is not TF binding.** Motif sequence and chromVAR scores provide candidates, not direct occupancy measurements.
3. **Association is not causation.** A peak and gene can covary because both reflect cell identity or another shared influence.
4. **Cells are not independent biological replicates.** Preserve the donor/sample structure when making population-level claims.

Try these replacements: “observed higher accessibility” instead of “activated the gene”; “candidate TF family” instead of “the TF is bound”; “candidate regulatory association” instead of “proven enhancer target”; and “cells from four donor pairs” instead of “thousands of biological replicates.”

### Guided lab — Inspect a small accessibility matrix

**Question:** What can a tiny matrix tell us about observed accessibility patterns, and what can it not tell us?

**Setup:** create a new folder named `scatac_course` and open an R project or set that folder as your working directory through your editor. Create `module01.R`. Paste the following code into that script, not into a shell. Run one section at a time and predict each output before running it. The code creates only a small `module01_output` folder relative to your working directory.

```r
# Module 1: synthetic teaching data, not experimental PBMC observations.
# Complete base-R script; no extra packages or downloads needed.

# 1. Choose an output folder without changing the working directory.
out_dir <- "module01_output"
dir.create(out_dir, showWarnings = FALSE)

# 2. Metadata: six fictional cells from one fictional donor.
cells <- data.frame(
  cell_id = c("A1", "A2", "A3", "B1", "B2", "B3"),
  group = c("A", "A", "A", "B", "B", "B"),
  donor = rep("demo_donor_1", 6),
  stringsAsFactors = FALSE
)
rownames(cells) <- cells$cell_id

# 3. Five fictional genomic regions by six cells.
# Names and counts were invented to make the concepts easy to inspect.
counts <- rbind(
  shared_region = c(4L, 2L, 3L, 4L, 2L, 3L),
  region_A = c(3L, 2L, 1L, 0L, 0L, 0L),
  region_B = c(0L, 0L, 0L, 3L, 1L, 2L),
  ambiguous_region = c(1L, 0L, 0L, 0L, 1L, 0L),
  unobserved_region = c(0L, 0L, 0L, 0L, 0L, 0L)
)
colnames(counts) <- cells$cell_id

# 4. Inspect and validate the matrix before summarizing it.
str(cells)
dim(counts)
print(counts)
stopifnot(identical(colnames(counts), rownames(cells)))
stopifnot(all(is.finite(counts)), all(counts >= 0))
stopifnot(all(counts == floor(counts)))

# 5. Count only what this five-region panel observed.
cells$panel_counts <- colSums(counts)
cells$detected_regions <- colSums(counts > 0)
zero_fraction <- mean(counts == 0)
regions <- data.frame(
  region = rownames(counts),
  total_count = rowSums(counts),
  detection_fraction = rowMeans(counts > 0),
  row.names = NULL
)
print(cells)
print(regions)
print(zero_fraction)

# 6. Make a simple plot with an honest axis label.
png(file.path(out_dir, "panel_counts.png"),
    width = 1000, height = 650, res = 120)
barplot(
  cells$panel_counts,
  names.arg = cells$cell_id,
  col = ifelse(cells$group == "A", "steelblue", "tan3"),
  xlab = "Fictional cell",
  ylab = "Observed counts across five teaching regions",
  main = "Synthetic accessibility panel"
)
legend("topright", legend = c("Group A", "Group B"),
       fill = c("steelblue", "tan3"), bty = "n")
dev.off()

# 7. Save results and the actual local R environment information.
write.csv(cells, file.path(out_dir, "cell_summary.csv"), row.names = FALSE)
write.csv(regions, file.path(out_dir, "region_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))
```

**What the important commands do:**

- `<-` assigns a value to an object. `out_dir` is an object; `"module01_output"` is a text string.
- `data.frame()` creates a table whose columns may have different types. `rbind()` stacks vectors as matrix rows; all count entries here are integers, indicated by `L`.
- `rownames()` and `colnames()` give explicit identifiers. `identical()` checks both identity and ordering, so a shuffled sample table cannot silently pass.
- `stopifnot()` stops on a failed contract. It is safer to investigate a mismatch than to proceed with wrongly attached labels.
- `counts > 0` creates a TRUE/FALSE matrix. R treats TRUE as 1 and FALSE as 0 in these summaries, so `rowMeans(counts > 0)` is the fraction of cells with a detected count in each region.
- `colSums()` adds down rows to give a value for each cell; `rowSums()` adds across cells to give a value for each region.
- `png()` sends a figure to a file; `dev.off()` closes that graphics device and finishes writing it.
- `file.path()` constructs platform-appropriate file paths. `sessionInfo()` records the environment actually used on your machine.

**Arithmetic checkpoints, not execution claims:**

| Check | Expected value from the supplied synthetic numbers |
|---|---|
| Matrix dimensions | 5 regions × 6 cells |
| Panel totals for A1, A2, A3, B1, B2, B3 | 8, 4, 4, 7, 4, 5 |
| Total counts across the panel | 32 |
| Number of zero entries | 16 of 30, approximately 53.3% |
| Detection fractions for shared, A, B, ambiguous, unobserved regions | 1, 0.5, 0.5, 1/3, 0 |
| Independent real donors in these data | None: all data are synthetic; the fictional design represents one donor |

If your values differ, check row construction, names and whether you changed the toy counts. These values are checks on arithmetic, not target biological outcomes for later PBMC analyses.

**Interpret the lab:** region_A has observed counts in the three A cells and no observed counts in the B cells. This is a descriptive pattern in the teaching panel. It does not prove that region_A regulates a gene, that group A is a real cell type, or that the pattern generalizes across donors. The all-zero region contributes no observed evidence under this toy measurement, and its label does not establish true biological closure.

### Three exercises

**Exercise 1 — Guided practice, approximately 20 minutes.**

Retrieve `counts["region_A", ]`. Fill in a command to calculate the fraction of the six cells with a positive value. Report the total counts for that row and explain the difference between “count fraction” and “fraction of cells with a detected count.” Write one supported observation using the word “observed.”

Hint: use `mean(... > 0)` for detection and `sum(...)` for the row total. Do not divide the row total by the number of cells and call that the detection fraction.

**Exercise 2 — Partly guided analysis, approximately 30 minutes.**

Compare observed counts after scaling each cell's five-region panel to 1,000 counts. This is a teaching exercise about denominators, **not the Signac TF-IDF method** and not a recommended replacement for ATAC normalization.

```r
# Run after the main lab.
stopifnot(all(cells$panel_counts > 0))
panel_relative <- sweep(counts, 2, cells$panel_counts, FUN = "/") * 1000
group_columns <- split(seq_len(ncol(counts)), cells$group)
group_means <- sapply(group_columns, function(j) {
  rowMeans(panel_relative[, j, drop = FALSE])
})
print(group_means)
```

Explain `sweep(..., 2, ...)`: it divides each column by its own panel total; margin 2 means columns. `split()` gives the column indices belonging to each group. `drop=FALSE` preserves matrix structure even for a one-cell group.

Calculate the group-A mean for region_A by hand and compare with R. Explain why multiplying all counts in one cell by two leaves its panel-relative profile unchanged. Explain why adding another region to the panel can change these relative values, even when the original counts are unchanged. This denominator dependence will matter later for normalization and FRiP.

**Exercise 3 — Independent interpretation, approximately 40 minutes.**

Write 200–300 words answering: “What evidence would I need before claiming that region_A is an enhancer driving a group-A-specific gene program?” Include three observations supported by the synthetic matrix, three unsupported claims, a proposed follow-up measurement or perturbation and its control, and the role of independent donors. Do not invent a gene target or a biological cell-type label for the toy groups.

There is no single correct enhancer or TF to identify. The exercise is about the chain of evidence required to support a regulatory explanation.

### Module 1 quiz

Answer without reading the key. Each question is worth two points: the answer and its reasoning.

1. What molecule does ATAC-seq sequence?
2. Is an accessible promoter proof of high RNA expression? Explain.
3. What does a zero in a sparse peak-by-cell matrix mean?
4. Are 1,000 cells from one donor 1,000 independent biological replicates? Identify the relevant unit.
5. Does a distal accessible region identify its target gene? What extra evidence could help?

### Troubleshooting and interpretation checklist

| Symptom or mistake | Likely explanation and next step |
|---|---|
| `object 'counts' not found` | The creation section has not run, or the object name differs in capitalization. Run the script from the top. |
| The console shows `+` and waits | A quote, parenthesis or other expression is unfinished. Cancel that incomplete command and rerun the intact section. |
| `subscript out of bounds` | A requested row/column name is missing. Inspect `rownames(counts)` and `colnames(counts)`. |
| File appears missing | Check `getwd()` and `list.files(out_dir)`. Output paths are relative to the current project. |
| Plot file is empty/incomplete | Ensure the plotting commands completed and `dev.off()` ran. Do not close an unrelated device blindly. |
| Columns and labels do not agree | Check the `identical()` assertion; explicitly reorder metadata by cell ID rather than guessing its order. |
| One cell has more panel counts | It has more observed counts in these five regions. That alone does not identify its quality or total whole-library depth. |
| All-zero region is called “closed” | Reword as “no counts observed.” The experiment has limited sampling and this teaching matrix is synthetic. |
| Group A is assigned a real cell type | The toy group labels contain no validated biological identity. Actual annotation requires independent marker evidence. |
| A result is called statistically significant | No inferential test has been justified or performed. The lab is descriptive and synthetic. |

### Assignment and completion criteria

Submit `module01.R`, `cell_summary.csv`, `region_summary.csv`, `panel_counts.png`, `sessionInfo.txt` and a 300–500-word interpretation. Include your Exercise 2 calculation and a short note from Exercise 3; its text can be incorporated into the interpretation rather than duplicated.

The report should answer four questions: What was counted? What does a zero mean? Which group-level patterns are observed? What cannot be concluded about expression, TF binding, causation or population-level differences?

**Rubric:** correct/replayable code 30; interpretation 25; clear denominators and inference limits 20; saved files/environment 15; readable plot/report 10. Completion requires at least 70/100, the numerical checks above, no confusion of regions with cells, and no claim that synthetic values are real observations. Aim for 8/10 on the quiz after reviewing mistakes. Finish by writing one sentence about what still feels unclear; that will guide the next module.

### Module 1 exercise solutions and quiz key

**Exercise 1:** `mean(counts["region_A", ] > 0)` gives 0.5 because three of six cells have positive counts. `sum(counts["region_A", ])` gives 6. Its share of all 32 panel counts is 6/32 = 0.1875; that is a different denominator and question. Supported wording: “Counts for region_A were observed in all three synthetic A cells and none of the synthetic B cells.”

**Exercise 2:** the three group-A scaled values for region_A are 3/8×1000 = 375, 2/4×1000 = 500 and 1/4×1000 = 250; their mean is 375. Multiplying a whole cell's vector by two doubles both numerator and denominator, leaving its proportions unchanged. Adding a new counted region changes the denominator, so relative values for the original regions can decrease even though their raw counts do not. This demonstrates compositional dependence, not that scaling has recovered a biological truth.

**Exercise 3:** acceptable observations include the shared region's detection in all six cells, region_A's observed group pattern and the differing panel totals. Unsupported claims include a named target gene, a particular TF binding and a population-wide cell-type difference. A stronger design would include independent donors, a defined candidate region and target, an appropriate perturbation with a negative control, and RNA or other relevant functional readouts. A perturbation's off-target effects and changes in cell viability/state would also need assessment. Orthogonal RNA/binding/contact measurements strengthen context but do not all establish causation by themselves.

**Quiz answers:**

1. DNA fragments, because the assay samples transposase-accessible chromatin rather than RNA directly.
2. No. Accessibility can be compatible with transcription without determining current transcript abundance.
3. No count was observed for that feature/cell under the measurement and processing rules; biological absence is only one possible explanation.
4. No. The population-level biological source is the donor; cells are nested observations. One donor is still one donor.
5. No. Distance and accessibility nominate possibilities. Paired RNA associations, chromatin-contact data, binding data and controlled perturbations can provide different kinds of additional evidence; causal conclusions require appropriate experimental support.

### Continue to Module 2

The next module reuses the saved table and adds R filtering, functions, sample-table joins, missing-value handling and statistical reasoning about donors. When requesting Module 2, include any error message or tell me which readiness task was hardest. Keep your files and decision log so that each new module extends the preceding work.
