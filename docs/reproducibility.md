# Reproducibility and extension guide

## What is recorded

Each mode writes a provenance JSON containing effective parameters, input SHA256 hashes and sizes, analysis-source hashes, the available Git commit/dirty state, Python/Snakemake versions and platform. R jobs write `sessionInfo`; Snakemake writes rule logs and benchmark TSVs. Random seeds are explicit. Multi-threaded numerical libraries are limited in biological jobs; different operating systems and library builds may still produce small floating-point differences.

`requirements.lock` records the exact installed Python packages. It is a version freeze, not a pip artifact-hash lock. Python 3.12.14 and the operating system remain part of the environment contract. Conda YAMLs pin top-level analysis versions. They become full platform locks only when all dependencies have been solved and the resulting explicit specs have been recorded. Consult `validation.md` for what was actually installed/tested.

For a working Conda environment, save its real platform-specific explicit specification alongside the YAML. Snakemake supports the naming convention below:

```bash
# Replace the path with the actual environment created for this YAML.
# This is a command template: the placeholder is not a real environment name.
conda list --explicit --prefix /absolute/path/to/solved/environment > workflow/envs/atac.linux-64.pin.txt
```

Review that it contains public package URLs without embedded credentials. Commit the solved pin file, environment recipe and actual session information. Do not create an `renv.lock` or explicit package lock by guessing transitive versions. The course also explains an alternative R/renv installation path.

Input hashes deliberately fail on changed data. Investigate a mismatch before accepting a new data version. To accept an intentional revision, keep the old lock in Git history, record why the source changed, and create a new checked lock. The pinning script refuses to overwrite an existing one. Checksums establish content identity, not publisher authenticity or data availability.

## Analyst decisions remain explicit

The atlas first generates QC metrics without retaining or removing cells by universal cutoffs. It then requires a written threshold justification; LSI output similarly precedes a documented dimension choice. The initial QC thresholds and doublet prior are teaching suggestions. Rerun sensitivity analyses when a decision materially changes the cell population or interpretation. A dimensional component can reflect both depth and real biology; no rule automatically removes component 1.

Doublets are scored using scDblFinder's ATAC feature aggregation and an explicit capture-level prior. A random cell subset does not have a proportionally smaller true doublet rate. Subsetting can still change class balance and detection sensitivity, and homotypic doublets remain difficult to detect. Review the full capture where feasible.

The paired model checks donor/sample attribution before testing. It does not fix confounding, invent replication, infer missing donor identities or automatically discard poorly supported pairs. A study design change should be explicit and versioned.

## Extend the course into new rules

1. Keep the complete [course](course.md) as the learning and interpretation guide. Deliver subsequent modules as lessons one at a time; their lesson snippets are not automatically all pipeline-ready.
2. Define the new rule's input artifact and data contract first: assembly, raw/normalized state, peak universe, cell IDs and independent samples.
3. Put its parameterized R/Python script in `workflow/scripts`, add a rule and separate environment where needed, and list every output consumed downstream.
4. Add tiny tests for meaningful failure modes or known quantitative behavior. Do not claim a smoke test validates an entire biological method.
5. Preserve true upstream inputs and provenance; emit actual figures and full statistics. Keep human interpretation separate from automatic cluster names.

Motifs/chromVAR, peak–gene associations, reference mapping, trajectories, ArchR and paired multiome WNN are planned extensions described in the course and capstones. FASTQ alignment/barcode/duplicate processing is intentionally a separate advanced lab, using a documented Cell Ranger ATAC installation and matching reference. No proprietary software, reference genome or multi-GB dataset is silently installed by the default example.

## GitHub checks

The CI workflow installs the versioned Python engine and R 4.5.1, tests input contracts and R/Python agreement, parses all R scripts, executes the synthetic workflow and checks that a second run has nothing pending. It does not download real data or claim to validate atlas results. GitHub Actions itself must be enabled by the repository/account; a committed workflow is not evidence that a remote run passed.

Official references: [Snakemake deployment and environment pin files](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html), [configuration validation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html), [Signac 1.16 workflow](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette), [scDblFinder ATAC methods](https://plger.github.io/scDblFinder/articles/scATAC.html), [edgeR](https://bioconductor.org/packages/3.22/bioc/html/edgeR.html).
