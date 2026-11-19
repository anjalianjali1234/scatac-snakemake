# Reproducibility and troubleshooting

[Start the tutorial](tutorial.md) · [Validation record](validation.md)

## What is recorded

Each run writes its effective parameters, actual input SHA256 hashes and sizes, analysis-source hashes, Python/Snakemake versions, platform and available Git commit/dirty state to `provenance.json`. R jobs write `sessionInfo`; Snakemake keeps rule logs and benchmark TSVs. Seeds are explicit, and numerical-library threads are limited for the biological steps.

A seed is one part of reproducibility. Input identity, feature definitions, annotation version, library builds and analyst decisions also matter. Numerical results can vary slightly across platforms even with the same seed.

## Environments

`requirements.lock` records the exact Python package versions from the installed engine environment. It is a version freeze, not a pip artifact-hash lock. The R recipes pin top-level versions; the supplied Linux explicit specifications record actual installed Conda package URLs and artifact identities. The `atac` environment includes the HDF5 reader required for the 10x H5 input.

Use Linux x86_64 for the reference path. Snakemake recognizes `atac.linux-64.pin.txt` alongside `atac.yaml` when deploying on that platform. Other platforms may solve the YAML recipe differently and require their own validation. Some Bioconductor data packages also download annotation resources in post-install steps; those resources depend on the package's recorded URLs and integrity checks.

An alternative to Snakemake-managed environments is to create and activate the explicit Linux environment yourself:

```bash
conda create --name scatac-atac --file workflow/envs/atac.linux-64.pin.txt
conda activate scatac-atac
```

Keep Snakemake available from the Python engine environment, for example by invoking its absolute executable path. With `Rscript` from the exact analysis environment on `PATH`, run the same Snakemake commands **without** `--sdm conda`. Do not accidentally run an unrelated system R.

The [validation record](validation.md) specifies which installation path was tested. A successful package solve alone is not evidence that the biological workflow executed.

## Input locks

For the provided public dataset, use the committed `config/checksums.atlas.json`. It records actual downloaded content identities. The workflow fails on mismatches instead of accepting a changed file silently.

When adapting the tutorial to new data, preserve the original manifest, create a new sample/config file and choose a new checksum path. Review source accession, release, assembly and file pairing before running:

```bash
# Command template: replace the config path with your real reviewed config.
python scripts/pin_inputs.py --config config/my_dataset.yaml
```

The script refuses to overwrite an existing lock. A checksum cannot recover a deleted file or establish the publisher's authenticity. Keep a durable copy of the approved inputs and the original source URLs.

## Efficient reruns

The atlas rules separate parameters by stage. QC depends on input files, the sample sheet, input lock, selected-cell count, doublet prior and seed. Filtering/LSI depends on the QC result, filter choices and requested rank. The final embedding depends on LSI output, selected dimensions, resolution and seed.

The full run manifest is audited, but its timestamp alone does not trigger another expensive QC scan when only a downstream review decision changes. Snakemake still tracks the explicit stage parameters and direct data inputs. Editing an analysis script or environment can legitimately rerun the affected work.

Use a dry-run to inspect the reason for any rerun:

```bash
snakemake --configfile config/atlas.yaml --cores 2 --sdm conda --dry-run
```

Do not use `--touch`, remove checksum checks, or fabricate output files to force a workflow to appear complete. Preserve failed-job logs before fixing the cause.

## Scientific decisions

The guided config requires a QC rationale and an LSI rationale. These are local scientific review checkpoints; they do not require a software account approval. Write observations from the figures, the choice you made, and an important alternative or limitation.

Doublet classification is imperfect, especially for similar cell types and subsets. QC thresholds can remove genuine populations. Depth-associated components can contain biology. Cluster resolution and cell selection affect granularity. Assess sensitivity where a choice changes the conclusion.

The gene-activity assay is called `ACTIVITY` because it is not RNA. Record multiple marker observations and contradictory evidence in your annotation table. Keep uncertain labels broad or unassigned. The single-donor tutorial does not establish a condition effect.

## Common problems

| Problem | Action |
|---|---|
| 404 on the repository | Sign into the owning or invited GitHub account; the repository is private |
| H5 import asks for `hdf5r` | Use the updated pinned atlas environment, not a previous/incomplete environment |
| Conda cannot create an environment | Save the solver log; verify platform, channels and disk space before changing versions |
| `Rscript` loads unexpected packages | Use `--sdm conda`, or explicitly check the activated environment and `sessionInfo()` |
| Interrupted download | Rerun `scripts/download_atlas.py`; it resumes `.partial` files |
| Hash mismatch | Check source/release and file completeness; do not regenerate hashes as a shortcut |
| QC/LSI review error | Review diagnostics and fill the relevant config rationale |
| Memory exhaustion | Reduce downstream cells or use a larger machine; full H5 import and fragments remain |
| A report image does not open | Keep `report.html` beside its generated PNGs/tables; use the PDF for a portable figure bundle |

## GitHub checks and extensions

The CI workflow installs the pinned Python engine and R, tests input contracts and R/Python agreement, parses R scripts, runs the synthetic installation check and verifies an unchanged rerun. It does not download several GB of public data on every push. The optional edgeR integration test executes when its R dependencies are available; skipped tests must not be reported as passed.

For a new method, first define its inputs, assembly, cell IDs, feature universe and statistical unit. Put the parameterized script in `workflow/scripts`, declare real dependencies/outputs, pin the environment and test meaningful failure modes. Preserve raw inputs and emit actual results rather than illustrative figures.

References: [Snakemake deployment](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html), [versioned Signac workflow](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette), [scDblFinder ATAC methods](https://plger.github.io/scDblFinder/articles/scATAC.html).
