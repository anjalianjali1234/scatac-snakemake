# Validation record

Checked on 2026-11-18 on Linux x86_64 with Python 3.12.14, Snakemake 9.27.2 and R 4.5.1.

## Executed checks

- The default Snakemake DAG parsed and completed; the output checker verified the expected synthetic counts, zero fraction, generated files and input provenance.
- A second dry-run after completion reported no pending jobs.
- Eleven tests passed, including R/Python numeric agreement; changed-input rejection; cell-ID/order and genome checks; complete donor pairing; sample/donor consistency; and rejection of cells masquerading as biological replicates.
- The paired integration test executed Snakemake, sparse aggregation, edgeR 4.8.2, diagnostics and HTML reporting on **synthetic** data: four donor pairs, two conditions, six cells per group and 1,000 regions. It checked group sizes, direction of a deliberately planted shift and valid FDR values. These are software tests, not biological discoveries or a calibration study.
- All analysis R scripts parsed under R 4.5.1.
- The pseudobulk Conda recipe was successfully installed. Its actual Linux package URLs and MD5 artifact identities were exported to `workflow/envs/pseudobulk.linux-64.pin.txt`. Some compiled Bioconductor packages report build-time R patch versions newer than 4.5.1; the tested pipeline completed in the recorded environment.

## Remaining validation boundaries

- The real PBMC atlas has not been executed. Its atlas environment recipe and API choices require installation and a real-data pilot before scientific use; no atlas output is claimed here.
- Public downloads, FASTQ processing, consensus-peak preparation, biological cell annotations, motif inference, trajectories and paired multiome analyses have not been performed by these tests.
- GitHub CI results must be inspected after the repository is published. A workflow file by itself is not evidence of a remote successful run.
- The tiny demonstration cannot validate realistic doublet sensitivity, QC thresholds, clustering stability, differential-accessibility calibration or causal regulatory claims.

The full course documents additional sensitivity analyses and biological checks needed for a scientific report.
