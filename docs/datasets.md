# Public data and input requirements

[Start the tutorial](tutorial.md)

## Main tutorial dataset

| Field | Value |
|---|---|
| Source | [10x Genomics: 10k Human PBMCs, ATAC v2, Chromium Controller](https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard) |
| Organism / material | Human peripheral blood mononuclear cells; nuclei from one healthy donor |
| Assay | Chromium single-cell ATAC v2 |
| Processing release | Cell Ranger ATAC 2.1.0 |
| Genome assembly | GRCh38/hg38 |
| Vendor-reported nuclei | 10,246 |
| Tutorial selection | 1,000 barcodes sampled with seed 42 before QC; retained count is an output |
| License | CC BY 4.0, as stated by 10x Genomics |
| Exact sources | [Download manifest](../config/downloads.atlas.tsv), consistent with [Signac 1.16](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette) |

This dataset is useful because the processed count matrix, barcode metadata and indexed fragments are available together. Starting from those files makes the tutorial accessible without a separate FASTQ alignment setup. One donor supports an exploratory example, not a replicated condition comparison.

The following sizes were measured from downloaded files on **2026-11-19**. GB/MB below are decimal units.

| Local filename under `data/pbmc10k/` | Bytes | Approximate size |
|---|---:|---:|
| `filtered_peak_bc_matrix.h5` | 169,866,643 | 170 MB |
| `singlecell.csv` | 32,821,340 | 33 MB |
| `fragments.tsv.gz` | 2,641,500,673 | 2.64 GB |
| `fragments.tsv.gz.tbi` | 1,152,837 | 1.2 MB |

The four inputs total approximately **2.85 GB**. The downloader saves them locally rather than committing them to Git. `config/checksums.atlas.json` contains their measured SHA256 hashes and the sample-sheet hash; these identify the downloaded bytes, not an independently authenticated publisher signature.

## File contracts

Use all four files from the same release. The index must be named `fragments.tsv.gz.tbi` beside the sorted, bgzip-compressed fragment file. H5 barcodes must match the metadata and fragment source. The metadata must contain `passed_filters` and `peak_region_fragments`.

Vendor FRiP is `peak_region_fragments / passed_filters`. Peak-count sums are not the same as all usable fragments. The blacklist fraction is calculated from the peak matrix, so state its denominator when comparing analyses. This workflow uses standard nuclear chromosomes and EnsDb.Hsapiens.v86 hg38 gene annotations; that annotation version need not be identical to the original alignment reference's gene annotation. No genome conversion or automatic barcode repair is attempted.

## Computing budget

The allowances below are planning estimates, excluding download speed and first environment installation. Consult [validation](validation.md) for any measured execution times.

| Route | RAM allowance | Free-disk allowance | Analysis-time allowance |
|---|---:|---:|---|
| Synthetic installation check | 1–2 GB | 1 GB | Seconds after installing the engine |
| 1,000-cell public-data tutorial | 16 GB recommended | 20 GB | About 15–90 minutes; fragment scans may take longer |
| Full 10,246-cell capture | 32 GB recommended | 30–50 GB | About 1–4 hours |
| Prepared paired pseudobulk model | Depends on the sparse input; 16 GB is a starting allowance | Input size plus several GB | Usually minutes for the model; preprocessing is separate |

`max_cells` reduces work after sparse H5 import; it does not reduce the raw download. A smaller subset may miss rare cell types and changes neighborhood structure and doublet sensitivity. If a machine cannot handle the import, use a capable machine to create a documented subset or use the synthetic installation check. Do not describe synthetic tests as biological results.

## Paired-donor extension

The optional `pseudobulk` workflow uses [GSE178428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178428): human SureCell scATAC, hg19. The supplied [sample sheet](../config/samples.paired.tsv) selects four matched donors at 6 h, control versus LPS, without Golgi inhibitor: Control6h GSM5390937–GSM5390940 and LPS6h GSM5390957–GSM5390960. Verify author metadata when preparing your own inputs. These hg19 files are not interchangeable with the tutorial's hg38 peaks.

This extension starts **after** per-capture QC/doublet handling, annotation, shared peak-set construction and recounting. Use the same documented peak universe for every retained cell. Separately called peak matrices must not simply be concatenated.

Required inputs:

- `counts.mtx.gz`: sparse raw integer peak-by-cell counts; no TF-IDF, log counts or integrated values.
- `features.tsv`: header `feature_id`, `chrom`, `start`, `end`, `genome`, in matrix row order. Use unique IDs and BED zero-based, half-open intervals.
- `barcodes.tsv`: globally unique cell IDs, one per line without a header, in matrix column order.
- `cells.tsv`: header `cell_id`, `sample`, `donor`, `condition`, `celltype`, in exactly barcode order. Preserve the mapping to original capture/barcode IDs.
- `config/samples.paired.tsv`: the sample-to-donor/condition/assembly mapping. It must agree with every cell's metadata.

After preparing and reviewing those inputs, edit `config/pseudobulk.yaml`, then run:

```bash
python scripts/pin_inputs.py --config config/pseudobulk.yaml
snakemake --configfile config/pseudobulk.yaml --cores 2 --sdm conda --dry-run
snakemake --configfile config/pseudobulk.yaml --cores 2 --sdm conda
```

The model sums raw counts within donor × condition × selected cell type, checks pair/cell support, filters low-count peaks, and uses edgeR with `~ donor + condition`. All tested peaks are retained in the statistics table with effect sizes and FDR. Missing donor pairs and conflicting sample metadata cause an error. Individual cells cannot create independent biological replication.

Additional covariates, unpaired designs and repeated time courses need a reviewed model change. TMM normalization can obscure global directional changes, and a completely confounded design cannot be repaired by adding terms. Real GSE178428 preprocessing and biological condition testing are not claimed as executed here; the extension has a separate synthetic integration test.
