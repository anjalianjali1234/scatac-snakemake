# Datasets, resources and contracts

No real biological dataset is bundled. Download links below were checked against the original course sources; no biological download or execution is implied.

| Use | Dataset | Organism / assay / assembly | Approximate scale and files |
|---|---|---|---|
| Default mechanics | Included synthetic panel | Invented; no genome | 6 cells × 5 regions; counts CSV and cells TSV |
| Atlas | [10x PBMC ATAC v2](https://www.10xgenomics.com/datasets/10k-human-pbmcs-atac-v2-chromium-controller-2-standard), exact downloads in [Signac 1.16](https://stuartlab.org/signac/1.16.0/articles/pbmc_vignette) | Human, scATAC, hg38; one donor | About 10,246 recovered nuclei; peak H5, `singlecell.csv`, fragments and index from Cell Ranger ATAC 2.1.0 |
| Paired condition inference | [GSE178428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178428) | Human, SureCell scATAC, hg19 | Full study: 42 libraries, about 67,581 cells; select 8 libraries/4 donor pairs. Full fragment archive about 6.7 GB. Author peak BED and cell metadata are also available; retained subset size must be measured |
| Multiome capstone | [10x paired PBMC](https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0), [Signac archived-file tutorial](https://stuartlab.org/signac/articles/pbmc_multiomic) | Human, paired RNA–ATAC, hg38; one donor | Use one processing version. Archived ARC 1.0.0 example: about 11,909 barcodes; feature H5 plus ATAC fragments/index |
| Trajectory / ArchR extension | [GSE129785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785), [ArchR subset](https://www.archrproject.com/articles/Articles/tutorial.html) | Human hematopoiesis, scATAC, hg19 | Prefer approximately 0.5 GB ArchR teaching download for mechanics; not a matched-donor treatment design |

## Computing budget

Planning estimates, **not measured performance**. Allow extra disk for environment installation and temporary copies. Download and installation time are excluded.

| Route | RAM | Free disk | Analysis time |
|---|---:|---:|---|
| Synthetic demo | 1–2 GB including engine | 1 GB including Python packages | Seconds after installation |
| 2,000-cell atlas | 8–16 GB; 16 GB preferred | 15–30 GB | 15–90 minutes, potentially longer for fragment scans and doublet calling |
| Full PBMC atlas | 16–32 GB | 25–50 GB | 1–4 hours |
| 8-library preparation and one-lineage pseudobulk | 32 GB preferred; 64 GB for comfortable joint preprocessing | 20–50 GB | 2–8 hours preparation, then minutes for the model |
| Full paired multiome capstone | 32–64 GB | 40–80 GB | 2–8 hours |
| Optional full FASTQ lab | 64–128 GB, 8–16 cores | 250–500 GB | 8–24+ hours |

Atlas `max_cells: 2000` reduces downstream work after sparse import; it does not shrink the input H5 or fragments. The smallest real-data mechanics option is Signac's bundled `atac_small` (100 cells, 323 features, hg19); do not feed it to the hg38 atlas config. No GPU is required.

## Atlas contract

Use the four exact files in `config/downloads.atlas.tsv`. The index must be named `fragments.tsv.gz.tbi` beside the position-sorted, bgzip-compressed fragment file. All barcodes in the H5 must occur in the matching metadata and fragment source. The metadata must contain `passed_filters` and `peak_region_fragments`.

Vendor FRiP is computed as `peak_region_fragments / passed_filters`. `nCount_peaks` is a peak-count sum; it is not interchangeable with all usable fragments. State the convention when comparing reports. The workflow uses hg38 standard chromosomes and EnsDb.Hsapiens.v86 annotations. No coordinate conversion or assembly repair is attempted.

## Paired pseudobulk contract

Prepare these after per-capture QC/doublet handling, joint peak-set construction and lineage annotation. Use a condition-blind, documented consensus/union procedure; remove inappropriate/blacklist intervals and recount every retained cell against that same interval set. The course covers this in the peak-set and donor-inference modules.

- `counts.mtx.gz`: Matrix Market sparse **raw integer** peak × cell counts. No TF-IDF, log counts, imputed values or batch-corrected representation.
- `features.tsv`: header `feature_id`, `chrom`, `start`, `end`, `genome`, in matrix row order. Use BED zero-based, half-open coordinates; identifiers must be unique. Any extra annotations are retained in the result table.
- `barcodes.tsv`: one globally unique cell ID per line, no header, in matrix column order. Prefix library IDs where original barcodes repeat.
- `cells.tsv`: header `cell_id`, `sample`, `donor`, `condition`, `celltype`, in exactly the barcode order. Retain an explicit mapping to the original capture/barcode for fragment access. Use reviewed lineage labels; this repository does not assert that GEO already provides the chosen annotation.
- `config/samples.paired.tsv`: every sample mapped to its independent donor, condition and assembly. Its donor/condition values must agree with every cell.

The example contrast uses Control6h GSM5390937–GSM5390940 and LPS6h GSM5390957–GSM5390960, with the same four donors in order and no Golgi inhibitor. Verify this mapping against GEO and author metadata when constructing your own files. Do not combine these hg19 files with the hg38 atlas peaks.

The supplied model aggregates within donor × condition × chosen lineage, requires all donor pairs and the configured minimum cell support, filters low-count peaks with edgeR `filterByExpr`, and fits `~ donor + condition`. All tested peaks remain in the result table; the separate feature filter records the testing universe. FDR adjustment covers the tested peaks for this contrast, not arbitrary additional lineages/contrasts. TMM can mask global directional changes; inspect normalization and consider external controls or a justified alternative where appropriate.

Hashes are generated from the actual local files by `scripts/pin_inputs.py`. None are invented for undownloaded data. Keep source URLs, accession mappings, data usage terms and preprocessing provenance alongside the files.
