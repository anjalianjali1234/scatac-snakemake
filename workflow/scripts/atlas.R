# Runnable analysis for the explicitly pinned Signac 1.16 / Seurat 5 API.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[1] %in% c("qc", "lsi", "embedding"))
suppressPackageStartupMessages({library(Signac); library(Seurat); library(jsonlite)})
versions <- c(Signac = "1.16.0", Seurat = "5.3.0", SeuratObject = "5.2.0")
for (p in names(versions)) {
  if (as.character(packageVersion(p)) != versions[[p]]) stop("Unexpected version of ", p)
}
cfg <- fromJSON(args[2])$parameters
a <- cfg$atlas
set.seed(cfg$seed)
out <- "results/atlas"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_session <- function(stage) writeLines(capture.output(sessionInfo()), file.path(out, paste0("session_", stage, ".txt")))

if (args[1] == "qc") {
  suppressPackageStartupMessages({
    library(GenomeInfoDb); library(GenomicRanges); library(EnsDb.Hsapiens.v86)
    library(SingleCellExperiment); library(scDblFinder); library(BiocParallel)
  })
  sm <- read.delim(cfg$samples, stringsAsFactors = FALSE)
  stopifnot(nrow(sm) == 1L, sm$genome == "hg38",
            sm$fragments_index == paste0(sm$fragments, ".tbi"))
  counts <- Read10X_h5(sm$counts_h5)
  if (is.list(counts)) stop("Expected an ATAC-only peak matrix, not a multiome feature list")
  meta <- read.csv(sm$metadata_csv, row.names = 1, check.names = FALSE)
  stopifnot(!anyDuplicated(colnames(counts)), !anyDuplicated(rownames(meta)),
            all(colnames(counts) %in% rownames(meta)),
            all(c("passed_filters", "peak_region_fragments") %in% colnames(meta)))
  original_n <- ncol(counts)
  # Subsampling does not reduce the raw download or initial H5 import memory.
  ids <- colnames(counts)
  if (a$max_cells > 0L && length(ids) > a$max_cells) ids <- sort(sample(ids, a$max_cells))
  counts <- counts[, ids, drop = FALSE]
  meta <- meta[ids, , drop = FALSE]
  stopifnot(all(meta$passed_filters > 0), all(meta$peak_region_fragments >= 0),
            all(meta$peak_region_fragments <= meta$passed_filters))
  ranges <- StringToGRanges(rownames(counts), sep = c(":", "-"))
  genome(ranges) <- "hg38"
  keep <- as.character(seqnames(ranges)) %in% paste0("chr", c(1:22, "X", "Y"))
  counts <- counts[keep, , drop = FALSE]; ranges <- ranges[keep]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  chrom <- CreateChromatinAssay(counts, ranges = ranges, min.cells = 1, min.features = 0,
                               fragments = sm$fragments, annotation = annotations)
  obj <- CreateSeuratObject(chrom, assay = "peaks", meta.data = meta)
  obj$sample <- sm$sample; obj$donor <- sm$donor; obj$condition <- sm$condition
  obj$original_capture_cells <- original_n
  obj <- NucleosomeSignal(obj)
  obj <- TSSEnrichment(obj, fast = TRUE)
  obj$frip_vendor <- obj$peak_region_fragments / obj$passed_filters
  obj$blacklist_fraction <- FractionCountsInRegion(obj, assay = "peaks", regions = blacklist_hg38_unified)
  # An explicit capture-level prior avoids deriving a lower rate from a cell subset.
  # This remains an imperfect doublet classifier; inspect its output and sensitivity.
  sce <- SingleCellExperiment(list(counts = GetAssayData(obj, assay = "peaks", layer = "counts")))
  sce <- scDblFinder(sce, aggregateFeatures = TRUE, dbr = a$expected_doublet_rate,
                     BPPARAM = SerialParam(RNGseed = cfg$seed))
  obj$doublet_score <- sce$scDblFinder.score
  obj$doublet_class <- as.character(sce$scDblFinder.class)
  m <- obj[[]]
  write.csv(data.frame(cell_id = rownames(m), m, row.names = NULL), file.path(out, "qc_metrics.csv"), row.names = FALSE)
  pdf(file.path(out, "qc.pdf"), width = 11, height = 8)
  par(mfrow = c(2, 3))
  for (metric in c("passed_filters", "TSS.enrichment", "nucleosome_signal", "frip_vendor", "blacklist_fraction", "doublet_score")) {
    x <- m[[metric]]; x <- x[is.finite(x)]
    if (!length(x)) stop("No finite values for QC metric ", metric)
    hist(x, breaks = 50, main = metric, xlab = metric, col = "steelblue")
  }
  dev.off()
  saveRDS(obj, file.path(out, "qc.rds")); save_session("qc")
}

if (args[1] == "lsi") {
  q <- a$qc
  if (!isTRUE(q$reviewed) || !nzchar(trimws(q$rationale))) stop("Review QC and record a rationale in config/atlas.yaml before filtering")
  obj <- readRDS(file.path(out, "qc.rds")); m <- obj[[]]
  checks <- data.frame(
    fragment_min = m$passed_filters >= q$min_fragments,
    fragment_max = m$passed_filters <= q$max_fragments,
    tss = m$TSS.enrichment >= q$min_tss,
    nucleosome = m$nucleosome_signal <= q$max_nucleosome,
    frip = m$frip_vendor >= q$min_frip,
    blacklist = m$blacklist_fraction <= q$max_blacklist_fraction,
    singlet = m$doublet_class == "singlet")
  checks[is.na(checks)] <- FALSE
  keep <- apply(checks, 1, all)
  write.csv(data.frame(cell_id = colnames(obj), checks, retained = keep), file.path(out, "cell_filtering.csv"), row.names = FALSE)
  if (sum(keep) <= a$lsi$n_components + 1) stop("Too few retained cells for the requested LSI rank")
  obj <- subset(obj, cells = colnames(obj)[keep])
  # Keep raw counts in the assay; TF-IDF is a separate representation.
  obj <- RunTFIDF(obj, method = 1)
  obj <- FindTopFeatures(obj, min.cutoff = "q0")
  if (length(VariableFeatures(obj)) <= a$lsi$n_components + 1) stop("Too few nonzero peaks for LSI")
  obj <- RunSVD(obj, n = a$lsi$n_components)
  embedding <- Embeddings(obj, "lsi")
  correlation <- apply(embedding, 2, function(x) cor(x, log10(obj$passed_filters), method = "spearman"))
  write.csv(data.frame(component = seq_along(correlation), spearman_log10_fragments = correlation), file.path(out, "depth_correlations.csv"), row.names = FALSE)
  pdf(file.path(out, "lsi_depth.pdf"), width = 9, height = 5)
  barplot(correlation, names.arg = seq_along(correlation), xlab = "LSI component", ylab = "Spearman correlation with log10 fragments")
  abline(h = 0); dev.off()
  saveRDS(obj, file.path(out, "lsi.rds")); save_session("lsi")
}

if (args[1] == "embedding") {
  lsi <- a$lsi
  if (!isTRUE(lsi$reviewed) || !nzchar(trimws(lsi$rationale))) stop("Review LSI depth associations and document dimensions before clustering")
  obj <- readRDS(file.path(out, "lsi.rds"))
  dims <- as.integer(lsi$dimensions)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = dims)
  obj <- FindClusters(obj, resolution = lsi$resolution, random.seed = cfg$seed, algorithm = 1)
  obj <- RunUMAP(obj, reduction = "lsi", dims = dims, seed.use = cfg$seed)
  activity <- GeneActivity(obj)
  obj[["ACTIVITY"]] <- CreateAssayObject(counts = activity)
  obj <- NormalizeData(obj, assay = "ACTIVITY", normalization.method = "LogNormalize", scale.factor = 10000)
  m <- obj[[]]; umap <- Embeddings(obj, "umap")
  write.csv(data.frame(cell_id = rownames(m), m, umap, row.names = NULL), file.path(out, "cells.csv"), row.names = FALSE)
  pdf(file.path(out, "embedding_and_activity.pdf"), width = 10, height = 7)
  print(DimPlot(obj, reduction = "umap", label = TRUE))
  print(FeaturePlot(obj, features = "passed_filters"))
  DefaultAssay(obj) <- "ACTIVITY"
  markers <- intersect(c("CD3D", "MS4A1", "LYZ", "NKG7"), rownames(obj[["ACTIVITY"]]))
  if (length(markers)) print(FeaturePlot(obj, features = markers))
  dev.off()
  DefaultAssay(obj) <- "peaks"
  saveRDS(obj, file.path(out, "atlas.rds")); save_session("embedding")
  write_json(list(synthetic = FALSE, cells_retained = ncol(obj), peaks = nrow(obj),
    clusters = length(unique(Idents(obj))), genome = "hg38", annotation = "EnsDb.Hsapiens.v86",
    interpretation = "Exploratory clusters and gene activity, not validated cell types or measured RNA. One donor cannot support a population-level condition test. Subsetting changes neighborhoods and doublet sensitivity.",
    artifacts = c("qc_metrics.csv", "qc.pdf", "cell_filtering.csv", "depth_correlations.csv", "lsi_depth.pdf", "embedding_and_activity.pdf", "cells.csv", "atlas.rds", "session_qc.txt", "session_lsi.txt", "session_embedding.txt")),
    file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
}
