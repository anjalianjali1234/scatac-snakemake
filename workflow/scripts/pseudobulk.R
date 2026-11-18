# Raw common-peak counts -> donor-condition sums -> paired edgeR QL model.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
suppressPackageStartupMessages({library(Matrix); library(edgeR); library(jsonlite)})
stopifnot(as.character(packageVersion("edgeR")) == "4.8.2")
cfg <- fromJSON(args[1])$parameters; p <- cfg$pseudobulk
set.seed(cfg$seed)
out <- "results/pseudobulk"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tab <- function(x, name, rows = FALSE) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = rows)
features <- read.delim(p$features, stringsAsFactors = FALSE, check.names = FALSE)
barcodes <- readLines(p$barcodes)
cells <- read.delim(p$cells, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("feature_id", "chrom", "start", "end", "genome") %in% names(features)),
          !anyDuplicated(features$feature_id), !anyDuplicated(barcodes),
          !anyDuplicated(cells$cell_id), identical(cells$cell_id, barcodes),
          all(features$genome == cfg$genome),
          all(is.finite(features$start)), all(is.finite(features$end)),
          all(features$start >= 0), all(features$end > features$start),
          all(features$start == floor(features$start)), all(features$end == floor(features$end)))
con <- if (grepl("[.]gz$", p$matrix)) gzfile(p$matrix) else file(p$matrix)
raw <- as(readMM(con), "CsparseMatrix")  # readMM closes the connection it consumes.
stopifnot(nrow(raw) == nrow(features), ncol(raw) == length(barcodes),
          all(is.finite(raw@x)), all(raw@x >= 0), all(raw@x == floor(raw@x)))
rownames(raw) <- features$feature_id; colnames(raw) <- barcodes
chosen <- cells$celltype == p$celltype
meta <- cells[chosen, , drop = FALSE]
groups <- unique(meta[, c("donor", "condition")])
groups$donor <- factor(groups$donor)
groups$condition <- factor(groups$condition, levels = c(p$control, p$treatment))
groups <- groups[order(groups$donor, groups$condition), , drop = FALSE]
stopifnot(!anyNA(groups), nlevels(groups$donor) >= p$min_pairs,
          nrow(groups) == 2L * nlevels(groups$donor), all(table(groups$donor, groups$condition) == 1))
key <- function(x) paste(x$donor, x$condition, sep = "\t")
index <- match(key(meta), key(groups))
stopifnot(!anyNA(index))
membership <- sparseMatrix(i = seq_len(nrow(meta)), j = index, x = 1,
                          dims = c(nrow(meta), nrow(groups)))
pb <- as.matrix(raw[, chosen, drop = FALSE] %*% membership)
groups$group_id <- paste0("pb", seq_len(nrow(groups)))
groups$n_cells <- as.integer(colSums(membership))
stopifnot(all(groups$n_cells >= p$min_cells), all(colSums(pb) > 0))
groups$samples <- vapply(seq_len(nrow(groups)), function(i) paste(sort(unique(meta$sample[index == i])), collapse = ";"), character(1))
colnames(pb) <- groups$group_id
design <- model.matrix(~ donor + condition, groups)
rownames(design) <- groups$group_id
stopifnot(qr(design)$rank == ncol(design), nrow(design) > ncol(design))
y <- DGEList(counts = pb)
tested <- filterByExpr(y, design = design)
if (sum(tested) < 2L) stop("Too few peaks pass the prespecified expression-count filter")
tab(data.frame(features, tested = tested), "feature_filter.tsv")
y <- y[tested, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = ncol(design))
result <- topTags(test, n = Inf, sort.by = "none")$table
annotations <- features[match(rownames(result), features$feature_id), , drop = FALSE]
tab(data.frame(annotations, result, row.names = NULL), "differential_accessibility.tsv")
groups$library_size_tested <- y$samples$lib.size
groups$TMM_norm_factor <- y$samples$norm.factors
tab(data.frame(feature_id = rownames(pb), pb, check.names = FALSE), "counts.tsv")
tab(groups, "groups.tsv")
tab(data.frame(group_id = rownames(design), design, check.names = FALSE), "design.tsv")
pdf(file.path(out, "diagnostics.pdf"), width = 9, height = 6)
plotMDS(y, labels = paste(groups$donor, groups$condition, sep = ":"), col = as.integer(groups$condition))
plotBCV(y)
plotQLDisp(fit)
plotMD(test, status = as.integer(result$FDR < 0.05), main = paste(p$treatment, "vs", p$control))
abline(h = c(-1, 1), lty = 2, col = "grey")
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
write_json(list(synthetic = isTRUE(p$synthetic), celltype = p$celltype,
  independent_donor_pairs = nlevels(groups$donor), pseudobulk_groups = nrow(groups),
  cells = nrow(meta), peaks_tested = nrow(result), peaks_FDR_below_0_05 = sum(result$FDR < 0.05),
  contrast = paste(p$treatment, "minus", p$control), design = "~ donor + condition",
  interpretation = "Positive logFC means higher relative accessibility in treatment. Donors are replicates; cells are subsamples. TMM assumes most features are not undergoing an unbalanced global change. Association is not causation. This fixed model does not adjust for arbitrary additional covariates.",
  artifacts = c("counts.tsv", "groups.tsv", "design.tsv", "feature_filter.tsv", "differential_accessibility.tsv", "diagnostics.pdf", "sessionInfo.txt")),
  file.path(out, "summary.json"), auto_unbox = TRUE, pretty = TRUE)
