# Fictional counts with a planted treatment shift, used only for an integration test.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
suppressPackageStartupMessages(library(Matrix))
set.seed(2026)
out <- args[1]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells <- expand.grid(cell = seq_len(6), condition = c("control", "treated"), donor = paste0("d", 1:4), stringsAsFactors = FALSE)
cells$sample <- paste(cells$donor, cells$condition, sep = "_")
cells$cell_id <- paste(cells$sample, cells$cell, sep = "_")
cells$celltype <- "synthetic_lineage"
features <- data.frame(feature_id = paste0("peak", 1:1000), chrom = "chr1", start = seq(0, by = 500, length.out = 1000), end = seq(200, by = 500, length.out = 1000), genome = "hg19")
mu <- matrix(rep(rlnorm(1000, log(15), 0.5), nrow(cells)), nrow = 1000)
donor_effect <- matrix(rlnorm(4000, 0, 0.2), nrow = 1000)
for (j in seq_len(nrow(cells))) mu[, j] <- mu[, j] * donor_effect[, match(cells$donor[j], paste0("d", 1:4))]
mu[1:20, cells$condition == "treated"] <- mu[1:20, cells$condition == "treated"] * 4
counts <- Matrix(matrix(rnbinom(length(mu), mu = as.vector(mu), size = 5), nrow = 1000), sparse = TRUE)
writeMM(counts, file.path(out, "counts.mtx"))
writeLines(cells$cell_id, file.path(out, "barcodes.tsv"))
write.table(features, file.path(out, "features.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cells[, c("cell_id", "sample", "donor", "condition", "celltype")], file.path(out, "cells.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
samples <- unique(cells[, c("sample", "donor", "condition")]); samples$genome <- "hg19"
write.table(samples, file.path(out, "samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
