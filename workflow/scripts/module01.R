# Same synthetic data and summaries as the course; parameterized for Snakemake.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
counts <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
cells <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(!anyDuplicated(cells$cell_id), !anyDuplicated(rownames(counts)))
stopifnot(identical(colnames(counts), cells$cell_id))
stopifnot(all(is.finite(counts)), all(counts >= 0), all(counts == floor(counts)))
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells$panel_counts <- colSums(counts)
cells$detected_regions <- colSums(counts > 0)
regions <- data.frame(region = rownames(counts), total_count = rowSums(counts),
                      detection_fraction = rowMeans(counts > 0), row.names = NULL)
write.csv(cells, file.path(out, "cell_summary.csv"), row.names = FALSE)
write.csv(regions, file.path(out, "region_summary.csv"), row.names = FALSE)
pdf(file.path(out, "panel_counts.pdf"), width = 8, height = 5)
barplot(cells$panel_counts, names.arg = cells$cell_id,
        col = ifelse(cells$group == "A", "steelblue", "tan3"),
        xlab = "Fictional cell", ylab = "Counts across five teaching regions",
        main = "Synthetic accessibility panel")
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
