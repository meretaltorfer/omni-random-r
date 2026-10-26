#!/usr/bin/env Rscript

# Random clustering baseline. The two cluster columns contain the same
# random assignment solely to match the existing metrics module's interface.
library(argparse)
library(jsonlite)
library(zellkonverter)

parser <- ArgumentParser(description = "Reproducible random clustering baseline")
parser$add_argument("--output_dir", "-o", required = TRUE)
parser$add_argument("--name", "-n", required = TRUE)
parser$add_argument("--data.h5ad", dest = "data_path", required = TRUE)
parser$add_argument("--data.clusters_truth_num", dest = "k_path", required = TRUE)
parser$add_argument("--seed", type = "integer", default = 42L)
parser$add_argument("--n_comp", type = "integer", default = 30L)
args <- parser$parse_args()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- function(suffix) {
  file.path(args$output_dir, paste0(args$name, suffix))
}

# Use the same H5AD reader as the existing, working R pipeline.
sce <- readH5AD(args$data_path, reader = "python")
cell_ids <- colnames(sce)
stopifnot(length(cell_ids) == ncol(sce), !anyNA(cell_ids),
          !anyDuplicated(cell_ids))

# Only the number of reference groups is used, never their cell memberships.
k <- as.integer(readLines(args$k_path, n = 1L, warn = FALSE))
stopifnot(length(k) == 1L, !is.na(k), k >= 2L, args$n_comp >= 1L)

# Real expression PCA is needed by the existing metrics module for its
# silhouette, Davies-Bouldin and Calinski-Harabasz scores. It does not
# determine the random clusters. Use counts, falling back to the main assay.
assay_names <- SummarizedExperiment::assayNames(sce)
stopifnot(length(assay_names) > 0L)
count_assay <- if ("counts" %in% assay_names) "counts" else assay_names[[1L]]
counts <- SummarizedExperiment::assay(sce, count_assay)
library_sizes <- colSums(counts)
keep <- is.finite(library_sizes) & library_sizes > 0
n_zero_library <- sum(!keep)
counts <- counts[, keep, drop = FALSE]
cell_ids <- cell_ids[keep]
n <- length(cell_ids)
stopifnot(n > k, nrow(counts) >= 2L)

# Build a fresh object so no reference annotations enter preprocessing.
sce_pca <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts)
)
start_pca <- proc.time()[["elapsed"]]
sce_pca <- scuttle::logNormCounts(sce_pca)
n_comp <- min(args$n_comp, n - 1L, nrow(counts) - 1L, 1499L)
set.seed(args$seed)
sce_pca <- scater::runPCA(
  sce_pca, ncomponents = n_comp, ntop = min(1500L, nrow(counts)),
  BSPARAM = BiocSingular::ExactParam()
)
pca <- SingleCellExperiment::reducedDim(sce_pca, "PCA")
pca_seconds <- proc.time()[["elapsed"]] - start_pca
stopifnot(nrow(pca) == n, ncol(pca) >= 1L, all(is.finite(pca)))

# Each cell is independently assigned to one of k groups with equal probability.
# Reset the seed here so PCA does not change the random assignments.
set.seed(args$seed)
start_random <- proc.time()[["elapsed"]]
random_labels <- sample.int(k, size = n, replace = TRUE)
random_seconds <- proc.time()[["elapsed"]] - start_random
if (length(unique(random_labels)) < 2L) {
  stop("Only one group was sampled; use another seed for structure metrics.")
}

clusters <- data.frame(
  cell_id = cell_ids,
  louvain = random_labels,
  leiden = random_labels
)
write.table(clusters, output_path(".clusters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pca <- as.data.frame(pca)
names(pca) <- paste0("PC", seq_len(ncol(pca)))
write.table(data.frame(cell_id = cell_ids, pca), output_path(".pca.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write_json(
  list(pca = pca_seconds, random_assignment = random_seconds),
  output_path(".timings.json"), auto_unbox = TRUE, pretty = TRUE
)
write_json(
  list(
    method = "random_uniform", seed = args$seed,
    requested_clusters = k, observed_clusters = length(unique(random_labels)),
    n_cells = n, excluded_zero_or_invalid_library_cells = n_zero_library,
    pca_components = ncol(pca),
    note = paste("The leiden and louvain columns contain identical random",
                 "labels for interface compatibility; neither algorithm was run.")
  ),
  output_path(".clustering_info.json"), auto_unbox = TRUE, pretty = TRUE
)
cat(sprintf("Wrote random baseline outputs for %d cells, seed %d, %d groups.\n",
            n, args$seed, length(unique(random_labels))))
