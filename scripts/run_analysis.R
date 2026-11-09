if (!file.exists("config.R")) stop("Run from the project root")
if (dir.exists(".Rlib")) .libPaths(c(normalizePath(".Rlib"), .libPaths()))
needed <- c("Seurat", "Matrix", "png", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Run Rscript scripts/setup.R first. Missing: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(Seurat))
source("config.R")
source("R/boundaries.R")
source("R/plotting.R")
source("R/report.R")
set.seed(config$seed)
options(future.globals.maxSize = 4 * 1024^3)
for (dir in c("tables", "figures", "objects")) {
  dir.create(file.path(config$out_dir, dir), recursive = TRUE, showWarnings = FALSE)
}
table_path <- function(name) file.path(config$out_dir, "tables", paste0(name, ".csv"))
save_table <- function(x, name) write.csv(x, table_path(name), row.names = FALSE)

message("1/7 Loading the original 10x spot-level counts and coordinates")
matrix_dir <- file.path(config$raw_dir, "filtered_feature_bc_matrix")
if (!dir.exists(matrix_dir)) stop("Run Rscript scripts/download_data.R first")
counts <- Seurat::Read10X(matrix_dir)
if (is.list(counts)) counts <- counts[["Gene Expression"]]
coords <- read.csv(file.path(config$raw_dir, "spatial/tissue_positions_list.csv"), header = FALSE,
  col.names = c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
coords <- coords[coords$in_tissue == 1 & coords$barcode %in% colnames(counts), ]
counts <- counts[, coords$barcode, drop = FALSE]
mt <- grepl("^mt-", rownames(counts), ignore.case = TRUE)
if (!any(mt)) stop("No mitochondrial genes found; check gene identifiers")
qc <- data.frame(barcode = colnames(counts), nCount = Matrix::colSums(counts),
  nFeature = Matrix::colSums(counts > 0),
  percent_mt = 100 * Matrix::colSums(counts[mt, , drop = FALSE]) / pmax(1, Matrix::colSums(counts)))
qc$keep <- qc$nFeature >= config$min_features & qc$nCount >= config$min_counts &
  qc$percent_mt <= config$max_percent_mt
save_table(qc, "spot_qc")
coords <- coords[qc$keep, ]; rownames(coords) <- coords$barcode
counts <- counts[, qc$keep, drop = FALSE]
if (ncol(counts) < 100L) stop("Fewer than 100 retained spots; inspect QC")
meta <- qc[qc$keep, c("nCount", "nFeature", "percent_mt")]
rownames(meta) <- coords$barcode
brain <- CreateSeuratObject(counts, assay = "Spatial", meta.data = meta,
  min.cells = config$min_gene_spots, project = "WhereDoesTheBrainEnd")
image <- Seurat::Read10X_Image(file.path(config$raw_dir, "spatial"),
  assay = "Spatial", slice = "anterior1", filter.matrix = TRUE)
brain[["anterior1"]] <- image[colnames(brain)]

message("2/7 Splitting genes before feature selection; normalizing and fitting PCA")
# Exclude mitochondrial/ribosomal genes from BOTH feature sets, but retain their QC metrics.
eligible <- rownames(brain)[!grepl("^(mt-|Rpl|Rps)", rownames(brain), ignore.case = TRUE)]
discovery <- sample(eligible, floor(length(eligible) * config$discovery_fraction))
validation <- setdiff(eligible, discovery)
save_table(data.frame(gene = eligible, set = ifelse(eligible %in% discovery, "discovery", "validation")), "gene_split")
brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
# Variable feature ranking is restricted to the discovery genes.
discovery_object <- subset(brain, features = discovery)
discovery_object <- FindVariableFeatures(discovery_object, selection.method = "vst",
  nfeatures = min(config$variable_features, length(discovery)), verbose = FALSE)
features <- VariableFeatures(discovery_object)
VariableFeatures(brain) <- features
rm(discovery_object)
brain <- ScaleData(brain, features = features, verbose = FALSE)
max_pc <- max(config$pcs)
brain <- RunPCA(brain, features = features, npcs = max_pc, seed.use = config$seed, verbose = FALSE)
embedding <- Embeddings(brain, "pca")
assert_spot_alignment(coords, embedding)
edges <- hex_edges(coords)

message("3/7 Comparing PC counts, clustering resolutions, and random seeds")
grid <- expand.grid(pcs = config$pcs, resolution = config$resolutions,
  seed = config$clustering_seeds, KEEP.OUT.ATTRS = FALSE)
labels <- matrix(NA_character_, nrow = ncol(brain), ncol = nrow(grid),
  dimnames = list(colnames(brain), NULL))
for (p in config$pcs) {
  brain <- FindNeighbors(brain, reduction = "pca", dims = seq_len(p),
    graph.name = c("discovery_nn", "discovery_snn"), verbose = FALSE)
  for (i in which(grid$pcs == p)) {
    brain <- FindClusters(brain, graph.name = "discovery_snn", resolution = grid$resolution[i],
      algorithm = 1, random.seed = grid$seed[i], verbose = FALSE)
    labels[, i] <- as.character(Idents(brain))
    message("  ", i, "/", nrow(grid), ": ", p, " PCs, resolution ", grid$resolution[i], ", seed ", grid$seed[i])
  }
}
colnames(labels) <- paste0("PC", grid$pcs, "_r", grid$resolution, "_s", grid$seed)
baseline_idx <- which(grid$pcs == config$baseline_pcs &
  grid$resolution == config$baseline_resolution & grid$seed == config$clustering_seeds[1])
stopifnot(length(baseline_idx) == 1L)
baseline <- labels[, baseline_idx]
grid$run <- colnames(labels)
grid$clusters <- apply(labels, 2, function(x) length(unique(x)))
grid$ari_to_baseline <- apply(labels, 2, function(x) adjusted_rand(baseline, x))
save_table(grid, "clustering_runs")
save_table(data.frame(barcode = rownames(labels), labels, check.names = FALSE), "cluster_assignments")
votes <- boundary_votes(labels, edges)
baseline_boundary <- votes[, baseline_idx]
grid$boundary_jaccard_to_baseline <- vapply(seq_len(ncol(votes)), function(k) {
  union <- sum(votes[, k] | baseline_boundary)
  if (union == 0) 1 else sum(votes[, k] & baseline_boundary) / union
}, numeric(1))
save_table(grid, "clustering_runs")
edges$stability <- rowMeans(votes)
edges$pc_stability_at_baseline_resolution <- rowMeans(votes[, grid$resolution == config$baseline_resolution, drop = FALSE])
edges$resolution_stability_at_baseline_pcs <- rowMeans(votes[, grid$pcs == config$baseline_pcs, drop = FALSE])
edges$baseline_crossing <- baseline[edges$i] != baseline[edges$j]
edges$stable_boundary <- edges$stability >= config$stable_threshold
edges$stable_interior <- edges$stability <= config$interior_threshold
edges$barcode_i <- coords$barcode[edges$i]; edges$barcode_j <- coords$barcode[edges$j]

message("4/7 Measuring expression changes in genes withheld from clustering")
normalized <- GetAssayData(brain, assay = "Spatial", layer = "data")
# Select validation genes by variance ONLY, without cluster labels or edge scores.
val <- normalized[validation, , drop = FALSE]
variance <- Matrix::rowMeans(val^2) - Matrix::rowMeans(val)^2
validation <- names(sort(variance, decreasing = TRUE))[seq_len(min(config$validation_features, sum(variance > 0)))]
z <- scale(t(as.matrix(normalized[validation, , drop = FALSE])))
z <- z[, apply(z, 2, function(x) all(is.finite(x))), drop = FALSE]
assert_spot_alignment(coords, z)
edges$heldout_rms <- edge_rms(z, edges)
edges$block <- spatial_blocks(coords, edges, config$block_size)
effect <- block_bootstrap_ratio(edges$heldout_rms, edges$stable_boundary,
  edges$stable_interior, edges$block, config$bootstrap_reps, config$seed)
save_table(effect, "boundary_contrast")
# Sensitivity to spatial block size matters as well as clustering parameters.
block_effects <- do.call(rbind, lapply(c(4L, 8L, 12L), function(size) {
  data.frame(block_size = size, block_bootstrap_ratio(edges$heldout_rms,
    edges$stable_boundary, edges$stable_interior, spatial_blocks(coords, edges, size),
    config$bootstrap_reps, config$seed))
}))
save_table(block_effects, "block_size_sensitivity")
gene_contrasts <- do.call(rbind, lapply(seq_len(ncol(z)), function(g) {
  delta <- abs(z[edges$i, g] - z[edges$j, g])
  stable <- mean(delta[edges$stable_boundary]); interior <- mean(delta[edges$stable_interior])
  data.frame(gene = colnames(z)[g], mean_abs_change_boundary = stable,
    mean_abs_change_interior = interior, ratio = stable / interior)
}))
gene_contrasts <- gene_contrasts[order(gene_contrasts$ratio, decreasing = TRUE), ]
save_table(gene_contrasts, "heldout_gene_contrasts")
# Technical discontinuities may coincide with expression discontinuities.
edges$abs_log_depth_change <- abs(log1p(meta$nCount[edges$i]) - log1p(meta$nCount[edges$j]))
edges$abs_mt_change <- abs(meta$percent_mt[edges$i] - meta$percent_mt[edges$j])
matching <- qc_matched_contrast(edges$heldout_rms, edges$stable_boundary,
  edges$stable_interior, data.frame(
    mean_log_depth = (log1p(meta$nCount[edges$i]) + log1p(meta$nCount[edges$j])) / 2,
    abs_log_depth_change = edges$abs_log_depth_change,
    mean_percent_mt = (meta$percent_mt[edges$i] + meta$percent_mt[edges$j]) / 2,
    abs_mt_change = edges$abs_mt_change))
save_table(matching$summary, "qc_matched_contrast")
save_table(matching$balance, "qc_matching_balance")
save_table(edges, "spatial_edges")
coords$baseline_cluster <- baseline
coords$boundary_score <- spot_boundary_score(nrow(coords), edges, edges$stability)
save_table(coords, "spot_boundaries")
qc_cor <- data.frame(metric = c("absolute_log_depth_change", "absolute_mitochondrial_change"),
  spearman_with_heldout_change = c(cor(edges$heldout_rms, edges$abs_log_depth_change, method = "spearman"),
    cor(edges$heldout_rms, edges$abs_mt_change, method = "spearman")))
save_table(qc_cor, "technical_associations")

message("5/7 Comparing step and gradient profiles with spatial block cross-validation")
profiles <- interface_profiles(coords, edges, baseline, z, max_hops = 5L, block_size = config$block_size)
if (!is.null(profiles$summary) && nrow(profiles$summary)) {
  save_table(profiles$summary, "interface_models")
  save_table(profiles$profiles, "interface_profiles")
}
message("6/7 Checking overall spatial organization against a position-shuffled null")
heldout_pc <- prcomp(z, center = FALSE, scale. = FALSE, rank. = 10)$x[, 1:10, drop = FALSE]
null <- spatial_randomization(heldout_pc, edges, config$permutation_reps, config$seed)
save_table(data.frame(observed = null$observed, null_mean = mean(null$null),
  p_lower = null$p_lower, permutations = length(null$null)), "global_spatial_null")
save_table(data.frame(iteration = seq_along(null$null), mean_neighbor_distance = null$null), "spatial_null_draws")
message("7/7 Saving figures, analysis objects, provenance, and the Markdown report")
plot_analysis(coords, edges, labels, baseline, z, profiles, null, qc, config$raw_dir,
  config$out_dir, gene_contrasts)
brain$baseline_cluster <- baseline
brain$boundary_score <- coords$boundary_score
saveRDS(brain, file.path(config$out_dir, "objects/brain.rds"))
saveRDS(list(config = config, coords = coords, edges = edges, runs = grid,
  profiles = profiles, effect = effect, null = null), file.path(config$out_dir, "objects/results.rds"))
writeLines(capture.output(sessionInfo()), file.path(config$out_dir, "session-info.txt"))
used_packages <- sort(unique(c(needed, loadedNamespaces())))
save_table(data.frame(package = used_packages, version = vapply(used_packages,
  function(p) as.character(packageVersion(p)), character(1))), "package_versions")
write_report(config, qc, coords, edges, grid, profiles, effect, null, qc_cor, gene_contrasts, matching)
message("Done. Open outputs/report.md and outputs/figures/.")
