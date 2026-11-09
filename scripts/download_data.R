source("config.R")
dir.create(config$raw_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = 1200)
prefix <- paste0("https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/",
  config$dataset, "/", config$dataset)
files <- c(matrix = "filtered_feature_bc_matrix.tar.gz", spatial = "spatial.tar.gz")
manifest <- data.frame(file = character(), url = character(), md5 = character())
for (suffix in files) {
  url <- paste0(prefix, "_", suffix)
  dest <- file.path(config$raw_dir, suffix)
  if (!file.exists(dest)) {
    message("Downloading ", suffix)
    partial <- paste0(dest, ".part")
    download.file(url, partial, mode = "wb", method = "libcurl")
    if (!file.rename(partial, dest)) stop("Could not move downloaded file: ", dest)
  }
  members <- utils::untar(dest, list = TRUE)
  if (any(grepl("(^/|(^|/)\\.\\.(/|$))", members))) stop("Unsafe archive paths")
  utils::untar(dest, exdir = config$raw_dir)
  manifest <- rbind(manifest, data.frame(file = suffix, url = url,
    md5 = unname(tools::md5sum(dest))))
}
manifest$retrieved_on <- as.character(Sys.Date())
write.csv(manifest, file.path(config$raw_dir, "download_manifest.csv"), row.names = FALSE)
required <- c("filtered_feature_bc_matrix/matrix.mtx.gz",
  "filtered_feature_bc_matrix/features.tsv.gz", "filtered_feature_bc_matrix/barcodes.tsv.gz",
  "spatial/tissue_positions_list.csv", "spatial/tissue_lowres_image.png",
  "spatial/scalefactors_json.json")
if (!all(file.exists(file.path(config$raw_dir, required)))) {
  stop("Archive layout differs from expected Space Ranger 1.1 files; inspect data/raw/anterior1.")
}
message("Data ready. See data/raw/anterior1/download_manifest.csv for source URLs and checksums.")
