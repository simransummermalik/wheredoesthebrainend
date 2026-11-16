# Install only into this project's ignored library, never the user's global library.
dir.create(".Rlib", showWarnings = FALSE)
.libPaths(c(normalizePath(".Rlib"), .libPaths()))
required <- c("Seurat", "Matrix", "png", "jsonlite", "knitr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
  install.packages(missing, lib = ".Rlib", Ncpus = 2L)
}
failed <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(failed)) stop("Installation incomplete: ", paste(failed, collapse = ", "))
message("R dependencies are ready in .Rlib/.")
