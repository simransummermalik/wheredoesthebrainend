# Rebuild the narrative from saved results without recomputing the analysis.
source("config.R")
source("R/report.R")
result_path <- file.path(config$out_dir, "objects/results.rds")
if (!file.exists(result_path)) stop("Run scripts/run_analysis.R first")
result <- readRDS(result_path)
config <- result$config
read_table <- function(name) read.csv(file.path(config$out_dir, "tables", paste0(name, ".csv")))
write_report(config, read_table("spot_qc"), result$coords, result$edges,
  result$runs, result$profiles, result$effect, result$null,
  read_table("technical_associations"), read_table("heldout_gene_contrasts"),
  list(summary = read_table("qc_matched_contrast"), balance = read_table("qc_matching_balance")))
message("Updated ", file.path(config$out_dir, "report.md"), " from saved results.")
