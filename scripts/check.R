# Meaningful checks for graph geometry, label invariance, uncertainty, and known transitions.
# Requires only base R; does not download data or install packages.
source("R/boundaries.R")
source("R/report.R")
grid <- expand.grid(array_row = 0:31, column = 0:31)
coords <- transform(grid, array_col = 2 * column + array_row %% 2,
  barcode = paste0("spot", seq_len(nrow(grid))))
edges <- hex_edges(coords)
stopifnot(all(edges$i < edges$j), !anyDuplicated(edges),
  max(tabulate(c(edges$i, edges$j), nrow(coords))) == 6L)
delta_row <- abs(coords$array_row[edges$i] - coords$array_row[edges$j])
delta_col <- abs(coords$array_col[edges$i] - coords$array_col[edges$j])
stopifnot(all((delta_row == 0 & delta_col == 2) | (delta_row == 1 & delta_col == 1)))
hole <- coords[-which(coords$array_row == 15), ]
hole_edges <- hex_edges(hole)
stopifnot(!any(abs(hole$array_row[hole_edges$i] - hole$array_row[hole_edges$j]) > 1))
baseline <- ifelse(coords$column < 16, "left", "right")
switched <- ifelse(baseline == "left", "B", "A")
votes <- boundary_votes(cbind(baseline, switched), edges)
stopifnot(identical(votes[, 1], votes[, 2]), abs(adjusted_rand(baseline, switched) - 1) < 1e-12)
scores <- spot_boundary_score(nrow(coords), edges, rep(1, nrow(edges)))
stopifnot(all(scores == 1))
set.seed(10)
step_signal <- as.numeric(baseline == "right")
gradient_signal <- coords$column / 31
noise <- matrix(rnorm(nrow(coords) * 40, sd = 0.12), nrow(coords), 40)
make_expression <- function(signal) {
  z <- scale(sweep(noise, 1, signal, "+"))
  rownames(z) <- coords$barcode
  z
}
sharp <- make_expression(step_signal)
gradual <- make_expression(gradient_signal)
assert_spot_alignment(coords, sharp)
bad <- sharp[rev(seq_len(nrow(sharp))), ]
stopifnot(inherits(try(assert_spot_alignment(coords, bad), silent = TRUE), "try-error"))
sharp_result <- interface_profiles(coords, edges, baseline, sharp)
gradual_result <- interface_profiles(coords, edges, baseline, gradual)
flat_result <- interface_profiles(coords, edges, baseline, make_expression(rep(0, nrow(coords))))
stopifnot(nrow(sharp_result$summary) == 1L, nrow(gradual_result$summary) == 1L,
  sharp_result$summary$model_preference == "step-like",
  gradual_result$summary$model_preference == "gradient-like",
  flat_result$summary$model_preference == "unresolved")
boundary <- baseline[edges$i] != baseline[edges$j]
blocks <- spatial_blocks(coords, edges, 4)
effect <- block_bootstrap_ratio(ifelse(boundary, 2, 1), boundary, !boundary, blocks, reps = 100)
stopifnot(abs(effect$ratio - 2) < 1e-12, abs(effect$lower - 2) < 1e-12,
  abs(effect$upper - 2) < 1e-12)
empty <- block_bootstrap_ratio(rep(1, nrow(edges)), rep(FALSE, nrow(edges)),
  rep(TRUE, nrow(edges)), blocks)
stopifnot(is.na(empty$ratio))
technical <- rep(c(1, 4), each = 40)
selected <- c(rep(TRUE, 10), rep(FALSE, 30), rep(TRUE, 30), rep(FALSE, 10))
matched <- qc_matched_contrast(technical, selected, !selected, data.frame(depth = technical))
stopifnot(mean(technical[selected]) / mean(technical[!selected]) > 1.5,
  abs(matched$summary$ratio - 1) < 1e-12, matched$summary$boundary_edges == 40L)
null <- spatial_randomization(cbind(gradient_signal, coords$array_row / 31), edges, reps = 39)
stopifnot(null$observed < min(null$null), null$p_lower == 1 / 40)
stopifnot(length(markdown_table(data.frame(a = 1.234, b = "x|y"))) == 3L)
for (file in c(list.files("R", pattern = "\\.R$", full.names = TRUE),
    list.files("scripts", pattern = "\\.R$", full.names = TRUE), "config.R")) parse(file)
cat("Passed: six-neighbor geometry; missing-spot gaps; cluster-label invariance;\n",
  "spot alignment; known sharp/gradual recovery; spatial block bootstrap;\n",
  "empty-boundary handling; QC matching; spatial randomization; report tables; R syntax.\n", sep = "")
