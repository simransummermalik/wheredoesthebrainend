# Functions use spot order explicitly: rows of coordinates / embeddings = graph vertex order.
assert_spot_alignment <- function(coords, embedding) {
  stopifnot(!anyDuplicated(coords$barcode),
    identical(as.character(coords$barcode), rownames(embedding)),
    all(is.finite(as.matrix(embedding))))
}

hex_edges <- function(coords) {
  stopifnot(all(c("barcode", "array_row", "array_col") %in% names(coords)),
    !anyDuplicated(coords$barcode))
  key <- paste(coords$array_row, coords$array_col, sep = ":")
  if (anyDuplicated(key)) stop("Duplicate Visium grid positions")
  # Exactly the three forward edges in the hexagonal Visium lattice.
  offsets <- rbind(c(0, 2), c(1, -1), c(1, 1))
  edges <- lapply(seq_len(nrow(offsets)), function(k) {
    neighbor <- paste(coords$array_row + offsets[k, 1],
      coords$array_col + offsets[k, 2], sep = ":")
    j <- match(neighbor, key)
    i <- which(!is.na(j))
    data.frame(i = pmin(i, j[i]), j = pmax(i, j[i]))
  })
  result <- unique(do.call(rbind, edges))
  if (!nrow(result)) stop("No adjacent retained spots; check coordinates and QC")
  result <- result[order(result$i, result$j), , drop = FALSE]
  rownames(result) <- NULL
  result
}

boundary_votes <- function(labels, edges) {
  labels <- as.matrix(labels)
  if (anyNA(labels)) stop("Missing clustering labels")
  votes <- vapply(seq_len(ncol(labels)), function(k) {
    labels[edges$i, k] != labels[edges$j, k]
  }, logical(nrow(edges)))
  matrix(votes, nrow = nrow(edges), ncol = ncol(labels),
    dimnames = list(NULL, colnames(labels)))
}

edge_rms <- function(z, edges) {
  # z: spots x held-out genes, standardized separately for each gene.
  # Chunking keeps the validation matrix bounded while vectorizing permutations.
  result <- numeric(nrow(edges))
  for (start in seq.int(1L, nrow(edges), by = 1000L)) {
    idx <- seq.int(start, min(nrow(edges), start + 999L))
    difference <- z[edges$i[idx], , drop = FALSE] - z[edges$j[idx], , drop = FALSE]
    result[idx] <- sqrt(rowMeans(difference^2))
  }
  result
}

spot_boundary_score <- function(n, edges, scores) {
  sums <- numeric(n)
  totals <- rowsum(c(scores, scores), c(edges$i, edges$j), reorder = FALSE)
  sums[as.integer(rownames(totals))] <- totals[, 1]
  degree <- tabulate(c(edges$i, edges$j), nbins = n)
  ifelse(degree > 0, sums / pmax(degree, 1L), NA_real_)
}

spatial_blocks <- function(coords, edges, size = 8) {
  # Midpoints assign each physical edge to exactly one block.
  x <- (coords$array_col[edges$i] + coords$array_col[edges$j]) / 4
  y <- (coords$array_row[edges$i] + coords$array_row[edges$j]) / 2
  paste(floor(x / size), floor(y / size), sep = ":")
}

block_bootstrap_ratio <- function(value, boundary, interior, blocks, reps = 500, seed = 925) {
  keep <- is.finite(value) & (boundary | interior)
  d <- data.frame(value = value[keep], boundary = boundary[keep],
    interior = interior[keep], block = blocks[keep])
  if (!any(d$boundary) || !any(d$interior)) {
    return(data.frame(ratio = NA_real_, lower = NA_real_, upper = NA_real_,
      blocks = length(unique(d$block)), successful_reps = 0L))
  }
  ratio <- mean(d$value[d$boundary]) / mean(d$value[d$interior])
  sums <- rowsum(cbind(bsum = d$value * d$boundary, bn = as.integer(d$boundary),
    isum = d$value * d$interior, inn = as.integer(d$interior)), d$block, reorder = FALSE)
  set.seed(seed)
  draws <- replicate(reps, {
    total <- colSums(sums[sample.int(nrow(sums), nrow(sums), replace = TRUE), , drop = FALSE])
    if (total["bn"] == 0 || total["inn"] == 0 || total["isum"] <= 0) NA_real_ else
      (total["bsum"] / total["bn"]) / (total["isum"] / total["inn"])
  })
  ok <- is.finite(draws)
  interval <- if (sum(ok) >= 20 && nrow(sums) >= 5) {
    unname(quantile(draws[ok], c(0.025, 0.975)))
  } else c(NA_real_, NA_real_)
  data.frame(ratio = ratio, lower = interval[1], upper = interval[2],
    blocks = nrow(sums), successful_reps = sum(ok))
}

spatial_randomization <- function(z, edges, reps = 499, seed = 925) {
  # Global null ONLY: profiles are exchangeable over tissue positions.
  # This destroys spatial autocorrelation; it cannot validate individual borders.
  observed <- mean(edge_rms(z, edges))
  set.seed(seed)
  null <- replicate(reps, {
    permutation <- sample.int(nrow(z))
    edge_rms(z, transform(edges, i = permutation[i], j = permutation[j])) |> mean()
  })
  list(observed = observed, null = null,
    p_lower = (1 + sum(null <= observed)) / (reps + 1))
}

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n2 <- choose2(sum(tab))
  if (n2 == 0) return(NA_real_)
  expected <- sum(choose2(rowSums(tab))) * sum(choose2(colSums(tab))) / n2
  upper <- (sum(choose2(rowSums(tab))) + sum(choose2(colSums(tab)))) / 2
  if (abs(upper - expected) < .Machine$double.eps) return(1)
  (sum(choose2(tab)) - expected) / (upper - expected)
}

qc_matched_contrast <- function(value, boundary, interior, covariates, bins = 3L, min_per_group = 3L) {
  stopifnot(nrow(covariates) == length(value), all(is.finite(as.matrix(covariates))))
  candidate <- boundary | interior
  strata <- lapply(covariates, function(x) {
    unique_values <- sort(unique(x[candidate]))
    if (length(unique_values) <= bins) return(match(x, unique_values))
    breaks <- unique(quantile(x[candidate], seq(0, 1, length.out = bins + 1L)))
    cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  })
  stratum <- do.call(paste, c(strata, sep = ":"))
  counts <- aggregate(cbind(boundary = as.integer(boundary), interior = as.integer(interior)),
    list(stratum = stratum), sum)
  support <- counts$stratum[counts$boundary >= min_per_group & counts$interior >= min_per_group]
  keep_b <- boundary & stratum %in% support
  keep_i <- interior & stratum %in% support
  if (!any(keep_b) || !any(keep_i)) return(list(
    summary = data.frame(ratio = NA_real_, boundary_edges = 0L, interior_edges = 0L,
      supported_strata = 0L, excluded_boundary_edges = sum(boundary)), balance = data.frame()))
  weights <- counts$boundary / pmax(counts$interior, 1)
  names(weights) <- counts$stratum
  w <- weights[stratum[keep_i]]
  bmean <- mean(value[keep_b])
  imean <- weighted.mean(value[keep_i], w)
  balance <- do.call(rbind, lapply(names(covariates), function(name) {
    x <- covariates[[name]]
    data.frame(covariate = name, original_boundary_mean = mean(x[boundary]),
      original_interior_mean = mean(x[interior]),
      common_support_boundary_mean = mean(x[keep_b]),
      weighted_interior_mean = weighted.mean(x[keep_i], w))
  }))
  list(summary = data.frame(ratio = bmean / imean, boundary_mean = bmean,
    matched_interior_mean = imean, boundary_edges = sum(keep_b), interior_edges = sum(keep_i),
    supported_strata = length(support), excluded_boundary_edges = sum(boundary & !keep_b)),
    balance = balance)
}

interface_model_cv <- function(coords, baseline, z, hops, pair, max_hops, block_size = 8) {
  keep <- which(baseline %in% pair & hops <= max_hops)
  blocks <- paste(floor(coords$array_col[keep] / (2 * block_size)),
    floor(coords$array_row[keep] / block_size), sep = ":")
  unique_blocks <- sort(unique(blocks))
  if (length(unique_blocks) < 5L) return(c(flat = NA, step = NA, linear = NA))
  # Whole spatial blocks stay together; no spot-wise random train/test split.
  set.seed(925)
  fold <- sample(rep(seq_len(5L), length.out = length(unique_blocks)))
  spot_fold <- fold[match(blocks, unique_blocks)]
  error <- matrix(NA_real_, nrow = length(keep), ncol = 3,
    dimnames = list(NULL, c("flat", "step", "linear")))
  for (k in seq_len(5L)) {
    train <- keep[spot_fold != k]; test <- keep[spot_fold == k]
    left <- train[baseline[train] == pair[1] & hops[train] >= 2]
    right <- train[baseline[train] == pair[2] & hops[train] >= 2]
    if (min(length(left), length(right)) < 5L) next
    origin <- colMeans(z[left, , drop = FALSE])
    direction <- colMeans(z[right, , drop = FALSE]) - origin
    denom <- sum(direction^2)
    if (denom <= 1e-10) next
    frame <- function(idx) data.frame(
      score = as.numeric(sweep(z[idx, , drop = FALSE], 2, origin, "-") %*% direction / denom),
      distance = (hops[idx] + 0.5) * ifelse(baseline[idx] == pair[1], -1, 1),
      side = as.numeric(baseline[idx] == pair[2]))
    train_data <- frame(train); test_data <- frame(test)
    fits <- list(flat = lm(score ~ 1, train_data), step = lm(score ~ side, train_data),
      linear = lm(score ~ distance, train_data))
    for (model in names(fits)) {
      error[spot_fold == k, model] <- (test_data$score - predict(fits[[model]], test_data))^2
    }
  }
  if (sum(complete.cases(error)) < 0.8 * nrow(error)) return(c(flat = NA, step = NA, linear = NA))
  colMeans(error, na.rm = TRUE)
}

interface_profiles <- function(coords, edges, baseline, z, max_hops = 5L, block_size = 8) {
  # Conditional descriptive profiles; domains are discovered from other genes.
  a <- as.character(baseline[edges$i]); b <- as.character(baseline[edges$j])
  crossing <- which(a != b)
  if (!length(crossing)) return(list(summary = data.frame(), profiles = data.frame()))
  keys <- paste(pmin(a[crossing], b[crossing]), pmax(a[crossing], b[crossing]), sep = "|")
  groups <- split(crossing, keys)
  groups <- groups[lengths(groups) >= 10L]
  summaries <- list(); profiles <- list()
  adjacency <- split(c(edges$j, edges$i), c(edges$i, edges$j))
  for (key in names(groups)) {
    pair <- strsplit(key, "|", fixed = TRUE)[[1]]
    active <- baseline %in% pair
    interface <- unique(c(edges$i[groups[[key]]], edges$j[groups[[key]]]))
    hops <- rep(Inf, nrow(coords)); hops[interface] <- 0
    for (h in seq_len(max_hops)) {
      frontier <- which(hops == h - 1L)
      neighbors <- unique(unlist(adjacency[as.character(frontier)], use.names = FALSE))
      neighbors <- neighbors[active[neighbors] & is.infinite(hops[neighbors])]
      hops[neighbors] <- h
    }
    left <- which(baseline == pair[1] & hops >= 2 & hops <= max_hops)
    right <- which(baseline == pair[2] & hops >= 2 & hops <= max_hops)
    if (min(length(left), length(right)) < 10) next
    # Orient a held-out expression axis using domain interiors; no p-values.
    center_left <- colMeans(z[left, , drop = FALSE])
    axis <- colMeans(z[right, , drop = FALSE]) - center_left
    denom <- sum(axis^2)
    if (denom <= 1e-10) next
    score <- as.numeric(sweep(z, 2, center_left, "-") %*% axis / denom)
    keep <- which(active & hops <= max_hops)
    distance <- (hops[keep] + 0.5) * ifelse(baseline[keep] == pair[1], -1, 1)
    df <- data.frame(interface = key, barcode = coords$barcode[keep],
      signed_hops = distance, expression_score = score[keep])
    bins <- aggregate(expression_score ~ signed_hops, df, mean)
    # Fraction of total profile variation concentrated at the immediate interface.
    # High values suggest a step at spot resolution, not a cellular discontinuity.
    jump <- abs(mean(score[baseline == pair[2] & hops == 0]) -
      mean(score[baseline == pair[1] & hops == 0]))
    total_change <- sum(abs(diff(bins$expression_score)))
    cv <- interface_model_cv(coords, baseline, z, hops, pair, max_hops, block_size)
    preference <- "unresolved"
    if (all(is.finite(cv)) && min(cv[c("step", "linear")]) < 0.95 * cv["flat"]) {
      if (cv["step"] < 0.95 * cv["linear"]) preference <- "step-like"
      if (cv["linear"] < 0.95 * cv["step"]) preference <- "gradient-like"
    }
    summaries[[key]] <- data.frame(interface = key, crossing_edges = length(groups[[key]]),
      spots = nrow(df), interface_jump = jump,
      mean_stability = if ("stability" %in% names(edges)) mean(edges$stability[groups[[key]]]) else NA_real_,
      mean_heldout_change = if ("heldout_rms" %in% names(edges)) mean(edges$heldout_rms[groups[[key]]]) else NA_real_,
      jump_fraction = if (total_change > 0) jump / total_change else NA_real_,
      cv_flat_mse = unname(cv["flat"]), cv_step_mse = unname(cv["step"]),
      cv_linear_mse = unname(cv["linear"]), model_preference = preference)
    profiles[[key]] <- df
  }
  list(summary = do.call(rbind, summaries), profiles = do.call(rbind, profiles))
}
