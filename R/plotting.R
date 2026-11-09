save_png <- function(path, draw, width = 1800, height = 1400) {
  grDevices::png(path, width = width, height = height, res = 170)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mar = c(3.5, 3.5, 3.5, 1), bg = "#faf9f6", fg = "#243137",
    col.axis = "#243137", col.lab = "#243137", col.main = "#243137")
  draw()
}

draw_tissue <- function(coords, raw_dir, title = "") {
  image <- png::readPNG(file.path(raw_dir, "spatial/tissue_lowres_image.png"))
  scales <- jsonlite::fromJSON(file.path(raw_dir, "spatial/scalefactors_json.json"))
  w <- dim(image)[2] / scales$tissue_lowres_scalef
  h <- dim(image)[1] / scales$tissue_lowres_scalef
  plot(NA, xlim = c(0, w), ylim = c(h, 0), asp = 1, xlab = "", ylab = "",
    axes = FALSE, main = title, xaxs = "i", yaxs = "i")
  rasterImage(image, 0, h, w, 0)
}

plot_analysis <- function(coords, edges, labels, baseline, z, profiles, null, qc,
                          raw_dir, out_dir, gene_contrasts) {
  fig <- file.path(out_dir, "figures")
  dir.create(fig, recursive = TRUE, showWarnings = FALSE)
  palette <- grDevices::hcl.colors(100, "Inferno")
  map_color <- function(x) palette[pmax(1, pmin(100, 1 + floor(99 * x)))]
  save_png(file.path(fig, "01_qc.png"), function() {
    par(mfrow = c(1, 3))
    for (v in c("nCount", "nFeature", "percent_mt")) {
      hist(qc[[v]], breaks = 40, col = "#3d8179", border = NA,
        main = v, xlab = "All tissue spots")
      rug(qc[[v]][!qc$keep], col = "#d35744")
    }
    mtext("Red ticks: spots removed by the prespecified QC rule", outer = TRUE, line = -1.5)
  }, width = 2200, height = 850)
  save_png(file.path(fig, "02_clusters.png"), function() {
    draw_tissue(coords, raw_dir, "Expression clusters (discovery genes only)")
    groups <- sort(unique(baseline))
    colors <- setNames(grDevices::hcl.colors(length(groups), "Dynamic"), groups)
    points(coords$pxl_col, coords$pxl_row, pch = 16, cex = 0.45,
      col = adjustcolor(colors[as.character(baseline)], 0.85))
    for (g in groups) {
      idx <- which(baseline == g)
      text(median(coords$pxl_col[idx]), median(coords$pxl_row[idx]), g,
        cex = 0.85, font = 2, col = "black")
    }
  })
  save_png(file.path(fig, "03_boundary_stability.png"), function() {
    draw_tissue(coords, raw_dir, "How often do neighboring spots fall in different clusters?")
    order <- order(edges$stability)
    segments(coords$pxl_col[edges$i[order]], coords$pxl_row[edges$i[order]],
      coords$pxl_col[edges$j[order]], coords$pxl_row[edges$j[order]],
      col = adjustcolor(map_color(edges$stability[order]), 0.9), lwd = 0.8)
    legend("bottomleft", legend = c("0: always together", "0.5: parameter-sensitive", "1: always separate"),
      col = map_color(c(0, 0.5, 1)), lwd = 3, bg = "white", cex = 0.65)
  })
  save_png(file.path(fig, "04_validation_change.png"), function() {
    draw_tissue(coords, raw_dir, "Expression change between neighbors: withheld genes")
    cut <- quantile(edges$heldout_rms, 0.98)
    value <- pmin(edges$heldout_rms / cut, 1)
    idx <- order(value)
    segments(coords$pxl_col[edges$i[idx]], coords$pxl_row[edges$i[idx]],
      coords$pxl_col[edges$j[idx]], coords$pxl_row[edges$j[idx]],
      col = map_color(value[idx]), lwd = 0.8)
    legend("bottomleft", legend = c("Small change", "Large change (98th percentile+)"),
      col = map_color(c(0, 1)), lwd = 3, bg = "white", cex = 0.65)
  })
  save_png(file.path(fig, "05_clustering_agreement.png"), function() {
    runs <- ncol(labels)
    ari <- outer(seq_len(runs), seq_len(runs), Vectorize(function(i, j) {
      adjusted_rand(labels[, i], labels[, j])
    }))
    par(mar = c(7, 7, 3, 1))
    image(seq_len(runs), seq_len(runs), ari, zlim = c(0, 1),
      col = grDevices::hcl.colors(100, "YlGnBu"), axes = FALSE,
      xlab = "", ylab = "", main = "Clustering agreement: adjusted Rand index")
    axis(1, at = seq_len(runs), labels = colnames(labels), las = 2, cex.axis = 0.45)
    axis(2, at = seq_len(runs), labels = colnames(labels), las = 2, cex.axis = 0.45)
  }, width = 1800, height = 1800)
  if (!is.null(profiles$profiles) && nrow(profiles$profiles)) {
    groups <- split(profiles$profiles, profiles$profiles$interface)
    selected <- names(sort(vapply(groups, nrow, integer(1)), decreasing = TRUE))[seq_len(min(6, length(groups)))]
    save_png(file.path(fig, "06_interface_profiles.png"), function() {
      par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
      for (key in selected) {
        d <- groups[[key]]
        bins <- aggregate(expression_score ~ signed_hops, d, mean)
        result <- profiles$summary[profiles$summary$interface == key, ]
        plot(bins$signed_hops, bins$expression_score, type = "b", pch = 16,
          col = "#28766f", xlab = "Signed distance (graph hops + 0.5)",
          ylab = "Held-out expression axis", main = paste("Clusters", key))
        mtext(result$model_preference, side = 3, line = 0.1, cex = 0.75)
        abline(v = 0, lty = 2, col = "#bd4e3e")
      }
    }, width = 2200, height = 1400)
  }
  save_png(file.path(fig, "07_spatial_null.png"), function() {
    hist(null$null, breaks = 30, col = "#a0bab6", border = "white",
      xlim = range(c(null$null, null$observed)),
      main = "Global organization: position-shuffled reference",
      xlab = "Mean neighboring distance in held-out PC space")
    abline(v = null$observed, lwd = 3, col = "#b34539")
    legend("topleft", "Observed tissue", col = "#b34539", lwd = 3, bty = "n")
  })
  save_png(file.path(fig, "08_heldout_genes.png"), function() {
    par(mfrow = c(2, 3), mar = c(1, 1, 3, 1))
    for (gene in head(gene_contrasts$gene, 6)) {
      draw_tissue(coords, raw_dir, gene)
      value <- z[, gene]
      limits <- quantile(value, c(0.02, 0.98))
      color <- map_color(pmax(0, pmin(1, (value - limits[1]) / max(diff(limits), 1e-10))))
      points(coords$pxl_col, coords$pxl_row, pch = 16, cex = 0.35, col = color)
      mtext("Within-gene expression; 2nd-98th percentile color scale", side = 3, line = 0.1, cex = 0.5)
    }
  }, width = 2200, height = 1500)
}
