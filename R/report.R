markdown_table <- function(x, digits = 3) {
  if (is.null(x) || !nrow(x)) return("No estimable results for this comparison.")
  for (name in names(x)) {
    if (is.numeric(x[[name]])) x[[name]] <- ifelse(is.finite(x[[name]]),
      format(round(x[[name]], digits), trim = TRUE, scientific = FALSE), "NA")
    x[[name]] <- gsub("|", "/", as.character(x[[name]]), fixed = TRUE)
  }
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    apply(x, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |")))
}

write_report <- function(config, qc, coords, edges, grid, profiles, effect, null,
                         qc_cor, gene_contrasts, matching) {
  fmt <- function(x, digits = 3) formatC(x, digits = digits, format = "f")
  n_stable <- sum(edges$stable_boundary)
  model_results <- profiles$summary
  model_counts <- if (!is.null(model_results) && nrow(model_results)) {
    as.data.frame(table(model_results$model_preference), responseName = "interfaces")
  } else data.frame()
  ratio_text <- if (is.finite(effect$ratio)) {
    paste0("Stable boundary edges had **", fmt(effect$ratio),
      " times** the average held-out expression change of stable interior edges.",
      if (all(is.finite(c(effect$lower, effect$upper)))) paste0(
        " The conditional spatial block bootstrap 95% interval was **[",
        fmt(effect$lower), ", ", fmt(effect$upper), "]**.") else
        " There were too few usable blocks to report an interval.")
  } else "The stable-boundary/interior contrast could not be estimated with these settings."
  abstract <- paste0("**Main results.** The analysis retained ", sum(qc$keep), " spots and audited ",
    nrow(grid), " clustering runs. It identified ", n_stable,
    " stable boundary edges. Among the sampled cluster-pair interfaces, ",
    if (!is.null(model_results)) sum(model_results$model_preference == "step-like") else 0,
    " favored a step, ",
    if (!is.null(model_results)) sum(model_results$model_preference == "gradient-like") else 0,
    " favored a gradient, and ",
    if (!is.null(model_results)) sum(model_results$model_preference == "unresolved") else 0,
    " were unresolved by the specified comparison. The boundary/interior held-out change ratio was ",
    fmt(effect$ratio), " before QC matching and ", fmt(matching$summary$ratio),
    " on the QC-matched common-support comparison. These comparisons use different edge weights/subsets; the matched result should be considered alongside the raw contrast.")
  gradient_example <- ""
  if (!is.null(model_results) && any(model_results$model_preference == "gradient-like")) {
    examples <- model_results[model_results$model_preference == "gradient-like", ]
    example <- examples[which.max(examples$mean_stability), ]
    gradient_example <- paste0("**A concrete example:** cluster pair ",
      gsub("|", "/", example$interface, fixed = TRUE), " had mean separation stability ",
      fmt(example$mean_stability), " but favored the linear gradient description under the specified rule. ",
      "A consistently separated pair of molecular domains can therefore still have a gradual expression profile across their interface in this analysis. This is a conditional, spot-level finding, not evidence of a new anatomical structure.")
  }
  gene_context <- character()
  if ("Snap25" %in% head(gene_contrasts$gene, 12)) {
    gene_context <- c(gene_context,
      "*Snap25* is among the highest-ranked held-out genes. Its protein participates in neurotransmitter release and vesicle fusion, according to the [NCBI mouse SNAP-25 protein record](https://pubchem.ncbi.nlm.nih.gov/protein/P60879). Its spatial pattern is therefore a biologically plausible lead for examining differences in neuronal contributions across the tissue, rather than evidence that SNAP-25 creates the interface.")
  }
  if ("Ndrg4" %in% head(gene_contrasts$gene, 12)) {
    gene_context <- c(gene_context,
      "*Ndrg4* also appears among the highest-ranked genes. The [NCBI mouse Ndrg4 record](https://www.ncbi.nlm.nih.gov/gene/234593) annotates nervous-system expression and roles including vesicle docking and endocytic recycling. These annotations provide biological context for follow-up. They do not establish a specific brain-region identity, a change within a single cell type, or a causal mechanism in this section.")
  }
  text <- c(
    "# Where Does the Brain End?",
    "## Stable molecular borders and gradual transitions in a mouse brain section", "",
    "**An exploratory analysis implemented entirely in R.** This report is generated from the downloaded biological data; numerical statements below are calculated by the workflow.", "",
    abstract, "",
    "### Research question", "",
    "Which molecular interfaces remain separated across reasonable clustering settings, and do expression profiles across those interfaces behave more like abrupt steps or gradual gradients at Visium spot resolution?", "",
    "The working hypothesis is that brain organization includes both persistent molecular separations and transitions that a discrete clustering map can make appear sharper than they are. A parameter-sensitive border is a candidate for further inspection; instability alone does not establish a biological transition.", "",
    "### Why this analysis adds something", "",
    "This project produces an interface audit rather than stopping at a cluster map. It combines parameter sensitivity, genes withheld from discovery, and spatially blocked comparisons of step and gradient profiles on a specific adult mouse section. The contribution is a reproducible assessment of the strength and limits of each molecular interface in this dataset.", "",
    "Spatial segmentation and molecular boundary detection already have a substantial literature. BANKSY incorporates neighborhood information for spatial clustering, while StarTrail studies localized expression gradients. Consequently, neither graph-based spatial analysis nor detecting brain boundaries is claimed as a new invention here. See [BANKSY (Nature Genetics, 2024)](https://www.nature.com/articles/s41588-024-01664-3) and [StarTrail (JASA, 2026)](https://doi.org/10.1080/01621459.2026.2654225).", "",
    "### Data and scope", "",
    paste0("Source: [10x Genomics, Mouse Brain Serial Section 1 (Sagittal-Anterior)](", config$dataset_url, ")."), "",
    "The source describes one sagittal-anterior section from an adult male C57BL/6 mouse, analyzed with Space Ranger 1.1.0. The download includes gene counts, spot coordinates, a tissue image, and image scale factors. The source dataset is CC BY 4.0; attribution is to 10x Genomics. Raw-file URLs and locally computed MD5 checksums are in `data/raw/anterior1/download_manifest.csv`.", "",
    "These are multicellular capture spots, not individual cells. Cluster numbers label expression-defined domains and are not automatically anatomical labels. One section from one donor cannot support population-level claims or estimates of cellular boundary widths.", "",
    "### Prespecified analysis choices", "",
    paste0("- Retain tissue spots with at least ", config$min_features, " detected genes, ",
      config$min_counts, " UMIs, and at most ", config$max_percent_mt, "% mitochondrial counts."),
    "- Exclude mitochondrial and ribosomal genes from discovery and validation features, while retaining their counts for QC.",
    "- Split eligible genes at random before feature selection. Choose variable discovery genes using Seurat's VST ranking, then log-normalize, scale, and fit PCA. Normalization uses total counts, including validation genes; this is feature withholding rather than complete statistical independence.",
    paste0("- Compare ", paste(config$pcs, collapse = "/"), " PCs, clustering resolutions ",
      paste(config$resolutions, collapse = "/"), ", and ", length(config$clustering_seeds),
      " random seeds (", nrow(grid), " total Louvain clustering runs)."),
    "- Connect only immediate neighbors on the original Visium hexagonal lattice. Missing or QC-removed spots are not bridged by artificial nearest-neighbor edges.",
    paste0("- Call edges with separation frequency >= ", config$stable_threshold,
      " stable boundaries and <= ", config$interior_threshold,
      " stable interiors. These thresholds are descriptive choices."),
    paste0("- Evaluate molecular differences using ", config$validation_features,
      " high-variance held-out genes. They are selected without cluster labels or boundary scores."), "",
    "The Seurat steps follow its [spatial analysis documentation](https://satijalab.org/seurat/articles/spatial_vignette.html), with log-normalization used here to keep the gene split explicit. The workflow uses R implementations and Louvain clustering; it does not require Python.", "",
    "### Quality control", "",
    paste0("Of **", nrow(qc), "** tissue spots, **", sum(qc$keep), "** were retained and **",
      sum(!qc$keep), "** removed. The retained spatial graph has **", nrow(edges), "** edges."), "",
    "![Quality control distributions](figures/01_qc.png)", "",
    "The thresholds are deliberately permissive because spatial depth differences can reflect tissue biology. `spot_qc.csv` contains every exclusion. This run does not establish robustness to alternative QC or normalization choices.", "",
    "### Result 1: cluster maps change with analysis settings", "",
    paste0("The number of clusters ranged from **", min(grid$clusters), " to ", max(grid$clusters),
      "**. The baseline used ", config$baseline_pcs, " PCs, resolution ",
      config$baseline_resolution, ", and seed ", config$clustering_seeds[1], "."), "",
    markdown_table(grid[, c("pcs", "resolution", "seed", "clusters", "ari_to_baseline", "boundary_jaccard_to_baseline")]), "",
    "![Discovery gene clusters on histology](figures/02_clusters.png)", "",
    "![Agreement between clustering runs](figures/05_clustering_agreement.png)", "",
    "Adjusted Rand index compares partitions without depending on arbitrary cluster numbering. Boundary Jaccard compares the sets of physical edges cut by each run and the baseline. A high value indicates similar assignments or interfaces, not proof of anatomical accuracy. The edge table also records PC stability at the baseline resolution and resolution stability at the baseline PC count, since changing resolution intentionally changes domain granularity.", "",
    "### Result 2: some edges are more persistent than others", "",
    paste0("**", n_stable, " of ", nrow(edges), " edges (", fmt(100 * n_stable / nrow(edges), 1),
      "%)** met the stable-boundary threshold; **", sum(edges$stable_interior),
      "** met the stable-interior threshold."), "",
    "![Edge separation stability on tissue](figures/03_boundary_stability.png)", "",
    "For an edge connecting neighboring spots i and j, stability is the fraction of specified runs assigning i and j to different clusters. The runs share the same counts and are correlated. This fraction measures sensitivity over this parameter grid; it is not a probability that a biological boundary exists.", "",
    "### Result 3: expression differences in withheld genes", "", ratio_text, "",
    "For each edge, the score is the root mean squared difference across gene-standardized, log-normalized held-out expression values. The contrast above divides its mean on stable boundaries by its mean on stable interiors. Intermediate edges are excluded from this contrast.", "",
    "The bootstrap resamples whole spatial blocks of edge midpoints, conditional on the observed clustering and feature split. It does not resample mice, refit the pipeline, or remove all long-range dependence. Block-size sensitivity is in `tables/block_size_sensitivity.csv`. This interval describes within-section sensitivity and is not a biological-replicate confidence interval.", "",
    "![Withheld gene changes on the tissue graph](figures/04_validation_change.png)", "",
    "Genes with the largest boundary/interior ratios (descriptive ranking, no significance claims):", "",
    markdown_table(head(gene_contrasts, 12)), "",
    "![Spatial expression of the highest-ranked withheld genes](figures/08_heldout_genes.png)", "",
    gene_context, "",
    "Genes covary in biological programs, so withholding genes reduces direct reuse but does not make validation independent. Gene rankings require checking expression levels, known biology, and independent data before interpreting mechanisms.", "",
    "### Result 4: step versus gradient descriptions", "",
    "For every baseline cluster pair sharing at least 10 spatial edges, graph distance is measured inward from the interface on both sides, up to five neighbor hops. An axis from the left to right domain interior summarizes held-out expression. The plotted profile uses all relevant spots and is descriptive.", "",
    "For model comparison, that axis is recomputed using training blocks only in each of five folds. A constant model, a two-level step at the discovered interface, and a linear trend across signed graph distance are compared by held-out spot mean squared error. Models need at least five spatial blocks and adequate interiors on both sides. A model is labeled step-like or gradient-like only if its MSE is at least 5% lower than both alternatives; otherwise the result is unresolved. The 5% rule is a heuristic, not a statistical test.", "",
    markdown_table(model_counts), "",
    markdown_table(if (!is.null(model_results) && nrow(model_results))
      model_results[, c("interface", "mean_stability", "jump_fraction", "cv_flat_mse",
        "cv_step_mse", "cv_linear_mse", "model_preference")] else data.frame()), "",
    gradient_example, "",
    if (!is.null(model_results) && nrow(model_results)) "![Expression profiles across the largest interfaces](figures/06_interface_profiles.png)" else "No interfaces had sufficient data for profile analysis.", "",
    "A linear preference supports a gradual description within the sampled window; it does not prove a continuous cellular state. Step-like means abrupt relative to the sampled graph spacing. Each cluster-pair profile can pool spatially disconnected interface segments. Branching geometry, mixed-cell spots, curved interfaces, and other nearby domains can affect both fits, so local inspection remains necessary. Graph hops are not calibrated micrometers, and no sub-spot transition widths are inferred.", "",
    "### Result 5: overall spatial organization", "",
    paste0("The observed mean neighboring distance in the first 10 held-out PCs was **", fmt(null$observed),
      "**, compared with a position-shuffled null mean of **", fmt(mean(null$null)),
      "**. The lower-tail Monte Carlo p-value was **", fmt(null$p_lower, 4),
      "**, using ", length(null$null), " permutations and the +1 correction."), "",
    "![Global spatial organization reference](figures/07_spatial_null.png)", "",
    "This is a test of the global exchangeability null: complete expression profiles can be assigned arbitrarily to tissue positions. It destroys spatial autocorrelation and tests overall organization only. It cannot test whether any individual border is significant, identify its cause, or establish anatomical agreement.", "",
    "### Technical associations", "", markdown_table(qc_cor), "",
    "These are descriptive Spearman correlations across graph edges. Edges share spots and are not independent replicates. A correlation with depth or mitochondrial differences flags a possible technical or biological association; it does not separate the two.", "",
    "**Exploratory QC-matched sensitivity analysis.** This check was added after observing the QC associations, rather than specified before the first data run. Edges are stratified by tertiles of four covariates: average log depth, absolute log-depth difference, average mitochondrial percentage, and absolute mitochondrial-percentage difference. Strata must contain at least three boundary and three interior edges. Interior edges are weighted to match the number of boundary edges within each supported stratum. The comparison is descriptive and has no separately estimated confidence interval.", "",
    markdown_table(matching$summary), "", markdown_table(matching$balance), "",
    "Compare the matched ratio with the original unadjusted ratio, and inspect covariate balance and excluded edges. Coarse matching cannot remove every technical influence or distinguish a biological depth difference from a technical one. It changes the comparison to the common-support subset; it is not a causal adjustment or proof that the molecular contrast is independent of QC.", "",
    "### What can be concluded", "",
    "The analysis quantifies how stable expression-defined interfaces are under specified clustering settings, checks whether different genes also change there, and compares simple step and gradient descriptions with spatial block cross-validation. The numerical results support a dataset-specific characterization; they do not determine a universal answer about where the brain or any named region ends.", "",
    "The H&E overlays allow visual inspection, but no independent anatomical labels are supplied by this pipeline. Therefore anatomical agreement is not measured here. A stronger follow-up would annotate histology while blinded to expression results, register an appropriate mouse-brain atlas, and repeat the analysis in independent sections or animals. Adjacent sections from one animal would assess technical/spatial reproducibility, not independent biological replication.", "",
    "### Reproducibility", "",
    "Run `Rscript scripts/setup.R`, `Rscript scripts/download_data.R`, and `Rscript scripts/run_analysis.R` from the repository root. The core validation checks run with `Rscript scripts/check.R`. All settings are in `config.R`; raw data, installed libraries, and large Seurat objects are excluded from Git. Package versions are in `tables/package_versions.csv`; the full environment is recorded in `session-info.txt`.", "",
    "R source files: `R/boundaries.R`, `R/plotting.R`, `R/report.R`, and `scripts/run_analysis.R`. The analysis saves tables, figures, the Seurat object, and the complete result object locally. See the project README for interpretation and rerunning instructions."
  )
  writeLines(text, file.path(config$out_dir, "report.md"))
}
