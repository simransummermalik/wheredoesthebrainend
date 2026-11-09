# Where Does the Brain End?

An R analysis of stable molecular interfaces and gradual expression transitions in the **10x mouse brain sagittal-anterior section 1**. This implements the selected Option 4 research project.

**Question:** Which neighboring tissue locations remain molecularly separated across Seurat settings, and which interfaces are better described by gradual expression trends?

The project combines an expression-only Seurat baseline with a physical hexagonal spot graph. It withholds genes from boundary discovery, audits 18 clustering runs, and compares constant, step, and gradient profiles using spatial block cross-validation. It also estimates within-section uncertainty with a spatial block bootstrap.

Spatial boundary analysis already exists; this project does not claim a new segmentation algorithm or a first discovery. Its contribution is a transparent, reproducible assessment of interface stability and transition shape in this dataset. [BANKSY](https://www.nature.com/articles/s41588-024-01664-3) and [StarTrail](https://doi.org/10.1080/01621459.2026.2654225) provide relevant prior work.

## Run

Open `wheredoesthebrainend.Rproj` in RStudio, or use a terminal from this directory:

```sh
Rscript scripts/setup.R
Rscript scripts/download_data.R
Rscript scripts/check.R
Rscript scripts/run_analysis.R
```

These commands invoke R; all downloading, computation, plotting, and report generation is written in R. There is no Python component. You can instead run each script with `source("scripts/setup.R")`, etc. in an R session, in the same order.

Use a recent R version compatible with Seurat (R 4.4+ recommended). Setup installs missing dependencies into the project's ignored `.Rlib/` directory. Initial setup needs internet access and may take several minutes, particularly when packages compile from source. The data download is approximately 64 MB; analysis and installed libraries require additional local disk space. Rerunning the analysis recomputes its results from the downloaded raw inputs.

Start reading at **[outputs/report.md](outputs/report.md)** after a successful run. Every result in that report comes from the real downloaded dataset. It is generated directly by R and can be read as Markdown without Pandoc or a document-rendering service.

After editing the report narrative in `R/report.R`, run `Rscript scripts/render_report.R` to regenerate the report from saved results without rerunning clustering.

## Outputs

| File or directory | Contents |
| --- | --- |
| `outputs/report.md` | Research methods, calculated findings, figures, and limitations |
| `outputs/figures/` | QC, clusters, edge stability, held-out expression, agreement, interface profiles, null distribution |
| `outputs/tables/clustering_runs.csv` | PC, resolution, seed, cluster count, adjusted Rand index |
| `outputs/tables/spatial_edges.csv` | Every adjacent spot pair and its stability, held-out change, and QC contrasts |
| `outputs/tables/interface_models.csv` | Spatial cross-validation errors for constant, step, and linear descriptions |
| `outputs/tables/boundary_contrast.csv` | Boundary/interior expression-change ratio and conditional block interval |
| `outputs/tables/block_size_sensitivity.csv` | Intervals using different spatial block sizes |
| `outputs/tables/qc_matched_contrast.csv` | Exploratory comparison after coarse depth/mitochondrial matching |
| `outputs/tables/qc_matching_balance.csv` | Covariate balance before and after matching |
| `outputs/tables/heldout_gene_contrasts.csv` | Descriptive gene ranking; no gene-level significance claims |
| `outputs/tables/spot_qc.csv` | QC metrics and inclusion decisions for every tissue spot |
| `outputs/tables/gene_split.csv` | Reproducible assignment to discovery and validation genes |
| `outputs/tables/package_versions.csv` | Versions used in the completed analysis |
| `outputs/objects/brain.rds` | Seurat object with baseline clusters and spot boundary scores; ignored by Git |
| `outputs/objects/results.rds` | Complete results for follow-up in R; ignored by Git |
| `outputs/session-info.txt` | R environment; ignored by Git |

If no interface has sufficient data, the interface tables/figure are omitted and the report says so. Change the analysis settings in `config.R` before rerunning. The tests use synthetic controls to check graph geometry, known sharp/gradual transitions, spatial resampling, and clustering-label invariance. Synthetic controls are not reported as mouse-brain findings.

## Data

[10x Genomics — Mouse Brain Serial Section 1 (Sagittal-Anterior)](https://www.10xgenomics.com/datasets/mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0), Space Ranger 1.1.0, CC BY 4.0. The original dataset has 2,695 tissue-covered spots from one adult male mouse section. The downloader records exact source URLs and MD5 checksums in `data/raw/anterior1/download_manifest.csv`; these locally calculated checksums record the files used, rather than verify against a publisher checksum.

## Interpretation

Spots contain multiple cells. Cluster labels are molecular domains, not verified brain-region names. Separation frequencies summarize the chosen parameter grid and are not biological-boundary probabilities. Graph hops do not establish micrometer transition widths. Withheld genes are correlated with discovery genes and therefore do not supply fully independent validation.

The bootstrap interval is conditional on this section and feature split. The global position-shuffle test assesses overall spatial structure, not individual boundary significance. No spot-wise t-tests, per-edge significance map, or claim of replication across animals is made. Histology overlays support inspection; anatomical agreement requires independent annotations or atlas registration and remains unmeasured here.

Seurat's [spatial workflow](https://satijalab.org/seurat/articles/spatial_vignette.html) documents the underlying count analysis. The report explains each estimator and its limits alongside its calculated results.
