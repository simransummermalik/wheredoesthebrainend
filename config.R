# Run scripts from the project root. All analysis is implemented in R.
config <- list(
  seed = 925L,
  dataset = "V1_Mouse_Brain_Sagittal_Anterior",
  dataset_url = paste0("https://www.10xgenomics.com/datasets/",
    "mouse-brain-serial-section-1-sagittal-anterior-1-standard-1-1-0"),
  raw_dir = "data/raw/anterior1",
  min_features = 200L,
  min_counts = 1000L,
  max_percent_mt = 30,
  min_gene_spots = 10L,
  discovery_fraction = 0.5,
  variable_features = 2000L,
  validation_features = 1000L,
  pcs = c(15L, 25L, 35L),
  resolutions = c(0.3, 0.6, 0.9),
  clustering_seeds = c(925L, 2026L),
  baseline_pcs = 25L,
  baseline_resolution = 0.6,
  stable_threshold = 0.8,
  interior_threshold = 0.2,
  block_size = 8L, # Visium array rows; x is array_col / 2
  bootstrap_reps = 500L,
  permutation_reps = 499L,
  transect_min_spots = 30L,
  out_dir = "outputs"
)
