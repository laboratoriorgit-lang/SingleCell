# Step 1 - Single-cell analysis (R)
# Auto-derived from the methods paper (Part 2). Run inside the
# psblab/scrnaseq container with the repository mounted at /workspace.
# Helper functions live in workflow/ and are sourced by the init block below.

# ==============================================================================
# Initialization
# ==============================================================================
# INITIALIZATION - Edit only this block before running the pipeline
# =============================================================================

# data_dir: root of the project (the cloned repository is mounted as /workspace).
#   Input paths are resolved relative to it, and all result files are written
#   to data_dir/results/<step>/.
# seed: random seed, fixed so that reruns give comparable results.
# init_pipeline() loads the packages and helper functions, sets the seed,
# moves R into data_dir and creates the results folders (dir_01 ... dir_08,
# dir_objects). It also defines DATA_DIR, PIPELINE_DIR, base_dir and output_dir.
source("/workspace/workflow/init.R")
init_pipeline(data_dir = "/workspace/.", seed = 1807)


# =============================================================================
# ########################  PART 1 - SINGLE-CELL ANALYSIS  ####################
# =============================================================================

# ==============================================================================
# Section 1 - Data loading and pre-filter QC
# ==============================================================================
# SECTION 1 - DATA LOADING AND PRE-FILTER QC
# =============================================================================
# Each sample is loaded from its input file and the mitochondrial read
# fraction is computed per cell to guide filtering thresholds in
# Section 2. Pre-filter violin plots are saved to 01_qc/.
#
# +- CHANGE FOR YOUR ORGANISM --------------------------------------------------
#   Arabidopsis : mt_pattern = "^ATMG"
#   Human       : mt_pattern = "^MT-"
#   Mouse       : mt_pattern = "^mt-"
# +-----------------------------------------------------------------------------
output_dir <- dir_01

# -- Sample manifest (CellRanger filtered_feature_bc_matrix) ------------------
# Add one entry per sample. Each entry needs:
#   file      - path to the filtered_feature_bc_matrix/ directory (relative to DATA_DIR)
#   label     - unique name for this sample (appears in all plots)
#   condition - experimental group this sample belongs to
samples <- list(
  list(
    file      = "ScWT/outs/filtered_feature_bc_matrix",
    label     = "WT",
    condition = "WT"
  ),
  list(
    file      = "Scpifq/outs/filtered_feature_bc_matrix",
    label     = "pifq",
    condition = "pifq"
  )
)


# -- Plot colors (one color per sample label) -----------------------------------
colors <- c(
  "WT"  = "#66c2a5",
  "pifq"  = "#fc8d62"
)


mt_pattern <- "^ATMG"  # Arabidopsis mitochondrial genes


# Load all samples using the helper function
seurat_list_raw <- load_seurat_samples(samples = samples,
                                       DATA_DIR = DATA_DIR,
                                       mt_pattern = mt_pattern)

plot_qc_batch(seurat_list_raw, colors, "qc_prefilter.pdf")


# =============================================================================

message("\nOK SECTION 1 COMPLETE: QC pre-filter plots saved")

# ==============================================================================
# Section 2 - Cell filtering and doublet detection
# ==============================================================================
output_dir <- dir_01

seurat_list <- filter_seurat_samples(
  seurat_list_raw, min_features = 200, max_mt = 5, run_doubletfinder = TRUE
)

plot_qc_batch(seurat_list, colors, "qc_postfilter.pdf")

saveRDS(seurat_list, file.path(dir_objects, "seurat_list_postfilter.rds"))

message("\nOK SECTION 2 COMPLETE: filtering and doublet detection complete")

# ==============================================================================
# Section 3 - Merge and initial preprocessing
# ==============================================================================
output_dir <- dir_01

ath_sc <- merge_and_preprocess(seurat_list, nfeatures = 2000, npcs = 30)

save_pdf(
  plot_umap(ath_sc, "orig.ident", colors),
  "umap_preharmony.pdf"
)

saveRDS(ath_sc, file.path(dir_objects, "ath_sc_preharmony.rds"))

message("\nOK SECTION 3 COMPLETE: merge and preprocessing complete")

# ==============================================================================
# Section 4 - Harmony integration
# ==============================================================================
output_dir <- dir_01

ath_sc <- run_harmony_umap(ath_sc, group_by = "orig.ident", dims = 1:30)

save_pdf(
  plot_umap(ath_sc, "orig.ident", colors),
  "umap_postharmony.pdf"
)

saveRDS(ath_sc, file.path(dir_objects, "ath_sc_postharmony.rds"))

message("\nOK SECTION 4 COMPLETE: Harmony integration complete")

# ==============================================================================
# Section 5 - Elbow plot (optional) and clustering resolution
# ==============================================================================
resolutions_test <- c(0.15, 0.35, 0.45, 0.55, 1.0)

output_dir <- dir_02

plot_resolution_elbow(ath_sc, output_dir = output_dir)

clu <- run_resolution_sweep(
  ath_sc, resolutions = resolutions_test, output_dir = output_dir
)

message("\nOK SECTION 5 COMPLETE: elbow plot and clustree saved")

# ==============================================================================
# Section 6 - Final clustering
# ==============================================================================
cluster_resolution <- 0.35
cluster_dims       <- 1:30   # default; fewer dimensions change the UMAP and clusters (see the elbow plot)
output_dir <- dir_02

ath_sc <- cluster_cells(
  ath_sc, resolution = cluster_resolution, dims = cluster_dims, k_param = 20
)

Idents(ath_sc) <- "seurat_clusters"
save_pdf(
  plot_umap(ath_sc, "seurat_clusters"),
  "umap_seuratclusters.pdf"
)

message("\nOK SECTION 6 COMPLETE: final clustering complete")

# ==============================================================================
# Section 7 - Cell-type annotation
# ==============================================================================
output_dir <- dir_03

biblio_marks_file <- file.path(DATA_DIR, "data", "biblio_marks.txt")
# Marker table with columns: celltype, symbol, ID (see standardize_marker_table()).
# The old two-column format (cell type, gene ID) is also accepted.
marker_table <- standardize_marker_table(biblio_marks_file)

markers <- find_markers(
  ath_sc,
  output_file = file.path(output_dir, "FindAllMarkers.tsv"),
  force = TRUE
)

ath_sc <- annotate_by_markers(
  ath_sc, markers,
  reference_file = biblio_marks_file
)

plot_marker_dotplot(
  ath_sc,
  marker_table,
  annot_col = "celltype",
  name_gene = TRUE,
  outfile   = file.path(
    output_dir, "dotplot_marker_table_annotation_biblio.pdf"
  ),
  width = 18, height = 18
)

save_pdf(
  plot_umap(ath_sc, "celltype"),
  "umap_annotation_biblio.pdf"
)

message("\nOK SECTION 7 COMPLETE: cell-type annotation complete")

# ==============================================================================
# Section 8 - Gene expression visualization
# ==============================================================================
gene <- "AT4G21750"   # one gene, or several: c("gen1", "gen2")

output_dir <- dir_04

ath_sc <- JoinLayers(ath_sc)
Idents(ath_sc) <- "celltype"

# Name each file after the plotted gene(s): AT4G21750  ->  ..._AT4G21750.pdf
# c("gen1","gen2") -> ..._gen1_gen2.pdf
gene_tag <- paste(gene, collapse = "_")

save_vln(
  plot_gene_violin(ath_sc, gene),
  paste0("vln_gene_", gene_tag, ".pdf")
)

save_pdf(
  plot_gene_umap(ath_sc, gene),
  paste0("feature_gene_", gene_tag, ".pdf")
)

message("\nOK SECTION 8 COMPLETE: gene expression visualization saved")

# ==============================================================================
# Section 9 - Export the annotated object
# ==============================================================================
# Checkpoint - restore with:
#   ath_sc <- readRDS(file.path(dir_objects, "ath_sc_annotated.rds"))
saveRDS(ath_sc, file.path(dir_objects, "ath_sc_annotated.rds"))

# Export the annotated object to AnnData h5ad format for Python-based
# trajectory and velocity analyses.
export_to_scanpy(
  ath_sc,
  file.path(dir_objects, "ath_sc_annotated.h5ad")
)

# To export a specific cell type only:
# export_to_scanpy(
#   subset(ath_sc, subset = celltype == "guard cell"),
#   file.path(dir_objects, "GuardCell.h5ad")
# )

message("\nOK SECTION 9 COMPLETE: annotated object saved (.rds + .h5ad)")
