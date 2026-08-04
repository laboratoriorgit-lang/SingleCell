# Step 1 - Single-cell analysis (R)
# Auto-derived from the methods paper (Part 2). Run inside the
# matigara/scrnaseq container with the repository mounted at /workspace.
# Helper functions live in workflow/ and are sourced by the init block below.

# ==============================================================================
# Configuration
# ==============================================================================
# CONFIGURATION - Edit only this block before running the pipeline
# =============================================================================

# Directory containing the pipeline helper scripts.
# The cloned repository itself is mounted as /workspace.
PIPELINE_DIR <- "/workspace/workflow"

# Root directory for your project data and results.
# All result files will be written to DATA_DIR/resultados/<step>/
DATA_DIR   <- "/workspace/."
base_dir   <- file.path(DATA_DIR, "resultados")

# =============================================================================

# ==============================================================================
# Initialization
# ==============================================================================
# INITIALIZATION
# =============================================================================

# -- Load helper scripts --------------------------------------------------------
# Each file is fully documented at https://github.com/cliford2001/SingleCell
source(file.path(PIPELINE_DIR, "load_libraries.R"))          # all R packages
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))# analysis functions

set.seed(1807)
options(Seurat.allow.s4 = FALSE)  # required for Seurat 5 compatibility
setwd(DATA_DIR)

# -- Create per-step output directories ----------------------------------------
list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)  # creates output folders and loads their paths as variables

# output_dir is the global variable used by save_pdf / save_qc / save_vln helpers.
# It is reassigned at the start of each section to the appropriate step directory.
output_dir <- base_dir


# =============================================================================
# ########################  PART 1 - SINGLE-CELL ANALYSIS  ####################
# =============================================================================

# ==============================================================================
# Section 1 - Data loading and pre-filter QC
# ==============================================================================
# SECTION 1 - DATA LOADING AND PRE-FILTER QC
# =============================================================================
# Each sample is loaded from its input file and mitochondrial / chloroplast
# read fractions are computed per cell to guide filtering thresholds in
# Section 2. Pre-filter violin plots are saved to 01_qc/.
#
# +- CHANGE FOR YOUR ORGANISM --------------------------------------------------
#   Arabidopsis : mt_pattern = "^ATMG"  |  cp_pattern = "^ATCG"
#   Human       : mt_pattern = "^MT-"   |  cp_pattern = NULL
#   Mouse       : mt_pattern = "^mt-"   |  cp_pattern = NULL
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
cp_pattern <- "^ATCG"  # Arabidopsis chloroplast genes


# Load all samples using the helper function
seurat_list_raw <- load_seurat_samples(samples = samples,
                                       DATA_DIR = DATA_DIR,
                                       mt_pattern = mt_pattern,
                                       cp_pattern = cp_pattern)

plot_qc_batch(seurat_list_raw, colors, "qc_prefilter.pdf")


# =============================================================================

message("\nOK SECTION 1 COMPLETE: QC pre-filter plots saved")

# ==============================================================================
# Section 2 - Cell filtering and doublet detection
# ==============================================================================
output_dir <- dir_01

seurat_list <- filter_seurat_samples(seurat_list_raw, min_features = 200, max_mt = 5)

plot_qc_batch(seurat_list, colors, "qc_postfilter.pdf")

saveRDS(seurat_list, file.path(dir_objects, "seurat_list_postfilter.rds"))

message("\nOK SECTION 2 COMPLETE: filtering and doublet detection complete")

# ==============================================================================
# Section 3 - Merge and initial preprocessing
# ==============================================================================
output_dir <- dir_01

pbmc_harmony <- reduce(seurat_list, merge) %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(
    selection.method = "vst", nfeatures = 2000, verbose = FALSE
  ) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

save_pdf(
  DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
  "umap_preharmony.pdf"
)

saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_preharmony.rds"))

message("\nOK SECTION 3 COMPLETE: merge and preprocessing complete")

# ==============================================================================
# Section 4 - Harmony batch correction
# ==============================================================================
output_dir <- dir_01

pbmc_harmony <- pbmc_harmony %>%
  RunHarmony("orig.ident", plot_convergence = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE)

save_pdf(
  DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
  "umap_postharmony.pdf"
)

saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_postharmony.rds"))

message("\nOK SECTION 4 COMPLETE: Harmony integration complete")

# ==============================================================================
# Section 5 - Resolution optimization
# ==============================================================================
resolutions_test <- c(0.15, 0.35, 0.45, 0.55, 1.0)

output_dir <- dir_02

k_range <- 1:31
pca_data <- Embeddings(pbmc_harmony, "pca")[, 1:30]
wss <- sapply(
  k_range,
  function(k) kmeans(pca_data, centers = k, nstart = 4)$tot.withinss
)

elbow_plot <- ggplot(data.frame(k = k_range, wss = wss), aes(k, wss)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Number of clusters (k)",
    y = "Within-cluster sum of squares"
  ) +
  theme_minimal()

save_pdf(elbow_plot, "elbow_plot.pdf", w = 18, h = 18)

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, k.param = 20, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(clustree(clu, prefix = "RNA_snn_res."), "clustree2.pdf", w = 18, h = 18)

message("\nOK SECTION 5 COMPLETE: elbow plot and clustree saved")

# ==============================================================================
# Section 6 - Final clustering
# ==============================================================================
cluster_resolution <- 0.35
output_dir <- dir_02

pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, k.param = 20, verbose = FALSE) %>%
  FindClusters(resolution = cluster_resolution, algorithm = 4, verbose = FALSE)

Idents(pbmc_harmony) <- "seurat_clusters"
save_pdf(
  DimPlot(pbmc_harmony, group.by = "seurat_clusters", label = TRUE),
  "umap_seuratclusters.pdf"
)

message("\nOK SECTION 6 COMPLETE: final clustering complete")

# ==============================================================================
# Section 7 - Cell-type annotation
# ==============================================================================
output_dir <- dir_03

biblio_marks_file <- file.path(DATA_DIR, "data", "biblio_marks.txt")
# Only the first two columns are used, by position: column 1 is the cell type
# and column 2 is the gene. The header row is skipped but its names are ignored,
# so the marker file can use any header labels.
marker_table <- read.table(
  biblio_marks_file, header = TRUE, sep = "\t", quote = ""
)
marker_table <- setNames(marker_table[, 1:2], c("cell.types", "gene"))

markers <- find_markers(
  pbmc_harmony,
  output_file = file.path(output_dir, "FindAllMarkers.tsv")
)

pbmc_harmony <- annotate_by_markers(
  pbmc_harmony, markers,
  reference_file = biblio_marks_file
)

plot_marker_dotplot(
  pbmc_harmony,
  marker_table,
  annot_col = "celltype",
  outfile   = file.path(
    output_dir, "dotplot_marker_table_annotation_biblio.pdf"
  ),
  width = 18, height = 18
)

save_pdf(
  DimPlot(
    pbmc_harmony, group.by = "celltype",
    label = TRUE, repel = TRUE, raster = FALSE
  ),
  "umap_annotation_biblio.pdf"
)

message("\nOK SECTION 7 COMPLETE: cell-type annotation complete")

# ==============================================================================
# Section 8 - Annotated clustree
# ==============================================================================
output_dir <- dir_03

Mode <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

stopifnot(exists("clu"))
stopifnot("celltype" %in% colnames(pbmc_harmony@meta.data))

celltype_label <- as.character(pbmc_harmony$celltype)
names(celltype_label) <- Cells(pbmc_harmony)
clu$celltype_label <- celltype_label[Cells(clu)]

print(table(clu$celltype_label, useNA = "ifany"))

save_pdf(
  clustree(
    clu, prefix = "RNA_snn_res.",
    node_label = "celltype_label", node_label_aggr = "Mode"
  ),
  "clustree_annotated.pdf", w = 14, h = 14
)

message("\nOK SECTION 8 COMPLETE: annotated clustree saved")

# ==============================================================================
# Section 9 - Gene expression visualization
# ==============================================================================
gene <- "AT5G26000"   # one gene, or several: c("gen1", "gen2")

output_dir <- dir_04

pbmc_harmony <- JoinLayers(pbmc_harmony)
Idents(pbmc_harmony) <- "celltype"

# Name each file after the plotted gene(s): AT5G26000  ->  ..._AT5G26000.pdf
# c("gen1","gen2") -> ..._gen1_gen2.pdf
gene_tag <- paste(gene, collapse = "_")

save_vln(
  VlnPlot(pbmc_harmony, features = gene),
  paste0("vln_gene_", gene_tag, ".pdf")
)

save_pdf(
  FeaturePlot(pbmc_harmony, features = gene),
  paste0("feature_gene_", gene_tag, ".pdf")
)

message("\nOK SECTION 9 COMPLETE: gene expression visualization saved")

# ==============================================================================
# Section 10 - Cell-type grouping [optional]
# ==============================================================================
grouping <- c()
# Example:
# grouping <- c(
#   "Epidermis Hypocotyl.1" = "Epidermis Hypocotyl",
#   "Epidermis Hypocotyl.2" = "Epidermis Hypocotyl"
# )

output_dir <- dir_05

if (length(grouping) > 0) {
  pbmc_harmony$celltype_grouped <- recode(pbmc_harmony$celltype, !!!grouping)
} else {
  pbmc_harmony$celltype_grouped <- pbmc_harmony$celltype
}

save_pdf(
  DimPlot(
    pbmc_harmony, group.by = "celltype_grouped",
    label = TRUE, repel = TRUE, raster = FALSE
  ),
  "umap_grouped.pdf"
)

message("\nOK SECTION 10 COMPLETE: cell-type grouping complete")

# ==============================================================================
# Section 11 - Subcluster marker inspection [optional]
# ==============================================================================
output_dir <- dir_05
Idents(pbmc_harmony) <- "celltype"

# 1. Inspect a cluster (saves the combined subcluster + marker figure)
inspect_subcluster_markers(
  pbmc_harmony, cluster_id = "2",
  marker_table = marker_table, output_dir = output_dir
)

# 2. Map each cluster's subclusters to labels ("others" = any ID not listed),
#    then curate them in one call. Example for two clusters:
reassign <- list(
  "2" = c("0" = "Identity A", "others" = "Unresolved"),
  "4" = c("1" = "Identity B", "others" = "Unresolved")
)
pbmc_harmony <- curate_clusters(
  pbmc_harmony, reassign,
  marker_table = marker_table, output_dir = output_dir
)

# 3. Save the curated annotation UMAP
save_pdf(
  DimPlot(pbmc_harmony, group.by = "celltype_curated",
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_curated.pdf"
)

message("\nOK SECTION 11 COMPLETE: curation complete")

# ==============================================================================
# Section 12 - Export the curated object
# ==============================================================================
# Checkpoint - restore with:
#   pbmc_harmony <- readRDS(file.path(dir_objects, "pbmc_harmony_curated.rds"))
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_curated.rds"))

# Export the curated object to AnnData h5ad format for Python-based
# trajectory and velocity analyses.
export_to_scanpy(
  pbmc_harmony,
  file.path(dir_objects, "pbmc_harmony_curated.h5ad")
)

# To export a specific cell type only:
# export_to_scanpy(
#   subset(pbmc_harmony, subset = celltype_curated == "guard cell"),
#   file.path(dir_objects, "GuardCell.h5ad")
# )

message("\nOK SECTION 12 COMPLETE: curated object saved (.rds + .h5ad)")
