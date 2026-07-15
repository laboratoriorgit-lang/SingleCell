# =============================================================================
# capitulo1_single_cell.R — Single-Cell RNA-seq Analysis (Part 1)
# =============================================================================
# Covers QC, integration, clustering, annotation, and export to h5ad.
# Output: pbmc_harmony_curated.rds  (used by capitulo2_pseudobulk_de.R)
#         pbmc_harmony_curated.h5ad (used by capitulo3_pseudotime.ipynb)
# =============================================================================

# =============================================================================
# Single-Cell RNA-seq Analysis Pipeline
# =============================================================================
# Protocol: Methods in Molecular Biology
#
# Description:
#   End-to-end scRNA-seq pipeline covering:
#     Part 1 — QC, integration, clustering, annotation, and export
#     Part 2 — Pseudobulk differential expression and GO enrichment
#
# Helper scripts (fully documented at https://github.com/cliford2001/SingleCell):
#   load_libraries.R           — loads all required R packages
#   ScRNA_Analysis_Functions.R — core functions (QC, clustering, annotation, DE, GO)
#   custom_seurat.R            — custom Seurat plot utilities (cluster bar charts)
#
# Usage:
#   Run sections sequentially. Section 11 (Cell-Type Curation) must be
#   run interactively — do NOT source the entire script with it active.
#
# Organism:
#   Default: Arabidopsis thaliana. Adapting to another organism requires
#   updating only the organelle patterns, GO database, and key type
#   (each marked clearly in the relevant section below).
#
# Reproducibility:
#   Fixed random seed (1807). Intermediate results written to disk.
# =============================================================================

# =============================================================================
# CONFIGURATION — Edit only this block before running the pipeline
# =============================================================================

# Directory containing the pipeline helper scripts.
# Inside the Docker container this is typically /workspace/SingleCell
PIPELINE_DIR <- "/workspace/SingleCell/workflow"

# Root directory for your project data and results.
# All result files will be written to DATA_DIR/resultados/<step>/
DATA_DIR   <- "/workspace/."
base_dir   <- file.path(DATA_DIR, "resultados")

# =============================================================================
# INITIALIZATION
# =============================================================================

# ── Load helper scripts ────────────────────────────────────────────────────────
# Each file is fully documented at https://github.com/cliford2001/SingleCell
source(file.path(PIPELINE_DIR, "load_libraries.R"))          # all R packages
source(file.path(PIPELINE_DIR, "custom_seurat.R"))           # plot_integrated_clusters()
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))# analysis functions

set.seed(1807)
options(Seurat.allow.s4 = FALSE)  # required for Seurat 5 compatibility
setwd(DATA_DIR)

# ── Create per-step output directories ────────────────────────────────────────
list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)  # creates output folders and loads their paths as variables

# output_dir is the global variable used by save_pdf / save_qc / save_vln helpers.
# It is reassigned at the start of each section to the appropriate step directory.
output_dir <- base_dir


# =============================================================================
# ████████████████████████  PART 1 — SINGLE-CELL ANALYSIS  ████████████████████
# =============================================================================


# SECTION 1 — DATA LOADING AND PRE-FILTER QC
# =============================================================================
# Each sample is loaded from its input file and mitochondrial / chloroplast
# read fractions are computed per cell to guide filtering thresholds in
# Section 2. Pre-filter violin plots are saved to 01_qc/.
#
# ┌─ CHANGE FOR YOUR ORGANISM ──────────────────────────────────────────────────
#   Arabidopsis : mt_pattern = "^ATMG"  |  cp_pattern = "^ATCG"
#   Human       : mt_pattern = "^MT-"   |  cp_pattern = NULL
#   Mouse       : mt_pattern = "^mt-"   |  cp_pattern = NULL
# └─────────────────────────────────────────────────────────────────────────────
output_dir <- dir_01

# ── Sample manifest (CellRanger filtered_feature_bc_matrix) ──────────────────
# Add one entry per sample. Each entry needs:
#   file      — path to the filtered_feature_bc_matrix/ directory (relative to DATA_DIR)
#   label     — unique name for this sample (appears in all plots)
#   condition — experimental group this sample belongs to
samples <- list(
  list(file = "cellranger/Sample_0N/outs/filtered_feature_bc_matrix",  label = "Sample_0N",  condition = "0N"),
  list(file = "cellranger/Sample_05N/outs/filtered_feature_bc_matrix", label = "Sample_05N", condition = "0.5N"),
  list(file = "cellranger/Sample_5N/outs/filtered_feature_bc_matrix",  label = "Sample_5N",  condition = "5N")
)


# ── Plot colors (one color per sample label) ───────────────────────────────────
colors <- c(
  "Sample_0N"  = "#66c2a5",
  "Sample_05N" = "#41ae76",
  "Sample_5N"  = "#fc8d62"
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

message("\n✓ SECTION 1 COMPLETE: QC pre-filter plots saved")
# SECTION 2 — CELL FILTERING AND DOUBLET DETECTION
# =============================================================================
# Thresholds are set based on the pre-filter violin plots (01_qc/qc_prefilter.pdf).
# Cells below min_features or above max_mt are removed. DoubletFinder is
# applied per sample to discard putative doublets.
#
# ┌─ ADJUST THRESHOLDS AFTER INSPECTING 01_qc/qc_prefilter.pdf ─────────────────
#   min_features : minimum number of detected genes per cell (default 200)
#   max_mt       : maximum mitochondrial read percentage  (default 5 %)
# └─────────────────────────────────────────────────────────────────────────────
output_dir <- dir_01   # both pre- and post-filter QC plots go to 01_qc/

seurat_list <- filter_seurat_samples(seurat_list_raw, min_features = 200, max_mt = 5)

plot_qc_batch(seurat_list, colors, "qc_postfilter.pdf")

# Checkpoint — restore with: seurat_list <- readRDS("resultados/objects/seurat_list_postfilter.rds")
saveRDS(seurat_list, file.path(dir_objects, "seurat_list_postfilter.rds"))


# =============================================================================

message("\n✓ SECTION 2 COMPLETE: Filtering and doublet detection complete")
# SECTION 3 — MERGE AND INITIAL PREPROCESSING
# =============================================================================
# Filtered samples are merged and preprocessed: log-normalization, variable
# feature selection (VST, 2,000 features), scaling, PCA (30 PCs), and UMAP.
# The resulting UMAP shows batch effects before integration.

output_dir <- dir_01

pbmc_harmony <- reduce(seurat_list, merge) %>%  # merge all samples into one object
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_preharmony.pdf")

# Checkpoint — restore with: pbmc_harmony <- readRDS("resultados/objects/pbmc_harmony_preharmony.rds")
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_preharmony.rds"))


# =============================================================================

message("\n✓ SECTION 3 COMPLETE: Merge and preprocessing complete")
# SECTION 4 — HARMONY BATCH CORRECTION
# =============================================================================
# Harmony adjusts cell embeddings to remove sample-level batch effects while
# preserving biological variation. All downstream steps use the "harmony"
# reduction instead of "pca".

output_dir <- dir_01

pbmc_harmony <- pbmc_harmony %>%
  RunHarmony("orig.ident", plot_convergence = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE)

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_postharmony.pdf")

# Checkpoint — restore with: pbmc_harmony <- readRDS("resultados/objects/pbmc_harmony_postharmony.rds")
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_postharmony.rds"))


# =============================================================================

message("\n✓ SECTION 4 COMPLETE: Harmony integration complete")
# SECTION 5 — RESOLUTION OPTIMIZATION
# =============================================================================
# Clustree tracks how Leiden communities split or merge across candidate
# resolutions. Choose the lowest resolution where clusters stabilize.
#
# ┌─ PARAMETERS ────────────────────────────────────────────────────────────────
#   resolutions_test : Leiden resolutions swept by clustree
#   → Inspect clustree.pdf before setting cluster_resolution in Section 6.
# └─────────────────────────────────────────────────────────────────────────────
resolutions_test <- c(0.15, 0.30, 0.50, 0.8, 1.0)

output_dir <- dir_02

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30,
                k.param = 30, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(clustree(clu, prefix = "RNA_snn_res."), "clustree.pdf", w = 18, h = 18)


# =============================================================================

message("\n✓ SECTION 5 COMPLETE: Clustree saved — inspect clustree.pdf before setting resolution")
# SECTION 6 — FINAL CLUSTERING
# =============================================================================
# Apply the selected resolution for the final cluster assignment.
# After clustering, a UMAP coloured by Seurat cluster and a bar chart of
# cells per sample are saved to 02_clustering/.
#
# ┌─ SET RESOLUTION AFTER INSPECTING elbow_plot.pdf AND clustree.pdf ──────────
#   cluster_resolution : Leiden resolution for final clustering (default 0.3)
# └─────────────────────────────────────────────────────────────────────────────
cluster_resolution <- 0.3
output_dir <- dir_02

# clu (Section 5) was temporary — re-run on pbmc_harmony to embed the final clusters
pbmc_harmony <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30,
                k.param = 30, verbose = FALSE) %>%
  FindClusters(resolution = cluster_resolution, algorithm = 4, verbose = FALSE)

Idents(pbmc_harmony) <- "seurat_clusters"
save_pdf(DimPlot(pbmc_harmony, group.by = "seurat_clusters", label = TRUE),
         "umap_seuratclusters.pdf")


# =============================================================================

message("\n✓ SECTION 6 COMPLETE: Final clustering complete")
# SECTION 7 — CELL-TYPE ANNOTATION
# =============================================================================
# Cell-type labels are assigned via bibliography-based markers
# (biblio_marks.txt): the marker table read above is crossed with
# cluster-level differential genes to assign initial labels automatically.
# Result stored in pbmc_harmony$celltype.
#
# A second strategy — reference transfer from a published Seurat object
# (e.g. an Arabidopsis leaf atlas) via FindTransferAnchors/TransferData — is
# not run by default. The full block, plus the downstream sections wired to
# it, is preserved in capitulo1_single_cell_reference_backup.R.

output_dir <- dir_03

biblio_marks_file <- file.path(DATA_DIR, "biblio_marks.txt")
marker_table      <- read.table(biblio_marks_file, header = TRUE, sep = "\t", quote = "")

# ── Bibliography-based annotation ─────────────────────────────────────────────
markers <- find_markers(pbmc_harmony,
                        output_file = file.path(output_dir, "FindAllMarkers.tsv"))

pbmc_harmony <- annotate_by_markers(pbmc_harmony, markers,
                                    reference_file = biblio_marks_file)
# Annotation stored in: pbmc_harmony$celltype

plot_marker_dotplot(
  pbmc_harmony,
  marker_table,
  annot_col = "celltype", # uses the newly assigned annotation column
  outfile   = file.path(output_dir, "dotplot_marker_table_annotation_biblio.pdf"),
  width = 18, height = 18
)

save_pdf(DimPlot(pbmc_harmony, group.by = "celltype",
                 label = TRUE, repel = TRUE, raster = FALSE),
         "umap_annotation_biblio.pdf")


# =============================================================================

message("\n✓ SECTION 7 COMPLETE: Cell-type annotation complete")
# SECTION 8 — ANNOTATED CLUSTREE
# =============================================================================
# Re-runs the resolution sweep with cell-type labels overlaid on each node.
# Confirms that the chosen resolution cleanly separates known cell types.

output_dir <- dir_03

# Mode: returns the most frequent value in a vector
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

clu <- pbmc_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30,
                k.param = 30, verbose = FALSE)

for (res in resolutions_test)
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)

save_pdf(
  clustree(clu, prefix = "RNA_snn_res.",
           node_label = "celltype", node_label_aggr = "Mode"),
  "clustree_annotated.pdf", w = 18, h = 18
)


# =============================================================================

message("\n✓ SECTION 8 COMPLETE: Annotated clustree saved")
# SECTION 9 — GENE EXPRESSION VISUALIZATION
# =============================================================================
# Violin and feature plots for individual genes or gene sets, generated both
# across all cell types and within a specific cell type of interest.
# JoinLayers() is required in Seurat 5 before subsetting after merge.
#
# ┌─ SET YOUR GENES AND CELL TYPE OF INTEREST ──────────────────────────────────
#   gene              : single gene to inspect
#   genes_of_interest : gene set to inspect together
#   celltype          : cell type to zoom into (must match celltype)
# └─────────────────────────────────────────────────────────────────────────────
gene              <- "AT5G26000"
genes_of_interest <- c("AT5G26000", "AT5G54250")
celltype          <- "guard cell"

output_dir <- dir_04

# JoinLayers is required in Seurat 5 before subsetting after merge
pbmc_harmony <- JoinLayers(pbmc_harmony)
Idents(pbmc_harmony) <- "celltype"

n_genes <- length(genes_of_interest)

# ── 10a. All cell types ───────────────────────────────────────────────────────
save_vln(VlnPlot(pbmc_harmony, features = gene),                  "vln_gene_all.pdf")
save_pdf(FeaturePlot(pbmc_harmony, features = gene),              "feature_gene_all.pdf")
save_vln(VlnPlot(pbmc_harmony, features = genes_of_interest),     "vln_geneset_all.pdf",     n = n_genes)
save_pdf(FeaturePlot(pbmc_harmony, features = genes_of_interest), "feature_geneset_all.pdf", w = 18, h = 18)

# ── 10b. Cell type of interest ────────────────────────────────────────────────
sub_obj <- subset(pbmc_harmony, idents = celltype)

save_vln(VlnPlot(sub_obj, features = gene),                  "vln_gene_celltype.pdf")
save_pdf(FeaturePlot(sub_obj, features = gene),              "feature_gene_celltype.pdf")
save_vln(VlnPlot(sub_obj, features = genes_of_interest),     "vln_geneset_celltype.pdf",     n = n_genes)
save_pdf(FeaturePlot(sub_obj, features = genes_of_interest), "feature_geneset_celltype.pdf", w = 18, h = 18)


# =============================================================================

message("\n✓ SECTION 9 COMPLETE: Expression visualization saved")
# SECTION 10 — CELL-TYPE GROUPING  [OPTIONAL]
# =============================================================================
# Fine-grained labels are collapsed into broader categories for downstream
# analyses. Cell types NOT listed in 'grouping' keep their original label.
# Skip this section if you do not need to merge cell types.
#
# ┌─ EDIT THIS MAP TO MATCH YOUR CELL TYPES ───────────────────────────────────
#   Left side  : original label (must match exactly)
#   Right side : new broader label to assign
# └─────────────────────────────────────────────────────────────────────────────
# Empty here: the bibliography annotation for this dataset (Section 7) didn't
# produce multiple fine-grained labels that share an obvious broader category,
# so celltype_grouped is a straight pass-through of celltype.
grouping <- c()

output_dir <- dir_05

# !!! unpacks the grouping vector as named arguments to recode()
# recode() errors on an empty replacement list, so an empty `grouping` (no
# merging needed) falls back to a straight pass-through instead.
if (length(grouping) > 0) {
  pbmc_harmony$celltype_grouped <- recode(pbmc_harmony$celltype, !!!grouping)
} else {
  pbmc_harmony$celltype_grouped <- pbmc_harmony$celltype
}

save_pdf(
  DimPlot(pbmc_harmony, group.by = "celltype_grouped",  # grouped broad cell types
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_grouped.pdf"
)


# =============================================================================

message("\n✓ SECTION 10 COMPLETE: Cell-type grouping complete")
# SECTION 11 — INTERACTIVE CELL-TYPE CURATION  [OPTIONAL]
# =============================================================================
# !!! WARNING: Run this section interactively, step by step.
# !!! Do NOT source the entire script with this section active.
# Skip this section if you are satisfied with the annotation from Section 7.
#
# Purpose: subcluster populations that appear heterogeneous in the UMAP,
# inspect them, and reassign cells to the correct cell type manually.
#
# Step 1 → subcluster the heterogeneous types
# Step 2 → save inspection figures; open them and decide on reassignments
# Step 3 → fill in the reassignment table below
# Step 4 → apply corrections to the global object

output_dir   <- dir_05
curation_col <- "celltype_grouped"   # starting annotation column for curation
Idents(pbmc_harmony) <- curation_col
table(pbmc_harmony[[curation_col]])

# ── Step 1. Subcluster ────────────────────────────────────────────
# "mesophyll" is the only population from Section 7 large enough to be worth
# inspecting here; the others (epidermis, root procambium.1/.2, endodermis,
# companion cell, guard cell) are left as-is.
mesophyll_umap <- subcluster_cell_type(pbmc_harmony, "mesophyll", annot_col = curation_col)

# Check how many subclusters this type produced:
# table(mesophyll_umap$cluster_subtipo)

# ── Step 2. Inspection figures ────────────────────────────────────
# Each call creates the DimPlot, saves it as PDF, and returns it for the composite
p_meso_dim <- plot_subcluster_umap(mesophyll_umap, "mesophyll", output_dir)

# Composite: each row = [ UMAP | marker genes ] for one cell type
save_subcluster_composite(
  subcluster_list = list(
    list(umap_plot = p_meso_dim, obj = mesophyll_umap)
  ),
  marker_table = marker_table,
  output_dir   = output_dir
)

# ── Step 3. Reassignment table ────────────────────────────────────
# Map subcluster IDs \u2192 final cell-type labels.
# "others" is a catch-all for any subcluster ID not listed.
# Names must match the variable names used in Step 1 exactly.
reassign <- list(
  mesophyll_umap = c(
    "0"      = "mesophyll",
    "1"      = "mesophyll",
    "2"      = "mesophyll",
    "others" = "mesophyll"
  )
)

# ── Step 4. Apply corrections ───────────────────────────────────────────────
subcluster_list <- list(
  mesophyll_umap = mesophyll_umap
)

pbmc_harmony <- apply_subcluster_reassignment(
  obj             = pbmc_harmony,
  subcluster_list = subcluster_list,
  reassign        = reassign,
  source_col      = curation_col,
  dest_col        = "celltype_curated"
)

save_pdf(
  DimPlot(pbmc_harmony, group.by = "celltype_curated",
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_curated.pdf"
)

# Checkpoint — restore with: pbmc_harmony <- readRDS("resultados/objects/pbmc_harmony_curated.rds")
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_curated.rds"))

# =============================================================================

message("\n✓ SECTION 11 COMPLETE: Curation complete — curated object saved")
# SECTION 12 — EXPORT TO H5AD (Scanpy / Python)
# =============================================================================
# Exports the curated object to AnnData h5ad format for Python-based
# trajectory and velocity analyses (Scanpy, scFates, Palantir — all
# pre-installed in the Docker image).

export_to_scanpy(pbmc_harmony,
                 file.path(dir_objects, "pbmc_harmony_curated.h5ad"))

# To export a specific cell type:
# export_to_scanpy(
#   subset(pbmc_harmony, subset = celltype_curated == "Guard Cell"),
#   file.path(dir_objects, "GuardCell.h5ad")
# )


# =============================================================================

message("\n✓ SECTION 12 COMPLETE: Export to h5ad complete")
