# =============================================================================
# Chapter 1 - WT vs pifq single-cell workflow
# Project: /home/mvergara/projects2/Sc_DB_test
# Output:  /home/mvergara/projects2/Sc_DB_test/resultados_wt
# =============================================================================

PIPELINE_DIR <- "/workspace/SingleCell/workflow"
DATA_DIR     <- "/workspace"
base_dir     <- file.path(DATA_DIR, "resultados_wt")

source(file.path(PIPELINE_DIR, "load_libraries.R"))
source(file.path(PIPELINE_DIR, "custom_seurat.R"))
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))

set.seed(1807)
options(Seurat.allow.s4 = FALSE)
setwd(DATA_DIR)

list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)
output_dir <- base_dir

log_msg <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

log_msg("Starting Chapter 1 WT/pifq run")

# Section 0 - Cell Ranger provenance ------------------------------------------
dir_00 <- file.path(base_dir, "00_cellranger")
dir.create(dir_00, recursive = TRUE, showWarnings = FALSE)

cellranger_inputs <- data.frame(
  sample = c("WT", "pifq"),
  cellranger_id = c("ScWT", "Scpifq"),
  matrix_dir = c(
    "ScWT/outs/filtered_feature_bc_matrix",
    "Scpifq/outs/filtered_feature_bc_matrix"
  ),
  metrics_summary = c(
    "ScWT/outs/metrics_summary.csv",
    "Scpifq/outs/metrics_summary.csv"
  ),
  condition = c("WT", "pifq"),
  stringsAsFactors = FALSE
)

write.table(
  cellranger_inputs,
  file.path(dir_00, "cellranger_inputs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

for (metrics in cellranger_inputs$metrics_summary) {
  if (!file.exists(metrics)) stop("Missing Cell Ranger metrics file: ", metrics)
  file.copy(metrics, file.path(dir_00, paste0(dirname(dirname(metrics)), "_metrics_summary.csv")),
            overwrite = TRUE)
}

writeLines(c(
  "Cell Ranger processing documented from local scripts:",
  "  cell_CRA010863_runs.sh",
  "  cell_CRA010863_ensure_ScWT.sh",
  "",
  "Observed command pattern:",
  "  cellranger count --localcores=80 --id=<ScWT|Scpifq> --fastqs=/home/mvergara/projects2/Sc_DB_test/CRA010863 --sample=<ScWT|Scpifq> --transcriptome=/home/mvergara/projects2/eleo/ftp1.cruk.cam.ac.uk/Seba/Ara --no-bam"
), file.path(dir_00, "cellranger_processing_notes.txt"))

# Section 1 - Data loading and pre-filter QC ----------------------------------
output_dir <- dir_01

samples <- list(
  list(file = "ScWT/outs/filtered_feature_bc_matrix",   label = "WT",   condition = "WT"),
  list(file = "Scpifq/outs/filtered_feature_bc_matrix", label = "pifq", condition = "pifq")
)

colors <- c(WT = "#66c2a5", pifq = "#fc8d62")
mt_pattern <- "^ATMG"
cp_pattern <- "^ATCG"

seurat_list_raw <- load_seurat_samples(
  samples = samples,
  DATA_DIR = DATA_DIR,
  mt_pattern = mt_pattern,
  cp_pattern = cp_pattern
)

plot_qc_batch(seurat_list_raw, colors, "qc_prefilter.pdf")
saveRDS(seurat_list_raw, file.path(dir_objects, "seurat_list_raw.rds"))
log_msg("Section 1 complete")

# Section 2 - Cell filtering and doublet detection ----------------------------
output_dir <- dir_01

seurat_list <- filter_seurat_samples(
  seurat_list_raw,
  min_features = 200,
  max_mt = 5
)

plot_qc_batch(seurat_list, colors, "qc_postfilter.pdf")
saveRDS(seurat_list, file.path(dir_objects, "seurat_list_postfilter.rds"))
log_msg("Section 2 complete")

# Section 3 - Merge and initial preprocessing ---------------------------------
output_dir <- dir_01

pbmc_harmony <- reduce(seurat_list, merge) |>
  NormalizeData(verbose = FALSE) |>
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |>
  RunPCA(npcs = 30, verbose = FALSE) |>
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_preharmony.pdf")
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_preharmony.rds"))
log_msg("Section 3 complete")

# Section 4 - Harmony batch correction ----------------------------------------
output_dir <- dir_01

pbmc_harmony <- pbmc_harmony |>
  RunHarmony("orig.ident", plot_convergence = FALSE) |>
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE)

save_pdf(DimPlot(pbmc_harmony, group.by = "orig.ident", cols = colors),
         "umap_postharmony.pdf")
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_postharmony.rds"))
log_msg("Section 4 complete")

# Section 5 - Resolution optimization -----------------------------------------
resolutions_test <- c(0.15, 0.35, 0.45, 0.55, 1.0)
output_dir <- dir_02

k_range <- 1:31
pca_data <- Embeddings(pbmc_harmony, "pca")[, 1:30]
wss <- sapply(k_range, function(k) kmeans(pca_data, centers = k, nstart = 4)$tot.withinss)

elbow_plot <- ggplot(data.frame(k = k_range, wss = wss), aes(k, wss)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of clusters (k)", y = "Within-cluster sum of squares") +
  theme_minimal()

save_pdf(elbow_plot, "elbow_plot.pdf", w = 18, h = 18)

clu <- pbmc_harmony |>
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, k.param = 20, verbose = FALSE)

for (res in resolutions_test) {
  clu <- FindClusters(clu, resolution = res, algorithm = 4, verbose = FALSE)
}

save_pdf(clustree(clu, prefix = "RNA_snn_res."), "clustree2.pdf", w = 18, h = 18)
saveRDS(clu, file.path(dir_objects, "resolution_sweep_clu.rds"))
log_msg("Section 5 complete")

# Section 6 - Final clustering -------------------------------------------------
cluster_resolution <- 0.35
output_dir <- dir_02

pbmc_harmony <- pbmc_harmony |>
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) |>
  FindNeighbors(reduction = "harmony", dims = 1:30, k.param = 20, verbose = FALSE) |>
  FindClusters(resolution = cluster_resolution, algorithm = 4, verbose = FALSE)

Idents(pbmc_harmony) <- "seurat_clusters"
save_pdf(DimPlot(pbmc_harmony, group.by = "seurat_clusters", label = TRUE),
         "umap_seuratclusters.pdf")
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_clusters.rds"))
log_msg("Section 6 complete")

# Section 7 - Bibliography marker annotation ----------------------------------
output_dir <- dir_03

biblio_marks_file <- file.path(DATA_DIR, "biblio_marks_custom.txt")
marker_table <- read.table(biblio_marks_file, header = TRUE, sep = "\t", quote = "")

markers <- find_markers(
  pbmc_harmony,
  output_file = file.path(output_dir, "FindAllMarkers.tsv")
)

pbmc_harmony <- annotate_by_markers(
  pbmc_harmony,
  markers,
  reference_file = biblio_marks_file
)

plot_marker_dotplot(
  pbmc_harmony,
  marker_table,
  annot_col = "celltype",
  outfile = file.path(output_dir, "dotplot_marker_table_annotation_biblio.pdf"),
  width = 18,
  height = 18
)

save_pdf(DimPlot(pbmc_harmony, group.by = "celltype", label = TRUE, repel = TRUE, raster = FALSE),
         "umap_annotation_biblio.pdf")
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_annotated.rds"))
log_msg("Section 7 complete")

# Section 8 - Annotated clustree ----------------------------------------------
output_dir <- dir_03

Mode <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

celltype_label <- as.character(pbmc_harmony$celltype)
names(celltype_label) <- Cells(pbmc_harmony)
clu$celltype_label <- celltype_label[Cells(clu)]

save_pdf(
  clustree(clu, prefix = "RNA_snn_res.",
           node_label = "celltype_label", node_label_aggr = "Mode"),
  "clustree_annotated.pdf", w = 14, h = 14
)
log_msg("Section 8 complete")

# Section 9 - Simple gene expression visualization ----------------------------
gene <- "AT5G26000"
output_dir <- dir_04

pbmc_harmony <- JoinLayers(pbmc_harmony)
Idents(pbmc_harmony) <- "celltype"

save_vln(VlnPlot(pbmc_harmony, features = gene), "vln_gene_all.pdf")
save_pdf(FeaturePlot(pbmc_harmony, features = gene), "feature_gene_all.pdf")
log_msg("Section 9 complete")

# Section 10 - Grouping pass-through ------------------------------------------
grouping <- c()
output_dir <- dir_05

if (length(grouping) > 0) {
  pbmc_harmony$celltype_grouped <- recode(pbmc_harmony$celltype, !!!grouping)
} else {
  pbmc_harmony$celltype_grouped <- pbmc_harmony$celltype
}

save_pdf(
  DimPlot(pbmc_harmony, group.by = "celltype_grouped", label = TRUE, repel = TRUE, raster = FALSE),
  "umap_grouped.pdf"
)

pbmc_harmony$celltype_curated <- pbmc_harmony$celltype_grouped
saveRDS(pbmc_harmony, file.path(dir_objects, "pbmc_harmony_curated.rds"))
log_msg("Section 10 complete")

# Section 11 - Optional subcluster inspection template ------------------------
output_dir <- dir_05
inspection_note <- file.path(output_dir, "section_11_subcluster_inspection_template.R")
writeLines(c(
  "output_dir <- dir_05",
  "curation_col <- \"celltype\"",
  "Idents(pbmc_harmony) <- curation_col",
  "print(table(pbmc_harmony[[curation_col]], useNA = \"ifany\"))",
  "types_to_inspect <- c(\"Epidermis Hypocotyl.1\")",
  "subcluster_resolution <- 0.3",
  "subcluster_dims <- 1:20",
  "subcluster_list <- list()",
  "inspection_entries <- list()",
  "for (ct in types_to_inspect) {",
  "  safe_name <- make.names(ct)",
  "  file_tag <- gsub(\"[^A-Za-z0-9]+\", \"_\", ct)",
  "  sub_obj <- subcluster_cell_type(pbmc_harmony, tipo = ct, annot_col = curation_col, resolution = subcluster_resolution, dims = subcluster_dims)",
  "  subcluster_list[[safe_name]] <- sub_obj",
  "  p_dim <- plot_subcluster_umap(sub_obj, ct, output_dir)",
  "  plot_marker_dotplot(sub_obj, marker_table, annot_col = \"cluster_subtipo\", outfile = file.path(output_dir, paste0(\"dotplot_markers_subclusters_\", file_tag, \".pdf\")), width = 18, height = 18)",
  "  inspection_entries[[safe_name]] <- list(umap_plot = p_dim, obj = sub_obj)",
  "}",
  "save_subcluster_composite(inspection_entries, marker_table, output_dir, filename = \"subclustering_marker_inspection.pdf\")"
), inspection_note)
log_msg("Section 11 template written")

# Section 12 - h5ad export -----------------------------------------------------
export_to_scanpy(
  pbmc_harmony,
  file.path(dir_objects, "pbmc_harmony_curated.h5ad")
)
log_msg("Section 12 complete")

writeLines(capture.output(sessionInfo()), file.path(base_dir, "sessionInfo.txt"))
log_msg("Chapter 1 WT/pifq run complete")
