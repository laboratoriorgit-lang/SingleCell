# =============================================================================
# capitulo2_wt_pifq.R — Chapter 2 pseudobulk DE for WT vs pifq
# =============================================================================
# Requires Chapter 1 output:
#   resultados_wt/objects/pbmc_harmony_annotated.rds
#
# This WT/pifq adaptation uses the direct bibliography annotation column
# `celltype`. It does not use celltype_grouped or celltype_curated.
# =============================================================================

PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", unset = "/workspace/workflow")
DATA_DIR     <- Sys.getenv("DATA_DIR",     unset = "/workspace")
base_dir     <- file.path(DATA_DIR, "resultados_wt")

source(file.path(PIPELINE_DIR, "load_libraries.R"))
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))

set.seed(1807)
setwd(DATA_DIR)
list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)

log_msg <- function(...) message("\n", paste0(...))

# =============================================================================
# SECTION 13 — LOAD CHAPTER 1 ANNOTATED OBJECT
# =============================================================================

pbmc_harmony <- readRDS(file.path(dir_objects, "pbmc_harmony_annotated.rds"))

pseudobulk_annot_col <- "celltype"
pseudobulk_conditions <- c("WT", "pifq")
comparison_tag <- "WT_vs_pifq"

stopifnot(pseudobulk_annot_col %in% colnames(pbmc_harmony@meta.data))
stopifnot("condition" %in% colnames(pbmc_harmony@meta.data))

output_dir <- dir_06

celltype_condition_counts <- as.data.frame.matrix(
  table(pbmc_harmony[[pseudobulk_annot_col]][, 1], pbmc_harmony$condition)
)
celltype_condition_counts$celltype <- rownames(celltype_condition_counts)
celltype_condition_counts <- celltype_condition_counts[, c("celltype", pseudobulk_conditions)]

write.table(
  celltype_condition_counts,
  file.path(dir_06, "celltype_condition_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

counts_long <- tidyr::pivot_longer(
  celltype_condition_counts,
  cols = all_of(pseudobulk_conditions),
  names_to = "condition",
  values_to = "cells"
)
counts_long$celltype <- factor(counts_long$celltype, levels = celltype_condition_counts$celltype)

p_cell_counts <- ggplot(counts_long, aes(celltype, cells, fill = condition)) +
  geom_col(position = "dodge", width = 0.75) +
  scale_y_continuous(trans = "log10") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = "Cells per cell type (log10 scale)", fill = NULL)

save_pdf(p_cell_counts, "celltype_condition_counts.pdf", w = 12, h = 8)

log_msg("SECTION 13 COMPLETE: Annotated object loaded and metadata checked")

# =============================================================================
# SECTION 14 — CELL-TYPE SUBSETS USING `celltype`
# =============================================================================

min_cells_per_condition <- 20

eligible_celltypes <- celltype_condition_counts$celltype[
  apply(celltype_condition_counts[, pseudobulk_conditions, drop = FALSE],
        1, function(x) all(x >= min_cells_per_condition))
]

writeLines(
  c(
    "Cell types retained for Chapter 2 pseudobulk DE:",
    eligible_celltypes
  ),
  file.path(dir_06, "celltypes_retained_for_pseudobulk.txt")
)

cell_type_subsets <- create_cell_type_subsets(
  pbmc_harmony,
  annot_col = pseudobulk_annot_col
)

eligible_subset_names <- gsub("[^[:alnum:]_]", "_", eligible_celltypes)
cell_type_subsets <- cell_type_subsets[names(cell_type_subsets) %in% eligible_subset_names]

log_msg("SECTION 14 COMPLETE: Cell-type subsets created from celltype")

# =============================================================================
# SECTION 15 — PSEUDO-REPLICATE ASSIGNMENT
# =============================================================================

n_pseudoreps <- 3

cell_type_subsets_replicates <- assign_pseudoreplicates_batch(
  cell_type_subsets,
  pseudobulk_conditions = pseudobulk_conditions,
  n_pseudoreps = n_pseudoreps
)

replicate_counts <- do.call(rbind, lapply(names(cell_type_subsets_replicates), function(ct) {
  tab <- table(cell_type_subsets_replicates[[ct]]$condition,
               cell_type_subsets_replicates[[ct]]$replicate)
  data.frame(
    celltype = ct,
    condition = rownames(tab)[row(tab)],
    replicate = colnames(tab)[col(tab)],
    cells = as.vector(tab),
    row.names = NULL
  )
}))

write.table(
  replicate_counts,
  file.path(dir_06, "pseudo_replicate_cell_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p_reps <- ggplot(replicate_counts, aes(replicate, cells, fill = condition)) +
  geom_col(width = 0.75) +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Cells per pseudo-replicate", fill = NULL)

save_pdf(p_reps, "pseudo_replicate_cell_counts.pdf", w = 16, h = 12)

log_msg("SECTION 15 COMPLETE: Pseudo-replicates assigned")

# =============================================================================
# SECTION 16 — PSEUDOBULK TABLES AND DESEQ2
# =============================================================================

comparisons <- list(
  list(conds = pseudobulk_conditions, tag = comparison_tag)
)

output_dir <- dir_06

deseq2_results <- run_pseudobulk_deseq2_analysis(
  cell_type_subsets_replicates = cell_type_subsets_replicates,
  comparisons = comparisons,
  output_dir = output_dir,
  cell_types = NULL,
  pseudobulk_dir = file.path(dir_objects, "pseudobulk_replicas_celltype")
)

log_msg("SECTION 16 COMPLETE: Pseudobulk tables and DESeq2 complete")

# =============================================================================
# SECTION 17 — VOLCANO PLOTS
# =============================================================================

volcano_tag <- comparison_tag
padj_cut <- 0.05
lfc_cut <- 1

render_volcano_plots(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag, "volcano"),
  pdf_name    = paste0("VolcanoPlots_", volcano_tag, ".pdf"),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut
)

log_msg("SECTION 17 COMPLETE: Volcano plots saved")

# =============================================================================
# SECTION 18 — DIFFERENTIAL TABLES AND LOG2FC HEATMAP
# =============================================================================

diff_prefix <- paste0("diff_table_", volcano_tag)

diff_tables <- build_differential_tables(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut,
  prefix      = diff_prefix
)

heatmap_limits <- c(-5, 5)

build_logfc_heatmap(
  logfc_table  = diff_tables$logfc,
  contrast_tag = volcano_tag,
  output_dir   = file.path(dir_06, volcano_tag),
  limits       = heatmap_limits
)

log_msg("SECTION 18 COMPLETE: Differential tables and heatmap saved")

# =============================================================================
# SECTION 19 — GO ENRICHMENT
# =============================================================================

go_space <- "BP"
go_orgdb <- org.At.tair.db
go_keytype <- "TAIR"
padj_cutoff <- 0.05

go_results <- run_go_enrichment_for_contrast(
  results_dir  = file.path(dir_06, volcano_tag),
  output_dir   = file.path(dir_07, volcano_tag),
  orgdb        = go_orgdb,
  keytype      = go_keytype,
  go_space     = go_space,
  padj_cutoff  = padj_cutoff,
  contrast_tag = volcano_tag
)

log_msg("SECTION 19 COMPLETE: GO enrichment complete")

writeLines(capture.output(sessionInfo()), file.path(base_dir, "sessionInfo_chapter2.txt"))
log_msg("Chapter 2 WT/pifq run complete")

