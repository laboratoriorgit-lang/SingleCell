# Step 2 - Pseudobulk differential expression and GO (R)
# Auto-derived from the methods paper (Part 3). Loads the curated object
# saved by Step 1, then runs DE, volcano, tables, GO and the TF network.

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
# Section 12 - Export the curated object
# ==============================================================================
pbmc_harmony <- readRDS(file.path(dir_objects, "pbmc_harmony_curated.rds"))

# Bibliography marker table (defined during Step 1 annotation; reloaded here so
# this script is self-contained). Used by the Section 19 heatmap row/column
# ordering and the Section 20 network.
biblio_marks_file <- file.path(DATA_DIR, "data", "biblio_marks.txt")
marker_table <- read.table(
  biblio_marks_file, header = TRUE, sep = "\t", quote = ""
)

# ==============================================================================
# Section 13 - Cell-type subsets
# ==============================================================================
pseudobulk_annot_col <- "celltype_curated"

# Inspect how many cells each cell type has before subsetting
table(pbmc_harmony[[pseudobulk_annot_col]])

cell_type_subsets <- create_cell_type_subsets(
  pbmc_harmony, annot_col = pseudobulk_annot_col
)

message("\nOK SECTION 13 COMPLETE: cell-type subsets created")

# ==============================================================================
# Section 14 - Pseudo-replicate assignment
# ==============================================================================
pseudobulk_conditions <- NULL
n_pseudoreps          <- 3

cell_type_subsets_replicates <- assign_pseudoreplicates_batch(
  cell_type_subsets,
  pseudobulk_conditions = pseudobulk_conditions,
  n_pseudoreps          = n_pseudoreps
)

# Inspect cells per pseudo-replicate for one cell type
table(cell_type_subsets_replicates[["Mesophyll"]]$replicate)

message("\nOK SECTION 14 COMPLETE: pseudo-replicates assigned")

# ==============================================================================
# Section 15 - Pseudobulk tables and DESeq2
# ==============================================================================
comparisons <- list(
  list(conds = c("WT", "pifq"), tag = "WT_vs_pifq")
)
cell_types_to_analyze <- NULL
output_dir <- dir_06

deseq2_results <- run_pseudobulk_deseq2_analysis(
  cell_type_subsets_replicates = cell_type_subsets_replicates,
  comparisons    = comparisons,
  output_dir     = output_dir,
  cell_types     = cell_types_to_analyze,
  pseudobulk_dir = file.path(dir_objects, "pseudobulk_replicas")
)

message("\nOK SECTION 15 COMPLETE: pseudobulk aggregation and DESeq2 complete")

# ==============================================================================
# Section 16 - Volcano plots
# ==============================================================================
volcano_tag <- "WT_vs_pifq"
padj_cut    <- 0.05
lfc_cut     <- 1

render_volcano_plots(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag, "volcano"),
  pdf_name    = paste0("VolcanoPlots_", volcano_tag, ".pdf"),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut
)

message("\nOK SECTION 16 COMPLETE: volcano plots saved")

# ==============================================================================
# Section 17 - Differential gene tables
# ==============================================================================
diff_prefix <- paste0("diff_table_", volcano_tag)

diff_tables <- build_differential_tables(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut,
  prefix      = diff_prefix
)

message("\nOK SECTION 17 COMPLETE: differential gene tables saved")

# ==============================================================================
# Section 18 - GO enrichment
# ==============================================================================
go_space    <- "BP"
go_orgdb    <- org.At.tair.db
go_keytype  <- "TAIR"
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

message("\nOK SECTION 18 COMPLETE: GO enrichment complete")

# ==============================================================================
# Section 19 - Log2FC heatmap
# ==============================================================================
heatmap_limits <- c(-5, 5)

build_logfc_heatmap(
  logfc_table  = diff_tables$logfc,
  contrast_tag = volcano_tag,
  output_dir   = file.path(dir_06, volcano_tag),
  limits       = heatmap_limits,
  cell_order   = unique(marker_table$cell.types)
)

message("\nOK SECTION 19 COMPLETE: log2FC heatmap saved")

# ==============================================================================
# Section 20 - Gene-gene coexpression and TF network
# ==============================================================================
run_tf_coexpression_network(
  seurat_obj      = pbmc_harmony,
  de_table_path   = file.path(dir_06, volcano_tag, "tabla_log2FC_fc1_padj_005.tsv"),
  tf_list_path    = file.path(DATA_DIR, "data", "AtTFDB_loci.txt"),
  output_dir      = dir_08,
  annot_col       = pseudobulk_annot_col,
  contrast_tag    = volcano_tag,
  n_metacells     = 50,     # metacells per cell type x sample
  min_module_size = 20,     # minimum genes per co-expression module
  deep_split      = 2,      # module-splitting sensitivity (0-4)
  soft_power      = NULL,   # NULL = auto-detect the soft power
  tom_threshold   = 0.2,    # minimum co-expression weight to keep an edge
  n_hub_label     = 15      # number of top hub TFs to label
)

message("\nOK SECTION 20 COMPLETE: TF coexpression network saved")
