# Methods Paper - Chapter 2 Skeleton

Working title: Pseudobulk differential expression and functional enrichment
from the Chapter 1 WT vs pifq object.

This is a planning skeleton for the next methods-paper section. It should be
converted into the final PDF/HTML format after the WT/pifq Chapter 2 workflow
has been tested end to end.

## Step 1 - Start From the Chapter 1 Object

Goal: load the bibliography-annotated Seurat object produced by Chapter 1,
before grouping or manual curation.

Expected input:

```text
resultados_wt/objects/pbmc_harmony_annotated.rds
```

Draft code:

```r
PIPELINE_DIR <- "/workspace/workflow"
DATA_DIR     <- "/workspace"
base_dir     <- file.path(DATA_DIR, "resultados_wt")

source(file.path(PIPELINE_DIR, "load_libraries.R"))
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))

set.seed(1807)
setwd(DATA_DIR)
list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)

pbmc_harmony <- readRDS(file.path(dir_objects, "pbmc_harmony_annotated.rds"))
```

## Step 2 - Check Metadata Required for Pseudobulk

Goal: verify the annotation and condition columns before any aggregation.

Draft code:

```r
table(pbmc_harmony$condition, useNA = "ifany")
table(pbmc_harmony$orig.ident, useNA = "ifany")
table(pbmc_harmony$celltype, useNA = "ifany")
```

Decision needed: confirm that WT and pifq have enough biological or technical
replicate structure for the intended inference. If only one library per
condition is available, document pseudo-replicates as a demonstration and avoid
overstating biological replication.

## Step 3 - Create Cell-Type Subsets

Goal: split the annotated object by direct bibliography cell-type annotation.

Draft code:

```r
pseudobulk_annot_col <- "celltype"

cell_type_subsets <- create_cell_type_subsets(
  pbmc_harmony,
  annot_col = pseudobulk_annot_col
)
```

Expected output: one subset per cell type.

## Step 4 - Assign Replicates

Goal: create replicate labels for pseudobulk aggregation.

Draft code:

```r
pseudobulk_conditions <- c("WT", "pifq")
n_pseudoreps <- 3

cell_type_subsets_replicates <- assign_pseudoreplicates_batch(
  cell_type_subsets,
  pseudobulk_conditions = pseudobulk_conditions,
  n_pseudoreps = n_pseudoreps
)
```

Documentation note: this step must explain whether these are true biological
replicates or pseudo-replicates.

## Step 5 - Run Pseudobulk DESeq2

Goal: aggregate counts and run WT vs pifq differential expression per cell type.

Draft code:

```r
comparisons <- list(
  list(conds = c("WT", "pifq"), tag = "WT_vs_pifq")
)

cell_types_to_analyze <- NULL
output_dir <- dir_06

deseq2_results <- run_pseudobulk_deseq2_analysis(
  cell_type_subsets_replicates = cell_type_subsets_replicates,
  comparisons = comparisons,
  output_dir = output_dir,
  cell_types = cell_types_to_analyze,
  pseudobulk_dir = file.path(dir_objects, "pseudobulk_replicas")
)
```

## Step 6 - Volcano Plots

Goal: make one volcano per cell type for the WT vs pifq contrast.

Draft code:

```r
volcano_tag <- "WT_vs_pifq"
padj_cut <- 0.05
lfc_cut <- 1

render_volcano_plots(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag, "volcano"),
  pdf_name    = paste0("VolcanoPlots_", volcano_tag, ".pdf"),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut
)
```

## Step 7 - Differential Gene Summary Tables

Goal: combine DESeq2 results across cell types and prepare matrices for plots.

Draft code:

```r
diff_prefix <- paste0("diff_table_", volcano_tag)

diff_tables <- build_differential_tables(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut,
  prefix      = diff_prefix
)
```

## Step 8 - GO Enrichment

Goal: run GO enrichment for significant genes per cell type.

Draft code:

```r
go_space <- "BP"
go_orgdb <- org.At.tair.db
go_keytype <- "TAIR"

go_results <- run_go_enrichment_for_contrast(
  results_dir  = file.path(dir_06, volcano_tag),
  output_dir   = file.path(dir_07, volcano_tag),
  orgdb        = go_orgdb,
  keytype      = go_keytype,
  go_space     = go_space,
  padj_cutoff  = padj_cut,
  contrast_tag = volcano_tag
)
```

## Step 9 - Log2FC Heatmap

Goal: summarize cell-type-specific log2FC patterns.

Draft code:

```r
heatmap_limits <- c(-5, 5)

build_logfc_heatmap(
  logfc_table  = diff_tables$logfc,
  contrast_tag = volcano_tag,
  output_dir   = file.path(dir_06, volcano_tag),
  limits       = heatmap_limits
)
```

## Step 10 - Optional Network Analysis

Goal: run hdWGCNA only if the DE table is large enough and the object has enough
cells per cell type/sample combination.

Draft code:

```r
wgcna_name <- "WT_pifq"
n_metacells <- 50
soft_power <- NULL
min_module_size <- 20
deep_split <- 2

run_unified_hdwgcna(
  seurat_obj    = pbmc_harmony,
  de_table_path = file.path(dir_06, volcano_tag, "tabla_log2FC_fc1_padj_005.tsv"),
  output_dir    = dir_08,
  annot_col     = pseudobulk_annot_col,
  sample_col    = "orig.ident",
  wgcna_name    = wgcna_name,
  n_metacells   = n_metacells,
  soft_power    = soft_power,
  min_module_size = min_module_size,
  deep_split    = deep_split
)
```

## Figures to Include in the Chapter 2 Methods Draft

Use small reference images, following the Chapter 1 style:

- Pseudobulk sample/condition check.
- One representative volcano plot.
- Combined volcano PDF reference.
- Log2FC heatmap.
- One representative GO bubble plot.
- Optional network figure if hdWGCNA is kept.

## Open Decisions Before Finalizing Chapter 2

- Confirm whether WT and pifq have true biological replicates.
- Decide how to describe pseudo-replicates if true replicates are absent.
- Decide whether Chapter 2 should stop at DESeq2/GO or include hdWGCNA.
- Test the full WT vs pifq run before rendering final figures.
- Keep the final PDF free of large tables; use text and compact figures.
