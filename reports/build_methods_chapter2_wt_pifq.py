from __future__ import annotations

import html
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "methods_paper_chapter2_wt_pifq.html"


def esc(value: object) -> str:
    return html.escape("" if value is None else str(value))


def code_block(language: str, text: str) -> str:
    return (
        '<div class="code-wrap">'
        f'<div class="code-label">{esc(language)}</div>'
        f"<pre><code>{esc(text.strip())}</code></pre>"
        "</div>"
    )


def img(filename: str, caption: str, cls: str = "") -> str:
    return (
        f'<figure class="{esc(cls)}">'
        f'<img src="assets/{esc(filename)}" alt="{esc(caption)}">'
        f"<figcaption>{esc(caption)}</figcaption>"
        "</figure>"
    )


def build() -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Chapter 2 WT pifq Methods Tutorial</title>
  <style>
    @page {{
      size: Letter;
      margin: 0.72in 0.78in;
      @bottom-center {{
        content: counter(page);
        color: #6b7280;
        font-size: 8.5pt;
      }}
    }}

    body {{
      color: #1f2933;
      font-family: Georgia, "Times New Roman", serif;
      font-size: 10.8pt;
      line-height: 1.43;
    }}

    .title-page {{
      margin-top: 0.55in;
    }}

    .title {{
      color: #123f5a;
      font-size: 23pt;
      font-weight: 700;
      line-height: 1.12;
      margin: 0 0 0.2em;
      text-align: center;
    }}

    .title em {{
      display: block;
      font-size: 20pt;
    }}

    .subtitle {{
      color: #2f6b49;
      font-size: 13.2pt;
      font-style: italic;
      margin: 1.25em 0 1.15em;
    }}

    .meta {{
      color: #374151;
      margin: 0.2em 0;
    }}

    .rule {{
      border-top: 2px solid #7ba7c0;
      margin: 1.05em 0 1.65em;
    }}

    h1, h2 {{
      color: #123f5a;
      font-family: Georgia, "Times New Roman", serif;
      font-weight: 700;
      line-height: 1.16;
      page-break-after: avoid;
    }}

    h1 {{
      border-bottom: 1.6px solid #a9c5d6;
      font-size: 17pt;
      margin: 1.2em 0 0.75em;
      padding-bottom: 0.22em;
    }}

    h2 {{
      border-bottom: 1px solid #b8ccd9;
      font-size: 13.3pt;
      margin: 1.05em 0 0.55em;
      padding-bottom: 0.15em;
    }}

    p {{
      margin: 0.42em 0 0.62em;
      text-align: justify;
    }}

    ul {{
      margin: 0.35em 0 0.75em 1.35em;
      padding: 0;
    }}

    li {{
      margin: 0.2em 0;
    }}

    code {{
      color: #174a67;
      font-family: "SFMono-Regular", Menlo, Consolas, monospace;
      font-size: 9.4pt;
    }}

    .code-wrap {{
      border: 1px solid #b7b7b7;
      margin: 0.62em 0 1.05em;
      page-break-inside: auto;
      position: relative;
    }}

    .code-label {{
      background: white;
      color: #2f6b49;
      font-size: 8.2pt;
      font-style: italic;
      left: 1.1em;
      padding: 0 0.35em;
      position: absolute;
      top: -0.68em;
    }}

    pre {{
      color: #174a67;
      font-family: "SFMono-Regular", Menlo, Consolas, monospace;
      font-size: 9.2pt;
      line-height: 1.33;
      margin: 0;
      overflow-wrap: break-word;
      padding: 0.72em 0.85em;
      white-space: pre-wrap;
    }}

    .page-section {{
      break-before: page;
    }}

    .small-note {{
      color: #4b5563;
      font-size: 9.4pt;
    }}

    figure {{
      margin: 0.45em auto 0.75em;
      page-break-inside: avoid;
      text-align: center;
    }}

    figure img {{
      display: block;
      height: auto;
      margin: 0 auto;
      max-height: 3.15in;
      max-width: 48%;
    }}

    figure.wide img {{
      max-height: 3.45in;
      max-width: 62%;
    }}

    figure.tall img {{
      max-height: 4.55in;
      max-width: 62%;
    }}

    .two-col {{
      display: grid;
      gap: 0.35in;
      grid-template-columns: 1fr 1fr;
      page-break-inside: avoid;
    }}

    .two-col figure img {{
      max-height: 2.85in;
      max-width: 94%;
    }}

    figcaption {{
      color: #4b5563;
      font-size: 8.9pt;
      margin-top: 0.35em;
    }}
  </style>
</head>
<body>
  <section class="title-page">
    <h1 class="title">Computational Methods for Single-Cell RNA-Seq Analysis of<br><em>Arabidopsis thaliana</em></h1>
    <p class="subtitle">Chapter 2 - WT vs pifq tutorial: pseudobulk differential expression and functional enrichment</p>
    <p class="meta">SingleCell Pipeline</p>
    <p class="meta">2026-07-15</p>
    <div class="rule"></div>

    <h1>Purpose</h1>
    <p>
      Chapter 2 starts from the annotated object produced in Chapter 1 and runs
      pseudobulk differential expression by cell type. For this WT/pifq draft,
      the analysis intentionally uses the direct bibliography annotation column
      <code>celltype</code>. It does not use <code>celltype_grouped</code> and
      does not use <code>celltype_curated</code>.
    </p>
    <p>
      The WT/pifq test currently has one Cell Ranger library per condition.
      Therefore the pseudo-replicates used here are a practical demonstration
      for method development and figure generation, not independent biological
      replicates. This limitation should remain visible in the final methods
      text.
    </p>
  </section>

  <section class="page-section">
    <h1>Step 1 - Start from Chapter 1 celltype annotation</h1>
    <p>
      The input is the Chapter 1 annotated Seurat object. The object is loaded
      before grouping or manual curation, so the active annotation column is
      <code>celltype</code>. This keeps numeric labels such as <code>2</code>,
      <code>4</code>, <code>5</code>, and <code>7</code> when no bibliography
      marker match was assigned.
    </p>
    {code_block("R", """
PIPELINE_DIR <- Sys.getenv("PIPELINE_DIR", unset = "/workspace/workflow")
DATA_DIR     <- Sys.getenv("DATA_DIR",     unset = "/workspace")
base_dir     <- file.path(DATA_DIR, "resultados_wt")

source(file.path(PIPELINE_DIR, "load_libraries.R"))
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))

set.seed(1807)
setwd(DATA_DIR)
list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)

pbmc_harmony <- readRDS(file.path(dir_objects,
                                  "pbmc_harmony_annotated.rds"))

pseudobulk_annot_col <- "celltype"
pseudobulk_conditions <- c("WT", "pifq")
comparison_tag <- "WT_vs_pifq"
""")}
  </section>

  <section class="page-section">
    <h1>Step 2 - Celltype and condition check</h1>
    <p>
      The first diagnostic counts cells per <code>celltype</code> and
      condition. Cell types with fewer than 20 cells in either WT or pifq are
      not used for DESeq2 in this draft. With the current data, 10 cell-type
      labels pass this filter.
    </p>
    {code_block("R", """
min_cells_per_condition <- 20

celltype_condition_counts <- as.data.frame.matrix(
  table(pbmc_harmony[[pseudobulk_annot_col]][, 1],
        pbmc_harmony$condition)
)
celltype_condition_counts$celltype <- rownames(celltype_condition_counts)
celltype_condition_counts <- celltype_condition_counts[
  , c("celltype", pseudobulk_conditions)
]

eligible_celltypes <- celltype_condition_counts$celltype[
  apply(celltype_condition_counts[, pseudobulk_conditions, drop = FALSE],
        1, function(x) all(x >= min_cells_per_condition))
]
""")}
    {img("cap2_celltype_condition_counts.png", "Cells per direct celltype annotation and condition.", "wide")}
  </section>

  <section class="page-section">
    <h1>Step 3 - Build celltype subsets</h1>
    <p>
      One Seurat subset is created per direct <code>celltype</code> label.
      Labels are only sanitized for object names and filenames; they are not
      grouped or manually curated.
    </p>
    {code_block("R", """
cell_type_subsets <- create_cell_type_subsets(
  pbmc_harmony,
  annot_col = pseudobulk_annot_col
)

eligible_subset_names <- gsub("[^[:alnum:]_]", "_", eligible_celltypes)
cell_type_subsets <- cell_type_subsets[
  names(cell_type_subsets) %in% eligible_subset_names
]
""")}
    <p class="small-note">
      The retained labels were: Epidermis Hypocotyl.1, Epidermis Cotyledon,
      Mesophyll, Epidermis Hypocotyl.2, PHLOEM, ENDO1, 2, 4, 5, and 7.
    </p>
  </section>

  <section class="page-section">
    <h1>Step 4 - Pseudo-replicate assignment</h1>
    <p>
      Cells inside each retained cell type are randomly assigned to three
      pseudo-replicates per condition. This creates replicate-level count
      columns for DESeq2. Because the dataset has one library per condition,
      these pseudo-replicates are used to exercise the workflow and should not
      be described as biological replicates.
    </p>
    {code_block("R", """
n_pseudoreps <- 3

cell_type_subsets_replicates <- assign_pseudoreplicates_batch(
  cell_type_subsets,
  pseudobulk_conditions = pseudobulk_conditions,
  n_pseudoreps = n_pseudoreps
)
""")}
    {img("cap2_pseudo_replicate_cell_counts.png", "Cells assigned to each pseudo-replicate for retained cell types.", "wide")}
  </section>

  <section class="page-section">
    <h1>Step 5 - Pseudobulk tables and DESeq2</h1>
    <p>
      Counts are aggregated by pseudo-replicate for each direct celltype label.
      DESeq2 is then run with <code>WT</code> as the reference and
      <code>pifq</code> as the treatment, so positive log2FC values indicate
      higher expression in pifq.
    </p>
    {code_block("R", """
comparisons <- list(
  list(conds = c("WT", "pifq"), tag = "WT_vs_pifq")
)

deseq2_results <- run_pseudobulk_deseq2_analysis(
  cell_type_subsets_replicates = cell_type_subsets_replicates,
  comparisons = comparisons,
  output_dir = dir_06,
  cell_types = NULL,
  pseudobulk_dir = file.path(dir_objects,
                             "pseudobulk_replicas_celltype")
)
""")}
    <p>
      This run produced one DESeq2 result file for each of the 10 retained
      direct celltype labels under
      <code>resultados_wt/06_de_results/WT_vs_pifq/</code>.
    </p>
  </section>

  <section class="page-section">
    <h1>Step 6 - Volcano plots</h1>
    <p>
      Volcano plots are rendered per cell type with adjusted p-value threshold
      0.05 and absolute log2FC threshold 1. The examples below show Mesophyll
      and Epidermis Hypocotyl.1.
    </p>
    {code_block("R", """
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
""")}
    <div class="two-col">
      {img("cap2_volcano_mesophyll_WT_vs_pifq.png", "Representative volcano plot: Mesophyll.")}
      {img("cap2_volcano_epidermis_hypocotyl_1_WT_vs_pifq.png", "Representative volcano plot: Epidermis Hypocotyl.1.")}
    </div>
  </section>

  <section class="page-section">
    <h1>Step 7 - Differential tables and log2FC heatmap</h1>
    <p>
      Per-cell-type DESeq2 outputs are combined into a differential-gene matrix
      and a log2FC matrix. The heatmap uses the same staircase logic as the
      helper function: genes are ordered by the cell type where the absolute
      log2FC is largest. The current WT/pifq draft includes 8,440 genes after
      filtering with adjusted p-value 0.05 and absolute log2FC 1.
    </p>
    {code_block("R", """
diff_tables <- build_differential_tables(
  results_dir = file.path(dir_06, volcano_tag),
  output_dir  = file.path(dir_06, volcano_tag),
  padj_cut    = padj_cut,
  lfc_cut     = lfc_cut,
  prefix      = paste0("diff_table_", volcano_tag)
)

build_logfc_heatmap(
  logfc_table  = diff_tables$logfc,
  contrast_tag = volcano_tag,
  output_dir   = file.path(dir_06, volcano_tag),
  limits       = c(-5, 5)
)
""")}
    {img("cap2_heatmap_WT_vs_pifq.png", "Combined log2FC heatmap across retained direct celltype labels.", "tall")}
  </section>

  <section class="page-section">
    <h1>Step 8 - GO enrichment</h1>
    <p>
      GO enrichment is run per retained cell type from significant genes in the
      WT vs pifq DESeq2 outputs. The Arabidopsis OrgDb annotation is used with
      TAIR gene identifiers and biological process ontology.
    </p>
    {code_block("R", """
go_results <- run_go_enrichment_for_contrast(
  results_dir  = file.path(dir_06, volcano_tag),
  output_dir   = file.path(dir_07, volcano_tag),
  orgdb        = org.At.tair.db,
  keytype      = "TAIR",
  go_space     = "BP",
  padj_cutoff  = 0.05,
  contrast_tag = volcano_tag
)
""")}
    {img("cap2_go_mesophyll_bp.png", "Representative GO biological-process bubble plot for Mesophyll.", "wide")}
  </section>

  <section class="page-section">
    <h1>Outputs and next decisions</h1>
    <p>
      The tested Chapter 2 script is
      <code>workflow/capitulo2_wt_pifq.R</code>. Main outputs are written under
      <code>resultados_wt/06_de_results/WT_vs_pifq/</code> and
      <code>resultados_wt/07_go/WT_vs_pifq/</code>. The next decision is
      whether the final methods paper should stop at DESeq2/GO for this dataset
      or continue into the optional hdWGCNA/network sections after the
      pseudo-replicate limitation is addressed in the text.
    </p>
  </section>
</body>
</html>
"""


if __name__ == "__main__":
    OUT.write_text(build(), encoding="utf-8")
    print(OUT)

