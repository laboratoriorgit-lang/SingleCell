# SingleCell — Arabidopsis scRNA-seq Pipeline

A complete single-cell RNA-seq pipeline for *Arabidopsis thaliana* (readily extensible to other organisms), from raw Cell Ranger processing through QC, integration, clustering, cell-type annotation, pseudobulk differential expression, GO enrichment, gene co-expression / regulatory network inference, and pseudotime trajectory analysis. The analytical logic lives in two documented function libraries (R + Python); the `workflow/stepN_*` scripts are thin, numbered drivers that call them in sequence and mirror the accompanying **Methods in Molecular Biology** protocol (`reports/arabidopsis_scRNAseq_pipeline.Rmd`/`.pdf`).

Worked example throughout: *Arabidopsis thaliana* root, **WT vs. *pifq*** mutant (public dataset CRA010863, Han et al. 2023, *Nature Plants*).

---

## Project Map

```text
┌──────────────────────────────────────────────────────────────────────────┐
│                 Raw data: FASTQ (CRA010863) · TAIR10.1 genome            │
│                          · Araport11 annotation                          │
└───────────────────────────────┬──────────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  Step 0 — Cell Ranger (bash)                 workflow/step0_cellranger.sh│
│                                                                          │
│  download FASTQ → cellranger mkref (TAIR10.1+Araport11)                 │
│    → cellranger count  →  filtered_feature_bc_matrix/ per sample        │
└───────────────────────────────┬──────────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────────┐
│  Part 1 — Single-cell analysis (R)  ·  Sections 1-12                    │
│                                          workflow/step1_singlecell.R    │
│                                                                          │
│   load + QC → filter → doublets → merge → Harmony integration           │
│     → resolution sweep → clustering → cell-type annotation              │
│     → subcluster curation → export (.rds + .h5ad)                       │
│                                                                          │
│   Helpers: ScRNA_Analysis_Functions.R                                   │
│   Output:  results/01_qc … 05_curation/, results/objects/               │
└───────────────────────────────┬──────────────────────────────────────────┘
                                 │  ath_sc_curated.rds / ath_sc_curated.h5ad
              ┌──────────────────┴───────────────────┐
              ▼                                       ▼
┌─────────────────────────────────────┐   ┌─────────────────────────────────┐
│ Part 2 — Pseudobulk DE / GO /        │   │ Part 3 — Pseudotime (Python)    │
│  co-expression networks (R)          │   │  Sections 23-29                 │
│  Sections 13-20                      │   │  workflow/step3_pseudotime.py   │
│  workflow/step2_degs.R               │   │                                 │
│                                       │   │  cell-type selection            │
│  pseudobulk per cell type            │   │    → Palantir diffusion maps    │
│    → DESeq2 (WT vs pifq)             │   │    → scFates PPT tree           │
│    → volcano / DE tables             │   │    → root + pseudotime          │
│    → GO enrichment (clusterProfiler) │   │    → gene trends (top-N +       │
│    → log2FC heatmap                  │   │      custom highlight)          │
│    → hdWGCNA / GENIE3 / WGCNA /      │   │                                 │
│      synergy networks → TF-DE network│   │  Helpers:                       │
│                                       │   │  ScRNA_Pseudotime_Functions.py  │
│  Helpers: ScRNA_Analysis_Functions.R │   │                                 │
│  Output: results/06_de_results …     │   │  Output: results/09_pseudotime/ │
│          08_networks/                │   │                                 │
└───────────────────────────────────────┘   └─────────────────────────────────┘
```

The **Rmd** in `reports/` is the literate, publication-ready version of the same three parts (plus Cell Ranger as Part 1) — it interleaves the identical code chunks from `step0`–`step3` with narrative text, figures and interpretation, and renders to `reports/arabidopsis_scRNAseq_pipeline.pdf`.

---

## Quick Start with Docker

`docker-compose.yml` is set up to pull the pre-built image by default:

```bash
docker compose pull      # fetch matigara/scrnaseq:latest from Docker Hub
docker compose build     # or build your own local version instead
                          #   (overwrites the same tag, this host only)
docker compose run --rm r          # open an R console with the mounted repo
docker compose run --rm r bash     # or a shell
```

`docker-compose.yml` mounts the current directory (`.`) into `/workspace` — the convention every `workflow/step*` script assumes.

**Manual, without docker-compose:**

```bash
docker pull matigara/scrnaseq:latest                     # from Docker Hub
# or: docker build -t matigara/scrnaseq:latest .          # build locally

docker run -it -v "$(pwd)":/workspace matigara/scrnaseq:latest bash
```

Once inside the container (`/workspace` = repo root):

```bash
R                                          # interactive R console
python3                                    # interactive Python console (venv on PATH)
Rscript workflow/step1_singlecell.R        # run Part 1 end-to-end
Rscript workflow/step2_degs.R              # run Part 2 end-to-end
python3 workflow/step3_pseudotime.py       # run Part 3 end-to-end
```

Cell Ranger itself is **not** bundled in the image (proprietary, license-gated download) — run `workflow/step0_cellranger.sh` in an environment where `cellranger` is on `PATH` before Part 1.

---

## Repository layout

```
SingleCell/
├── Dockerfile                          # rocker/r-ver:4.5 + full R/Python/Bioc stack
├── docker-compose.yml                  # pulls matigara/scrnaseq:latest, mounts repo at /workspace
├── .Rprofile                           # sets RETICULATE_PYTHON for the container's venv
├── README.md
├── data/
│   ├── Arabidopsis_thaliana_TAIR10.1_Araport11_genome.fasta.gz
│   ├── Arabidopsis_thaliana_TAIR10.1_Araport11_annotation.gtf.gz
│   ├── biblio_marks.txt                # cell-type marker table (cell.types <tab> gene)
│   └── AtTFDB_loci.txt                 # Arabidopsis transcription-factor loci, one per line
├── workflow/
│   ├── step0_cellranger.sh             # Part 0: FASTQ download + mkref + cellranger count
│   ├── step1_singlecell.R              # Part 1: Sections 1-12 (QC → curated object)
│   ├── step2_degs.R                    # Part 2: Sections 13-20 (pseudobulk DE → networks)
│   ├── step3_pseudotime.py             # Part 3: Sections 23-29 (pseudotime trajectory)
│   ├── load_libraries.R                # loads every R package the pipeline needs
│   ├── load_libraries_python.py        # loads every Python package Part 3 needs
│   ├── ScRNA_Analysis_Functions.R      # 75 R functions backing Parts 1-2 (documented below)
│   └── ScRNA_Pseudotime_Functions.py   # Python functions backing Part 3 (documented below)
└── reports/
    ├── arabidopsis_scRNAseq_pipeline.Rmd   # literate Methods-paper source (WT vs pifq)
    ├── arabidopsis_scRNAseq_pipeline.pdf   # rendered report
    ├── lang_box.lua / latex_header.tex     # Pandoc/LaTeX helpers for the PDF render
    └── assets/                             # figures embedded in the report
```

`workflow/step1_singlecell.R` and `step2_degs.R` both call `create_pipeline_dirs()` on `results/`, which creates:

```
results/
├── 01_qc/            02_clustering/     03_annotation/     04_expression/
├── 05_curation/       06_de_results/     07_go/             08_networks/
├── 09_pseudotime/     (created directly by step3_pseudotime.py)
└── objects/           # .rds / .h5ad checkpoints (ath_sc_curated.rds, ath_sc_curated.h5ad, ...)
```

`results/` itself is not tracked in git — it is generated fresh by running the pipeline against your own data.

---

## Pipeline walkthrough

### Part 0 — Cell Ranger (`workflow/step0_cellranger.sh`)

Reference script (not meant to run unattended): install Cell Ranger 9.0.1, download FASTQ for accession CRA010863 (WT = CRR775298, *pifq* = CRR775297) from the BIG Data Center, build the Cell Ranger reference from the bundled TAIR10.1 genome + Araport11 GTF (`cellranger mkref`), then run `cellranger count` per sample. Requires a signed, account-tied 10x download URL that cannot be hard-coded — copy it from the 10x Genomics downloads page.

### Part 1 — Single-cell analysis (`workflow/step1_singlecell.R`, Sections 1–12)

| Section | What happens |
|---|---|
| 1 | Load each CellRanger `filtered_feature_bc_matrix/` (`load_seurat_samples`), compute `percent.mt`/`percent.cp`, save pre-filter QC violins (`plot_qc_batch`) |
| 2 | Filter cells (`filter_seurat_samples`/`filter_sample`), save post-filter QC, checkpoint `.rds` |
| 3 | Merge samples, `NormalizeData` → `FindVariableFeatures` → `ScaleData` → `RunPCA` → `RunUMAP` (pre-Harmony) |
| 4 | `RunHarmony("orig.ident")` batch correction → post-Harmony UMAP |
| 5 | Elbow plot (k-means WSS over PCs) + `clustree` across 5 candidate resolutions |
| 6 | Final `FindNeighbors`/`FindClusters` (Leiden, `algorithm = 4`) at the chosen resolution |
| 7 | Marker-based annotation: `find_markers` (`FindAllMarkers`, TSV-cached) + `annotate_by_markers` against `data/biblio_marks.txt` + `plot_marker_dotplot` |
| 8 | Annotated `clustree` (cluster tree labeled by majority cell type) |
| 9 | Gene expression violin + feature plots for genes of interest |
| 10 | *[optional]* Merge duplicate cluster labels into `celltype_grouped` |
| 11 | *[optional]* Subcluster inspection and curation (`inspect_subcluster_markers`, `curate_clusters`) → `celltype_curated` |
| 12 | Export the curated object as `.rds` and, via `export_to_scanpy`, as `.h5ad` for Part 3 |

### Part 2 — Pseudobulk DE, GO enrichment, networks (`workflow/step2_degs.R`, Sections 13–20)

| Section | What happens |
|---|---|
| 13 | `create_cell_type_subsets` — one Seurat subset per curated cell type |
| 14 | `assign_pseudoreplicates_batch` — random pseudo-replicate groups per condition (no biological replicates needed) |
| 15 | `run_pseudobulk_deseq2_analysis` — pseudobulk count tables + `run_deseq2` per cell type, per comparison (WT vs *pifq*) |
| 16 | `render_volcano_plots` — per-cell-type volcano plots + combined PDF |
| 17 | `build_differential_tables` — combined classification (up/down/none) and log2FC tables across cell types |
| 18 | `run_go_enrichment_for_contrast` — `clusterProfiler::enrichGO` per cell type on significant genes |
| 19 | `build_logfc_heatmap` — staircase-ordered log2FC heatmap across cell types |
| 20 | `run_tf_coexpression_network` — hdWGCNA co-expression network on DE genes, filtered to transcription factors (`data/AtTFDB_loci.txt`), colored by up/down direction |

### Part 3 — Pseudotime trajectory (`workflow/step3_pseudotime.py`, Sections 23–29)

| Section | What happens |
|---|---|
| 23 | Setup: load `load_libraries_python.py` + `ScRNA_Pseudotime_Functions.py`, create `results/09_pseudotime/` |
| 24 | `load_curated_object` — read the Part 1 `.h5ad`, fix Seurat→scanpy `obsm` key names, plot overview UMAP |
| 25 | `preview_trajectory_selection` — subset to the cell type(s) that form the trajectory (e.g. `"Epidermis"`) |
| 26 | `run_trajectory_runs` → `build_pseudotime_trajectory` — Palantir diffusion maps → scFates Principal Polynomial Tree (PPT) → automatic root-cell selection → pseudotime assignment; supports sweeping multiple parameter sets in one call |
| 27 | Gene expression plotted directly on the force-directed trajectory graph |
| 28 | `run_step29_gene_trends` — `scFates.tl.test_association` + `.fit`, top-N and custom-gene trend plots along pseudotime |
| 29 | Export the final trajectory object as `results/objects/pbmc_pseudotime_final.h5ad` |

---

## Function reference

### 1. QC and Visualization

#### `plot_qc_batch(seurat_list, colors, file)`
Iterates a named list of Seurat objects, builds one QC violin panel per sample (`plot_qc_violin_grid`), stacks them and saves a single PDF via `save_qc`.

#### `plot_qc_violin_grid(obj1, label, color)`
Violin grid for `nFeature_RNA`, `nCount_RNA`, `percent.mt`, and `percent.cp` if present, for one Seurat object; title includes cell count.

#### `summarize_nfeature_plot(obj_list, labels = NULL, colores = NULL)`
Boxplot + jitter of `nFeature_RNA` across a list of Seurat objects, combined with printed quartile/quintile summary tables (`cowplot::plot_grid`). `labels` auto-generates `"Group1"`, `"Group2"`, … if `NULL`; `colores` defaults to a 5-color palette.

---

### 2. Doublet Detection and Sample Filtering

#### `preprocess_and_doubletfinder(seurat_obj, pcs = 1:20, expected_doublet_rate = 0.075, project_id = "sample")`
Minimal preprocessing (Normalize → FindVariableFeatures → Scale → PCA) then a single DoubletFinder parameter sweep (`paramSweep`/`find.pK`) to classify doublets. Adds `pANN_*`/`DF.classifications_*` metadata columns.

#### `doubletfinder_pipeline(obj, etiqueta = "Sample", PCs = 1:20, resolution = 0.5, return_singlets = TRUE, sct = FALSE)`
Full doublet pipeline: normalize, scale, PCA, `FindNeighbors`/`FindClusters`, `BCreal`-based parameter sweep, `doubletFinder`, adds `doublet_class`; optionally subsets to singlets.

#### `filter_seurat_samples(seurat_list, ...)`
Applies `filter_sample` to every element of a named Seurat list via `lapply`, preserving names; `...` is forwarded.

#### `filter_sample(obj, min_features = 200, max_features = Inf, min_counts = 0, max_counts = Inf, max_mt = 5, max_cp = 100, run_doubletfinder = TRUE)`
QC-threshold `subset()` on nFeature/nCount/%mt/%cp, then optionally runs `doubletfinder_pipeline` on the filtered cells.

---

### 3. Bulk / Pseudobulk Utilities

#### `normalize_bulk_pseudobulk(pseudobulk_counts, bulk_counts)`
DESeq2 size-factor normalization + log2(x+1) for a pseudobulk/bulk count-vector pair, restricted to their common genes (`stop()`s if fewer than 10 genes overlap).

#### `classify_residuals(df, umbral = 5)`
Fits `lm(bulk ~ pseudobulk)` and labels each gene `"Upregulated"` / `"Downregulated"` / `"Consistent"` by residual magnitude vs. `umbral`.

#### `generate_pseudobulk(seurat_obj, group_by = "orig.ident", merge_replicates = TRUE)`
Aggregates raw counts by a metadata grouping column (merging multi-layer count matrices if present). Returns a matrix, or `list(by_sample, by_condition)` if `merge_replicates = TRUE`.

#### `plot_replicate_correlation(pseudobulk_mat, main = "Replicate Correlation")`
Pairwise Pearson correlation across pseudobulk matrix columns, rendered as a `pheatmap`.

---

### 4. Seurat Utilities

#### `unify_names(obj)`
Strips numeric suffixes (`.1`, `_2`, …) from the active Idents' factor levels via `RenameIdents` — cleans up duplicate labels created by merging objects.

#### `show_annotation_table(filtered_vec, reference_vec, title = "Annotations")`
Side-by-side cell-type count comparison table (with totals row) drawn via `gridExtra::tableGrob`.

#### `export_to_scanpy(seurat_obj, outfile, assay_name = "RNA", use_reduc = c("pca","umap","harmony"), X_name = "logcounts", overwrite = TRUE)`
Converts Seurat → `SingleCellExperiment` (stripping non-serializable slots and duplicate metadata) and writes `.h5ad`, trying `anndataR` → `zellkonverter` → `SeuratDisk` in that order.

#### `safe_vln(obj, feature, colors)`
Thin `VlnPlot` wrapper grouped by `orig.ident` with a custom fill palette — safe for non-interactive RMarkdown rendering.

#### `join_layers_counts(obj, capas)`
Merges named RNA-assay count layers into one sparse matrix (`RowMergeSparseMatrices`); returns the single layer unchanged if only one is given.

---

### 5. Cell-Type Annotation

#### `find_markers(seurat_obj, output_file = "results/FindAllMarkers.tsv", only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0.25, force = FALSE)`
Runs `FindAllMarkers` after `JoinLayers`; caches to `output_file` and reloads from cache unless `force = TRUE`.

#### `annotate_by_markers(seurat_obj, markers, reference_file = NULL)`
Cross-references `find_markers()` output with a two-column reference table (col 1 = cell type, col 2 = gene, read **by position**, header ignored) to pick the best match per cluster by lowest `p_val_adj`; adds `.1`/`.2` suffixes when multiple clusters share a cell type. `reference_file = NULL` opens an interactive `file.choose()`. Sets metadata `celltype`.

#### `annotate_by_reference(seurat_obj, reference_obj = NULL, reference_col = NULL, dims = 1:30)`
Label transfer via `FindTransferAnchors`/`TransferData` (CCA) from a reference Seurat object. `reference_obj = NULL` opens `file.choose()`; `reference_col = NULL` prompts interactively (`readline`). Sets metadata `celltype_reference`.

---

### 6. Advanced Subclustering and Curation

#### `plot_subcluster_umap(obj, label, output_dir)`
Labeled `DimPlot` of a subclustered object (grouped by `cluster_subtipo`); saves `subcluster_<label>.pdf` (18×12 in, 300 dpi).

#### `subcluster_cell_type(obj, ctype, annot_col = "celltype_grouped", resolution = 0.3, dims = 1:20)`
Subsets to one or more cell types and recomputes PCA/UMAP/neighbors/clusters from scratch (`npcs` auto-capped to `min(ncells-1, 30)` for small subsets). Adds `cluster_subtipo`.

#### `plot_markers_for_subset(subset_obj, marker_table)`
`FeaturePlot` per marker gene present in the object (silently skips genes not found).

#### `save_subcluster_composite(subcluster_list, marker_table, output_dir, filename = "subclustering_inspection.pdf", n_marker_cols = 4)`
Multi-page PDF: page 1 = all subcluster UMAPs side by side; following pages = marker `FeaturePlot` grids, one page per cell type.

#### `inspect_subcluster_markers(obj, cluster_id, marker_table, annot_col = "celltype", resolution = 0.3, dims = 1:20, output_dir, filename = NULL, n_marker_cols = 6)`
Subclusters a single cluster/cell type, combines UMAP + marker dotplot + FeaturePlot grid into one `patchwork` figure, saves as PDF. Returns `list(plot, sub_obj)`.

#### `apply_subcluster_reassignment(obj, subcluster_list, reassign, source_col = "celltype_grouped", dest_col = "celltype_curated")`
Maps subcluster IDs from subclustered objects back onto the global object's cell-type labels, per a user-defined `reassign` table (`"others"` as catch-all).

#### `curate_clusters(obj, reassign, marker_table, output_dir, source_col = "celltype", dest_col = "celltype_curated", resolution = 0.3, dims = 1:20)`
Convenience wrapper: for every cluster named in `reassign`, subclusters + saves an inspection figure (`inspect_subcluster_markers`), then reassigns labels (`apply_subcluster_reassignment`) — all in one call.

---

### 7. Pseudobulk, DESeq2, Volcano, Heatmap (Differential Expression)

#### `assign_pseudo_replicates(obj, conditions = NULL, n_reps = 3, seed = 1807)`
Randomly assigns cells within each `condition` to `n_reps` pseudo-replicate groups. Returns `NULL` if fewer than 2 conditions are present.

#### `run_pseudobulk(obj)`
Sums counts per `replicate` group (handles multi-layer assays); returns a genes × replicate-group data.frame, columns alphabetically sorted. `stop()`s if `replicate` is missing.

#### `save_pseudobulk_tables(obj_list, output_dir, prefix = "Pseudobulk_Reps_", min_cells = 10, min_replicates = 2)`
Builds and writes a pseudobulk CSV per cell type in `obj_list` (skipping types below `min_cells`/`min_replicates`).

#### `run_deseq2(counts_mat, comparisons, output_dir, ctype = NULL)`
Builds a `DESeqDataSet`, auto-detects condition levels from column names (stripping `_repN`), runs DESeq2 per comparison (skipping any with <2 replicates per condition), writes result CSVs to `output_dir/<tag>/`.

#### `plot_volcano(file, padj_cut = 0.05, lfc_cut = 1)`
Reads a DESeq2 CSV, classifies genes by significance/direction, returns a colored volcano `ggplot` (does not save to disk itself).

#### `process_deseq2_result(file_path, output_dir, padj_cut = 0.05, lfc_cut = 1)`
Reads a DESeq2 CSV, classifies genes up(1)/down(-1)/unchanged(0), extracts log2FC for significant genes, writes a filtered CSV. Returns `list(class, logfc)`.

#### `render_volcano_plots(results_dir, output_dir, pdf_name = "VolcanoPlots.pdf", padj_cut = 0.05, lfc_cut = 1)`
Reads every `DESeq2_*.csv` in `results_dir`, saves one PNG per result plus a combined multi-panel PDF.

#### `build_differential_tables(results_dir, output_dir, padj_cut = 0.05, lfc_cut = 1, prefix = "diff_table")`
Processes every DESeq2 CSV in a contrast directory (`process_deseq2_result`), assembles combined classification and log2FC matrices across cell types (`full_join` on `gene_id`), filters to genes with ≥1 nonzero change, writes combined TSVs.

#### `plot_heatmap(mat, min_genes = 1, deepSplit_val = 0, breaks = c(-5, 5))`
Rows (genes) clustered by Euclidean distance + `cutreeDynamic`; columns (conditions) by PCA-based distance (components explaining ≥90% variance); renders a `pheatmap` with per-gene cluster color annotation.

#### `plot_marker_dotplot(seurat_obj, marks, annot_col = "celltype_grouped", cell_order = NULL, clusters_remove = NULL, rename_map = NULL, outfile = NULL, width = 18, height = 18, dot_scale = 12, base_size = 18)`
Coordinate-flipped `DotPlot` in "staircase" order — marker genes grouped under the cell type they mark, cell types ordered to match — so significant dots fall along a near-diagonal. Supports renaming, exclusion, and manual `cell_order`.

---

### 8. GO Enrichment

#### `compute_go_enrichment(tbl, universo, go_space, orgdb = org.At.tair.db, keytype = "TAIR", qvalueCutoff = 0.05, pvalueCutoff = 0.05, simplificar = FALSE, umbral_simply = 0.7, output_dir = "results/Enrichment")`
Iterates columns of a binary genes × comparisons matrix, runs `clusterProfiler::enrichGO` per set of genes marked `1`, writes raw + symbol-readable (`setReadable`) result tables, optionally simplifies redundant terms.

#### `prune_go(go_result, nivel, go_space, qvalueCutoff, simplificar = FALSE, output_dir = "results/Enrichment")`
Applies `gofilter` to each `enrichResult`, keeping only terms at or below a given GO level; writes filtered tables.

#### `plot_go_bubbles(go_result)`
Balloon/bubble chart (`ggpubr::ggballoonplot`) — bubble size = fold enrichment, fill = -log10(q-value). `stop()`s if no results are available.

#### `run_go_enrichment_suite(diff_table, output_dir, orgdb, keytype, go_space = "BP", qvalue_cutoff = 0.05, pvalue_cutoff = 0.05, simplify_cutoff = 0.7, go_level = 6, pdf_name = "GO_enrichment.pdf")`
Auto-derives the background gene universe from the OrgDb/keytype pair, runs full + simplified `enrichGO`, optional level-pruning, and exports a multi-page balloon-plot PDF. Returns `list(tbl, universo, total, simple, total_podado, simple_podado)`.

#### `run_go_for_gene_clusters(assignments, cluster_col, output_dir, orgdb, keytype, go_space = "BP", qvalue_cutoff = 0.05, pvalue_cutoff = 0.05, simplify_cutoff = 0.7, go_level = 6, pdf_name = "GO_clusters.pdf")`
Converts a gene→cluster/module assignment table into a binary membership matrix, then runs `run_go_enrichment_suite` on it.

#### `run_simple_go_enrichment(diff_table, output_dir, orgdb, keytype = "TAIR", go_space = "BP", padj_cutoff = 0.05, cell_type = NULL, contrast_tag = NULL)`
Single-pass `enrichGO` (no up/down split, `readable = TRUE`) on all genes in a diff table's first column; writes a TSV and a dotplot PDF.

#### `run_go_enrichment_for_contrast(results_dir, output_dir, orgdb, keytype = "TAIR", go_space = "BP", padj_cutoff = 0.05, contrast_tag)`
Finds every DESeq2 CSV for a contrast, filters genes by `padj_cutoff` per cell type, runs `run_simple_go_enrichment` for each cell type with significant genes.

---

### 9. Co-expression / TOM Modules (generic WGCNA, no hdWGCNA)

#### `run_coexpression_cluster_suite(diff_table, output_dir, selected_cols = NULL, min_genes = 1, deepSplit_val = 0, breaks = c(-5, 5), network_power = 6, network_type = c("signed","unsigned"), cor_method = "spearman", go_orgdb, go_keytype, go_space = "BP", go_qvalue_cutoff = 0.05, go_pvalue_cutoff = 0.05, go_simplify_cutoff = 0.7, go_level = 6, heatmap_pdf = "coexpression_heatmap.pdf", tom_pdf = "tom_heatmap.pdf", go_pdf = "GO_clusters.pdf")`
All-in-one: builds a log2FC heatmap with dynamic tree cut, computes a rank-based co-expression network (correlation → adjacency → TOM via `WGCNA::TOMsimilarity`), clusters TOM modules, and runs GO enrichment per positive module.

#### `prepare_coexpression_matrix(diff_table, selected_cols = NULL)`
Loads a combined log2FC table, selects columns, returns a numeric matrix with gene IDs as rownames.

#### `build_heatmap_clusters(Mz, output_dir, min_genes = 1, deepSplit_val = 0, breaks = c(-5, 5), heatmap_pdf = "coexpression_heatmap.pdf")`
Hierarchical heatmap of a log2FC matrix (same logic as `plot_heatmap`, self-contained with PDF/TSV export). Returns cluster assignments.

#### `build_coexpression_modules(Mz, output_dir, min_genes = 1, deepSplit_val = 0, network_power = 6, network_type = c("signed","unsigned"), cor_method = "spearman", tom_pdf = "tom_heatmap.pdf")`
Rank-based gene-gene correlation → adjacency → TOM → TOM modules (`cutreeDynamic`) from a log2FC matrix; exports TOM heatmap + module assignment table.
> **Known bug:** the returned `tom_heatmap` element references an undefined `tom_plot` object — the `pheatmap()` call inside the `pdf()`/`dev.off()` block is never captured with `<-`.

#### `run_go_for_gene_clusters(...)` — see section 8 (used here to enrich the modules produced above).

---

### 10. hdWGCNA (Hierarchical Co-expression Networks)

#### `run_unified_hdwgcna(seurat_obj, de_table_path, output_dir, annot_col, sample_col = "orig.ident", wgcna_name = "unified", n_metacells = 25, soft_power = NULL, min_module_size = 30, deep_split = 2, make_plots = TRUE)`
Builds one hdWGCNA network restricted to genes present in a log2FC/DE table. Full pipeline: `SetupForWGCNA` → `MetacellsByGroups` → `NormalizeMetacells` → `SetDatExpr` → `TestSoftPowers` (auto soft power at `SFT.R.sq ≥ 0.8`) → `ConstructNetwork` → `ModuleEigengenes` → `ModuleConnectivity`. Saves modules, hub genes, summary, and the hdWGCNA object; optionally 3 plots (module UMAP, eigengene heatmap, dendrogram).

#### `run_hdwgcna(seurat_obj, annot_col = "celltype_reference", output_dir, cell_types = NULL, n_metacells = 25, soft_power = NULL, max_modules = 8, gene_select = "fraction", fraction = 0.05, de_genes_per_ct = NULL)`
Same pipeline run **per cell type**. Trims to the `max_modules` largest modules per type (rest relabeled `"grey"`); supports resume (skips a type if its `.rds` already exists) and restriction to per-cell-type DE gene lists.

#### `plot_hdwgcna_network(hdwgcna_dir, output_dir = hdwgcna_dir, tom_threshold = 0.1, cell_types = NULL, n_hub_label = 5, max_modules = NULL, make_plot = TRUE)`
Reads saved hdWGCNA `.rds` objects, pulls the TOM matrix directly from the on-disk `.rda` files (bypasses `GetTOM`), filters edges by `tom_threshold`, exports edge/node lists and (optionally) a `ggraph` network plot per cell type.

#### `plot_hdwgcna_network_tf(network_dir, tf_list_file, wgcna_name = "unified", output_dir = network_dir, n_hub_label = 10)`
Restricts an exported hdWGCNA network to TF–TF and TF–target edges (drops target–target) using a flat TF locus list; plots TFs as triangles, targets as circles.

#### `filter_hdwgcna_by_de(hdwgcna_dir, de_dirs, output_dir = hdwgcna_dir, padj_cut = 0.05, lfc_cut = 1, n_hub_label = 10, tom_threshold = 0.1, max_modules = NULL)`
Filters an hdWGCNA network to DE genes only (union across all `de_dirs` contrasts, keeping the largest |log2FC| record per gene); exports DE-filtered edges/nodes, a network plot (node size = |log2FC|), and an eigengene heatmap restricted to modules with ≥3 DE genes.

---

### 11. GENIE3 / WGCNA per Cluster and Synergy Networks

#### `get_tfs_from_orgdb(orgdb, keytype = "TAIR")`
Extracts transcription-factor gene IDs from any Bioconductor OrgDb via GO:0003700 ("DNA-binding transcription factor activity") — organism-agnostic.

#### `load_pseudobulk_matrix(pseudobulk_dir, normalize = TRUE)`
Loads and merges pseudobulk replicate CSVs (`Pseudobulk_Reps_<celltype>.csv`) from a directory; optional CPM + log2 normalization.

#### `run_genie3_per_cluster(cluster_assignments, pseudobulk_dir, output_dir, orgdb, keytype = "TAIR", custom_tfs = NULL, n_top_clusters = 3, cor_min = 0.90, genie3_ntrees = 100, n_cores = 4, min_var_filter = 0.01)`
Directed (TF → target) regulatory network inference via `GENIE3`, per gene cluster; filters edges by minimum absolute Pearson correlation.

#### `run_wgcna_per_cluster(cluster_assignments, pseudobulk_dir, output_dir, n_top_clusters = 3, soft_power = 6, network_type = "signed", tom_threshold = 0.10, min_var_filter = 0.01)`
Undirected TOM-based co-expression network (`WGCNA::adjacency` → `TOMsimilarity`) per gene cluster, filtered by TOM threshold. Clusters with <10 genes after variance filtering are skipped.

#### `run_synergistic_network(cluster_assignments, pseudobulk_dir, output_dir, orgdb, keytype = "TAIR", custom_tfs = NULL, n_top_clusters = 3, soft_power = 6, network_type = "signed", genie3_ntrees = 100, n_cores = 4, min_var_filter = 0.01, cor_min = 0.90, tom_min = 0.15)`
Combines GENIE3 directionality with WGCNA TOM robustness via geometric mean of rank-normalized scores: `score_synergy = sqrt(rank(weight_genie3) × rank(TOM))`. Keeps edges only if they pass **both** the Pearson and TOM thresholds (not a simple intersection).

#### `generate_network_pdf(results, output_dir, method_name, weight_col, directed = FALSE, n_top_hubs = 10, max_nodes = 80, edge_color = "#1f78b4")`
Reusable PDF report for any of the three methods above: cover page, cross-cluster summary table, then per cluster: network plot (`igraph`), top-hub bar chart, top-15-edges table + interpretation text.

#### `visualize_network_per_cluster(network_results, cluster_assignments, output_dir, method_name = "GENIE3")`
Force-directed (Fruchterman–Reingold) network plot per cluster with weight normalization, node color/size by degree.

#### `generate_cluster_profile_report(cluster_assignments, pseudobulk_matrix, output_dir, method_name = "MIXED")`
Scaled-expression (z-score, clipped to [-3, 3]) heatmap via `ComplexHeatmap` plus per-cluster summary statistics (mean/sd/min/max).

#### `run_network_inference_pipeline(heatmap_results, pseudobulk_dir, output_base_dir, methods = c("GENIE3","WGCNA","SYNERGY"), orgdb = org.At.tair.db, keytype = "TAIR", custom_tfs = NULL, cor_min = 0.90, genie3_ntrees = 100, n_cores = 4, soft_power = 6, network_type = "signed", tom_threshold = 0.15, n_top_clusters = 3, min_var_filter = 0.01)`
Orchestrator: runs the selected methods and generates reports for each via `generate_network_pdf`.
> **Known fragility:** the GENIE3 branch calls `generate_network_pdf(results[[1]], ...)` instead of `results$GENIE3` — correct only while GENIE3 runs first in `methods`.

#### `test_network_thresholds(heatmap_results, pseudobulk_dir, output_dir, method = "SYNERGY", orgdb = org.At.tair.db, keytype = "TAIR", custom_tfs = NULL, genie3_ntrees = 100, n_cores = 4, soft_power = 6, network_type = "signed", n_top_clusters = 3, min_var_filter = 0.01)`
Sweeps 5 predefined threshold combinations (Strict/Moderate/Exploratory/Permissive/Very Lax), summarizes cluster/edge coverage per combination, generates a comparison PDF with a recommendation (highest cluster coverage).

---

### 12. Transcription-Factor / Differential-Expression Networks

#### `build_tf_network(edges, tfs, de_mat)`
Filters an edge list to pairs involving ≥1 TF, classifies each gene's direction ("up"/"down"/"mixed") from a multi-contrast log2FC matrix, builds an undirected `igraph` object with node attributes.

#### `plot_tf_de_network(net, output_dir, layout = "stress", n_hub_label = 15, contrast_tag = "condition_1_vs_condition_2", output_pdf = "network_tf_DE_direction.pdf", output_width = 12, output_height = 10)`
Draws the TF graph via `ggraph`, colored by DE direction (red = up, blue = down, grey = mixed), shaped by role (triangle = TF, circle = co-expression partner), labels the top hub TFs. Supports `"fr"`, `"lgl"`, or any native `ggraph` layout (e.g. `"stress"`, `"graphopt"`).

#### `run_tf_coexpression_network(seurat_obj, de_table_path, tf_list_path, output_dir, annot_col, contrast_tag, sample_col = "orig.ident", n_metacells = 50, min_module_size = 20, deep_split = 2, soft_power = NULL, tom_threshold = 0.2, n_hub_label = 15)`
Single-call wrapper chaining: `run_unified_hdwgcna` (gene-gene network on significant DE genes) → `plot_hdwgcna_network` (export edges only) → filter to TF-containing edges → `build_tf_network` + `plot_tf_de_network`. This is what `step2_degs.R` Section 20 calls.

---

### 13. Chapter-Level Pipelines

#### `run_pseudobulk_pipeline(obj, comparisons, orgdb, keytype, annot_col = "celltype_grouped", n_reps = 3, padj_cut = 0.05, lfc_cut = 1, go_space = "BP", qval = 0.05, nivel_poda = 6, dir_pseudobulk, dir_deseq2, dir_volcano, dir_heatmaps, dir_go)`
Runs Sections 14–20 in one call: replicate correlation → per-cell-type subsets + pseudo-replicates → DESeq2 → volcano → DE tables + heatmaps → GO enrichment (full + simplified, raw + level-pruned).

#### `run_pseudobulk_deseq2_analysis(cell_type_subsets_replicates, comparisons, output_dir, cell_types = NULL, pseudobulk_dir = NULL)`
Intermediate orchestrator: writes pseudobulk count tables (`save_pseudobulk_tables`), filters to requested cell types, runs `run_deseq2` for each with progress separators. This is what `step2_degs.R` Section 15 calls.

#### `run_network_inference_pipeline(...)` — see section 11 (GENIE3/WGCNA/SYNERGY orchestrator).

#### `create_pipeline_dirs(base_dir)`
Creates the standard output folder structure under `base_dir` — `01_qc` … `08_networks` + `objects` — and returns a named list of paths (`list2env()`'d into the global environment by `step1_singlecell.R`/`step2_degs.R`).

---

### 14. Data Loading and Subsetting

#### `load_seurat_samples(samples, DATA_DIR, mt_pattern, cp_pattern)`
Loads multiple CellRanger samples (`Read10X`), creates Seurat objects (`min.cells = 3`, `min.features = 200`), tags condition, computes `percent.mt`/`percent.cp` via `PercentageFeatureSet`.

#### `create_cell_type_subsets(seurat_obj, annot_col = "celltype_grouped")`
One Seurat subset per unique (non-NA) cell type in `annot_col`, with sanitized list names.

#### `assign_pseudoreplicates_batch(cell_type_subsets, pseudobulk_conditions = NULL, n_pseudoreps = 3)`
Batch-applies `assign_pseudo_replicates` across cell types, dropping any that return `NULL` (<2 conditions present). **Requires a global `set.seed()` before calling** for reproducibility.

---

### 15. DE-Specific log2FC Heatmap

#### `build_logfc_heatmap(logfc_table, contrast_tag, output_dir, limits = c(-5, 5), cell_order = NULL)`
`ComplexHeatmap`-based log2FC heatmap in staircase gene order (each gene grouped under the cell type where it peaks in |log2FC|); columns follow `cell_order` (or natural order) instead of a column dendrogram, so up/down blocks fall on the diagonal.

---

### 16. Save Helpers

Three lightweight `ggsave` wrappers. **All three read `output_dir` from the calling global environment — it is not a function parameter.**

#### `save_pdf(plot, file, w = 18, h = 18)`
Saves any ggplot at 300 dpi, `limitsize = FALSE`.

#### `save_vln(plot, file, n = 1)`
One-line wrapper around `save_pdf` (fixed 18×18). `n` is documented but unused in the body.

#### `save_qc(plot_list, file)`
Stacks a list of plots into one column (`patchwork::wrap_plots(ncol = 1)`), saves at 300 dpi with a white background.

---

## Python function reference (`ScRNA_Pseudotime_Functions.py`)

| Function | Purpose |
|---|---|
| `save_plot_18x18(fig, path)` | Saves a matplotlib figure at 18×18 in, 200 dpi, white background |
| `ensure_expression_in_X(adata, layer="logcounts")` | Copies `adata.layers[layer]` into `.X` if `.X` is empty — required before any scanpy gene-color plot |
| `load_curated_object(input_h5ad, dir_pseudotime, annotation_col, n_jobs=4)` | Step 25: loads the Part 1 `.h5ad`, remaps Seurat-style `obsm` keys (`UMAP`→`X_umap`, `PCA`→`X_pca`, `HARMONY`→`X_harmony`), plots overview UMAP with a side legend, prints per-cell-type counts |
| `preview_trajectory_selection(adata, clusters, annotation_col, dir_pseudotime)` | Step 26: subsets to selected cell types, plots a preview UMAP; `raise`s if a requested cluster name doesn't exist |
| `trajectory_run(nodes=50, sigma=0.2, lambda_value=60, eigs=20, seed=3)` | Builds a named parameter-set dict for `run_trajectory_runs`; the name encodes every parameter so identical configs share a folder |
| `build_pseudotime_trajectory(adata, clusters, root_cluster, annotation_col, nodes=150, sigma=0.2, ppt_lambda=60, n_components=50, n_eigs=20, n_neighbors=50, seed=3)` | Core Step 27 pipeline: subset → Palantir diffusion maps → neighbors graph → `scFates.tl.tree` (PPT) → root selection (most-connected cell in `root_cluster`) → `scFates.tl.pseudotime` |
| `build_dendrogram(adata)` | Adds a dendrogram to a tree-fitted `adata` (`scFates.tl.dendrogram`) |
| `draw_graph_basis(adata)` | Returns whichever `X_draw_graph_{fa,fr}` embedding scanpy actually computed (`fa` = ForceAtlas2 if `fa2-modified` is installed, else falls back to `fr`) |
| `plot_trajectory_graphs(adata, name, output_dir, annotation_col, show_inline=False, param_label=None)` | Saves the force-directed tree graph colored by cell type |
| `export_pseudotime_table(adata, name, output_dir, annotation_col=None)` | Writes a per-cell pseudotime TSV (`t` + optional annotation column) |
| `plot_pseudotime_trajectory(adata, name, output_dir, annotation_col, show_inline=False, param_label=None)` | Overlays the principal graph + pseudotime coloring on the force-directed layout |
| `plot_root_cell(adata, name, output_dir, show_inline=False, param_label=None)` | QC plot: highlights the auto-selected root cell in red on the FA1/FA2 embedding |
| `run_trajectory_runs(adata, clusters, root_cluster, annotation_col, output_base_dir, runs, selected_run=None, show_inline=True)` | Step 27 entry point: runs one or more parameter sets from `trajectory_run()`, saves each to its own folder (h5ad, plots, pseudotime table), returns the selected run's `adata` + folder + all runs |
| `strip_segment_frames(ax)` | Removes stray unfilled rectangle artifacts scFates leaves behind on trend plots |
| `save_trends_plot(axes_list, pdf_path)` | Saves gene-trend axes at 24×18 in (wider than the standard 18×18 — the dendrogram panel needs the extra width) |
| `run_step29_gene_trends(adata, run_dir, name="trajectory", custom_genes=None, top_n=50, ordering="max", n_jobs=4, spline_df=5, A_cut=0.3, p_val_cut=0.001)` | Step 28/29: `scFates.tl.test_association` + `.fit`, then two trend heatmaps — top-N most variable genes and a custom highlighted gene list (force-included even if not significant) |

`load_libraries_python.py` additionally auto-detects Jupyter vs. script execution (`IPKernelApp` in the IPython shell config) to pick the inline vs. `Agg` matplotlib backend, and silences warnings globally.

---

## Data files (`data/`)

| File | Format | Used by |
|---|---|---|
| `Arabidopsis_thaliana_TAIR10.1_Araport11_genome.fasta.gz` | FASTA, gzipped | `cellranger mkref` (Part 0) |
| `Arabidopsis_thaliana_TAIR10.1_Araport11_annotation.gtf.gz` | GTF, gzipped | `cellranger mkref` (Part 0) |
| `biblio_marks.txt` | TSV, columns `cell types`/`gene` (read **by position**, not by header name) | `annotate_by_markers`, `plot_marker_dotplot`, `build_logfc_heatmap` (Parts 1–2) |
| `AtTFDB_loci.txt` | plain text, one AGI locus per line | `run_tf_coexpression_network` / `plot_hdwgcna_network_tf` (Part 2, Section 20) |

---

## Versions and dependencies

| Component | Version / source |
|---|---|
| Base image | `rocker/r-ver:4.5` (R 4.5) |
| Cell Ranger | 9.0.1 (not bundled — proprietary, licensed download) |
| Reference genome | TAIR10.1 + Araport11 annotation (bundled, gzipped, under `data/`) |
| CRAN / Bioconductor packages | **not version-pinned** in the Dockerfile — installed at build time from `packagemanager.posit.co` (latest for R 4.5) and Bioconductor's release matching R 4.5 |
| Python packages | **not version-pinned** — `pip install` inside a venv at build time (`scanpy`, `scFates`, `palantir`, `fa2-modified`, `rpy2`, `pandas`, `numpy`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`) |
| GitHub-only R packages | `SeuratDisk` (mojaveazure), `DoubletFinder` (chris-mcginnis-ucsf), `SeuratWrappers` (satijalab, needs `GITHUB_PAT` build arg), `monocle3` (cole-trapnell-lab) |
| Seurat ecosystem | `Seurat`, `SeuratObject` (CRAN, latest), `Signac` |
| Key Bioconductor packages | `DESeq2`, `clusterProfiler`, `org.At.tair.db`, `GENIE3`, `ComplexHeatmap`, `zellkonverter`, `basilisk`, `SingleCellExperiment`, `scuttle`, `scater` |
| Other clustering / network packages | `WGCNA` (+ `impute`, `preprocessCore` from Bioconductor), `dynamicTreeCut`, `igraph`, `ggraph`, `tidygraph`, `leidenbase` (Leiden clustering, `algorithm = 4`) |
| Docker image (published) | `matigara/scrnaseq:latest` on Docker Hub |

> **Known gap:** `workflow/load_libraries.R` calls `library(hdWGCNA)`, and several R functions (`run_unified_hdwgcna`, `run_hdwgcna`, `plot_hdwgcna_network*`, `filter_hdwgcna_by_de`, `run_tf_coexpression_network`) depend on it, but the Dockerfile never installs `hdWGCNA` (typically `remotes::install_github("smorabit/hdWGCNA")`). Install it manually before running Part 2, Section 20.

> `.Rprofile` sets `RETICULATE_PYTHON` to a hard-coded path (`/home/mvergara/projects3/app/miniconda/envs/scrna_seba/bin/python`) from the original development host — this only matters if you call Python from within an R session via `reticulate`; it does not affect running `workflow/step3_pseudotime.py` directly with the container's own `python3`. Edit it if you use `reticulate` from R on a different machine.

---

## Adapting to other organisms

Swap the organism-specific parameters below when applying this pipeline to species other than *Arabidopsis thaliana*. `mt_pattern`/`cp_pattern` go to `load_seurat_samples()`; `orgdb`/`keytype` go to every GO/network function (`run_go_enrichment_for_contrast`, `run_go_enrichment_suite`, `run_network_inference_pipeline`, `test_network_thresholds`, `get_tfs_from_orgdb`, …).

| Organism | `mt_pattern` | `cp_pattern` | `orgdb` | `keytype` |
|---|---|---|---|---|
| *Arabidopsis thaliana* | `"^ATMG"` | `"^ATCG"` | `org.At.tair.db` | `"TAIR"` |
| *Homo sapiens* | `"^MT-"` | `NULL` | `org.Hs.eg.db` | `"ENSEMBL"` |
| *Mus musculus* | `"^mt-"` | `NULL` | `org.Mm.eg.db` | `"ENSEMBL"` |
| *Oryza sativa* | organism-specific | `NULL` | `org.Os.eg.db` | `"GID"` |

For any new organism: install the matching `org.*.eg.db` package, update `mt_pattern`/`cp_pattern` to your genome annotation's organelle gene-naming convention, replace `data/biblio_marks.txt` with your own marker table, and update the `universo`/background gene vector used by the GO enrichment functions to reflect your target organism's full gene set. `get_tfs_from_orgdb()` works unmodified for any OrgDb with GO annotations; for the TF-network functions (`run_tf_coexpression_network`, `plot_hdwgcna_network_tf`), replace `data/AtTFDB_loci.txt` with an equivalent flat TF locus list for your organism.
