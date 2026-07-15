# Methods Paper Working Memory

Date: 2026-07-15

This file is the working memory for the methods-paper draft. It records the
current decisions, file layout, and next steps so Chapter 2 can continue with
the same style and assumptions used for Chapter 1.

## Current Scope

The active methods-paper draft is focused on a tutorial-style workflow for the
WT vs pifq Chapter 1 analysis. It starts from repository cloning, raw data
download, reference staging, Cell Ranger count generation, Docker startup, and
then proceeds through the R/Seurat workflow section by section.

The current rendered files are:

- `methods_paper_chapter1_wt_pifq.pdf`
- `methods_paper_chapter1_wt_pifq.html`

The active server copy is stored under:

- `/home/mvergara/projects2/Sc_DB_test/resultados_wt/`

The user-facing MacBook copy is stored on the Desktop.

## Style Rules Established for the PDF

Keep the document as a practical tutorial, not a dense report.

- Start from `git clone`, then use `PROJECT_DIR="$(pwd)"`.
- Keep commands copy/paste friendly.
- Avoid automated wrapper scripts for the main tutorial flow.
- Do not manually create `resultados_wt/`; the R pipeline creates it.
- Use small reference images, not full-size report panels.
- Start each major step on a new page.
- Avoid large tables in the PDF body.
- Use prose or compact lists instead of tables when the information is simple.
- Keep software versions at the end as a compact appendix.
- Do not include base R package tables.
- Keep code blocks visually clean and avoid unnecessary comments in rendered
  code blocks when the surrounding text already explains the purpose.

## Chapter 1 Workflow State

Chapter 1 now documents the following sequence:

1. Clone the `SingleCell` repository.
2. Define `PROJECT_DIR` as the cloned repository root.
3. Download raw CRA010863 FASTQ files for WT and pifq.
4. Download the hosted Arabidopsis Cell Ranger reference from VirtualPlant.
5. Run Cell Ranger outside Docker using the repository launcher.
6. Start the Docker image `matigara/scrnaseq:latest`.
7. Run R sections manually and inspect outputs after each step.

The R sections currently documented are:

- Section 1: data loading and pre-filter QC.
- Section 2: filtering and doublet detection.
- Section 3: merge and initial preprocessing.
- Section 4: Harmony batch correction.
- Section 5: resolution optimization with elbow plot and clustree.
- Section 6: final clustering.
- Section 7: bibliography marker annotation.
- Section 8: annotated clustree.
- Section 9: simple gene expression visualization.
- Section 10: optional grouping, default pass-through.
- Section 11: optional interactive cell-type inspection.
- Section 12: export to h5ad.

## Important Chapter 1 Decisions

Bibliography annotation is the active annotation path. Reference-transfer
annotation was intentionally left out.

The marker table controls the marker dotplot order:

- Y axis follows the marker-file row order.
- X axis follows the matching cell-type order implied by the marker file.
- Numeric suffixes such as `Epidermis Hypocotyl.1` and
  `Epidermis Hypocotyl.2` stay adjacent to their parent marker category.

The grouping section stays inactive by default:

```r
# grouping <- c(
#   "Epidermis Hypocotyl.1" = "Epidermis Hypocotyl",
#   "Epidermis Hypocotyl.2" = "Epidermis Hypocotyl"
# )

grouping <- c()
```

The Cell Ranger executable itself is not bundled in the Docker image. The repo
contains a launcher in:

- `tools/cellranger/cellranger`

Users can place an official Cell Ranger installation under `tools/cellranger/`
or rely on the server environment, then use:

```bash
export PATH="$PROJECT_DIR/tools/cellranger:$PATH"
cellranger --version
```

## Section 11 Interactive Inspection

Section 11 is intentionally not an automated curation step.

It now works as an interactive inspection example:

- The tested example is `Epidermis Hypocotyl.1`.
- The function calculates subclusters inside that subset.
- The visual UMAP uses the original global UMAP coordinates, subsetted to the
  inspected cells, so the global geometry is preserved.
- The function saves a bibliography-marker dotplot across subclusters.
- The function saves a small FeaturePlot panel for the bibliography markers.
- It does not modify `pbmc_harmony`.
- It does not create or overwrite `celltype_curated`.
- Manual reassignment should only be added later if the user decides it is
  biologically justified after inspecting the figures.

The helper added for this is:

- `inspect_subcluster_markers()` in `workflow/ScRNA_Analysis_Functions.R`

The relevant Git commits are:

- `ebbe0b8` - `Simplify interactive subcluster inspection`
- `cd3dffa` - `Plot subcluster inspection on original UMAP`

## Chapter 1 Output Assets Used in the PDF

The methods PDF uses small PNG versions of the analysis outputs:

- `qc_prefilter.png`
- `qc_postfilter.png`
- `umap_preharmony.png`
- `umap_postharmony.png`
- `elbow_plot.png`
- `clustree2.png`
- `umap_seuratclusters.png`
- `umap_annotation_biblio.png`
- `dotplot_marker_table_annotation_biblio.png`
- `clustree_annotated.png`
- `feature_gene_all.png`
- `vln_gene_all.png`
- `umap_grouped.png`
- `umap_original_subclusters_Epidermis_Hypocotyl_1.png`
- `dotplot_bibliomarks_subclusters_Epidermis_Hypocotyl_1.png`
- `featureplots_bibliomarks_Epidermis_Hypocotyl_1.png`

The source server assets are under:

- `/home/mvergara/projects2/Sc_DB_test/resultados_wt/report_assets/`

## Chapter 2 Direction

Chapter 2 should continue from the Chapter 1 bibliography-annotated object.
For the WT/pifq draft, use the direct `celltype` column and do not use
`celltype_grouped` or `celltype_curated`.

The existing pipeline entry point is:

- `workflow/capitulo2_pseudobulk_de.R`

Before drafting the Chapter 2 tutorial, adapt it to the WT/pifq layout:

- Use `base_dir <- file.path(DATA_DIR, "resultados_wt")`.
- Use the object produced before grouping/curation:
  `resultados_wt/objects/pbmc_harmony_annotated.rds`.
- Remove `custom_seurat.R` unless a tested function truly requires it.
- Set the comparison to WT vs pifq using the `condition` column.
- Confirm that the sample and condition metadata are suitable for pseudobulk.
- Decide whether pseudo-replicates are acceptable for this dataset or whether
  the text should clearly state this as a demonstration when biological
  replicates are limited.

Chapter 2 should preserve the base section order from
`workflow/capitulo2_pseudobulk_de.R`:

- Load annotated object from Part 1
- Section 13 - Cell-type subsets
- Section 14 - Pseudo-replicate assignment
- Section 15 - Pseudobulk tables and DESeq2
- Section 16 - Volcano plots
- Section 17 - Differential gene tables
- Section 18 - GO enrichment
- Section 19 - Log2FC heatmap
- Section 20 - Coexpression network (hdWGCNA) on significant genes
- Section 21 - Network export and visualization
- Section 21b - TF coexpression network

## Regeneration Notes

The current Chapter 1 PDF was generated from a local HTML/PDF builder during
drafting. Future cleanup should either:

- move that builder into `reports/` as a reproducible script, or
- convert the final tutorial into a tracked Rmd/Quarto document.

For Chapter 2, prefer starting directly from a tracked document template so the
next PDF can be regenerated from the repository without relying on temporary
working files.
