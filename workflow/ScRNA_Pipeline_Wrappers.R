# =============================================================================
# ScRNA_Pipeline_Wrappers.R  -  "version_RG" (condensed) pipeline
# =============================================================================
# High-level functions that each run several consecutive sections of the
# step-by-step pipeline (step1_singlecell.R / step2_degs.R) with the same
# calls, saving the same PDFs and objects. Load them AFTER init_pipeline():
#
#   source("/workspace/workflow/init.R")
#   init_pipeline(data_dir = "/workspace/.", seed = 1807)
#   source("/workspace/workflow/ScRNA_Pipeline_Wrappers.R")
#
# They use the globals created by init_pipeline() (DATA_DIR, dir_01..dir_08,
# dir_objects) and the helper functions of ScRNA_Analysis_Functions.R.
# =============================================================================

# -- internal: point save_pdf()/save_qc() to a folder ------------------------
.set_output_dir <- function(d) assign("output_dir", d, envir = globalenv())

#' Sections 1-4: load, QC, filter, merge, Harmony
#' @return list(ath_sc, seurat_list)
run_preprocessing <- function(samples, colors,
                              mt_pattern        = "^ATMG",
                              min_features      = 200,
                              max_mt            = 5,
                              run_doubletfinder = TRUE,
                              nfeatures         = 2000,
                              npcs              = 30,
                              dims              = 1:30) {
  .set_output_dir(dir_01)

  # Section 1 - loading and pre-filter QC
  seurat_list_raw <- load_seurat_samples(samples = samples, DATA_DIR = DATA_DIR,
                                         mt_pattern = mt_pattern)
  plot_qc_batch(seurat_list_raw, colors, "qc_prefilter.pdf")
  message("\nOK SECTION 1 COMPLETE: QC pre-filter plots saved")

  # Section 2 - filtering and doublets
  seurat_list <- filter_seurat_samples(
    seurat_list_raw, min_features = min_features, max_mt = max_mt,
    run_doubletfinder = run_doubletfinder
  )
  plot_qc_batch(seurat_list, colors, "qc_postfilter.pdf")
  saveRDS(seurat_list, file.path(dir_objects, "seurat_list_postfilter.rds"))
  message("\nOK SECTION 2 COMPLETE: filtering and doublet detection complete")

  # Section 3 - merge and preprocessing
  ath_sc <- merge_and_preprocess(seurat_list, nfeatures = nfeatures, npcs = npcs)
  save_pdf(plot_umap(ath_sc, "orig.ident", colors), "umap_preharmony.pdf")
  saveRDS(ath_sc, file.path(dir_objects, "ath_sc_preharmony.rds"))
  message("\nOK SECTION 3 COMPLETE: merge and preprocessing complete")

  # Section 4 - Harmony
  ath_sc <- run_harmony_umap(ath_sc, group_by = "orig.ident", dims = dims)
  save_pdf(plot_umap(ath_sc, "orig.ident", colors), "umap_postharmony.pdf")
  saveRDS(ath_sc, file.path(dir_objects, "ath_sc_postharmony.rds"))
  message("\nOK SECTION 4 COMPLETE: Harmony integration complete")

  list(ath_sc = ath_sc, seurat_list = seurat_list)
}

#' Sections 5-6: elbow, resolution sweep (clustree) and final clustering
#' @return the Seurat object with seurat_clusters
run_clustering <- function(ath_sc,
                           resolutions_test   = c(0.15, 0.35, 0.45, 0.55, 1.0),
                           cluster_resolution = 0.35,
                           dims               = 1:30,
                           k_param            = 20) {
  # Section 5 - elbow plot and clustree
  plot_resolution_elbow(ath_sc, output_dir = dir_02)
  run_resolution_sweep(ath_sc, resolutions = resolutions_test, output_dir = dir_02)
  message("\nOK SECTION 5 COMPLETE: elbow plot and clustree saved")

  # Section 6 - final clustering
  .set_output_dir(dir_02)
  ath_sc <- cluster_cells(ath_sc, resolution = cluster_resolution,
                          dims = dims, k_param = k_param)
  Idents(ath_sc) <- "seurat_clusters"
  save_pdf(plot_umap(ath_sc, "seurat_clusters"), "umap_seuratclusters.pdf")
  message("\nOK SECTION 6 COMPLETE: final clustering complete")
  ath_sc
}

#' Sections 7-8: markers, annotation, expression plots
#' @return the annotated Seurat object
run_annotation <- function(ath_sc, marker_file, gene = NULL) {
  marker_file  <- file.path(DATA_DIR, marker_file)
  marker_table <- standardize_marker_table(marker_file)

  # Section 7 - annotation
  .set_output_dir(dir_03)
  markers <- find_markers(ath_sc, output_file = file.path(dir_03, "FindAllMarkers.tsv"),
                          force = TRUE)
  ath_sc <- annotate_by_markers(ath_sc, markers, reference_file = marker_file)
  plot_marker_dotplot(ath_sc, marker_table, annot_col = "celltype", name_gene = TRUE,
                      outfile = file.path(dir_03, "dotplot_marker_table_annotation_biblio.pdf"),
                      width = 18, height = 18)
  save_pdf(plot_umap(ath_sc, "celltype"), "umap_annotation_biblio.pdf")
  message("\nOK SECTION 7 COMPLETE: cell-type annotation complete")

  # Section 8 - gene expression (optional)
  if (!is.null(gene)) {
    .set_output_dir(dir_04)
    ath_sc <- JoinLayers(ath_sc)
    Idents(ath_sc) <- "celltype"
    gene_tag <- paste(gene, collapse = "_")
    save_vln(plot_gene_violin(ath_sc, gene), paste0("vln_gene_", gene_tag, ".pdf"))
    save_pdf(plot_gene_umap(ath_sc, gene), paste0("feature_gene_", gene_tag, ".pdf"))
    message("\nOK SECTION 8 COMPLETE: gene expression visualization saved")
  }

  ath_sc
}

#' Sections 10-16: subsets, pseudo-replicates, DESeq2, volcano, tables, GO, heatmap
#' @return list(deseq2, diff_tables, go)
run_pseudobulk_de <- function(ath_sc, marker_file,
                              comparisons,
                              annot_col             = "celltype",
                              n_pseudoreps          = 3,
                              pseudobulk_conditions = NULL,
                              cell_types            = NULL,
                              padj_cut              = 0.05,
                              lfc_cut               = 1,
                              go_space              = "BP",
                              go_orgdb              = org.At.tair.db,
                              go_keytype            = "TAIR",
                              heatmap_limits        = c(-5, 5)) {
  marker_table <- standardize_marker_table(file.path(DATA_DIR, marker_file))

  # Sections 10-11 - subsets and pseudo-replicates
  subsets <- create_cell_type_subsets(ath_sc, annot_col = annot_col)
  reps <- assign_pseudoreplicates_batch(subsets,
            pseudobulk_conditions = pseudobulk_conditions,
            n_pseudoreps = n_pseudoreps)
  message("\nOK SECTIONS 10-11 COMPLETE: subsets and pseudo-replicates")

  # Section 12 - pseudobulk + DESeq2 (all comparisons)
  deseq2_results <- run_pseudobulk_deseq2_analysis(
    cell_type_subsets_replicates = reps, comparisons = comparisons,
    output_dir = dir_06, cell_types = cell_types,
    pseudobulk_dir = file.path(dir_objects, "pseudobulk_replicates"))
  message("\nOK SECTION 12 COMPLETE: pseudobulk aggregation and DESeq2 complete")

  # Sections 13-16 - one pass per comparison
  diff_tables <- list(); go_results <- list()
  for (cmp in comparisons) {
    tag <- cmp$tag
    render_volcano_plots(
      results_dir = file.path(dir_06, tag),
      output_dir  = file.path(dir_06, tag, "volcano"),
      pdf_name    = paste0("VolcanoPlots_", tag, ".pdf"),
      padj_cut = padj_cut, lfc_cut = lfc_cut)
    diff_tables[[tag]] <- build_differential_tables(
      results_dir = file.path(dir_06, tag), output_dir = file.path(dir_06, tag),
      padj_cut = padj_cut, lfc_cut = lfc_cut, prefix = paste0("diff_table_", tag))
    go_results[[tag]] <- run_go_enrichment_for_contrast(
      results_dir = file.path(dir_06, tag), output_dir = file.path(dir_07, tag),
      orgdb = go_orgdb, keytype = go_keytype, go_space = go_space,
      padj_cutoff = padj_cut, contrast_tag = tag)
    build_logfc_heatmap(
      logfc_table = diff_tables[[tag]]$logfc, contrast_tag = tag,
      output_dir = file.path(dir_06, tag), limits = heatmap_limits,
      cell_order = unique(marker_table$cell.types))
  }
  message("\nOK SECTIONS 13-16 COMPLETE: volcano, tables, GO and heatmap saved")
  invisible(list(deseq2 = deseq2_results, diff_tables = diff_tables, go = go_results))
}

#' Section 17: hdWGCNA co-expression and TF network
run_tf_network <- function(ath_sc, contrast_tag, annot_col = "celltype",
                           tf_list = "data/AtTFDB_loci.txt",
                           n_metacells = 50, min_module_size = 20, deep_split = 2,
                           soft_power = NULL, tom_threshold = 0.2, n_hub_label = 15) {
  run_hdwgcna_tf_network(
    seurat_obj    = ath_sc,
    de_table_path = file.path(dir_06, contrast_tag, "log2FC_table_fc1_padj_005.tsv"),
    tf_list_path  = file.path(DATA_DIR, tf_list),
    output_dir    = dir_08, annot_col = annot_col, contrast_tag = contrast_tag,
    n_metacells = n_metacells, min_module_size = min_module_size,
    deep_split = deep_split, soft_power = soft_power,
    tom_threshold = tom_threshold, n_hub_label = n_hub_label)
  message("\nOK SECTION 17 COMPLETE: TF coexpression network saved")
}
