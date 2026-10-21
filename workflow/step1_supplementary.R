# Step 1 (supplementary) - optional steps of the single-cell analysis (R)
# Run after Section 7 of step1_singlecell.R (needs `ath_sc` annotated, `marker_table`,
# the directories created in the init block and, for S1, `clu` from Section 5).
# These steps are not needed by the rest of the workflow.

# ==============================================================================
# Supplementary Section S1 - Annotated clustree
# ==============================================================================
output_dir <- dir_03

plot_annotated_clustree(
  ath_sc, clu, annot_col = "celltype", output_dir = output_dir
)

message("\nOK SUPPLEMENTARY S1 COMPLETE: annotated clustree saved")


# ==============================================================================
# Supplementary Section S2 - Cell-type grouping
# ==============================================================================
# Merge the two mesophyll subclusters into a single label.
# Set grouping <- c() to skip this step (celltype_grouped == celltype).
grouping <- c(
  "Mesophyll.1" = "Mesophyll",
  "Mesophyll.2" = "Mesophyll"
)

output_dir <- dir_05

if (length(grouping) > 0) {
  current <- as.character(ath_sc$celltype)
  ath_sc$celltype_grouped <- ifelse(current %in% names(grouping), grouping[current], current)
} else {
  ath_sc$celltype_grouped <- ath_sc$celltype
}

save_pdf(
  DimPlot(
    ath_sc, group.by = "celltype_grouped",
    label = TRUE, repel = TRUE, raster = FALSE
  ),
  "umap_grouped.pdf"
)

message("\nOK SUPPLEMENTARY S2 COMPLETE: cell-type grouping complete")


# ==============================================================================
# Supplementary Section S3 - Subcluster marker inspection
# ==============================================================================
output_dir <- dir_05
Idents(ath_sc) <- "celltype"

# 1. Inspect a cluster (saves the combined subcluster + marker figure)
inspect_subcluster_markers(
  ath_sc, cluster_id = "Mesophyll",
  marker_table = marker_table, output_dir = output_dir
)

# 2. Map each cluster's subclusters to labels ("others" = any ID not listed),
#    then curate them in one call. Here Mesophyll's subclusters all map back
#    to "Mesophyll" (no relabeling) -- use this pattern when inspection
#    confirms the original label was already correct:
reassign <- list(
  "Mesophyll" = c("others" = "Mesophyll")
)
ath_sc <- curate_clusters(
  ath_sc, reassign,
  marker_table = marker_table, output_dir = output_dir
)

# 3. Save the curated annotation UMAP
save_pdf(
  DimPlot(ath_sc, group.by = "celltype_curated",
          label = TRUE, repel = TRUE, raster = FALSE),
  "umap_curated.pdf"
)

message("\nOK SUPPLEMENTARY S3 COMPLETE: curation complete")


