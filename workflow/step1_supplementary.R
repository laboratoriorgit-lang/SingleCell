# Step 1 (supplementary) - optional steps of the single-cell analysis (R)
# Run after Step 2.5 (needs `ath_sc` annotated and the directories created by init_pipeline()).
# These steps are not needed by the rest of the workflow.

# ==============================================================================
# Supplementary Section S1 - Cell-type grouping
# ==============================================================================
# Merge the two mesophyll subclusters into a single label.
# Set grouping <- c() to skip this step (celltype_grouped == celltype).
grouping <- c(
  "Mesophyll.1" = "Mesophyll",
  "Mesophyll.2" = "Mesophyll"
)

output_dir <- dir_05

current <- as.character(ath_sc$celltype)
ath_sc$celltype_grouped <- ifelse(
  current %in% names(grouping), grouping[current], current
)

save_pdf(
  plot_umap(ath_sc, "celltype_grouped"),
  "umap_grouped.pdf"
)

message("\nOK S1 COMPLETE: Cell-type grouping")


# ==============================================================================
# Supplementary Section S2 - Subcluster marker inspection
# ==============================================================================
output_dir <- dir_05
Idents(ath_sc) <- "celltype"

# 1. Inspect a cluster (saves the combined subcluster + marker figure)
inspect_subcluster_markers(
  ath_sc, cluster_id = "Mesophyll.1", output_dir = output_dir
)

# 2. Map each cluster's subclusters to labels ("others" = any ID not listed),
#    then curate them in one call. Here Mesophyll's subclusters all map back
#    to "Mesophyll" (no relabeling) -- use this pattern when inspection
#    confirms the original label was already correct.
#    The keys of `reassign` are labels of celltype_grouped (created in
#    Supplementary Section S1), so annot_col and source_col point to that column:
reassign <- list(
  "Mesophyll" = c("others" = "Mesophyll")
)
ath_sc <- curate_clusters(
  ath_sc, reassign, output_dir = output_dir,
  annot_col = "celltype_grouped", source_col = "celltype_grouped"
)

# 3. Save the curated annotation UMAP
save_pdf(
  plot_umap(ath_sc, "celltype_curated"),
  "umap_curated.pdf"
)

message("\nOK S2 COMPLETE: Subcluster inspection and manual curation")

# ==============================================================================
# Supplementary Section S3 - Gene expression violin and feature plots
# ==============================================================================
# Violin plot and FeaturePlot of any gene (Step 2.5 only draws the gene given to
# run_annotation()). Use one AGI locus code, or a vector of them (one panel per gene).
gene <- "AT4G21750"
gene_tag <- paste(gene, collapse = "_")

output_dir <- dir_04
ath_sc <- JoinLayers(ath_sc)
Idents(ath_sc) <- "celltype"

save_vln(
  plot_gene_violin(ath_sc, gene),
  paste0("vln_gene_", gene_tag, ".pdf")
)
save_pdf(
  plot_gene_umap(ath_sc, gene),
  paste0("feature_gene_", gene_tag, ".pdf")
)

message("\nOK S3 COMPLETE: Gene expression violin and feature plots")
