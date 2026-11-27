# Step 3 - Pseudotime trajectory analysis (Python)
# Auto-derived from the methods paper (Part 4). scanpy / Palantir / scFates
# are Python-only, so this step is .py. Run with the container python.

# ==============================================================================
# Section 18 - Setup
# ==============================================================================
# data_dir: root of the project (the cloned repository is mounted as /workspace).
# init_pipeline() loads the Python packages and helper functions and creates
# results/09_pseudotime/. It defines DATA_DIR, PIPELINE_DIR, base_dir and dir_pseudotime.
exec(open("/workspace/workflow/init.py").read(), globals())
init_pipeline(data_dir="/workspace/.")

# ==============================================================================
# Section 19 - Load curated object
# ==============================================================================
INPUT_H5AD     = "/workspace/results/objects/ath_sc_annotated.h5ad"
# "celltype" (Section 7) is always present. Switch to "celltype_grouped" if
# Supplementary Section S1 was used, or "celltype_curated" if Supplementary Section S2 was used.
ANNOTATION_COL = "celltype"
N_JOBS         = 4

adata, N_JOBS = load_curated_object(
    input_h5ad     = INPUT_H5AD,
    dir_pseudotime = dir_pseudotime,
    annotation_col = ANNOTATION_COL,
    n_jobs         = N_JOBS,
)

# ==============================================================================
# Section 20 - Cell type selection
# ==============================================================================
TRAJECTORY_CLUSTERS = ["Epidermis Cotyledon"]

adata_sub = preview_trajectory_selection(
    adata          = adata,
    clusters       = TRAJECTORY_CLUSTERS,
    annotation_col = ANNOTATION_COL,
    dir_pseudotime = dir_pseudotime,
)


# ==============================================================================
# Section 21 - Trajectory inference
# ==============================================================================
# ROOT_CLUSTER is the cell type used to choose the root when the trajectory has two or more
# cell types. root_method of run_trajectory_runs() selects how the root is chosen: "auto"
# (default, the tip with the highest mean number of detected genes), "score", "most_connected",
# "medoid" or "tip"; see Supplementary Section S4.
ROOT_CLUSTER = "Epidermis Cotyledon"

TRAJECTORY_RUNS = [
    trajectory_run(nodes=100, sigma=0.1, lambda_value=100, eigs=3, seed=3),
]

adata_traj, selected_trajectory_dir, trajectory_runs = run_trajectory_runs(
    adata           = adata,
    clusters        = TRAJECTORY_CLUSTERS,
    root_cluster    = ROOT_CLUSTER,
    annotation_col  = ANNOTATION_COL,
    output_base_dir = dir_pseudotime,
    runs            = TRAJECTORY_RUNS,
)

# ==============================================================================
# Section 22 - Plot genes on trajectory
# ==============================================================================
# ATML1, the gene of interest, with a readable panel title
plot_genes_on_trajectory(
    adata_traj,
    genes    = ["AT4G21750"],                         # ATML1
    titles   = ["ATML1 (AT4G21750)"],
    run_dir  = selected_trajectory_dir,
)

# ==============================================================================
# Section 23 - Gene trends
# ==============================================================================
SPLINE_DF = 3

# Step 1 - compute (slow, ~15 min). The fitted trends are saved in
# <run>/gene_trends/gene_trends_fitted.h5ad and reused if you run it again.
adata_fitted = compute_gene_trends(
    adata_traj, selected_trajectory_dir,
    n_jobs    = N_JOBS,
    spline_df = SPLINE_DF,
)

TOP_N           = 10
HIGHLIGHT_GENES = ["AT4G21750"]
ORDERING        = "max"

# Step 2 - select genes and plot (seconds). Change TOP_N, HIGHLIGHT_GENES or
# ORDERING and run only this block again; nothing is recomputed. A highlighted
# gene that is not in the fit is fitted on its own and added to it.
adata_fitted = plot_gene_trends(
    adata_fitted, selected_trajectory_dir,
    custom_genes = HIGHLIGHT_GENES,
    top_n        = TOP_N,
    ordering     = ORDERING,
    adata        = adata_traj,
    n_jobs       = N_JOBS,
)


# ==============================================================================
# Section 24 - Final export
# ==============================================================================
# Save the fully processed pseudotime object (trajectory tree, pseudotime `t`,
# milestones, segments and the force-directed layout) as the pipeline's final
# Chapter 3 output, next to the curated Chapter 1 object.
export_final_object(adata_traj, name="ath_sc_pseudotime_final")
