# Step 3 - Pseudotime trajectory analysis (Python)
# Auto-derived from the methods paper (Part 4). scanpy / Palantir / scFates
# are Python-only, so this step is .py. Run with the container python.

# ==============================================================================
# Section 18 - Setup
# ==============================================================================
import os

PIPELINE_DIR = "/workspace/workflow"
DATA_DIR     = "/workspace/."
base_dir     = os.path.join(DATA_DIR, "results")

exec(open(os.path.join(PIPELINE_DIR, "load_libraries_python.py")).read(), globals())
exec(open(os.path.join(PIPELINE_DIR, "ScRNA_Pseudotime_Functions.py")).read(), globals())

dir_pseudotime = os.path.join(base_dir, "09_pseudotime")
os.makedirs(dir_pseudotime, exist_ok=True)

print("SECTION 18 COMPLETE: setup done")

# ==============================================================================
# Section 19 - Load curated object
# ==============================================================================
INPUT_H5AD     = "/workspace/results/objects/ath_sc_annotated.h5ad"
# "celltype" (Section 7) is always present. Switch to "celltype_grouped" if
# Supplementary Section S2 was used, or "celltype_curated" if Supplementary Section S3 was used.
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

print("SECTION 20 COMPLETE: cell type selection saved")

# ==============================================================================
# Section 21 - Trajectory inference
# ==============================================================================
# ROOT_CLUSTER is the root cell type when the trajectory has two or more cell types
# (the root is its most-connected cell). With a single cell type, as here, the root
# is set automatically at the tip of the tree with the most detected genes (see
# "Root selection for pseudotime" in the README); root_score, root_lowest, root_tip
# and root_method of run_trajectory_runs() change it.
ROOT_CLUSTER = "Epidermis Cotyledon"

TRAJECTORY_RUNS = [
    trajectory_run(nodes=100, sigma=0.1, lambda_value=20, eigs=3, seed=3),
]

adata_traj, selected_trajectory_dir, trajectory_runs = run_trajectory_runs(
    adata           = adata,
    clusters        = TRAJECTORY_CLUSTERS,
    root_cluster    = ROOT_CLUSTER,
    annotation_col  = ANNOTATION_COL,
    output_base_dir = dir_pseudotime,
    runs            = TRAJECTORY_RUNS,
)

print("SECTION 21 COMPLETE: trajectory inference saved")

# ==============================================================================
# Section 22 - Plot genes on trajectory
# ==============================================================================
# ATML1, the gene of interest, with a readable panel title
plot_genes_on_trajectory(
    adata_traj,
    genes    = ["AT4G21750"],                         # ATML1
    titles   = ["ATML1 (AT4G21750)"],
    run_dir  = selected_trajectory_dir,
    filename = "trajectory_genes_ATML1.pdf",
)

print("SECTION 22 COMPLETE: gene trajectory plots saved")

# ==============================================================================
# Section 23 - Gene trends
# ==============================================================================
TOP_N           = 10
HIGHLIGHT_GENES = ["AT4G21750"]
ORDERING        = "max"
SPLINE_DF       = 3

name = os.path.basename(selected_trajectory_dir)

# Step 1 - compute (slow, ~15 min). The fitted curves are saved in
# <run>/gene_trends/gene_trends_fitted.h5ad and reused if you run it again.
adata_fitted = compute_gene_trends(
    adata        = adata_traj,
    run_dir      = selected_trajectory_dir,
    name         = name,
    custom_genes = HIGHLIGHT_GENES,
    n_jobs       = N_JOBS,
    spline_df    = SPLINE_DF,
)

# Step 2 - plot (seconds). Change TOP_N, HIGHLIGHT_GENES or ORDERING and run only
# this call again; nothing is recomputed.
top_path, highlight_path = plot_gene_trends(
    adata_fitted,
    run_dir      = selected_trajectory_dir,
    name         = name,
    custom_genes = HIGHLIGHT_GENES,
    top_n        = TOP_N,
    ordering     = ORDERING,
)

print("SECTION 23 COMPLETE: gene trend plots saved")

# ==============================================================================
# Section 24 - Final export
# ==============================================================================
# Save the fully processed pseudotime object (trajectory tree, pseudotime `t`,
# milestones, segments and the force-directed layout) as the pipeline's final
# Chapter 3 output, next to the curated Chapter 1 object.
final_h5ad = os.path.join(base_dir, "objects", "ath_sc_pseudotime_final.h5ad")
os.makedirs(os.path.dirname(final_h5ad), exist_ok=True)
adata_traj.write_h5ad(final_h5ad)

print(f"SECTION 24 COMPLETE: final pseudotime object saved -> {final_h5ad}")
