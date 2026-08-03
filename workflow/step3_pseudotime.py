# Step 3 - Pseudotime trajectory analysis (Python)
# Auto-derived from the methods paper (Part 4). scanpy / Palantir / scFates
# are Python-only, so this step is .py. Run with the container python.

# ==============================================================================
# Section 23 - Setup
# ==============================================================================
import os

PIPELINE_DIR = "/workspace/workflow"
DATA_DIR     = "/workspace/."
base_dir     = os.path.join(DATA_DIR, "resultados")

exec(open(os.path.join(PIPELINE_DIR, "load_libraries_python.py")).read(), globals())
exec(open(os.path.join(PIPELINE_DIR, "ScRNA_Pseudotime_Functions.py")).read(), globals())

dir_pseudotime = os.path.join(base_dir, "09_pseudotime")
os.makedirs(dir_pseudotime, exist_ok=True)

print("SECTION 23 COMPLETE: setup done")

# ==============================================================================
# Section 24 - Load curated object
# ==============================================================================
INPUT_H5AD     = "/workspace/resultados/objects/pbmc_harmony_curated.h5ad"
ANNOTATION_COL = "celltype_curated"
N_JOBS         = 4

adata, N_JOBS = load_curated_object(
    input_h5ad     = INPUT_H5AD,
    dir_pseudotime = dir_pseudotime,
    annotation_col = ANNOTATION_COL,
    n_jobs         = N_JOBS,
)

# ==============================================================================
# Section 25 - Cell type selection
# ==============================================================================
TRAJECTORY_CLUSTERS = ["Epidermis"]

adata_sub = preview_trajectory_selection(
    adata          = adata,
    clusters       = TRAJECTORY_CLUSTERS,
    annotation_col = ANNOTATION_COL,
    dir_pseudotime = dir_pseudotime,
)

# ==============================================================================
# Section 26 - Trajectory inference
# ==============================================================================
ROOT_CLUSTER = "Epidermis"

TRAJECTORY_RUNS = [
    trajectory_run(nodes=50, sigma=0.1, lambda_value=200, eigs=8, seed=3),
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
# Section 27 - Plot genes on trajectory
# ==============================================================================
GENES = ["AT3G24140", "AT4G21750", "AT4G20260"]
name  = os.path.basename(selected_trajectory_dir)

ensure_expression_in_X(adata_traj)

gene_plots_dir = os.path.join(selected_trajectory_dir, "gene_plots")
os.makedirs(gene_plots_dir, exist_ok=True)

fig = sc.pl.draw_graph(
    adata_traj, color=GENES, title=GENES, show=False, return_fig=True
)
save_plot_18x18(fig, os.path.join(gene_plots_dir, f"{name}_gene_plots.pdf"))
plt.close(fig)

print("SECTION 27 COMPLETE: gene trajectory plots saved")

# ==============================================================================
# Section 28 - Gene trends
# ==============================================================================
TOP_N           = 10
HIGHLIGHT_GENES = ["AT3G24140", "AT4G21750", "AT4G20260"]
ORDERING        = "max"
SPLINE_DF       = 3

top_path, highlight_path = run_step29_gene_trends(
    adata        = adata_traj,
    run_dir      = selected_trajectory_dir,
    name         = os.path.basename(selected_trajectory_dir),
    custom_genes = HIGHLIGHT_GENES,
    top_n        = TOP_N,
    ordering     = ORDERING,
    n_jobs       = N_JOBS,
    spline_df    = SPLINE_DF,
)

print("SECTION 28 COMPLETE: gene trend plots saved")

# ==============================================================================
# Section 29 - Final export
# ==============================================================================
# Save the fully processed pseudotime object (trajectory tree, pseudotime `t`,
# milestones, segments and the force-directed layout) as the pipeline's final
# Chapter 3 output, next to the curated Chapter 1 object.
final_h5ad = os.path.join(base_dir, "objects", "pbmc_pseudotime_final.h5ad")
os.makedirs(os.path.dirname(final_h5ad), exist_ok=True)
adata_traj.write_h5ad(final_h5ad)

print(f"SECTION 29 COMPLETE: final pseudotime object saved -> {final_h5ad}")
