# =============================================================================
# ScRNA_Pseudotime_Wrappers.py  -  "version_RG" (condensed) pipeline
# =============================================================================
# One function that runs Steps 3.1-3.7 of step3_pseudotime.py with the same
# calls. Load it AFTER init_pipeline():
#
#   exec(open("/workspace/workflow/init.py").read(), globals())
#   init_pipeline(data_dir="/workspace/.")
#   exec(open("/workspace/workflow/ScRNA_Pseudotime_Wrappers.py").read(), globals())
# =============================================================================
import os


def run_pseudotime(input_h5ad, clusters, root_cluster, runs,
                   annotation_col="celltype", n_jobs=4,
                   plot_genes=None, plot_titles=None, plot_filename="trajectory_genes.pdf",
                   spline_df=3, top_n=10, highlight_genes=None, ordering="max"):
    """Steps 3.2-3.6: load, select, trajectory, gene plots, gene trends.

    Returns (adata_traj, selected_trajectory_dir, adata_fitted).
    """
    # 3.2 - load annotated object
    adata, n_jobs = load_curated_object(
        input_h5ad=input_h5ad, dir_pseudotime=dir_pseudotime,
        annotation_col=annotation_col, n_jobs=n_jobs)

    # 3.3 - cell type selection
    preview_trajectory_selection(
        adata=adata, clusters=clusters, annotation_col=annotation_col,
        dir_pseudotime=dir_pseudotime)
    print("[3.2] cell type selection saved")

    # 3.4 - trajectory inference
    adata_traj, selected_trajectory_dir, trajectory_runs = run_trajectory_runs(
        adata=adata, clusters=clusters, root_cluster=root_cluster,
        annotation_col=annotation_col, output_base_dir=dir_pseudotime, runs=runs)
    print("[3.2] trajectory inference saved")

    # 3.5 - genes on the trajectory (optional)
    if plot_genes:
        plot_genes_on_trajectory(
            adata_traj, genes=plot_genes, titles=plot_titles or plot_genes,
            run_dir=selected_trajectory_dir, filename=plot_filename)
        print("[3.2] gene trajectory plots saved")

    # 3.6 - gene trends (compute once, then select and plot)
    adata_fitted = compute_gene_trends(
        adata_traj, selected_trajectory_dir, n_jobs=n_jobs, spline_df=spline_df)
    adata_fitted = plot_gene_trends(
        adata_fitted, selected_trajectory_dir,
        custom_genes=highlight_genes, top_n=top_n, ordering=ordering,
        adata=adata_traj, n_jobs=n_jobs)
    print("OK 3.2 COMPLETE: trajectory inference and gene trends complete")

    return adata_traj, selected_trajectory_dir, adata_fitted
