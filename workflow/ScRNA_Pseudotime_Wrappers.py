# =============================================================================
# ScRNA_Pseudotime_Wrappers.py  -  "version_RG" (condensed) pipeline
# =============================================================================
# One function that runs Sections 18-24 of step3_pseudotime.py with the same
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
    """Sections 19-24: load, select, trajectory, gene plots, gene trends, export.

    Returns (adata_traj, selected_trajectory_dir, adata_fitted).
    """
    # Section 19 - load annotated object
    adata, n_jobs = load_curated_object(
        input_h5ad=input_h5ad, dir_pseudotime=dir_pseudotime,
        annotation_col=annotation_col, n_jobs=n_jobs)

    # Section 20 - cell type selection
    preview_trajectory_selection(
        adata=adata, clusters=clusters, annotation_col=annotation_col,
        dir_pseudotime=dir_pseudotime)
    print("SECTION 20 COMPLETE: cell type selection saved")

    # Section 21 - trajectory inference
    adata_traj, selected_trajectory_dir, trajectory_runs = run_trajectory_runs(
        adata=adata, clusters=clusters, root_cluster=root_cluster,
        annotation_col=annotation_col, output_base_dir=dir_pseudotime, runs=runs)
    print("SECTION 21 COMPLETE: trajectory inference saved")

    # Section 22 - genes on the trajectory (optional)
    if plot_genes:
        plot_genes_on_trajectory(
            adata_traj, genes=plot_genes, titles=plot_titles or plot_genes,
            run_dir=selected_trajectory_dir, filename=plot_filename)
        print("SECTION 22 COMPLETE: gene trajectory plots saved")

    # Section 23 - gene trends (compute once, then select and plot)
    adata_fitted = compute_gene_trends(
        adata_traj, selected_trajectory_dir, n_jobs=n_jobs, spline_df=spline_df)
    adata_fitted = plot_gene_trends(
        adata_fitted, selected_trajectory_dir,
        custom_genes=highlight_genes, top_n=top_n, ordering=ordering,
        adata=adata_traj, n_jobs=n_jobs)
    print("SECTION 23 COMPLETE: gene trends saved")

    # Section 24 - final export
    final_h5ad = os.path.join(base_dir, "objects", "ath_sc_pseudotime_final.h5ad")
    os.makedirs(os.path.dirname(final_h5ad), exist_ok=True)
    adata_traj.write_h5ad(final_h5ad)
    print(f"SECTION 24 COMPLETE: final pseudotime object saved -> {final_h5ad}")
    return adata_traj, selected_trajectory_dir, adata_fitted
