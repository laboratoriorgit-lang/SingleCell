# =============================================================================
# ScRNA_Pseudotime_Functions.py — simple helpers for Capitulo 3
# =============================================================================
# Companion to capitulo3_pseudotime.ipynb.
# This file intentionally contains only the functions used by the notebook:
#   STEP 25 — load data
#   STEP 26 — choose cell types
#   STEP 27 — build pseudotime trajectory
#   STEP 29 — top-N and custom-highlight gene trends
# =============================================================================

# STEP 25 — load_curated_object
# =============================================================================
PLOT_FIGSIZE = (18, 18)
PLOT_DPI = 200


def save_plot_18x18(fig, path):
    fig.set_size_inches(*PLOT_FIGSIZE, forward=True)
    fig.savefig(path, dpi=PLOT_DPI, facecolor="white")
    return path


# scanpy/scFates gene-colored plots read expression from .X, but trajectory
# objects often carry the normalized matrix in a layer ("logcounts") with .X
# left empty. Call this before any gene plot so it never fails with a
# "'NoneType' object is not subscriptable" error.
def ensure_expression_in_X(adata, layer="logcounts"):
    if adata.X is None and layer in adata.layers:
        adata.X = adata.layers[layer]
    return adata


# The R workflow can export up to three annotation columns depending on which
# optional sections were run: "celltype" (Section 7, always present), plus
# "celltype_grouped" (Supplementary Section S2) and/or "celltype_curated" (Supplementary Section S3) if
# the curation steps were used. Checked here so a missing ANNOTATION_COL fails
# with a clear message instead of a raw KeyError deep inside sc.pl.umap.
KNOWN_ANNOTATION_COLS = ["celltype", "celltype_grouped", "celltype_curated"]


# Notebook wrapper for Step 25. Loads the full curated object, fixes Seurat-style
# coordinate names, plots the overview UMAP, and prints available cell types.
def load_curated_object(input_h5ad, dir_pseudotime, annotation_col, n_jobs=4):

    print(f"Loading object: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    # Convert obsm keys from R/Seurat format to scanpy format
    # Handles both DataFrame (from older anndata) and ndarray (from newer anndata)
    for key_from, key_to in [("UMAP", "X_umap"), ("PCA", "X_pca"), ("HARMONY", "X_harmony")]:
        if key_from in adata.obsm and key_to not in adata.obsm:
            val = adata.obsm[key_from]
            adata.obsm[key_to] = val.values if hasattr(val, "values") else val

    available_annots = [c for c in KNOWN_ANNOTATION_COLS if c in adata.obs.columns]
    print(f"\nAvailable annotation columns in this object: {available_annots}")
    if annotation_col not in available_annots:
        raise ValueError(
            f"ANNOTATION_COL = '{annotation_col}' not found in this object. "
            f"Available annotation columns: {available_annots}. Set ANNOTATION_COL "
            "to one of these -- 'celltype_grouped'/'celltype_curated' are only "
            "present if Supplementary Section S2/11 were run in Step 1."
        )

    sc.settings.figdir = dir_pseudotime
    sc.set_figure_params(figsize=PLOT_FIGSIZE, dpi=100, dpi_save=PLOT_DPI)

    fig = sc.pl.umap(
        adata,
        color              = annotation_col,
        legend_loc         = "on data",
        legend_fontsize    = 6,
        legend_fontoutline = 2,
        frameon            = False,
        size               = 10,
        show               = False,
        return_fig         = True,
    )

    # Add a side legend with one color swatch per cell type, in addition to
    # the on-data labels above (scanpy only supports one legend_loc at a time).
    ax  = fig.axes[0]
    cats   = adata.obs[annotation_col].cat.categories
    colors = adata.uns[f"{annotation_col}_colors"]
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8)
        for c in colors
    ]
    ax.legend(
        handles, cats,
        loc="center left", bbox_to_anchor=(1.0, 0.5),
        fontsize=6, frameon=False,
    )

    fig.subplots_adjust(right=0.78)
    save_plot_18x18(fig, os.path.join(dir_pseudotime, "umap_overview.pdf"))
    plt.close(fig)

    counts = adata.obs[annotation_col].value_counts()
    print(f"\nCell types in '{annotation_col}':\n")
    for i, ct in enumerate(sorted(adata.obs[annotation_col].unique()), 1):
        print(f"  {i:2d}.  {ct:<30s}  {counts[ct]:>5d} cells")

    print("\nSTEP 25 COMPLETE: full object loaded")
    return adata, n_jobs


# =============================================================================
# STEP 26 — preview_trajectory_selection
# =============================================================================
# Notebook wrapper for Step 26. Subsets the selected cell types and displays a
# UMAP preview so the user can decide whether the selection makes biological sense.
def preview_trajectory_selection(adata, clusters, annotation_col, dir_pseudotime):
    missing = [x for x in clusters if x not in set(adata.obs[annotation_col])]
    if missing:
        raise ValueError(
            "These cell types were not found in "
            f"'{annotation_col}': {missing}. Copy names exactly from Step 25."
        )

    adata_sub = adata[adata.obs[annotation_col].isin(clusters)].copy()
    print(f"Selected {len(adata_sub)} cells: {clusters}")

    fig = sc.pl.umap(
        adata_sub,
        color              = annotation_col,
        legend_loc         = "on data",
        legend_fontsize    = 6,
        legend_fontoutline = 2,
        frameon            = False,
        size               = 10,
        show               = False,
        return_fig         = True,
    )

    # Side legend with one color swatch per cell type, in addition to the
    # on-data labels above (scanpy only supports one legend_loc at a time).
    ax     = fig.axes[0]
    cats   = adata_sub.obs[annotation_col].cat.categories
    colors = adata_sub.uns[f"{annotation_col}_colors"]
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=8)
        for c in colors
    ]
    ax.legend(
        handles, cats,
        loc="center left", bbox_to_anchor=(1.0, 0.5),
        fontsize=6, frameon=False,
    )

    fig.subplots_adjust(right=0.78)
    save_plot_18x18(fig, os.path.join(dir_pseudotime, "umap_selection.pdf"))
    plt.close(fig)

    print("\nSTEP 26 COMPLETE: continue only if this UMAP looks biologically sensible")
    return adata_sub


# =============================================================================
# STEP 27 — trajectory_run
# =============================================================================
# Small helper to describe one trajectory parameter set as a dict. Used to
# build the TRAJECTORY_RUNS list in the notebook's Step 27 cell.
# The run's name/folder is generated from the parameters themselves (not a
# free-text label), so two runs are named the same only if every parameter
# matches — in that case they are the same run and should share a folder.
# seed controls the stochastic steps inside build_pseudotime_trajectory
# (neighbors, ForceAtlas2 layout, tree fitting) — same seed + same params
# always reproduces the same result.
def trajectory_run(nodes=50, sigma=0.2, lambda_value=60, eigs=20, seed=3):
    name = f"n{nodes}_s{sigma}_l{lambda_value}_e{eigs}_seed{seed}"
    return {
        "name": name,
        "nodes": nodes,
        "sigma": sigma,
        "lambda": lambda_value,
        "eigs": eigs,
        "seed": seed,
    }


# Hidden legacy helper, not called by the main workflow unless root_method is
# "auto" with several cell types or "most_connected": the most-connected cell of the
# root_cluster population (largest summed connectivity within that population).
# With a single cell type this is the center of the cloud and it changes with the
# graph parameters, which is why the single-type case uses the tip rule instead.
def _most_connected_cell(adata, annotation_col, root_cluster):
    mask          = (adata.obs[annotation_col] == root_cluster).to_numpy()
    cluster_cells = adata.obs_names[mask]
    sub_conn      = adata.obsp["connectivities"][mask][:, mask]
    degrees       = np.array(sub_conn.sum(axis=1)).flatten()
    return cluster_cells[degrees.argmax()]


# =============================================================================
# STEP 27 — build_pseudotime_trajectory
# =============================================================================
# Full pipeline: subset cells → Palantir diffusion maps → scFates tree →
# root selection → pseudotime assignment. Called by run_trajectory_runs,
# which is what the notebook's Step 27 cell actually calls.
def build_pseudotime_trajectory(
    adata,
    clusters,           # list of cell types to include (values from annotation_col)
    root_cluster,       # cell type where pseudotime = 0 (the biological progenitor)
    annotation_col,     # adata.obs column containing cell type labels
    nodes       = 150,  # tree nodes: more = finer branches, slower to compute
    sigma       = 0.2,  # smoothing: lower = tree follows cells more tightly
    ppt_lambda  = 60,   # complexity: higher = simpler tree with fewer branches
    n_components= 50,   # total diffusion map dimensions computed by Palantir
    n_eigs      = 20,   # dimensions retained (must be strictly < n_components)
    n_neighbors = 50,   # neighbors for final graph (higher = smoother layout)
    seed        = 3,    # random seed for reproducibility across runs
    root_method = "auto",           # "auto" (default), "score" or "tip"; see the root selection notes below
    root_score  = "nFeature_RNA",   # obs column or gene ID that anchors the root (single cell type)
    root_lowest = False,            # True = the root is where root_score is lowest
    root_tip    = None,             # "tip": node number of the tip (see the tree figure)
):
    # Subset to selected cell types
    if clusters:
        adata = adata[adata.obs[annotation_col].isin(clusters)].copy()

    # Palantir diffusion maps: captures global developmental geometry
    pca_proj = pd.DataFrame(adata.obsm["X_pca"], index=adata.obs_names)
    dm_res   = palantir.utils.run_diffusion_maps(pca_proj, n_components=n_components)
    ms_data  = palantir.utils.determine_multiscale_space(
        dm_res, n_eigs=min(n_eigs, n_components - 1)
    )
    val = ms_data.values if hasattr(ms_data, "values") else ms_data
    adata.obsm["X_palantir"] = val

    # Graph in Palantir space for tree construction
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_palantir", method="umap")
    adata.obsm["X_pca2d"] = adata.obsm["X_pca"][:, :2]
    sc.tl.draw_graph(adata, init_pos="X_pca2d")

    # Build principal tree (PPT = Principal Polynomial Tree)
    scf.tl.tree(
        adata,
        method     = "ppt",
        Nodes      = nodes,
        use_rep    = "palantir",
        plot       = False,
        device     = "cpu",
        seed       = seed,
        ppt_lambda = ppt_lambda,
        ppt_nsteps = 200,
        ppt_sigma  = sigma,
    )

    # Root selection (see "Root selection for pseudotime" in the README).
    #   "auto" (default):
    #       - more than one cell type in the trajectory: the root is the most-connected
    #         cell of the selected root_cluster (the root cell type chosen by the user).
    #       - a single cell type: the root is placed automatically at the TIP of the
    #         tree where the mean of root_score is highest (lowest if root_lowest).
    #         root_score defaults to the number of detected genes per cell, a proxy of
    #         developmental potential (Gulati et al. 2020). It is an assumption, so a
    #         message is printed; the direction of pseudotime is conventional.
    #   "score": always use the tip rule above, even with several cell types.
    #   "tip": the root is the tip with node number root_tip.
    #   "most_connected": hidden legacy option, the most-connected cell of root_cluster
    #         (see _most_connected_cell); not used by the main workflow.
    n_types = adata.obs[annotation_col].nunique()
    method  = root_method
    if method == "auto":
        method = "most_connected" if n_types > 1 else "score"
        label  = "selected root cell type" if n_types > 1 else "automatic tip rule (single cell type)"
    else:
        label  = method

    if method == "most_connected":
        root_cell = _most_connected_cell(adata, annotation_col, root_cluster)
        adata.obs["is_root"] = adata.obs_names == root_cell
        scf.tl.root(adata, root="is_root")
        scf.tl.pseudotime(adata)
    else:
        if method == "tip":
            if root_tip is None:
                raise ValueError("root_method = 'tip' needs root_tip (a tip node number).")
            scf.tl.root(adata, root=int(root_tip))
        elif method == "score":
            ensure_expression_in_X(adata)
            if root_score not in adata.obs.columns and root_score not in adata.var_names:
                # number of detected genes per cell, computed from the expression matrix
                adata.obs[root_score] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
            tip_means = None
            if root_score in adata.obs.columns:
                # mean of the score at each tip node (soft cell-to-node assignment, as scFates does)
                R = adata.obsm["X_R"]
                R = R.toarray() if hasattr(R, "toarray") else np.asarray(R)
                node_means = (R.T @ adata.obs[root_score].values.astype(float)) / np.maximum(R.sum(0), 1e-9)
                tip_means  = {int(t): float(node_means[int(t)]) for t in adata.uns["graph"]["tips"]}
            scf.tl.root(adata, root=root_score, tips_only=True, min_val=root_lowest)
        else:
            raise ValueError(f"Unknown root_method '{root_method}'.")
        scf.tl.pseudotime(adata)
        # root cell = the cell with the smallest pseudotime (closest to the root tip)
        root_cell = adata.obs["t"].idxmin()
        adata.obs["is_root"] = adata.obs_names == root_cell
    print(f"✓ Root: {label} | root node {int(adata.uns['graph']['root'])}")
    if method == "score":
        if tip_means:
            print(f"  mean {root_score} per tip: " + ", ".join(f"tip {t} = {m:.1f}" for t, m in tip_means.items()))
        print("  NOTE: with a single cell type the data do not say where the trajectory starts. "
              "The root was set by an automatic rule (an assumption), so the direction of pseudotime "
              "is conventional. Use root_score, root_tip or the tree figure to set it from your biology.")

    print(f"✓ Trajectory built — root cell: {root_cell}")
    return adata


# =============================================================================
# STEP 27 — output helpers for run_trajectory_runs
# =============================================================================
# Adds a dendrogram to an adata object that already has a scFates tree.
def build_dendrogram(adata):
    scf.tl.dendrogram(adata)
    return adata


# Returns the draw_graph embedding key actually present in adata. scanpy's
# sc.tl.draw_graph defaults to the ForceAtlas2 ("fa") layout, but silently falls
# back to Fruchterman-Reingold ("fr") when the optional 'fa2-modified' package is
# missing. We pick whichever embedding exists (fa preferred) so downstream plots
# work regardless of that package being installed.
def draw_graph_basis(adata):
    for layout in ("fa", "fr"):
        if f"X_draw_graph_{layout}" in adata.obsm:
            return f"draw_graph_{layout}"
    raise KeyError(
        "No X_draw_graph_* embedding found in adata.obsm; run sc.tl.draw_graph first."
    )


# Saves the force-directed tree graph colored by cell type.
def plot_trajectory_graphs(
    adata,
    name,
    output_dir,
    annotation_col,
    show_inline=False,
    param_label=None,
):
    os.makedirs(output_dir, exist_ok=True)
    title_suffix = f"\n{param_label}" if param_label else ""

    fig, ax = plt.subplots(figsize=PLOT_FIGSIZE)
    scf.pl.graph(adata, color_cells=annotation_col, ax=ax, show=False)
    ax.set_title(f"{ax.get_title()}{title_suffix}")
    fig.savefig(os.path.join(output_dir, f"{name}_annotation.pdf"), dpi=PLOT_DPI, facecolor="white")

    plt.close(fig)

    print(f"✓ Trajectory graph saved to {output_dir}")


def export_pseudotime_table(adata, name, output_dir, annotation_col=None):
    os.makedirs(output_dir, exist_ok=True)
    cell_time = adata.obs[["t"]].copy()
    cell_time.index.name = "cell_id"
    if annotation_col is not None and annotation_col in adata.obs.columns:
        cell_time[annotation_col] = adata.obs[annotation_col]
    out_file = os.path.join(output_dir, f"{name}_pseudotime_by_cell.tsv")
    cell_time.to_csv(out_file, sep="\t")
    print(f"Pseudotime table saved: {out_file}")
    return out_file


# Overlay the principal graph and pseudotime trajectory on the FA layout.
# Optional diagnostic (not called by the main workflow): mean of numeric obs columns
# and of selected genes along pseudotime, in n_bins equal-size bins, plus the share of
# each group of group_col (for example "condition"). Use it to check what the pseudotime
# axis mixes with, and how well supported the chosen root is.
def summarize_along_pseudotime(adata, columns=("nFeature_RNA",), genes=None, group_col=None, n_bins=10):
    t   = adata.obs["t"].values
    out = pd.DataFrame({"bin": pd.qcut(t, n_bins, labels=False)}, index=adata.obs_names)
    for c in columns:
        out[c] = adata.obs[c].values
    layer = adata.layers["logcounts"] if "logcounts" in adata.layers else adata.X
    for name, gid in (genes or {}).items():
        col = layer[:, list(adata.var_names).index(gid)]
        out[name] = np.asarray(col.todense() if hasattr(col, "todense") else col).ravel()
    if group_col is not None:
        for grp in adata.obs[group_col].astype(str).unique():
            out[f"{group_col}={grp}"] = (adata.obs[group_col].astype(str).values == grp).astype(float)
    res = out.groupby("bin").mean().round(2)
    res.insert(0, "n", out.groupby("bin").size())
    return res


# Marks the root cell with a red dot, 30% larger than the cells already drawn on the
# axis, on a plot that uses a draw_graph basis.
def mark_root_cell(ax, adata, basis):
    if "is_root" not in adata.obs.columns:
        return
    idx = np.where(adata.obs["is_root"].values)[0]
    key = f"X_{basis}"
    if len(idx) == 0 or key not in adata.obsm:
        return
    xy = np.asarray(adata.obsm[key])
    sizes = [np.median(c.get_sizes()) for c in ax.collections
             if hasattr(c, "get_sizes") and len(c.get_sizes()) > 0]
    base  = float(np.median(sizes)) if sizes else 26
    ax.scatter(xy[idx, 0], xy[idx, 1], s=1.3 * base, marker="o", c="red",
               edgecolor="black", linewidth=0.5, zorder=10)


def plot_pseudotime_trajectory(
    adata,
    name,
    output_dir,
    annotation_col,
    show_inline=False,
    param_label=None,
):
    os.makedirs(output_dir, exist_ok=True)

    basis = draw_graph_basis(adata)
    fig, ax = plt.subplots(figsize=PLOT_FIGSIZE)
    scf.pl.graph(
        adata,
        basis       = basis,
        color_cells = annotation_col,
        ax          = ax,
        show        = False,
    )
    scf.pl.trajectory(
        adata,
        color_seg  = "t",
        basis      = basis,
        frameon    = False,
        s          = 50,
        scale_path = 0.6,
        ax         = ax,
        show       = False,
    )
    mark_root_cell(ax, adata, basis)
    title = f"{name} — pseudotime trajectory (red dot = root cell)"
    if param_label:
        title += f"\n{param_label}"
    ax.set_title(title)
    save_plot_18x18(fig, os.path.join(output_dir, f"{name}_pseudotime_trajectory.pdf"))

    if show_inline:
        pass
    plt.close(fig)
    print(f"Pseudotime trajectory saved to {output_dir}")


# QC plot: show the automatically selected root cell as a red point on FA1/FA2.
def plot_root_cell(adata, name, output_dir, show_inline=False, param_label=None):
    os.makedirs(output_dir, exist_ok=True)
    if "is_root" not in adata.obs.columns:
        print("Root-cell plot skipped: is_root not found.")
        return None

    root_cells = adata.obs_names[adata.obs["is_root"]]
    if len(root_cells) == 0:
        print("Root-cell plot skipped: no root cell marked.")
        return None

    root_cell = root_cells[0]
    plot_data = adata.copy()
    plot_data.obs["root_cell"] = np.where(plot_data.obs["is_root"], "Root cell", "Other cells")
    plot_data = plot_data[plot_data.obs.sort_values("root_cell").index, :]
    plot_data.uns["root_cell_colors"] = ["lightgray", "red"]

    title = f"Root cell: {root_cell}"
    if param_label:
        title += f"\n{param_label}"

    fig = sc.pl.embedding(
        plot_data,
        basis              = draw_graph_basis(plot_data),
        color              = "root_cell",
        legend_loc         = "right margin",
        title              = title,
        frameon            = False,
        size               = 26,
        alpha              = 0.9,
        show               = False,
        return_fig         = True,
    )
    mark_root_cell(fig.axes[0], plot_data, draw_graph_basis(plot_data))
    out_file = os.path.join(output_dir, f"{name}_root_cell.pdf")
    save_plot_18x18(fig, out_file)

    if show_inline:
        pass
    plt.close(fig)
    print(f"Root-cell plot saved: {out_file}")
    return root_cell


# =============================================================================
# STEP 27 — run_trajectory_runs (function the notebook calls)
# =============================================================================
# Runs one or more trajectory parameter sets, saves each result in a separate
# folder, and returns the selected trajectory for downstream steps.
def run_trajectory_runs(
    adata,
    clusters,
    root_cluster,
    annotation_col,
    output_base_dir,
    runs,
    selected_run=None,
    show_inline=True,
    root_method="auto",
    root_score="nFeature_RNA",
    root_lowest=False,
    root_tip=None,
):
    trajectory_runs = {}

    for run in runs:
        run_name = run.get("name", run.get("id"))
        if not run_name:
            raise ValueError("Each trajectory run needs a 'name'.")

        run_dir = os.path.join(output_base_dir, "trajectory", run_name)
        os.makedirs(run_dir, exist_ok=True)

        nodes       = run.get("nodes", 50)
        sigma       = run.get("sigma", 0.2)
        ppt_lambda  = run.get("lambda", run.get("ppt_lambda", 60))
        n_eigs      = run.get("eigs", run.get("n_eigs", 20))
        n_neighbors = run.get("n_neighbors", 50)
        seed        = run.get("seed", 3)
        r_method    = run.get("root_method", root_method)
        r_score     = run.get("root_score", root_score)
        r_lowest    = run.get("root_lowest", root_lowest)
        r_tip       = run.get("root_tip", root_tip)

        param_label = (
            f"nodes={nodes}, sigma={sigma}, lambda={ppt_lambda}, eigs={n_eigs}"
        )

        print(f"\nRunning: {run_name}")
        print(f"Results folder: {run_dir}")

        adata_run = build_pseudotime_trajectory(
            adata          = adata,
            clusters       = clusters,
            root_cluster   = root_cluster,
            annotation_col = annotation_col,
            nodes          = nodes,
            sigma          = sigma,
            ppt_lambda     = ppt_lambda,
            n_eigs         = n_eigs,
            n_components   = 50,
            n_neighbors    = n_neighbors,
            seed           = seed,
            root_method    = r_method,
            root_score     = r_score,
            root_lowest    = r_lowest,
            root_tip       = r_tip,
        )

        params = dict(run)
        params.update({"root_cluster": root_cluster, "root_method": r_method,
                       "root_score": r_score, "root_lowest": r_lowest})
        adata_run.write_h5ad(os.path.join(run_dir, f"{run_name}_trajectory.h5ad"))
        pd.Series(params).to_csv(os.path.join(run_dir, "parameters.tsv"), sep="\t", header=False)

        plot_trajectory_graphs(
            adata          = adata_run,
            name           = run_name,
            output_dir     = run_dir,
            annotation_col = annotation_col,
            show_inline    = show_inline,
            param_label    = param_label,
        )

        export_pseudotime_table(
            adata          = adata_run,
            name           = run_name,
            output_dir     = run_dir,
            annotation_col = annotation_col,
        )

        plot_pseudotime_trajectory(
            adata          = adata_run,
            name           = run_name,
            output_dir     = run_dir,
            annotation_col = annotation_col,
            show_inline    = show_inline,
            param_label    = param_label,
        )

        plot_root_cell(
            adata       = adata_run,
            name        = run_name,
            output_dir  = run_dir,
            show_inline = show_inline,
            param_label = param_label,
        )

        trajectory_runs[run_name] = {
            "adata": adata_run,
            "output_dir": run_dir,
            "params": params,
        }

    if selected_run is None:
        selected_run = runs[-1].get("name", runs[-1].get("id"))

    if selected_run not in trajectory_runs:
        available = list(trajectory_runs)
        raise ValueError(f"selected_run '{selected_run}' not found. Available: {available}")

    selected = trajectory_runs[selected_run]
    print(f"\nSTEP 27 COMPLETE")
    print(f"Selected run: {selected_run}")
    print(f"Selected folder: {selected['output_dir']}")

    return selected["adata"], selected["output_dir"], trajectory_runs


# =============================================================================
# STEP 28 — Plot genes on trajectory
# =============================================================================
# Draws one or more genes on the force-directed trajectory graph and saves them in
# <run_dir>/gene_plots/trajectory_gene_<AGI>.pdf (18x18 in), where run_dir is the
# trajectory run folder (selected_trajectory_dir). The file is named after the
# plotted genes, as in Section 8: one gene -> trajectory_gene_AT4G21750.pdf, several
# -> trajectory_gene_AT3G24140_AT4G21750.pdf. Puts the expression in .X if needed,
# creates the folder, closes the figure, and returns the saved path.
def plot_genes_on_trajectory(adata, genes, run_dir, titles=None):
    ensure_expression_in_X(adata)
    out_dir = os.path.join(run_dir, "gene_plots")
    os.makedirs(out_dir, exist_ok=True)
    fig = sc.pl.draw_graph(adata, color=genes, title=titles if titles is not None else genes,
                           show=False, return_fig=True)
    path = os.path.join(out_dir, f"trajectory_gene_{'_'.join(genes)}.pdf")
    save_plot_18x18(fig, path)
    plt.close(fig)
    print(f"✓ Gene plots saved: {path}")
    return path


def strip_segment_frames(ax):
    xlim = ax.get_xlim()
    full_width = xlim[1] - xlim[0]
    for p in list(ax.patches):
        if isinstance(p, mpatches.Rectangle):
            ec          = p.get_edgecolor()
            is_black    = ec[:3] == (0.0, 0.0, 0.0)
            is_unfilled = (not p.get_fill()) or p.get_facecolor()[3] == 0
            if is_black and is_unfilled and p.get_width() < full_width * 0.95:
                p.remove()


# Trends plots get their own canvas (24x18, wider than tall) instead of the
# generic PLOT_FIGSIZE (18x18) used by other plot types — the dendrogram
# panel looks cramped at a square 18x18 once the tree is a single segment.
def save_trends_plot(axes_list, pdf_path):
    for ax in axes_list:
        strip_segment_frames(ax)
    fig = axes_list[0].get_figure()
    fig.set_size_inches(24, 18, forward=True)
    fig.subplots_adjust(left=0.04, right=0.85, bottom=0.06, top=0.92)
    fig.savefig(pdf_path, dpi=PLOT_DPI, facecolor="white")
    plt.close(fig)


def run_step29_gene_trends(
    adata,
    run_dir,
    name        = "trajectory",
    custom_genes = None,
    top_n        = 50,
    ordering     = "max",
    n_jobs       = 4,
    spline_df    = 5,
    A_cut        = 0.3,
    p_val_cut    = 0.001,
):
    output_dir = os.path.join(run_dir, "gene_trends")
    os.makedirs(output_dir, exist_ok=True)

    ensure_expression_in_X(adata)

    # spline_df: degrees of freedom of the per-segment GAM smoothing spline
    # used by test_association. Lower it (e.g. 3) if you hit "A term has
    # fewer unique covariate combinations than specified maximum degrees of
    # freedom" — that means some tree segment has too few distinct pseudotime
    # values for the default of 5.
    #
    # A_cut / p_val_cut: significance is redefined from scFates's own
    # combined fdr/A/st cutoffs (very conservative on small datasets) to a
    # lower amplitude cutoff plus the raw (non-FDR-corrected) p-value — the
    # global/whole-tree version of the technique used per-branch in pse.py.
    scf.tl.dendrogram(adata)
    scf.tl.test_association(adata, n_jobs=n_jobs, spline_df=spline_df, A_cut=A_cut)
    adata.var["signi"] = adata.var["p_val"] < p_val_cut
    significant_genes = adata.var_names[adata.var["signi"]].tolist()

    # Custom genes are force-included in the fit even if not statistically
    # significant, so a gene you explicitly asked for never silently
    # disappears from either heatmap.
    custom_genes  = custom_genes or []
    genes_present = [g for g in custom_genes if g in adata.var_names]
    missing       = sorted(set(custom_genes) - set(genes_present))
    if missing:
        print(f"custom selection: genes ignored (not found in adata): {', '.join(missing)}")
    fit_features = list(dict.fromkeys(significant_genes + genes_present))

    # Fit on a copy so the caller's adata is never subsetted. Both heatmaps
    # below reuse this single fit — rows = significant + custom genes, only
    # the highlighted/labeled subset differs between the two plots.
    adata_fitted = adata.copy()
    scf.tl.fit(adata_fitted, features=fit_features, n_jobs=n_jobs)

    top_genes = adata_fitted.var.sort_values("A", ascending=False).head(top_n).index.tolist()
    axes_list = scf.pl.trends(
        adata_fitted,
        highlight_features = top_genes,
        style               = "italic",
        add_outline         = True,
        basis               = "dendro",
        show_segs           = False,
        fontsize            = 18,
        ordering            = ordering,
        show                = False,
        title               = f"{name} — top {top_n} variable genes",
    )
    top_pdf = os.path.join(output_dir, f"step29_gene_trends_top{top_n}_FINAL_24x18.pdf")
    save_trends_plot(axes_list, top_pdf)
    print(f"✓ Top-{top_n} gene trends saved: {top_pdf}")

    custom_pdf = None
    if genes_present:
        axes_list = scf.pl.trends(
            adata_fitted,
            highlight_features = genes_present,
            style               = "italic",
            add_outline         = True,
            basis               = "dendro",
            show_segs           = False,
            fontsize            = 18,
            ordering            = ordering,
            show                = False,
            title               = f"{name} — custom selection",
        )
        custom_pdf = os.path.join(output_dir, "step29_gene_trends_highlight_list_FINAL_24x18.pdf")
        save_trends_plot(axes_list, custom_pdf)
        print(f"✓ Custom selection gene trends saved: {custom_pdf}")

    print("STEP 29 COMPLETE: gene trends saved")
    return top_pdf, custom_pdf
