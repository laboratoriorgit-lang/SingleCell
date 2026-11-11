#!/usr/bin/env python3
"""Make copies of the step scripts that read their values from params.yaml.

The scripts in workflow/ (and the report) contain plain values. This tool writes,
for each of them, a copy in which every editable value is wrapped as
    opt("section.name", <the value written in the script>)
so that params.yaml replaces it when a parameter file is used (PARAMS_FILE), and
the value of the script is kept otherwise. The copies are generated, never edited.

Usage:  python3 workflow/make_param_steps.py [OUTPUT_DIR]     (default: results/params_steps)
It stops with an error if a rule does not match, e.g. after changing the layout of a
line in a step script; update the rule in this file in that case.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "results", "params_steps")
ML = re.M
MLS = re.M | re.S


def ra(var, key):                       # R:  var <- value   [# comment]
    return (rf'^({var}\s*<-\s*)(.*?)(\s+#.*)?$', ML,
            lambda m: f'{m[1]}opt("{key}", {m[2]}){m[3] or ""}')


def pa(var, key):                       # Python:  VAR = value   [# comment]
    return (rf'^({var}\s*=\s*)(.*?)(\s+#.*)?$', ML,
            lambda m: f'{m[1]}opt("{key}", {m[2]}){m[3] or ""}')


def after(block_start, line):           # keep a multi-line block, add a line after it
    return (block_start, MLS, lambda m: m[0] + "\n" + line)


R_STEP1 = [
    after(r'^samples <- list\(.*?\n\)(?=\n)', 'samples <- opt("samples", samples)'),
    after(r'^colors <- c\(.*?\n\)(?=\n)', 'colors <- unlist(opt("colors", colors))'),
    ra("mt_pattern", "qc.mt_pattern"),
    (r'^  seurat_list_raw, min_features = (\S+), max_mt = (\S+), run_doubletfinder = (\w+)$', ML,
     lambda m: ('  seurat_list_raw,\n'
                f'  min_features      = opt("qc.min_features", {m[1]}),\n'
                f'  max_mt            = opt("qc.max_mt", {m[2]}),\n'
                f'  run_doubletfinder = opt("qc.run_doubletfinder", {m[3]})')),
    (r'^ath_sc <- merge_and_preprocess\(seurat_list, nfeatures = (\d+), npcs = (\d+)\)$', ML,
     lambda m: ('ath_sc <- merge_and_preprocess(\n'
                f'  seurat_list, nfeatures = opt("preprocessing.nfeatures", {m[1]}),\n'
                f'  npcs = opt("preprocessing.npcs", {m[2]})\n)')),
    (r'^ath_sc <- run_harmony_umap\(ath_sc, group_by = ("[^"]*"), dims = 1:(\d+)\)$', ML,
     lambda m: ('ath_sc <- run_harmony_umap(\n'
                f'  ath_sc, group_by = opt("integration.group_by", {m[1]}),\n'
                f'  dims = 1:opt("integration.n_dims", {m[2]})\n)')),
    ra("resolutions_test", "clustering.resolutions_test"),
    ra("cluster_resolution", "clustering.cluster_resolution"),
    (r'^(cluster_dims\s*<-\s*1:)(\d+)(.*)$', ML,
     lambda m: f'{m[1]}opt("clustering.n_dims", {m[2]}){m[3]}'),
    (r'^  ath_sc, resolution = cluster_resolution, dims = cluster_dims, k_param = (\d+)$', ML,
     lambda m: ('  ath_sc, resolution = cluster_resolution, dims = cluster_dims,\n'
                f'  k_param = opt("clustering.k_param", {m[1]})')),
    (r'^(biblio_marks_file <- file\.path\(DATA_DIR, )"([^"]+)", "([^"]+)"\)$', ML,
     lambda m: f'{m[1]}opt("annotation.marker_file", "{m[2]}/{m[3]}"))'),
    ra("gene", "expression.gene"),
    (r'^saveRDS\(ath_sc, file\.path\(dir_objects, "(ath_sc_annotated)\.rds"\)\)$', ML,
     lambda m: (f'annotated_object <- opt("outputs.annotated_object", "{m[1]}")   # file name, no extension\n'
                'saveRDS(ath_sc, file.path(dir_objects, paste0(annotated_object, ".rds")))')),
    (r'^  file\.path\(dir_objects, "ath_sc_annotated\.h5ad"\)$', ML,
     lambda m: '  file.path(dir_objects, paste0(annotated_object, ".h5ad"))'),
]

R_STEP2 = [
    (r'^(biblio_marks_file <- file\.path\(DATA_DIR, )"([^"]+)", "([^"]+)"\)$', ML,
     lambda m: f'{m[1]}opt("annotation.marker_file", "{m[2]}/{m[3]}"))'),
    (r'^ath_sc <- readRDS\(file\.path\(dir_objects, "(\w+)\.rds"\)\)$', ML,
     lambda m: f'ath_sc <- readRDS(file.path(\n  dir_objects, paste0(opt("outputs.annotated_object", "{m[1]}"), ".rds")\n))'),
    ra("pseudobulk_annot_col", "pseudobulk.annot_col"),
    ra("pseudobulk_conditions", "pseudobulk.conditions"),
    ra("n_pseudoreps", "pseudobulk.n_pseudoreps"),
    after(r'^comparisons <- list\(.*?\n\)(?=\n)', 'comparisons <- opt("de.comparisons", comparisons)'),
    ra("cell_types_to_analyze", "de.cell_types"),
    ra("volcano_tag", "de.volcano_tag"),
    ra("padj_cut", "de.padj_cut"),
    ra("lfc_cut", "de.lfc_cut"),
    ra("go_space", "go.space"),
    (r'^(go_orgdb\s*<-\s*)([\w.]+)$', ML, lambda m: f'{m[1]}get(opt("go.orgdb", "{m[2]}"))'),
    ra("go_keytype", "go.keytype"),
    ra("padj_cutoff", "go.padj_cutoff"),
    ra("heatmap_limits", "heatmap.limits"),
    (r'^(  tf_list_path    = file\.path\(DATA_DIR, )"([^"]+)", "([^"]+)"\),$', ML,
     lambda m: f'{m[1]}opt("network.tf_list", "{m[2]}/{m[3]}")),'),
] + [
    (rf'^(  {n}\s*=\s*)([^,\s#]+)(.*)$', ML, (lambda n: lambda m: f'{m[1]}opt("network.{n}", {m[2]}){m[3]}')(n))
    for n in ("n_metacells", "min_module_size", "deep_split", "soft_power", "tom_threshold", "n_hub_label")
]

R_SUPP = [
    after(r'^grouping <- c\(.*?\n\)(?=\n)', 'grouping <- unlist(opt("supplementary.grouping", grouping))'),
    (r'^  ath_sc, cluster_id = ("[^"]*"),$', ML,
     lambda m: f'  ath_sc, cluster_id = opt("supplementary.inspect_cluster", {m[1]}),'),
    after(r'^reassign <- list\(.*?\n\)(?=\n)', 'reassign <- lapply(opt("supplementary.reassign", reassign), unlist)'),
]

PY_STEP3 = [
    (r'^INPUT_H5AD     = "[^"]*/(\w+)\.h5ad"$', ML,
     lambda m: f'INPUT_H5AD     = os.path.join(base_dir, "objects", opt("outputs.annotated_object", "{m[1]}") + ".h5ad")'),
    pa("ANNOTATION_COL", "pseudotime.annotation_col"),
    pa("N_JOBS", "pseudotime.n_jobs"),
    pa("TRAJECTORY_CLUSTERS", "pseudotime.trajectory_clusters"),
    pa("ROOT_CLUSTER", "pseudotime.root_cluster"),
    after(r'^TRAJECTORY_RUNS = \[\n.*?\n\](?=\n)',
          '# params.yaml replaces the list of runs only if it is used\n'
          'TRAJECTORY_RUNS = [trajectory_run(**r) for r in opt("pseudotime.runs", [])] or TRAJECTORY_RUNS'),
    (r'^(    genes    = )(\[[^\]]*\])(,.*)$', ML, lambda m: f'{m[1]}opt("pseudotime.plot_genes", {m[2]}){m[3]}'),
    (r'^(    titles   = )(\[[^\]]*\])(,.*)$', ML, lambda m: f'{m[1]}opt("pseudotime.plot_gene_titles", {m[2]}){m[3]}'),
    (r'^(    filename = )("[^"]*")(,.*)$', ML, lambda m: f'{m[1]}opt("pseudotime.plot_filename", {m[2]}){m[3]}'),
    pa("SPLINE_DF", "pseudotime.spline_df"),
    pa("TOP_N", "pseudotime.top_n"),
    pa("HIGHLIGHT_GENES", "pseudotime.highlight_genes"),
    pa("ORDERING", "pseudotime.ordering"),
    (r'^final_h5ad = os\.path\.join\(base_dir, "objects", "(\w+)\.h5ad"\)$', ML,
     lambda m: f'final_h5ad = os.path.join(base_dir, "objects", opt("outputs.pseudotime_object", "{m[1]}") + ".h5ad")'),
]

JOBS = [("step1_singlecell.R", R_STEP1), ("step1_supplementary.R", R_SUPP),
        ("step2_degs.R", R_STEP2), ("step3_pseudotime.py", PY_STEP3)]


def main():
    os.makedirs(OUT, exist_ok=True)
    failed = False
    for name, rules in JOBS:
        text = open(os.path.join(HERE, name), encoding="utf-8").read()
        for pat, flags, repl in rules:
            text, n = re.subn(pat, repl, text, flags=flags)
            if n != 1:
                failed = True
                print(f"ERROR {name}: rule matched {n} times (expected 1): {pat[:70]}", file=sys.stderr)
        header = (f"# GENERATED by workflow/make_param_steps.py from workflow/{name} - do not edit.\n"
                  "# Every editable value is opt(\"section.name\", value): params.yaml replaces it when PARAMS_FILE is set.\n")
        with open(os.path.join(OUT, name), "w", encoding="utf-8") as fh:
            fh.write(header + text)
    if failed:
        sys.exit(1)
    print(f"Parameter-driven copies written to {os.path.abspath(OUT)}")


if __name__ == "__main__":
    main()
