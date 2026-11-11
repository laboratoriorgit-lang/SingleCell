# =============================================================================
# init.R - one-call initialization of the R session
# =============================================================================
# Usage:
#   source("/workspace/workflow/init.R")
#   init_pipeline(data_dir = "/workspace/.", seed = 1807)
#
# Optional: PARAMS_FILE=/workspace/params.yaml (environment variable, or the
# params_file argument) is used by the parameter-driven copies of the steps made by
# make_param_steps.py, in which values are written opt("section.name", value).
# The step scripts of workflow/ do not call opt().
# =============================================================================

#' Value from params.yaml, or the default
#'
#' @param key     "section.name", e.g. "clustering.cluster_resolution".
#' @param default Value used when no params file is loaded or the key is
#'                missing or null.
opt <- function(key, default = NULL) {
  v <- tryCatch(
    Reduce(function(x, k) x[[k]], strsplit(key, ".", fixed = TRUE)[[1]],
           get0("P", envir = .GlobalEnv, ifnotfound = list())),
    error = function(e) NULL)
  if (is.null(v)) default else v
}

#' Initialize the R Session for the scRNA-seq Pipeline
#'
#' Loads the R packages and helper functions, fixes the random seed, moves R
#' into the project folder and creates the results folders.
#'
#' Creates these variables in the global environment: P (parameters read from
#' params.yaml; an empty list if no file is used), DATA_DIR, PIPELINE_DIR,
#' base_dir, output_dir and dir_01 ... dir_08 + dir_objects.
#'
#' @param data_dir     Project root. Inputs are read relative to it and results
#'                     are written to data_dir/results/<step>/.
#' @param seed         Random seed (default 1807).
#' @param params_file  Optional params.yaml (default: the PARAMS_FILE variable).
#'                     Its project$data_dir and project$seed replace the two
#'                     arguments above.
#' @param pipeline_dir Folder with the workflow scripts (default data_dir/workflow).
#' @return Invisible NULL.
init_pipeline <- function(data_dir     = "/workspace/.",
                          seed         = 1807,
                          params_file  = Sys.getenv("PARAMS_FILE", ""),
                          pipeline_dir = NULL) {

  P <- list()
  if (nzchar(params_file)) {
    if (!file.exists(params_file))
      stop("params file not found: ", params_file)
    P <- yaml::read_yaml(params_file)
    message("Using parameters from ", params_file)
  }
  assign("P", P, envir = .GlobalEnv)

  data_dir <- opt("project.data_dir", data_dir)
  seed     <- opt("project.seed", seed)
  if (is.null(pipeline_dir)) pipeline_dir <- file.path(data_dir, "workflow")

  if (!dir.exists(data_dir))
    stop("data_dir does not exist: ", data_dir)
  if (!file.exists(file.path(pipeline_dir, "ScRNA_Analysis_Functions.R")))
    stop("ScRNA_Analysis_Functions.R not found in pipeline_dir: ", pipeline_dir)

  # Packages and helper functions (loaded into the global environment)
  source(file.path(pipeline_dir, "load_libraries.R"))
  source(file.path(pipeline_dir, "ScRNA_Analysis_Functions.R"))

  set.seed(seed)
  options(Seurat.allow.s4 = FALSE)            # required for Seurat 5 compatibility
  ComplexHeatmap::ht_opt$message <- FALSE     # silence ComplexHeatmap raster messages
  setwd(data_dir)

  base_dir <- file.path(data_dir, opt("project.results_dir", "results"))
  list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)  # dir_01 ... dir_objects

  assign("DATA_DIR",     data_dir,     envir = .GlobalEnv)
  assign("PIPELINE_DIR", pipeline_dir, envir = .GlobalEnv)
  assign("base_dir",     base_dir,     envir = .GlobalEnv)
  assign("output_dir",   base_dir,     envir = .GlobalEnv)  # reassigned in each section

  invisible(NULL)
}
