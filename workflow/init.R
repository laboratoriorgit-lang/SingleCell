# =============================================================================
# init.R - one-call initialization of the R session
# =============================================================================
# Usage:
#   source("/workspace/workflow/init.R")
#   init_pipeline(data_dir = "/workspace/.", seed = 1807)
# =============================================================================

#' Initialize the R Session for the scRNA-seq Pipeline
#'
#' Loads the R packages and helper functions, fixes the random seed, moves R
#' into the project folder and creates the results folders.
#'
#' Creates these variables in the global environment: DATA_DIR, PIPELINE_DIR,
#' base_dir, output_dir and dir_01 ... dir_08 + dir_objects.
#'
#' @param data_dir     Project root. Inputs are read relative to it and results
#'                     are written to data_dir/results/<step>/.
#' @param seed         Random seed (default 1807).
#' @param pipeline_dir Folder with the workflow scripts (default data_dir/workflow).
#' @return Invisible NULL.
init_pipeline <- function(data_dir     = "/workspace/.",
                          seed         = 1807,
                          pipeline_dir = file.path(data_dir, "workflow")) {

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

  base_dir <- file.path(data_dir, "results")
  list2env(create_pipeline_dirs(base_dir), envir = .GlobalEnv)  # dir_01 ... dir_objects

  assign("DATA_DIR",     data_dir,     envir = .GlobalEnv)
  assign("PIPELINE_DIR", pipeline_dir, envir = .GlobalEnv)
  assign("base_dir",     base_dir,     envir = .GlobalEnv)
  assign("output_dir",   base_dir,     envir = .GlobalEnv)  # reassigned in each section

  invisible(NULL)
}
