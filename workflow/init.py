# =============================================================================
# init.py - one-call initialization of the Python session
# =============================================================================
# Usage:
#   exec(open("/workspace/workflow/init.py").read(), globals())
#   init_pipeline(data_dir="/workspace/.")
# =============================================================================
import os


def init_pipeline(data_dir="/workspace/.", pipeline_dir=None):
    """Load the Python packages and helper functions and create results/09_pseudotime.

    Creates these variables in the caller's global namespace: DATA_DIR,
    PIPELINE_DIR, base_dir and dir_pseudotime.
    """
    if pipeline_dir is None:
        pipeline_dir = os.path.join(data_dir, "workflow")
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"data_dir does not exist: {data_dir}")

    g = globals()  # this function is defined via exec(..., globals()), so these are the caller's globals
    exec(open(os.path.join(pipeline_dir, "load_libraries_python.py")).read(), g)
    exec(open(os.path.join(pipeline_dir, "ScRNA_Pseudotime_Functions.py")).read(), g)

    base_dir = os.path.join(data_dir, "results")
    dir_pseudotime = os.path.join(base_dir, "09_pseudotime")
    os.makedirs(dir_pseudotime, exist_ok=True)

    g.update(DATA_DIR=data_dir, PIPELINE_DIR=pipeline_dir,
             base_dir=base_dir, dir_pseudotime=dir_pseudotime)
