# =============================================================================
# init.py - one-call initialization of the Python session
# =============================================================================
# Usage:
#   exec(open("/workspace/workflow/init.py").read(), globals())
#   init_pipeline(data_dir="/workspace/.")
#
# Optional: PARAMS_FILE=/workspace/params.yaml (environment variable, or the
# params_file argument) is used by the parameter-driven copies of the steps made by
# make_param_steps.py, in which values are written opt("section.name", value).
# The step scripts of workflow/ do not call opt().
# =============================================================================
import os

P = {}   # parameters read from params.yaml (empty if no file is used)


def opt(key, default=None):
    """Value from params.yaml ("section.name"), or the default if there is no
    params file or the key is missing or null."""
    v = P
    for k in key.split("."):
        if not isinstance(v, dict) or k not in v:
            return default
        v = v[k]
    return default if v is None else v


def init_pipeline(data_dir="/workspace/.", params_file=None, pipeline_dir=None):
    """Load the Python packages and helper functions and create results/09_pseudotime.

    Creates these variables in the caller's global namespace: P, DATA_DIR,
    PIPELINE_DIR, base_dir and dir_pseudotime.
    """
    g = globals()  # this function is defined via exec(..., globals()), so these are the caller's globals

    params_file = params_file or os.environ.get("PARAMS_FILE", "")
    P = {}
    if params_file:
        try:
            import yaml
        except ImportError as e:
            raise ImportError("PyYAML is required to read params.yaml: pip install pyyaml") from e
        if not os.path.isfile(params_file):
            raise FileNotFoundError(f"params file not found: {params_file}")
        with open(params_file) as fh:
            P = yaml.safe_load(fh) or {}
        print(f"Using parameters from {params_file}")
    g["P"] = P

    data_dir = opt("project.data_dir", data_dir)
    if pipeline_dir is None:
        pipeline_dir = os.path.join(data_dir, "workflow")
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"data_dir does not exist: {data_dir}")

    exec(open(os.path.join(pipeline_dir, "load_libraries_python.py")).read(), g)
    exec(open(os.path.join(pipeline_dir, "ScRNA_Pseudotime_Functions.py")).read(), g)

    base_dir = os.path.join(data_dir, opt("project.results_dir", "results"))
    dir_pseudotime = os.path.join(base_dir, "09_pseudotime")
    os.makedirs(dir_pseudotime, exist_ok=True)

    g.update(DATA_DIR=data_dir, PIPELINE_DIR=pipeline_dir,
             base_dir=base_dir, dir_pseudotime=dir_pseudotime)
