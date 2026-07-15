# Cell Ranger launcher

This folder contains a lightweight `cellranger` launcher so users can run:

```bash
export PATH="$PROJECT_DIR/tools/cellranger:$PATH"
cellranger --version
```

The actual 10x Genomics Cell Ranger binary is not committed to the repository.
The launcher searches for Cell Ranger in this order:

1. `$CELLRANGER_BIN`
2. `$CELLRANGER_HOME/cellranger`
3. `tools/cellranger/cellranger-*/cellranger`
4. the current server path `/home/projects3/jm/sc/cellranger-7.1.0/cellranger`
5. another `cellranger` already available in `PATH`

For a new machine, download Cell Ranger from 10x Genomics, extract it, and place
its folder here, for example:

```text
tools/cellranger/cellranger-7.1.0/cellranger
```

Then add this folder to `PATH` before running the Cell Ranger step.
