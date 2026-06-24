---
title: "Computational Methods for Single-Cell RNA-Seq Analysis of *Arabidopsis thaliana*"
subtitle: "Part 1 — Computational Environment and Quick Start"
author: "SingleCell Pipeline"
date: "2026-06-24"
geometry: margin=1in
fontsize: 11pt
toc: false
numbersections: true
---

# Materials and Methods

## Materials

The pipeline source code is available at
`https://github.com/cliford2001/SingleCell`:

```bash
git clone https://github.com/cliford2001/SingleCell.git
cd SingleCell
```

The repository is meant to sit inside a working directory that also holds
the raw data and the pipeline's output folder — for instance
`~/your-workdir/SingleCell`, with `~/your-workdir/resultados` alongside it. That
parent directory, not the repository alone, is what gets mounted into the
container, so that code, input matrices, and results are all visible
together at `/workspace`.

Docker packages R, Python, and every pipeline dependency into one portable
**image** (`matigara/scrnaseq:latest`), so the analysis runs identically on
any machine with Docker installed. Two equivalent ways to start it:

### Option 1 — single command

```bash
docker run -it --rm \
  -v ~/your-workdir:/workspace \
  matigara/scrnaseq:latest /bin/bash
```

Pulls the image if missing and opens a shell. Replace `~/your-workdir` with the
actual path to that working directory if it lives somewhere else.

### Option 2 — the repository's `docker-compose.yml`

Run as-is from the repository root:

```bash
docker compose pull
docker compose run --rm r bash
```

`pull` only reads that file and downloads the image it names — unrelated to
git, which just placed the file on disk. Same `~/your-workdir` → `/workspace`
mount as Option 1.

Either way, `/workspace` **is** `~/your-workdir` on the host: files written
there persist after the container exits.

Once inside the shell, Chapter 1-2 commands are typed into R, and Chapter 3
commands into Python:

```r
R
```

```python
python3
```

<div class="pagebreak"></div>

## Methods

### Chapter 1 — Single-cell processing

Running Option 1 opens a self-contained terminal inside the container,
separate from the host shell it was launched from:

```bash
docker run -it --rm \
  -v ~/your-workdir:/workspace \
  matigara/scrnaseq:latest /bin/bash
```

`R` is then started from that terminal, on its own:

```r
R
```

Chapter 1 (`capitulo1_single_cell.R`) runs entirely from this session.

#### Initialization

Paths are set once at the top of the script: `PIPELINE_DIR` points to the
helper scripts, and `DATA_DIR`/`base_dir` define where input data are read
from and where every result file is written. Three helper scripts are then
sourced — package loading, custom Seurat plotting utilities, and the core
analysis functions used throughout the chapter — before fixing the random
seed (1807) for reproducibility, disabling `Seurat.allow.s4` for Seurat 5
compatibility, and setting the working directory to `DATA_DIR`.

```r
PIPELINE_DIR <- "/workspace/ScRNASeq-Docker/workflow"
DATA_DIR     <- "/workspace/."
base_dir     <- file.path(DATA_DIR, "resultados")

source(file.path(PIPELINE_DIR, "load_libraries.R"))
source(file.path(PIPELINE_DIR, "custom_seurat.R"))
source(file.path(PIPELINE_DIR, "ScRNA_Analysis_Functions.R"))

set.seed(1807)
options(Seurat.allow.s4 = FALSE)
setwd(DATA_DIR)
```

#### Sample definition and data loading

Three samples spanning a nitrogen-availability gradient were analyzed —
`Sample_0N`, `Sample_05N`, `Sample_5N` (conditions 0N, 0.5N, 5N) — each
loaded from its CellRanger `filtered_feature_bc_matrix` output. For every
cell, the mitochondrial (`^ATMG`) and chloroplast (`^ATCG`) read fraction
was computed alongside the standard `nFeature_RNA`/`nCount_RNA` metrics, to
guide the filtering thresholds applied in the next step.

```r
samples <- list(
  list(file = "cellranger/Sample_0N/outs/filtered_feature_bc_matrix",
       label = "Sample_0N",  condition = "0N"),
  list(file = "cellranger/Sample_05N/outs/filtered_feature_bc_matrix",
       label = "Sample_05N", condition = "0.5N"),
  list(file = "cellranger/Sample_5N/outs/filtered_feature_bc_matrix",
       label = "Sample_5N",  condition = "5N")
)

colors <- c("Sample_0N" = "#66c2a5", "Sample_05N" = "#41ae76", "Sample_5N" = "#fc8d62")

seurat_list_raw <- load_seurat_samples(samples, DATA_DIR,
                                        mt_pattern = "^ATMG",
                                        cp_pattern = "^ATCG")
```
