# scRNA-seq Downstream Workflow

Code-only downstream single-cell RNA-seq workflows for public repository sharing.

This repository keeps analysis code, workflow entry points, lightweight configuration, and expected output structure, while excluding datasets, large intermediate objects, and generated results.

## Overview

Included workflow tracks:

- `pbmc/`: compact Scanpy-based PBMC downstream analysis
- `crc/`: colorectal cancer integration, annotation, and communication analysis workspace

Core downstream tasks covered across the repository:

- quality control and filtering
- normalization and dimensionality reduction
- clustering and marker discovery
- heuristic cell-type annotation
- subtype-focused follow-up analyses

## Repository Layout

```text
.
├── README.md
├── LICENSE
├── .gitignore
├── pbmc/
└── crc/
```

## Workflow Modules

### `pbmc/`

This module is intended as a minimal, runnable scRNA-seq downstream example with:

- QC
- normalization
- PCA / UMAP / t-SNE
- Leiden clustering
- marker extraction
- heuristic cell-type annotation
- optional single-gene summaries

### `crc/`

This module is a larger analysis workspace for colorectal cancer datasets, including:

- multi-dataset integration
- major cell-type annotation
- CD8 T-cell refinement
- macrophage or myeloid refinement
- CellChat-ready downstream communication analysis

## Quick Start

Run the workflow from the relevant module directory:

```bash
cd pbmc
./workflows/run_pbmc_analysis.sh
```

or:

```bash
cd crc
./workflows/run_crc_pipeline.sh
```

Refer to each module-level README for environment setup and workflow-specific parameters.

## Repository Policy

This repository is intentionally prepared for code and workflow sharing only.

Included:

- code
- workflow wrappers
- configuration
- documentation
- expected output skeletons

Excluded:

- raw or processed data
- downloaded public datasets
- AnnData or Seurat objects
- generated figures and tables
- local virtual environments
- local R package libraries
- generated `results/` or `analysis_outputs/`

## Expected Results

Each module contains an `expected_results/` directory that documents the intended output structure without shipping real outputs:

- `pbmc/expected_results/`
- `crc/expected_results/`

These folders are designed to make the repository reproducible and reviewable while keeping it lightweight for public distribution.
