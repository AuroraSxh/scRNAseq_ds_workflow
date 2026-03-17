# scRNAseq Downstream Workflow

This repository contains downstream single-cell RNA-seq workflows without bundled data.

It currently includes two downstream analysis tracks:

- `pbmc/` for a compact PBMC workflow built around Scanpy
- `crc/` for a larger colorectal cancer integration and subtype analysis workspace

The repository intentionally excludes:
- raw data
- processed AnnData or Seurat objects
- generated figures and tables
- local virtual environments
- local R package libraries

It keeps only:
- code
- workflow entry points
- config
- README files
- expected output structure

## Repository Layout

```text
.
├── README.md
├── pbmc/
└── crc/
```

## PBMC Workflow

Directory:
- `pbmc/`

Contents:
- `code/` for analysis scripts
- `workflows/` for runnable wrappers
- `expected_results/` for result structure only

Use this workflow when you want a minimal downstream scRNA-seq example with:
- QC
- normalization
- PCA / UMAP / tSNE
- Leiden clustering
- marker extraction
- heuristic cell-type annotation
- optional single-gene summary

## CRC Workflow

Directory:
- `crc/`

Contents:
- `config/`
- `src/`
- `workflows/`
- `requirements.txt`
- `expected_results/`

Use this workflow when you want a more complete downstream workspace with:
- multi-dataset integration
- major cell-type annotation
- CD8 T-cell refinement
- macrophage / myeloid refinement
- CellChat-ready downstream communication analysis

## No Data Policy

This repository is prepared for code and workflow sharing only.

Excluded on purpose:
- `data/`
- public dataset downloads
- `.venv/`
- `.r_libs/`
- generated `results/`
- generated `analysis_outputs/`

## Expected Outputs

Each workflow has its own output skeleton:

- `pbmc/expected_results/`
- `crc/expected_results/`

These describe the intended directory layout and major output artifact types without shipping any real result files.
