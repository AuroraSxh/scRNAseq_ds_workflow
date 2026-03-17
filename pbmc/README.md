# PBMC Downstream Workflow

This is a compact downstream scRNA-seq workflow for 10x PBMC matrices using Python and Scanpy.

## Included Files

```text
pbmc/
├── README.md
├── code/
│   ├── analyze_pbmc.py
│   └── gene_expression_summary.py
├── workflows/
│   └── run_pbmc_analysis.sh
└── expected_results/
```

## What The Workflow Does

`analyze_pbmc.py` performs:
- matrix loading
- QC metric calculation
- filtering
- normalization and HVG selection
- PCA
- Leiden clustering
- UMAP / tSNE
- marker extraction
- heuristic cell-type annotation

`gene_expression_summary.py` performs:
- per-gene expression summary over annotated cell types
- nonparametric comparison
- UMAP / tSNE gene-level visualization

## Expected Input

By default the main script expects a 10x matrix directory like:

```text
filtered_gene_bc_matrices/hg19/
├── barcodes.tsv
├── genes.tsv or features.tsv
└── matrix.mtx
```

You can override the input path with script arguments.

## How To Run

```bash
./workflows/run_pbmc_analysis.sh
```

Optional gene-summary example:

```bash
python code/gene_expression_summary.py PDCD1 IL7R
```

## Expected Outputs

See:
- `expected_results/README.md`
