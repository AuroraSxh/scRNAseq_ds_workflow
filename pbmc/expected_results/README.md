# Expected Results: PBMC Workflow

This workflow writes outputs under an `analysis_outputs/` directory.

Expected structure:

```text
analysis_outputs/
├── analysis_summary.tsv
├── qc_metrics_before_filtering.tsv
├── qc_metrics_after_filtering.tsv
├── qc_before_filtering.pdf
├── qc_after_filtering.pdf
├── umap_leiden_clusters.pdf
├── tsne_leiden_clusters.pdf
├── umap_cell_type_annotations.pdf
├── tsne_cell_type_annotations.pdf
├── canonical_marker_dotplot.pdf
├── cluster_annotations.tsv
├── cluster_cell_counts.tsv
├── cluster_cell_counts.pdf
├── cluster_markers_top20.tsv
├── cluster_markers_top10_filtered.tsv
├── cell_metadata.tsv
├── pbmc_processed.h5ad
└── Gene_expression/
```

## Gene Expression Subdirectory

For each requested gene, the expected structure is:

```text
analysis_outputs/Gene_expression/<GENE>/
├── expression_by_cell_type.tsv
├── expression_by_cell_type.pdf
├── umap_expression.pdf
├── tsne_expression.pdf
└── kruskal_wallis.txt
```

This repository ships only the structure description, not any real outputs.
