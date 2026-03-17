#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

CANONICAL_MARKERS = {
    "T_cells": ["CD3D", "CD3E", "IL7R", "LTB", "IL32"],
    "NK_cells": ["NKG7", "CCL5", "PRF1", "FGFBP2", "GZMB"],
    "B_cells": ["MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA"],
    "Monocytes": ["LYZ", "S100A9", "FCER1G", "TYROBP", "LST1"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a standard scRNA-seq workflow for 10x PBMC test data."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("filtered_gene_bc_matrices/hg19"),
        help="Directory containing barcodes.tsv, genes.tsv/features.tsv and matrix.mtx.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_outputs"),
        help="Directory for figures and analysis tables.",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=200,
        help="Minimum detected genes per cell.",
    )
    parser.add_argument(
        "--max-mt-pct",
        type=float,
        default=20.0,
        help="Maximum mitochondrial percentage per cell.",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.5,
        help="Leiden clustering resolution.",
    )
    return parser.parse_args()


def save_qc_pdf(
    qc_df: pd.DataFrame,
    output_path: Path,
    title_prefix: str,
) -> None:
    sns.set_theme(style="whitegrid")
    with PdfPages(output_path) as pdf:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        for ax, column, title in zip(
            axes,
            ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
            ["UMI counts", "Genes per cell", "Mito %"],
            strict=True,
        ):
            sns.violinplot(y=qc_df[column], ax=ax, inner="quartile", color="#5B8E7D")
            ax.set_title(f"{title_prefix}: {title}")
            ax.set_xlabel("")
            ax.set_ylabel(column)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        sns.scatterplot(
            data=qc_df,
            x="total_counts",
            y="n_genes_by_counts",
            hue="pct_counts_mt",
            palette="viridis",
            s=10,
            linewidth=0,
            ax=axes[0],
        )
        axes[0].set_title(f"{title_prefix}: counts vs genes")
        axes[0].legend(loc="best", title="Mito %", fontsize=8)

        sns.histplot(qc_df["pct_counts_mt"], bins=40, color="#C8553D", ax=axes[1])
        axes[1].set_title(f"{title_prefix}: mito percentage")
        axes[1].set_xlabel("pct_counts_mt")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def save_embedding_pdf(adata: sc.AnnData, basis: str, color: str, output_path: Path) -> None:
    fig = sc.pl.embedding(
        adata,
        basis=basis,
        color=[color],
        show=False,
        return_fig=True,
        frameon=False,
        legend_loc="on data" if color == "leiden" else "right margin",
    )
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def extract_marker_table(adata: sc.AnnData, groupby: str, top_n: int) -> pd.DataFrame:
    marker_frames: list[pd.DataFrame] = []
    for cluster in adata.obs[groupby].cat.categories:
        marker_df = sc.get.rank_genes_groups_df(adata, group=cluster).head(top_n).copy()
        marker_df.insert(0, "cluster", cluster)
        marker_frames.append(marker_df)
    return pd.concat(marker_frames, ignore_index=True)


def filter_marker_table(marker_table: pd.DataFrame, top_n: int) -> pd.DataFrame:
    filtered = marker_table[
        ~marker_table["names"].str.startswith(("RPS", "RPL", "MT-"), na=False)
    ].copy()
    filtered = filtered[filtered["logfoldchanges"] > 0].copy()
    return filtered.groupby("cluster", group_keys=False).head(top_n).reset_index(drop=True)


def annotate_clusters_by_markers(adata: sc.AnnData) -> tuple[pd.DataFrame, dict[str, str]]:
    annotations: list[dict[str, object]] = []
    annotation_map: dict[str, str] = {}

    for cluster in adata.obs["leiden"].cat.categories:
        subset = adata[adata.obs["leiden"] == cluster]
        scores: dict[str, float] = {}
        evidence: dict[str, str] = {}
        for cell_type, markers in CANONICAL_MARKERS.items():
            present = [gene for gene in markers if gene in adata.raw.var_names]
            expr = subset.raw[:, present].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray()
            expr = np.asarray(expr)
            mean_expr = expr.mean(axis=0)
            scores[cell_type] = float(mean_expr.mean())
            evidence[cell_type] = ",".join(
                gene for gene, value in zip(present, mean_expr, strict=True) if value > 1.0
            )

        best_cell_type = max(scores, key=scores.get)
        annotation_map[cluster] = best_cell_type
        annotations.append(
            {
                "cluster": cluster,
                "cell_type": best_cell_type,
                "cell_count": int(subset.n_obs),
                "supporting_markers": evidence[best_cell_type],
                **{f"{cell_type}_score": round(score, 4) for cell_type, score in scores.items()},
            }
        )

    return pd.DataFrame(annotations), annotation_map


def save_marker_dotplot(adata: sc.AnnData, output_path: Path) -> None:
    marker_genes = []
    for genes in CANONICAL_MARKERS.values():
        for gene in genes:
            if gene in adata.raw.var_names and gene not in marker_genes:
                marker_genes.append(gene)
    fig = sc.pl.dotplot(
        adata,
        var_names=marker_genes,
        groupby="cell_type",
        use_raw=True,
        show=False,
        return_fig=True,
        standard_scale="var",
    )
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close("all")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.verbosity = 2
    sc.set_figure_params(dpi=120, facecolor="white")

    adata = sc.read_10x_mtx(
        args.input_dir,
        var_names="gene_symbols",
        make_unique=True,
    )
    adata.var_names_make_unique()
    input_cells = adata.n_obs
    input_genes = adata.n_vars

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=[20],
        log1p=False,
        inplace=True,
    )

    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    qc_before = adata.obs[qc_cols].copy()
    qc_before.to_csv(args.output_dir / "qc_metrics_before_filtering.tsv", sep="\t")
    save_qc_pdf(qc_before, args.output_dir / "qc_before_filtering.pdf", "Before filtering")

    max_genes = int(np.ceil(adata.obs["n_genes_by_counts"].quantile(0.99)))
    max_counts = int(np.ceil(adata.obs["total_counts"].quantile(0.99)))

    keep_mask = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["total_counts"] <= max_counts)
        & (adata.obs["pct_counts_mt"] <= args.max_mt_pct)
    )
    adata = adata[keep_mask].copy()

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=[20],
        log1p=False,
        inplace=True,
    )
    qc_after = adata.obs[qc_cols].copy()
    qc_after.to_csv(args.output_dir / "qc_metrics_after_filtering.tsv", sep="\t")
    save_qc_pdf(qc_after, args.output_dir / "qc_after_filtering.pdf", "After filtering")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)

    n_comps = max(2, min(50, adata.n_obs - 1, adata.n_vars - 1))
    n_pcs = max(2, min(30, n_comps))
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata, random_state=0)
    sc.tl.tsne(adata, use_rep="X_pca", random_state=0)
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added="leiden",
        flavor="igraph",
        directed=False,
        n_iterations=2,
    )

    save_embedding_pdf(adata, "tsne", "leiden", args.output_dir / "tsne_leiden_clusters.pdf")
    save_embedding_pdf(adata, "umap", "leiden", args.output_dir / "umap_leiden_clusters.pdf")
    save_embedding_pdf(
        adata, "umap", "pct_counts_mt", args.output_dir / "umap_mito_percentage.pdf"
    )

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)
    marker_table = extract_marker_table(adata, groupby="leiden", top_n=20)
    marker_table.to_csv(args.output_dir / "cluster_markers_top20.tsv", sep="\t", index=False)
    filter_marker_table(marker_table, top_n=10).to_csv(
        args.output_dir / "cluster_markers_top10_filtered.tsv", sep="\t", index=False
    )

    annotation_df, annotation_map = annotate_clusters_by_markers(adata)
    annotation_df.to_csv(args.output_dir / "cluster_annotations.tsv", sep="\t", index=False)
    adata.obs["cell_type"] = adata.obs["leiden"].map(annotation_map).astype("category")

    cluster_counts = adata.obs["leiden"].value_counts().sort_index()
    cluster_df = cluster_counts.rename_axis("cluster").reset_index(name="cell_count")
    cluster_df.to_csv(args.output_dir / "cluster_cell_counts.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.barplot(data=cluster_df, x="cluster", y="cell_count", color="#4C78A8", ax=ax)
    ax.set_title("Cell counts per Leiden cluster")
    fig.tight_layout()
    fig.savefig(args.output_dir / "cluster_cell_counts.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

    save_embedding_pdf(
        adata, "umap", "cell_type", args.output_dir / "umap_cell_type_annotations.pdf"
    )
    save_embedding_pdf(
        adata, "tsne", "cell_type", args.output_dir / "tsne_cell_type_annotations.pdf"
    )
    save_marker_dotplot(adata, args.output_dir / "canonical_marker_dotplot.pdf")

    obs_export = adata.obs[
        ["total_counts", "n_genes_by_counts", "pct_counts_mt", "leiden", "cell_type"]
    ].copy()
    obs_export["barcode"] = obs_export.index
    obs_export = obs_export[
        ["barcode", "total_counts", "n_genes_by_counts", "pct_counts_mt", "leiden", "cell_type"]
    ]
    obs_export.to_csv(args.output_dir / "cell_metadata.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            ("input_cells", int(input_cells)),
            ("input_genes", int(input_genes)),
            ("filtered_cells", int(adata.n_obs)),
            ("highly_variable_genes", int(adata.n_vars)),
            ("min_genes_threshold", int(args.min_genes)),
            ("max_genes_threshold", int(max_genes)),
            ("max_counts_threshold", int(max_counts)),
            ("max_mt_pct_threshold", float(args.max_mt_pct)),
            ("leiden_resolution", float(args.resolution)),
            ("cluster_count", int(cluster_df.shape[0])),
            ("annotated_cell_types", int(annotation_df["cell_type"].nunique())),
        ],
        columns=["metric", "value"],
    )
    summary.to_csv(args.output_dir / "analysis_summary.tsv", sep="\t", index=False)

    adata.write_h5ad(args.output_dir / "pbmc_processed.h5ad")


if __name__ == "__main__":
    main()
