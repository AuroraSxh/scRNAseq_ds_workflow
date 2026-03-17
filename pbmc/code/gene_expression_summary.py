#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

GENE_ALIASES = {
    "PD-1": "PDCD1",
    "PD1": "PDCD1",
    "PROGRAMMED CELL DEATH 1": "PDCD1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize single-gene expression across annotated PBMC cell types."
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=Path("analysis_outputs/pbmc_processed.h5ad"),
        help="Processed AnnData file with cell_type annotations.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("analysis_outputs/Gene_expression"),
        help="Root directory for per-gene outputs.",
    )
    parser.add_argument(
        "genes",
        nargs="+",
        help="Gene symbols or aliases such as METTL16 or PD-1.",
    )
    return parser.parse_args()


def resolve_gene_name(requested_gene: str, adata: sc.AnnData) -> tuple[str, str]:
    canonical = GENE_ALIASES.get(requested_gene.upper(), requested_gene)
    if adata.raw is None or canonical not in adata.raw.var_names:
        raise ValueError(f"{requested_gene} resolved to {canonical}, but the gene is absent.")
    return requested_gene, canonical


def sanitize_dirname(name: str) -> str:
    return "".join(char if char.isalnum() or char in {"-", "_"} else "_" for char in name)


def extract_expression(adata: sc.AnnData, gene: str) -> np.ndarray:
    expr = adata.raw[:, gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    return np.asarray(expr).ravel()


def build_summary_table(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    summary = (
        df.groupby(["cell_type", "leiden"], observed=True)
        .agg(
            cell_count=(gene, "size"),
            mean_expr=(gene, "mean"),
            median_expr=(gene, "median"),
            positive_cells=(gene, lambda s: int((s > 0).sum())),
            positive_fraction=(gene, lambda s: float((s > 0).mean())),
        )
        .reset_index()
    )
    return summary.sort_values(["mean_expr", "positive_fraction"], ascending=False).reset_index(
        drop=True
    )


def run_kruskal(df: pd.DataFrame, gene: str) -> tuple[float, float]:
    groups = [group[gene].to_numpy() for _, group in df.groupby("cell_type", observed=True)]
    return stats.kruskal(*groups)


def save_plot(df: pd.DataFrame, summary: pd.DataFrame, gene: str, output_pdf: Path) -> None:
    order = summary["cell_type"].drop_duplicates().tolist()
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    sns.violinplot(
        data=df,
        x="cell_type",
        y=gene,
        order=order,
        inner="quartile",
        cut=0,
        color="#5B8E7D",
        ax=axes[0],
    )
    axes[0].set_title(f"{gene} expression by cell type")
    axes[0].tick_params(axis="x", rotation=25)

    sns.barplot(
        data=summary,
        x="cell_type",
        y="positive_fraction",
        order=order,
        color="#4C78A8",
        ax=axes[1],
    )
    axes[1].set_title(f"{gene} positive fraction")
    axes[1].set_ylabel("Fraction > 0")
    axes[1].tick_params(axis="x", rotation=25)

    fig.tight_layout()
    fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


def save_embedding_plot(
    adata: sc.AnnData,
    basis: str,
    gene: str,
    title: str,
    output_pdf: Path,
) -> None:
    fig = sc.pl.embedding(
        adata,
        basis=basis,
        color=[gene],
        use_raw=True,
        show=False,
        return_fig=True,
        frameon=False,
        color_map="viridis",
        title=title,
    )
    fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


def process_gene(
    adata: sc.AnnData, requested_gene: str, output_root: Path
) -> tuple[str, str, pd.DataFrame, float, float]:
    requested_gene, canonical_gene = resolve_gene_name(requested_gene, adata)
    gene_dir = output_root / sanitize_dirname(requested_gene)
    gene_dir.mkdir(parents=True, exist_ok=True)

    df = adata.obs[["cell_type", "leiden"]].copy()
    df[canonical_gene] = extract_expression(adata, canonical_gene)

    summary = build_summary_table(df, canonical_gene)
    summary.to_csv(gene_dir / "expression_by_cell_type.tsv", sep="\t", index=False)

    stat, pvalue = run_kruskal(df, canonical_gene)
    with (gene_dir / "kruskal_wallis.txt").open("w") as handle:
        handle.write(f"requested_gene\t{requested_gene}\n")
        handle.write(f"canonical_gene\t{canonical_gene}\n")
        handle.write(f"kruskal_wallis_stat\t{stat}\n")
        handle.write(f"kruskal_wallis_pvalue\t{pvalue}\n")

    save_plot(df, summary, canonical_gene, gene_dir / "expression_by_cell_type.pdf")
    save_embedding_plot(
        adata,
        "umap",
        canonical_gene,
        f"{requested_gene} ({canonical_gene}) on UMAP",
        gene_dir / "umap_expression.pdf",
    )
    save_embedding_plot(
        adata,
        "tsne",
        canonical_gene,
        f"{requested_gene} ({canonical_gene}) on tSNE",
        gene_dir / "tsne_expression.pdf",
    )
    return requested_gene, canonical_gene, summary, stat, pvalue


def main() -> None:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.h5ad)
    if "cell_type" not in adata.obs:
        raise ValueError("cell_type annotation is missing from the AnnData object.")

    for requested_gene in args.genes:
        original, canonical, summary, stat, pvalue = process_gene(
            adata, requested_gene, args.output_root
        )
        print(f"requested_gene={original}")
        print(f"canonical_gene={canonical}")
        print(summary.to_string(index=False))
        print(f"kruskal_wallis_stat={stat:.4f}")
        print(f"kruskal_wallis_pvalue={pvalue:.4e}")
        print()


if __name__ == "__main__":
    main()
