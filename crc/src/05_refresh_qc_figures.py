from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc

from crc_sc_integration.utils import save_qc_violin, setup_logging


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data/raw"
GSE178341_DIR = RAW_DIR / "GSE178341"
GSE146771_DIR = RAW_DIR / "GSE146771"


def _is_invalid_metrics(frame: pd.DataFrame) -> bool:
    return frame.empty or frame.isna().all().all()


def main() -> None:
    logger = setup_logging(ROOT / "logs/prepare/refresh_qc_figures.log")
    qc_dir = ROOT / "results/qc"
    qc_metrics_dir = qc_dir / "tables"
    qc_metrics_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading or reconstructing QC metric tables")
    after_1_path = qc_metrics_dir / "gse178341_after_qc_metrics.csv.gz"
    after_2_path = qc_metrics_dir / "gse146771_after_qc_metrics.csv.gz"
    if after_1_path.exists():
        after_1 = pd.read_csv(after_1_path, index_col=0)
        if _is_invalid_metrics(after_1):
            after_1 = sc.read_h5ad(ROOT / "data/interim/gse178341_qc.h5ad").obs[
                ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
            ].copy()
            after_1.to_csv(after_1_path, index=True, compression="gzip")
    else:
        adata_1 = sc.read_h5ad(ROOT / "data/interim/gse178341_qc.h5ad")
        after_1 = adata_1.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
        after_1.to_csv(after_1_path, index=True, compression="gzip")

    if after_2_path.exists():
        after_2 = pd.read_csv(after_2_path, index_col=0)
        if _is_invalid_metrics(after_2):
            after_2 = sc.read_h5ad(ROOT / "data/interim/gse146771_qc.h5ad").obs[
                ["total_counts", "n_genes_by_counts", "pct_counts_ribo"]
            ].copy()
            after_2.to_csv(after_2_path, index=True, compression="gzip")
    else:
        adata_2 = sc.read_h5ad(ROOT / "data/interim/gse146771_qc.h5ad")
        after_2 = adata_2.obs[["total_counts", "n_genes_by_counts", "pct_counts_ribo"]].copy()
        after_2.to_csv(after_2_path, index=True, compression="gzip")

    before_1_path = qc_metrics_dir / "gse178341_before_qc_metrics.csv.gz"
    before_2_path = qc_metrics_dir / "gse146771_before_qc_metrics.csv.gz"
    if before_1_path.exists():
        before_1 = pd.read_csv(before_1_path, index_col=0)
        if _is_invalid_metrics(before_1):
            adata_1_raw = sc.read_10x_h5(GSE178341_DIR / "GSE178341_crc10x_full_c295v4_submit.h5")
            adata_1_raw.var["mt"] = adata_1_raw.var_names.str.upper().str.startswith("MT-")
            sc.pp.calculate_qc_metrics(adata_1_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
            before_1 = adata_1_raw.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
            before_1.to_csv(before_1_path, index=True, compression="gzip")
    else:
        adata_1_raw = sc.read_10x_h5(GSE178341_DIR / "GSE178341_crc10x_full_c295v4_submit.h5")
        adata_1_raw.var["mt"] = adata_1_raw.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata_1_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        before_1 = adata_1_raw.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
        before_1.to_csv(before_1_path, index=True, compression="gzip")

    if before_2_path.exists():
        before_2 = pd.read_csv(before_2_path, index_col=0)
        if _is_invalid_metrics(before_2):
            metadata = pd.read_csv(GSE146771_DIR / "GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz", sep="\t")
            ribo = pd.to_numeric(metadata["ribo.per"], errors="coerce")
            before_2 = pd.DataFrame(
                {
                    "total_counts": pd.to_numeric(metadata["filter.nUMI"], errors="coerce").to_numpy(),
                    "n_genes_by_counts": pd.to_numeric(metadata["filter.nGene"], errors="coerce").to_numpy(),
                    "pct_counts_ribo": ribo.where(ribo > 1.5, ribo * 100.0).to_numpy(),
                },
                index=metadata["CellName"].astype(str),
            )
            before_2.to_csv(before_2_path, index=True, compression="gzip")
    else:
        metadata = pd.read_csv(GSE146771_DIR / "GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz", sep="\t")
        ribo = pd.to_numeric(metadata["ribo.per"], errors="coerce")
        before_2 = pd.DataFrame(
            {
                "total_counts": pd.to_numeric(metadata["filter.nUMI"], errors="coerce").to_numpy(),
                "n_genes_by_counts": pd.to_numeric(metadata["filter.nGene"], errors="coerce").to_numpy(),
                "pct_counts_ribo": ribo.where(ribo > 1.5, ribo * 100.0).to_numpy(),
            },
            index=metadata["CellName"].astype(str),
        )
        before_2.to_csv(before_2_path, index=True, compression="gzip")

    logger.info("Rebuilding QC violin plots with improved styling")
    save_qc_violin(
        before_1,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        "GSE178341 before QC",
        qc_dir / "figures/gse178341_before_qc_violin.png",
    )
    save_qc_violin(
        after_1,
        ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        "GSE178341 after QC",
        qc_dir / "figures/gse178341_after_qc_violin.png",
    )
    save_qc_violin(
        before_2,
        ["total_counts", "n_genes_by_counts", "pct_counts_ribo"],
        "GSE146771 before QC",
        qc_dir / "figures/gse146771_before_qc_violin.png",
    )
    save_qc_violin(
        after_2,
        ["total_counts", "n_genes_by_counts", "pct_counts_ribo"],
        "GSE146771 after QC",
        qc_dir / "figures/gse146771_after_qc_violin.png",
    )
    logger.info("QC figure refresh finished")


if __name__ == "__main__":
    main()
