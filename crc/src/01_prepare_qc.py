from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from crc_sc_integration.utils import (
    align_obs_metadata,
    coerce_percent,
    load_yaml,
    read_space_quoted_matrix,
    save_qc_violin,
    setup_logging,
    standardize_var_names,
)


ROOT = Path(__file__).resolve().parents[1]
PARAMS = load_yaml(ROOT / "config/params.yaml")
RAW_DIR = ROOT / "data/raw"
GSE178341_DIR = RAW_DIR / "GSE178341"
GSE146771_DIR = RAW_DIR / "GSE146771"


def prepare_gse178341(logger) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, dict]:
    params = PARAMS["gse178341"]
    logger.info("Reading GSE178341 10x HDF5")
    adata = sc.read_10x_h5(GSE178341_DIR / "GSE178341_crc10x_full_c295v4_submit.h5")
    if "gene_symbols" in adata.var.columns:
        adata.var_names = adata.var["gene_symbols"].astype(str)
    standardize_var_names(adata)
    metadata = pd.read_csv(GSE178341_DIR / "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz")
    metadata["cellID"] = metadata["cellID"].astype(str)
    align_obs_metadata(adata, metadata, "cellID")
    adata.obs["dataset"] = "GSE178341"
    adata.obs["sample_id"] = adata.obs["PatientTypeID"].astype(str)
    adata.obs["integration_batch"] = adata.obs["dataset"] + ":" + adata.obs["sample_id"]
    adata.obs["tissue"] = adata.obs["SPECIMEN_TYPE"].astype(str)
    adata.obs["platform"] = adata.obs["SINGLECELL_TYPE"].astype(str)

    sc.pp.filter_genes(adata, min_cells=params["min_cells_per_gene"])
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    before_qc = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()

    before = adata.n_obs
    keep = (
        adata.obs["n_genes_by_counts"].between(params["min_genes"], params["max_genes"])
        & adata.obs["total_counts"].between(params["min_counts"], params["max_counts"])
        & (adata.obs["pct_counts_mt"] <= params["max_mt_pct"])
    )
    adata = adata[keep].copy()
    after = adata.n_obs
    logger.info("GSE178341 retained %s / %s cells after QC", after, before)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=PARAMS["general"]["target_sum"])
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    after_qc = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()

    summary = {"dataset": "GSE178341", "cells_before_qc": before, "cells_after_qc": after, "genes": adata.n_vars}
    return adata, before_qc, after_qc, summary


def prepare_gse146771(logger) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, dict]:
    params = PARAMS["gse146771"]
    logger.info("Reading GSE146771 quoted-space expression matrix")
    metadata = pd.read_csv(GSE146771_DIR / "GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz", sep="\t")
    metadata["CellName"] = metadata["CellName"].astype(str)
    adata = read_space_quoted_matrix(
        GSE146771_DIR / "GSE146771_CRC.Leukocyte.10x.TPM.txt.gz",
        chunk_size=params["loader_chunk_size"],
        logger=logger,
    )
    standardize_var_names(adata)
    align_obs_metadata(adata, metadata, "CellName")

    adata.obs["dataset"] = "GSE146771"
    adata.obs["sample_id"] = adata.obs["Sample"].astype(str)
    adata.obs["integration_batch"] = adata.obs["dataset"] + ":" + adata.obs["sample_id"]
    adata.obs["tissue"] = adata.obs["Tissue"].astype(str)
    adata.obs["platform"] = adata.obs["Platform"].astype(str)
    adata.obs["source_global_cluster"] = adata.obs["Global_Cluster"].astype(str)
    adata.obs["source_subcluster"] = adata.obs["Sub_Cluster"].astype(str)
    adata.obs["n_genes_by_counts"] = pd.to_numeric(adata.obs["filter.nGene"], errors="coerce")
    adata.obs["total_counts"] = pd.to_numeric(adata.obs["filter.nUMI"], errors="coerce")
    adata.obs["pct_counts_ribo"] = coerce_percent(adata.obs["ribo.per"])
    before_qc = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_ribo"]].copy()

    before = adata.n_obs
    keep = (
        adata.obs["n_genes_by_counts"].between(params["min_genes"], params["max_genes"])
        & adata.obs["total_counts"].between(params["min_counts_proxy"], params["max_counts_proxy"])
        & (adata.obs["pct_counts_ribo"] <= params["max_ribo_pct"])
    )
    adata = adata[keep].copy()
    after = adata.n_obs
    logger.info("GSE146771 retained %s / %s cells after QC", after, before)

    if params["matrix_logged"]:
        adata.layers["lognorm"] = adata.X.copy()
    else:
        sc.pp.normalize_total(adata, target_sum=PARAMS["general"]["target_sum"])
        sc.pp.log1p(adata)
        adata.layers["lognorm"] = adata.X.copy()
    after_qc = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_ribo"]].copy()

    summary = {"dataset": "GSE146771", "cells_before_qc": before, "cells_after_qc": after, "genes": adata.n_vars}
    return adata, before_qc, after_qc, summary


def main() -> None:
    logger = setup_logging(ROOT / "logs/prepare/prepare_qc.log")
    qc_dir = ROOT / "results/qc"
    interim_dir = ROOT / "data/interim"
    qc_metrics_dir = qc_dir / "tables"
    interim_dir.mkdir(parents=True, exist_ok=True)
    qc_metrics_dir.mkdir(parents=True, exist_ok=True)

    adata_1, before_1, after_1, summary_1 = prepare_gse178341(logger)
    adata_2, before_2, after_2, summary_2 = prepare_gse146771(logger)

    logger.info("Writing interim h5ad files")
    adata_1.write_h5ad(interim_dir / "gse178341_qc.h5ad", compression="gzip")
    adata_2.write_h5ad(interim_dir / "gse146771_qc.h5ad", compression="gzip")
    before_1.to_csv(qc_metrics_dir / "gse178341_before_qc_metrics.csv.gz", index=True, compression="gzip")
    after_1.to_csv(qc_metrics_dir / "gse178341_after_qc_metrics.csv.gz", index=True, compression="gzip")
    before_2.to_csv(qc_metrics_dir / "gse146771_before_qc_metrics.csv.gz", index=True, compression="gzip")
    after_2.to_csv(qc_metrics_dir / "gse146771_after_qc_metrics.csv.gz", index=True, compression="gzip")

    pd.DataFrame([summary_1, summary_2]).to_csv(qc_dir / "qc_summary.csv", index=False)
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
    logger.info("Prepare/QC step finished")


if __name__ == "__main__":
    main()
