from __future__ import annotations

import gc
from pathlib import Path

import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.io import mmwrite

from crc_sc_integration.utils import load_full_common_atlas, setup_logging


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "results/communication/cellchat"
INPUT_DIR = OUTPUT_DIR / "input"
DB_GENE_PATH = OUTPUT_DIR / "cellchat_db_genes.tsv"
CD8_ORDER = ["Tpex_like", "Early_Tex", "Effector_memory", "Terminal_Tex", "Stress_program"]
MACROPHAGE_ORDER = ["C1QC_TAM", "FOLR2_TAM", "SPP1_TAM", "Inflammatory_Macro"]


def main() -> None:
    logger = setup_logging(ROOT / "logs/annotation/export_cellchat_inputs.log")
    INPUT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Loading full shared atlas and subgroup annotations")
    full_atlas = load_full_common_atlas(ROOT)
    cd8 = sc.read_h5ad(ROOT / "data/processed/cd8_merged.h5ad")
    macrophage = sc.read_h5ad(ROOT / "data/processed/macrophage_subclusters.h5ad")

    metadata = full_atlas.obs.loc[:, [col for col in ["dataset", "tissue"] if col in full_atlas.obs.columns]].copy()
    metadata["cell_id"] = metadata.index.astype(str)
    metadata["lineage"] = pd.NA
    metadata["cellchat_group"] = pd.NA

    metadata.loc[cd8.obs_names, "lineage"] = "CD8"
    metadata.loc[cd8.obs_names, "cellchat_group"] = cd8.obs["cd8_merged_subgroup"].astype(str).values
    metadata.loc[macrophage.obs_names, "lineage"] = "Myeloid"
    metadata.loc[macrophage.obs_names, "cellchat_group"] = macrophage.obs["macrophage_state"].astype(str).values

    selected_mask = metadata["cellchat_group"].notna()
    selected_meta = metadata.loc[selected_mask].copy()
    selected_meta["group_order"] = selected_meta["cellchat_group"].map(
        {group: idx for idx, group in enumerate(MACROPHAGE_ORDER + CD8_ORDER)}
    )
    selected_meta = selected_meta.sort_values(["group_order", "cell_id"], kind="stable").drop(columns="group_order")

    logger.info("Selected %s cells across %s groups", selected_meta.shape[0], selected_meta["cellchat_group"].nunique())
    subset = full_atlas[selected_meta.index, :].copy()
    if DB_GENE_PATH.exists():
        db_genes = pd.read_csv(DB_GENE_PATH, sep="\t", header=None)[0].astype(str)
        keep_genes = subset.var_names.isin(set(db_genes))
        logger.info("Restricting export to %s CellChat database genes", int(keep_genes.sum()))
        subset = subset[:, keep_genes].copy()
    matrix = subset.X
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()

    gene_mask = (matrix > 0).sum(axis=0)
    gene_mask = pd.Series(gene_mask.A1, index=subset.var_names) > 0
    subset = subset[:, gene_mask.to_numpy()].copy()
    matrix = subset.X
    if not sp.issparse(matrix):
        matrix = sp.csc_matrix(matrix)
    else:
        matrix = matrix.tocsc()

    logger.info("Writing sparse CellChat matrix with %s genes", subset.n_vars)
    mmwrite(INPUT_DIR / "expression.mtx", matrix.transpose())
    pd.Series(subset.var_names.astype(str), name="gene_symbol").to_csv(
        INPUT_DIR / "genes.tsv",
        sep="\t",
        index=False,
        header=False,
    )
    selected_meta.to_csv(INPUT_DIR / "cells.tsv", sep="\t", index=False)
    (
        selected_meta.groupby(["lineage", "cellchat_group"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .to_csv(OUTPUT_DIR / "cellchat_group_counts.csv", index=False)
    )

    logger.info("CellChat input export finished")
    del full_atlas, cd8, macrophage, subset, matrix
    gc.collect()


if __name__ == "__main__":
    main()
