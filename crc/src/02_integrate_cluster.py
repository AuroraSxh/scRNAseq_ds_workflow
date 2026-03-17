from __future__ import annotations

from pathlib import Path

import anndata as ad
import bbknn
import pandas as pd
import scanpy as sc

from crc_sc_integration.utils import load_yaml, save_umap, setup_logging


ROOT = Path(__file__).resolve().parents[1]
PARAMS = load_yaml(ROOT / "config/params.yaml")


def main() -> None:
    logger = setup_logging(ROOT / "logs/integrate/integrate_cluster.log")
    processed_dir = ROOT / "data/processed"
    results_dir = ROOT / "results/integration"
    processed_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading QC-passed objects")
    adata_1 = sc.read_h5ad(ROOT / "data/interim/gse178341_qc.h5ad")
    adata_2 = sc.read_h5ad(ROOT / "data/interim/gse146771_qc.h5ad")
    adata_1.layers.clear()
    adata_2.layers.clear()
    for col in ["source_global_cluster", "source_subcluster"]:
        if col not in adata_1.obs:
            adata_1.obs[col] = pd.NA
        if col not in adata_2.obs:
            adata_2.obs[col] = pd.NA

    common_genes = adata_1.var_names.intersection(adata_2.var_names)
    logger.info("Common genes retained for integration: %s", len(common_genes))
    adata_1 = adata_1[:, common_genes].copy()
    adata_2 = adata_2[:, common_genes].copy()

    logger.info("Concatenating datasets")
    adata = ad.concat(
        [adata_1, adata_2],
        join="inner",
        label="dataset",
        keys=["GSE178341", "GSE146771"],
        index_unique="-",
        merge="same",
    )

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=PARAMS["general"]["n_top_genes"],
        batch_key="dataset",
        flavor="seurat",
    )
    adata = adata[:, adata.var["highly_variable"]].copy()
    logger.info("Highly variable genes selected: %s", adata.n_vars)

    sc.pp.pca(
        adata,
        n_comps=PARAMS["general"]["n_pcs"],
        zero_center=False,
        svd_solver="randomized",
        random_state=PARAMS["general"]["random_state"],
    )

    batch_key = PARAMS["general"]["batch_key"]
    logger.info("Running BBKNN with batch_key=%s", batch_key)
    bbknn.bbknn(
        adata,
        batch_key=batch_key,
        neighbors_within_batch=PARAMS["general"]["neighbors_within_batch"],
        n_pcs=PARAMS["general"]["n_pcs"],
    )

    sc.tl.umap(
        adata,
        random_state=PARAMS["general"]["random_state"],
        min_dist=PARAMS["general"]["umap_min_dist"],
        spread=PARAMS["general"]["umap_spread"],
        negative_sample_rate=PARAMS["general"]["umap_negative_sample_rate"],
        gamma=PARAMS["general"]["umap_repulsion_strength"],
        init_pos=PARAMS["general"]["umap_init_pos"],
    )
    sc.tl.leiden(
        adata,
        resolution=PARAMS["general"]["leiden_resolution"],
        key_added="cluster",
    )
    sc.tl.rank_genes_groups(adata, groupby="cluster", method="wilcoxon")

    logger.info("Writing integration outputs")
    adata.write_h5ad(processed_dir / "crc_integrated_clusters.h5ad", compression="gzip")
    marker_df = sc.get.rank_genes_groups_df(adata, group=None)
    marker_df.to_csv(results_dir / "cluster_markers.csv", index=False)
    adata.obs.groupby(["cluster", "dataset"], observed=True).size().rename("n_cells").reset_index().to_csv(
        results_dir / "cluster_dataset_composition.csv", index=False
    )

    save_umap(
        adata,
        ["dataset", "cluster", "tissue"],
        results_dir / "figures/integrated_umap.png",
        "Integrated CRC single-cell atlas",
    )
    logger.info("Integration step finished")


if __name__ == "__main__":
    main()
