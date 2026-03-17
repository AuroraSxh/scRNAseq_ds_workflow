from __future__ import annotations

from pathlib import Path

import scanpy as sc

from crc_sc_integration.utils import (
    load_yaml,
    rerun_umap,
    save_integrated_overview,
    save_umap,
    set_categorical_colors,
    setup_logging,
)


ROOT = Path(__file__).resolve().parents[1]
PARAMS = load_yaml(ROOT / "config/params.yaml")
MAJOR_ORDER = [
    "Epithelial",
    "Fibroblast",
    "Mast",
    "Myeloid",
    "B",
    "Plasma",
    "NK",
    "CD4/Treg",
    "CD8 T",
    "Cycling",
]
MAJOR_PALETTE = {
    "Epithelial": "#5AA469",
    "Fibroblast": "#9CCB86",
    "Mast": "#D62828",
    "Myeloid": "#F4A9A8",
    "B": "#1F77B4",
    "Plasma": "#B39BC8",
    "NK": "#8C5FBF",
    "CD4/Treg": "#AEC7E8",
    "CD8 T": "#F28E2B",
    "Cycling": "#F6BD60",
}
DATASET_ORDER = ["GSE178341", "GSE146771"]
DATASET_PALETTE = {
    "GSE178341": "#1D5F9A",
    "GSE146771": "#E76F51",
}
TISSUE_ORDER = ["T", "N", "P"]
TISSUE_PALETTE = {
    "T": "#F28E2B",
    "N": "#4E79A7",
    "P": "#BFC7D5",
}


def apply_atlas_colors(adata: sc.AnnData) -> None:
    if "major_cell_type" in adata.obs:
        set_categorical_colors(adata, "major_cell_type", order=MAJOR_ORDER, palette=MAJOR_PALETTE)
    if "dataset" in adata.obs:
        set_categorical_colors(adata, "dataset", order=DATASET_ORDER, palette=DATASET_PALETTE)
    if "tissue" in adata.obs:
        set_categorical_colors(adata, "tissue", order=TISSUE_ORDER, palette=TISSUE_PALETTE)


def main() -> None:
    logger = setup_logging(ROOT / "logs/integrate/refresh_umaps.log")
    random_state = PARAMS["general"]["random_state"]

    logger.info("Refreshing integrated UMAP embedding and figures")
    integrated = sc.read_h5ad(ROOT / "data/processed/crc_integrated_clusters.h5ad")
    rerun_umap(
        integrated,
        random_state=random_state,
        min_dist=PARAMS["general"]["umap_min_dist"],
        spread=PARAMS["general"]["umap_spread"],
        negative_sample_rate=PARAMS["general"]["umap_negative_sample_rate"],
        gamma=PARAMS["general"]["umap_repulsion_strength"],
        init_pos=PARAMS["general"]["umap_init_pos"],
    )
    integrated.write_h5ad(ROOT / "data/processed/crc_integrated_clusters.h5ad", compression="gzip")
    save_umap(
        integrated,
        ["dataset", "cluster", "tissue"],
        ROOT / "results/integration/figures/integrated_umap.png",
        "Integrated CRC single-cell atlas",
    )

    logger.info("Refreshing annotated UMAP figure")
    annotated = sc.read_h5ad(ROOT / "data/processed/crc_integrated_annotated.h5ad")
    annotated.obsm["X_umap"] = integrated.obsm["X_umap"].copy()
    annotated.uns["umap"] = dict(integrated.uns.get("umap", {}))
    apply_atlas_colors(annotated)
    annotated.write_h5ad(ROOT / "data/processed/crc_integrated_annotated.h5ad", compression="gzip")
    save_umap(
        annotated,
        ["major_cell_type", "dataset"],
        ROOT / "results/annotation/figures/major_annotation_umap.png",
        "Major cell-type annotation",
    )
    save_integrated_overview(
        annotated,
        major_key="major_cell_type",
        side_keys=["dataset", "tissue"],
        output_path=ROOT / "results/integration/figures/integrated_umap.png",
        title="Integrated CRC single-cell atlas",
    )

    logger.info("Refreshing CD8 UMAP embedding and figure")
    cd8 = sc.read_h5ad(ROOT / "data/processed/cd8_subclusters.h5ad")
    rerun_umap(
        cd8,
        random_state=random_state,
        min_dist=PARAMS["cd8"]["umap_min_dist"],
        spread=PARAMS["cd8"]["umap_spread"],
        negative_sample_rate=PARAMS["cd8"]["umap_negative_sample_rate"],
        gamma=PARAMS["cd8"]["umap_repulsion_strength"],
        init_pos=PARAMS["cd8"]["umap_init_pos"],
    )
    cd8.write_h5ad(ROOT / "data/processed/cd8_subclusters.h5ad", compression="gzip")
    set_categorical_colors(cd8, "dataset", order=DATASET_ORDER, palette=DATASET_PALETTE)
    save_umap(
        cd8,
        ["cd8_cluster_label", "cd8_state", "dataset"],
        ROOT / "results/cd8/figures/cd8_umap.png",
        "CD8+ T Harmony UMAP",
    )
    logger.info("UMAP refresh finished")


if __name__ == "__main__":
    main()
