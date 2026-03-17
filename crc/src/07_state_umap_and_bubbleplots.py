from __future__ import annotations

from pathlib import Path

import scanpy as sc

from crc_sc_integration.utils import load_full_common_atlas, load_yaml, save_embedding, save_marker_bubbleplot, setup_logging


ROOT = Path(__file__).resolve().parents[1]
MARKERS = load_yaml(ROOT / "config/marker_panels.yaml")

CD8_STATE_ORDER = [
    "Tpex",
    "Early_Tex",
    "Teff_like",
    "Tem",
    "Term_Tex",
]

CD8_BUBBLE_MARKERS = {state: MARKERS["cd8_states"][state] for state in CD8_STATE_ORDER}

MAJOR_BUBBLE_MARKERS = {
    "Epithelial": ["KRT8", "KRT19"],
    "Endothelial": ["PECAM1", "VWF", "KDR"],
    "CD4/Treg": ["IL7R", "LTB", "FOXP3"],
    "Myeloid": ["LST1", "FCER1G", "C1QC"],
    "Plasma": ["MZB1", "JCHAIN"],
    "B": ["MS4A1", "CD79A", "CD74"],
    "CD8 T": ["CD8A", "CCL5"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "Cycling": ["MKI67", "TOP2A", "STMN1"],
    "Fibroblast": ["COL1A1", "DCN"],
    "Mast": ["TPSAB1", "KIT", "CPA3"],
}

MAJOR_ORDER = [
    "Epithelial",
    "Endothelial",
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


def main() -> None:
    logger = setup_logging(ROOT / "logs/annotation/state_umap_and_bubbleplots.log")

    logger.info("Loading annotated total object and CD8 object")
    annotated = sc.read_h5ad(ROOT / "data/processed/crc_integrated_annotated.h5ad")
    cd8 = sc.read_h5ad(ROOT / "data/processed/cd8_subclusters.h5ad")
    full_atlas = load_full_common_atlas(ROOT)
    full_atlas.obs = full_atlas.obs.join(annotated.obs[["major_cell_type"]], how="left")
    cluster_order = (
        cd8.obs["cd8_cluster_label"].cat.categories.tolist()
        if hasattr(cd8.obs["cd8_cluster_label"], "cat")
        else sorted(cd8.obs["cd8_cluster_label"].astype(str).unique().tolist())
    )

    logger.info("Saving CD8 subgroup UMAP")
    save_embedding(
        cd8,
        ["cd8_cluster_label"],
        ROOT / "results/cd8/figures/cd8_state_guided_umap.png",
        "CD8+ T subgroups",
        basis="umap",
    )

    logger.info("Saving CD8 marker bubble plot")
    save_marker_bubbleplot(
        cd8,
        groupby="cd8_cluster_label",
        marker_dict=CD8_BUBBLE_MARKERS,
        output_path=ROOT / "results/cd8/figures/cd8_marker_bubbleplot.png",
        title="CD8 subgroup marker expression",
        groups_order=cluster_order,
        marker_group_order=CD8_STATE_ORDER,
    )

    logger.info("Saving CD8 state-only marker bubble plot")
    save_marker_bubbleplot(
        cd8,
        groupby="cd8_state_score_only",
        marker_dict=CD8_BUBBLE_MARKERS,
        output_path=ROOT / "results/cd8/figures/cd8_state_marker_bubbleplot.png",
        title="CD8 state marker expression",
        groups_order=CD8_STATE_ORDER,
        marker_group_order=CD8_STATE_ORDER,
    )

    logger.info("Saving major cell-type marker bubble plot")
    save_marker_bubbleplot(
        full_atlas,
        groupby="major_cell_type",
        marker_dict=MAJOR_BUBBLE_MARKERS,
        output_path=ROOT / "results/annotation/figures/major_marker_bubbleplot.png",
        title="Major cell-type marker expression",
        groups_order=MAJOR_ORDER,
        marker_group_order=MAJOR_ORDER,
    )
    logger.info("CD8 subgroup UMAP and bubble plots finished")


if __name__ == "__main__":
    main()
