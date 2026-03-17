from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from crc_sc_integration.utils import (
    load_full_common_atlas,
    load_yaml,
    save_heatmap,
    set_categorical_colors,
    save_umap,
    score_gene_panels,
    setup_logging,
)


ROOT = Path(__file__).resolve().parents[1]
PARAMS = load_yaml(ROOT / "config/params.yaml")
MARKERS = load_yaml(ROOT / "config/marker_panels.yaml")

CD8_STATE_ORDER = ["Tpex", "Early_Tex", "Teff_like", "Tem", "Term_Tex"]
CD8_SPECIAL_PANELS = {
    "IFN_response": ["ISG15", "IFI6", "LY6E", "MX1", "OAS1", "IFI44L", "IFIT1", "IFIT3"],
    "Cycling": ["MKI67", "TOP2A", "TYMS", "STMN1", "PCNA"],
}
CD8_STATE_PALETTE = {
    "Tpex": "#3A7D44",
    "Early_Tex": "#D95F02",
    "Teff_like": "#C03A2B",
    "Tem": "#2C7FB8",
    "Term_Tex": "#7B3294",
}
CD8_SUBGROUP_CLASS_ORDER = CD8_STATE_ORDER + ["IFN_response", "Cycling"]


def mean_expression_by_group(
    adata: sc.AnnData,
    groupby: str,
    marker_dict: dict[str, list[str]],
) -> pd.DataFrame:
    group_labels = adata.obs[groupby].astype(str)
    rows = {}
    for label, genes in marker_dict.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if not present:
            rows[label] = pd.Series(np.nan, index=sorted(group_labels.unique(), key=lambda value: int(value)))
            continue
        matrix = adata[:, present].X
        if sp.issparse(matrix):
            matrix = matrix.toarray()
        expr = pd.DataFrame(matrix, index=adata.obs_names, columns=present)
        rows[label] = expr.groupby(group_labels, observed=True).mean().mean(axis=1)
    frame = pd.DataFrame(rows)
    return frame.loc[sorted(frame.index, key=lambda value: int(value))]


def zscore_columns(frame: pd.DataFrame) -> pd.DataFrame:
    centered = frame - frame.mean(axis=0)
    scaled = frame.std(axis=0, ddof=0).replace(0, 1.0)
    return centered.divide(scaled, axis=1).fillna(0.0)


def pick_cd8_cells(adata: sc.AnnData) -> sc.AnnData:
    score_gene_panels(adata, MARKERS["lineage_scores"], prefix="lineage")
    major_cols = score_gene_panels(adata, MARKERS["major_cell_types"], prefix="major")
    adata.obs["major_score_winner"] = adata.obs[major_cols].idxmax(axis=1).str.removeprefix("major_")

    platelet_genes = [gene for gene in ["PF4", "PPBP", "SDPR", "CLU", "NRGN"] if gene in adata.var_names]
    if len(platelet_genes) >= 2:
        sc.tl.score_genes(adata, platelet_genes, score_name="lineage_platelet", use_raw=False)
    else:
        adata.obs["lineage_platelet"] = 0.0

    ref_labels = adata.obs["source_global_cluster"].astype(str).replace("nan", "")
    cd8_ref = ref_labels.str.contains("CD8", case=False)
    mask = (
        (adata.obs["lineage_t_cell"] > adata.obs[["lineage_nk_cell", "lineage_treg"]].max(axis=1) + 0.05)
        & (
            adata.obs["lineage_t_cell"]
            > adata.obs[
                [
                    "major_myeloid",
                    "major_b",
                    "major_plasma",
                    "major_mast",
                    "major_fibroblast",
                    "major_endothelial",
                    "major_epithelial",
                ]
            ].max(axis=1)
        )
        & (adata.obs["lineage_platelet"] < adata.obs["lineage_t_cell"])
        & ((adata.obs["lineage_cd8_core"] >= adata.obs["lineage_cd4_core"] - 0.02) | cd8_ref)
        & (~adata.obs["major_score_winner"].isin(["myeloid", "b", "plasma", "mast", "fibroblast", "endothelial", "epithelial"]))
    )
    return adata[mask].copy()


def annotate_cd8_clusters(cd8_full: sc.AnnData, marker_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    panel_means = mean_expression_by_group(cd8_full, "cd8_cluster", MARKERS["cd8_states"])
    panel_z = zscore_columns(panel_means)
    panel_state = panel_z.idxmax(axis=1)

    special_means = mean_expression_by_group(cd8_full, "cd8_cluster", CD8_SPECIAL_PANELS)
    special_z = zscore_columns(special_means)
    subgroup_class = panel_state.copy()
    for cluster in subgroup_class.index:
        if special_z.empty:
            continue
        best_special = special_z.loc[cluster].idxmax()
        if float(special_z.loc[cluster, best_special]) >= max(1.0, float(panel_z.loc[cluster].max()) + 0.2):
            subgroup_class.loc[cluster] = best_special

    top_markers = (
        marker_df.groupby("group", observed=True)["names"]
        .apply(lambda series: ", ".join(series.astype(str).head(5)))
        .rename("top_markers")
    )
    cluster_sizes = cd8_full.obs["cd8_cluster"].value_counts().rename("n_cells")
    cluster_order = sorted(panel_means.index, key=lambda value: int(value))

    annotation = pd.DataFrame(
        {
            "cd8_cluster": cluster_order,
            "cd8_state": [panel_state.loc[cluster] for cluster in cluster_order],
            "cd8_subgroup": [subgroup_class.loc[cluster] for cluster in cluster_order],
            "n_cells": [int(cluster_sizes.loc[cluster]) for cluster in cluster_order],
            "top_markers": [top_markers.get(cluster, "") for cluster in cluster_order],
        }
    )
    annotation["cd8_cluster_label"] = annotation.apply(
        lambda row: f"{row['cd8_subgroup']}_c{row['cd8_cluster']}",
        axis=1,
    )
    return annotation, panel_z.loc[cluster_order, CD8_STATE_ORDER]


def main() -> None:
    logger = setup_logging(ROOT / "logs/cd8/cd8_subtype.log")
    results_dir = ROOT / "results/cd8"
    annotated = sc.read_h5ad(ROOT / "data/processed/crc_integrated_annotated.h5ad")
    label_df = annotated.obs[["major_cell_type", "source_global_cluster", "source_subcluster"]].copy()
    adata = load_full_common_atlas(ROOT)
    adata.obs = adata.obs.join(label_df, how="left")

    logger.info("Selecting candidate CD8+ T cells")
    cd8_full = pick_cd8_cells(adata)
    if cd8_full.n_obs < PARAMS["cd8"]["min_cells"]:
        raise RuntimeError(f"Too few CD8-like cells retained: {cd8_full.n_obs}")
    logger.info("CD8 candidate cells retained: %s", cd8_full.n_obs)
    cd8_full.raw = cd8_full.copy()

    cd8 = cd8_full.copy()
    sc.pp.highly_variable_genes(
        cd8,
        n_top_genes=PARAMS["cd8"]["n_top_genes"],
        batch_key="dataset",
        flavor="seurat",
    )
    cd8 = cd8[:, cd8.var["highly_variable"]].copy()
    sc.pp.scale(cd8, max_value=10)
    sc.pp.pca(
        cd8,
        n_comps=PARAMS["cd8"]["n_pcs"],
        svd_solver="arpack",
        random_state=PARAMS["general"]["random_state"],
    )
    harmony = hm.run_harmony(
        cd8.obsm["X_pca"],
        cd8.obs,
        PARAMS["general"]["batch_key"],
        max_iter_harmony=PARAMS["cd8"]["harmony_max_iter"],
        verbose=False,
        random_state=PARAMS["general"]["random_state"],
    )
    cd8.obsm["X_pca_harmony"] = harmony.Z_corr
    sc.pp.neighbors(
        cd8,
        n_neighbors=PARAMS["cd8"]["n_neighbors"],
        use_rep="X_pca_harmony",
        metric="cosine",
        random_state=PARAMS["general"]["random_state"],
    )
    sc.tl.umap(
        cd8,
        random_state=PARAMS["general"]["random_state"],
        min_dist=PARAMS["cd8"]["umap_min_dist"],
        spread=PARAMS["cd8"]["umap_spread"],
        negative_sample_rate=PARAMS["cd8"]["umap_negative_sample_rate"],
        gamma=PARAMS["cd8"]["umap_repulsion_strength"],
        init_pos=PARAMS["cd8"]["umap_init_pos"],
    )
    sc.tl.leiden(
        cd8,
        resolution=PARAMS["cd8"]["leiden_resolution"],
        key_added="cd8_cluster",
        flavor="igraph",
        directed=False,
        n_iterations=-1,
    )
    sc.tl.rank_genes_groups(cd8, groupby="cd8_cluster", method="wilcoxon")
    cd8_full.obs["cd8_cluster"] = cd8.obs["cd8_cluster"].astype(str).values
    cd8_full.obsm["X_pca"] = cd8.obsm["X_pca"].copy()
    cd8_full.obsm["X_pca_harmony"] = cd8.obsm["X_pca_harmony"].copy()
    cd8_full.obsm["X_umap"] = cd8.obsm["X_umap"].copy()
    cd8_full.obsp["connectivities"] = cd8.obsp["connectivities"].copy()
    cd8_full.obsp["distances"] = cd8.obsp["distances"].copy()
    cd8_full.uns["neighbors"] = deepcopy(cd8.uns.get("neighbors", {}))
    cd8_full.uns["umap"] = dict(cd8.uns.get("umap", {}))

    state_score_cols = score_gene_panels(cd8_full, MARKERS["cd8_states"], prefix="state", use_raw=False)
    state_to_label = {
        "state_tpex": "Tpex",
        "state_early_tex": "Early_Tex",
        "state_teff_like": "Teff_like",
        "state_tem": "Tem",
        "state_term_tex": "Term_Tex",
    }
    marker_df = sc.get.rank_genes_groups_df(cd8, group=None)
    state_df, heatmap = annotate_cd8_clusters(cd8_full, marker_df)

    state_map = state_df.set_index("cd8_cluster")["cd8_state"]
    subgroup_map = state_df.set_index("cd8_cluster")["cd8_subgroup"]
    cluster_label_map = state_df.set_index("cd8_cluster")["cd8_cluster_label"]

    cd8_full.obs["cd8_state_score_only"] = cd8_full.obs[state_score_cols].idxmax(axis=1).map(state_to_label)
    cd8_full.obs["cd8_state"] = cd8_full.obs["cd8_cluster"].map(state_map)
    cd8_full.obs["cd8_subgroup"] = cd8_full.obs["cd8_cluster"].map(subgroup_map)
    cd8_full.obs["cd8_cluster_label"] = cd8_full.obs["cd8_cluster"].map(cluster_label_map)

    cluster_label_order = (
        state_df.assign(
            subgroup_rank=state_df["cd8_subgroup"].map({name: idx for idx, name in enumerate(CD8_SUBGROUP_CLASS_ORDER)}),
            cluster_numeric=state_df["cd8_cluster"].astype(int),
        )
        .sort_values(["subgroup_rank", "cluster_numeric"], kind="stable")["cd8_cluster_label"]
        .tolist()
    )
    cluster_order = [str(cluster) for cluster in state_df["cd8_cluster"].tolist()]
    set_categorical_colors(cd8_full, "cd8_state", order=CD8_STATE_ORDER, palette=CD8_STATE_PALETTE)
    set_categorical_colors(cd8_full, "cd8_subgroup", order=CD8_SUBGROUP_CLASS_ORDER)
    set_categorical_colors(cd8_full, "cd8_cluster_label", order=cluster_label_order)
    set_categorical_colors(cd8_full, "cd8_cluster", order=cluster_order, palette=cd8_full.uns["cd8_cluster_label_colors"])
    set_categorical_colors(cd8_full, "dataset")

    logger.info("Writing CD8 state outputs")
    state_df.to_csv(results_dir / "cd8_cluster_annotation.csv", index=False)
    (
        cd8_full.obs["cd8_state"]
        .value_counts()
        .rename_axis("cd8_state")
        .reset_index(name="n_cells")
        .to_csv(results_dir / "cd8_state_counts.csv", index=False)
    )
    (
        cd8_full.obs["cd8_cluster_label"]
        .value_counts()
        .rename_axis("cd8_cluster_label")
        .reset_index(name="n_cells")
        .to_csv(results_dir / "cd8_subgroup_counts.csv", index=False)
    )
    marker_df["cd8_cluster_label"] = marker_df["group"].astype(str).map(cluster_label_map)
    marker_df["cd8_state"] = marker_df["group"].astype(str).map(state_map)
    marker_df["cd8_subgroup"] = marker_df["group"].astype(str).map(subgroup_map)
    marker_df.to_csv(results_dir / "cd8_cluster_markers.csv", index=False)
    cd8_full.write_h5ad(ROOT / "data/processed/cd8_subclusters.h5ad", compression="gzip")
    heatmap.index = [cluster_label_map.loc[str(cluster)] for cluster in heatmap.index.astype(str)]
    save_heatmap(heatmap, results_dir / "figures/cd8_state_score_heatmap.png", "CD8 subgroup state enrichment")
    save_umap(
        cd8_full,
        ["cd8_cluster_label", "cd8_state", "dataset"],
        results_dir / "figures/cd8_umap.png",
        "CD8+ T Harmony UMAP",
    )
    logger.info("CD8 subtype step finished")


if __name__ == "__main__":
    main()
