from __future__ import annotations

import gc
import os
from copy import deepcopy
from pathlib import Path

import harmonypy as hm
import matplotlib
os.environ.setdefault("MPLBACKEND", "Agg")
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns

from crc_sc_integration.utils import (
    cluster_mean_scores,
    finalize_axes,
    load_yaml,
    save_embedding,
    save_marker_bubbleplot,
    score_gene_panels,
    set_beautiful_style,
    set_categorical_colors,
    setup_logging,
)


ROOT = Path(__file__).resolve().parents[1]
PARAMS = load_yaml(ROOT / "config/params.yaml")
MARKERS = load_yaml(ROOT / "config/marker_panels.yaml")

CD8_MERGED_ORDER = ["Tpex_like", "Early_Tex", "Effector_memory", "Terminal_Tex", "Stress_program"]
CD8_MERGED_PALETTE = {
    "Tpex_like": "#2B7A78",
    "Early_Tex": "#E07A5F",
    "Effector_memory": "#3D5A80",
    "Terminal_Tex": "#7D4E57",
    "Stress_program": "#8D99AE",
}
MACROPHAGE_ORDER = ["C1QC_TAM", "FOLR2_TAM", "SPP1_TAM", "Inflammatory_Macro"]
MACROPHAGE_PALETTE = {
    "C1QC_TAM": "#4C956C",
    "FOLR2_TAM": "#2E86AB",
    "SPP1_TAM": "#C06C84",
    "Inflammatory_Macro": "#E76F51",
}
MYELOID_PANEL_ORDER = MACROPHAGE_ORDER + ["Monocyte_CD14", "Monocyte_CD16", "Dendritic"]
MACROPHAGE_STATES = set(MACROPHAGE_ORDER)
CANONICAL_MYELOID_LABELS = {name.lower(): name for name in MYELOID_PANEL_ORDER}

LIGAND_RECEPTOR_PAIRS = [
    {"ligand": "CXCL9", "receptor": "CXCR3", "pathway": "CXCL9-CXCR3"},
    {"ligand": "CXCL10", "receptor": "CXCR3", "pathway": "CXCL10-CXCR3"},
    {"ligand": "CXCL11", "receptor": "CXCR3", "pathway": "CXCL11-CXCR3"},
    {"ligand": "CXCL16", "receptor": "CXCR6", "pathway": "CXCL16-CXCR6"},
    {"ligand": "CD274", "receptor": "PDCD1", "pathway": "PD-L1"},
    {"ligand": "PDCD1LG2", "receptor": "PDCD1", "pathway": "PD-L2"},
    {"ligand": "LGALS9", "receptor": "HAVCR2", "pathway": "LGALS9-HAVCR2"},
    {"ligand": "TNF", "receptor": "TNFRSF1B", "pathway": "TNF-TNFRSF1B"},
    {"ligand": "TNF", "receptor": "TNFRSF1A", "pathway": "TNF-TNFRSF1A"},
    {"ligand": "TGFB1", "receptor": "TGFBR1", "pathway": "TGFB1-TGFBR1"},
    {"ligand": "TGFB1", "receptor": "TGFBR2", "pathway": "TGFB1-TGFBR2"},
    {"ligand": "SPP1", "receptor": "CD44", "pathway": "SPP1-CD44"},
    {"ligand": "SPP1", "receptor": "ITGAV", "pathway": "SPP1-ITGAV"},
    {"ligand": "SPP1", "receptor": "ITGB1", "pathway": "SPP1-ITGB1"},
    {"ligand": "ICAM1", "receptor": "ITGAL", "pathway": "ICAM1-ITGAL"},
    {"ligand": "ICAM1", "receptor": "ITGB2", "pathway": "ICAM1-ITGB2"},
    {"ligand": "CD80", "receptor": "CD28", "pathway": "CD80-CD28"},
    {"ligand": "CD86", "receptor": "CD28", "pathway": "CD86-CD28"},
    {"ligand": "CD86", "receptor": "CTLA4", "pathway": "CD86-CTLA4"},
    {"ligand": "TNFSF9", "receptor": "TNFRSF9", "pathway": "TNFSF9-TNFRSF9"},
    {"ligand": "IFNG", "receptor": "IFNGR1", "pathway": "IFNG-IFNGR1"},
    {"ligand": "IFNG", "receptor": "IFNGR2", "pathway": "IFNG-IFNGR2"},
    {"ligand": "CCL5", "receptor": "CCR5", "pathway": "CCL5-CCR5"},
    {"ligand": "HLA-E", "receptor": "KLRC1", "pathway": "HLAE-KLRC1"},
    {"ligand": "HLA-E", "receptor": "KLRC2", "pathway": "HLAE-KLRC2"},
]


def mean_expression_by_group(
    adata: sc.AnnData,
    groupby: str,
    marker_dict: dict[str, list[str]],
) -> pd.DataFrame:
    group_labels = adata.obs[groupby].astype(str)
    rows = {}
    group_order = sorted(group_labels.unique(), key=lambda value: int(value))
    for label, genes in marker_dict.items():
        present = [gene for gene in genes if gene in adata.var_names]
        if not present:
            rows[label] = pd.Series(np.nan, index=group_order)
            continue
        matrix = adata[:, present].X
        if sp.issparse(matrix):
            matrix = matrix.toarray()
        expr = pd.DataFrame(matrix, index=adata.obs_names, columns=present)
        rows[label] = expr.groupby(group_labels, observed=True).mean().mean(axis=1)
    frame = pd.DataFrame(rows)
    return frame.loc[group_order]


def zscore_columns(frame: pd.DataFrame) -> pd.DataFrame:
    centered = frame - frame.mean(axis=0)
    scaled = frame.std(axis=0, ddof=0).replace(0, 1.0)
    return centered.divide(scaled, axis=1).fillna(0.0)


def coarse_reference_myeloid_label(label: str | float | None) -> str:
    if label is None or pd.isna(label):
        return "Unknown"
    label = str(label)
    if "TAM-C1QC" in label:
        return "C1QC_TAM"
    if "TAM-SPP1" in label:
        return "SPP1_TAM"
    if "Macro-PLTP" in label:
        return "FOLR2_TAM"
    if "Macro-NLRP3" in label or "Macro-IL1B" in label:
        return "Inflammatory_Macro"
    if "Mono-CD14" in label or "Monolike-FCN1" in label:
        return "Monocyte_CD14"
    if "Mono-CD16" in label:
        return "Monocyte_CD16"
    if "cDC" in label or "pDC" in label:
        return "Dendritic"
    if "Mast" in label:
        return "Mast"
    return "Unknown"


def canonical_myeloid_label(label: str | float | None) -> str:
    if label is None or pd.isna(label):
        return "Unknown"
    label_str = str(label)
    return CANONICAL_MYELOID_LABELS.get(label_str.lower(), label_str)


def top_reference_label(frame: pd.DataFrame, cluster_key: str, ref_key: str) -> pd.DataFrame:
    ref = (
        frame.groupby([cluster_key, ref_key], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = ref.groupby(cluster_key, observed=True)["n_cells"].sum().rename("cluster_total")
    ref = ref.merge(totals, on=cluster_key, how="left")
    ref["fraction"] = ref["n_cells"] / ref["cluster_total"]
    ref = ref.sort_values([cluster_key, "fraction"], ascending=[True, False], kind="stable")
    return ref.drop_duplicates(cluster_key).rename(
        columns={ref_key: "top_reference_label", "fraction": "top_reference_fraction"}
    )


def contamination_label(top_markers: list[str]) -> str | None:
    marker_set = set(top_markers)
    if {"PHGR1", "EPCAM", "KRT8", "KRT18", "KRT19"} & marker_set:
        return "Epithelial_contam"
    if {"CD3D", "CD3E", "TRAC", "IL32", "CCL5"} & marker_set:
        return "Tcell_contam"
    if {"TPSAB1", "TPSB2", "CPA3"} & marker_set:
        return "Mast"
    if {"MKI67", "TOP2A", "TYMS", "STMN1"} & marker_set:
        return "Cycling"
    if {"FCGR3B", "CXCL8", "G0S2"} & marker_set:
        return "Neutrophil_like"
    return None


def merged_cd8_label(subgroup: str) -> str:
    if subgroup == "Tpex":
        return "Tpex_like"
    if subgroup == "Early_Tex":
        return "Early_Tex"
    if subgroup in {"Tem", "Teff_like"}:
        return "Effector_memory"
    if subgroup == "Term_Tex":
        return "Terminal_Tex"
    return "Stress_program"


def annotate_cd8_merged(adata: sc.AnnData) -> sc.AnnData:
    cd8_data = adata.copy()
    cd8_data.obs["cd8_merged_subgroup"] = cd8_data.obs["cd8_subgroup"].astype(str).map(merged_cd8_label)
    set_categorical_colors(
        cd8_data,
        "cd8_merged_subgroup",
        order=CD8_MERGED_ORDER,
        palette=CD8_MERGED_PALETTE,
    )
    return cd8_data


def save_cd8_outputs(cd8_plot: sc.AnnData) -> pd.DataFrame:
    output_dir = ROOT / "results/cd8"
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    annotation = (
        cd8_plot.obs.groupby(["cd8_cluster_label", "cd8_subgroup", "cd8_merged_subgroup"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    annotation["merge_reason"] = annotation["cd8_merged_subgroup"].map(
        {
            "Tpex_like": "TCF7/IL7R-like progenitor-memory program",
            "Early_Tex": "PDCD1/TOX low-to-mid exhaustion transition",
            "Effector_memory": "GZMK/GZMA effector-memory continuum",
            "Terminal_Tex": "HAVCR2/CXCL13 terminal exhaustion program",
            "Stress_program": "Interferon or cycling overlay program",
        }
    )
    annotation.to_csv(output_dir / "cd8_merged_annotation.csv", index=False)
    (
        cd8_plot.obs["cd8_merged_subgroup"]
        .value_counts()
        .rename_axis("cd8_merged_subgroup")
        .reset_index(name="n_cells")
        .to_csv(output_dir / "cd8_merged_counts.csv", index=False)
    )
    save_marker_bubbleplot(
        cd8_plot,
        groupby="cd8_merged_subgroup",
        marker_dict=MARKERS["cd8_merged_states"],
        output_path=figure_dir / "cd8_merged_marker_bubbleplot.png",
        title="Merged CD8 subgroup markers",
        groups_order=CD8_MERGED_ORDER,
        marker_group_order=CD8_MERGED_ORDER,
        show_group_labels=True,
    )
    save_embedding(
        cd8_plot,
        ["cd8_merged_subgroup"],
        figure_dir / "cd8_merged_umap.png",
        "Merged CD8 subgroups",
        basis="umap",
    )
    cd8_plot.write_h5ad(ROOT / "data/processed/cd8_merged.h5ad", compression="gzip")
    return annotation


def run_myeloid_reclustering(annotated: sc.AnnData, logger) -> tuple[sc.AnnData, pd.DataFrame]:
    myeloid_full = annotated[annotated.obs["major_cell_type"].astype(str) == "Myeloid"].copy()

    logger.info("Reclustering myeloid compartment for macrophage refinement")
    myeloid = myeloid_full.copy()
    sc.pp.highly_variable_genes(
        myeloid,
        n_top_genes=PARAMS["myeloid"]["n_top_genes"],
        batch_key="dataset",
        flavor="seurat",
    )
    myeloid = myeloid[:, myeloid.var["highly_variable"]].copy()
    sc.pp.pca(
        myeloid,
        n_comps=PARAMS["myeloid"]["n_pcs"],
        zero_center=False,
        svd_solver="randomized",
        random_state=PARAMS["general"]["random_state"],
    )
    harmony = hm.run_harmony(
        myeloid.obsm["X_pca"],
        myeloid.obs,
        PARAMS["general"]["batch_key"],
        max_iter_harmony=PARAMS["myeloid"]["harmony_max_iter"],
        verbose=False,
        random_state=PARAMS["general"]["random_state"],
    )
    myeloid.obsm["X_pca_harmony"] = harmony.Z_corr
    sc.pp.neighbors(
        myeloid,
        n_neighbors=PARAMS["myeloid"]["n_neighbors"],
        use_rep="X_pca_harmony",
        metric="cosine",
        random_state=PARAMS["general"]["random_state"],
    )
    sc.tl.umap(
        myeloid,
        random_state=PARAMS["general"]["random_state"],
        min_dist=PARAMS["myeloid"]["umap_min_dist"],
        spread=PARAMS["myeloid"]["umap_spread"],
        negative_sample_rate=PARAMS["myeloid"]["umap_negative_sample_rate"],
        gamma=PARAMS["myeloid"]["umap_repulsion_strength"],
        init_pos=PARAMS["myeloid"]["umap_init_pos"],
    )
    sc.tl.leiden(
        myeloid,
        resolution=PARAMS["myeloid"]["leiden_resolution"],
        key_added="myeloid_cluster",
        flavor="igraph",
        directed=False,
        n_iterations=-1,
    )
    sc.tl.rank_genes_groups(myeloid, groupby="myeloid_cluster", method="wilcoxon")

    myeloid_full.obs["myeloid_cluster"] = myeloid.obs["myeloid_cluster"].astype(str).values
    myeloid_full.obsm["X_pca"] = myeloid.obsm["X_pca"].copy()
    myeloid_full.obsm["X_pca_harmony"] = myeloid.obsm["X_pca_harmony"].copy()
    myeloid_full.obsm["X_umap"] = myeloid.obsm["X_umap"].copy()
    myeloid_full.obsp["connectivities"] = myeloid.obsp["connectivities"].copy()
    myeloid_full.obsp["distances"] = myeloid.obsp["distances"].copy()
    myeloid_full.uns["neighbors"] = deepcopy(myeloid.uns.get("neighbors", {}))
    myeloid_full.uns["umap"] = dict(myeloid.uns.get("umap", {}))

    marker_df = sc.get.rank_genes_groups_df(myeloid, group=None)
    return myeloid_full, marker_df


def annotate_myeloid_clusters(myeloid_full: sc.AnnData, marker_df: pd.DataFrame) -> tuple[sc.AnnData, pd.DataFrame]:
    panel_cols = score_gene_panels(myeloid_full, MARKERS["macrophage_states"], prefix="macro", use_raw=False)
    cluster_scores = cluster_mean_scores(myeloid_full, "myeloid_cluster", panel_cols)
    cluster_scores.columns = [canonical_myeloid_label(col.removeprefix("macro_")) for col in cluster_scores.columns]
    panel_z = zscore_columns(cluster_scores)
    panel_label = panel_z.idxmax(axis=1).map(canonical_myeloid_label)

    myeloid_full.obs["reference_myeloid"] = myeloid_full.obs["source_subcluster"].map(coarse_reference_myeloid_label)
    ref_df = top_reference_label(myeloid_full.obs, "myeloid_cluster", "reference_myeloid")

    top_markers = (
        marker_df.groupby("group", observed=True)["names"]
        .apply(lambda series: [str(value) for value in series.head(12)])
        .rename("top_markers")
    )
    cluster_sizes = myeloid_full.obs["myeloid_cluster"].value_counts().rename("n_cells")

    rows = []
    for cluster in sorted(cluster_sizes.index, key=lambda value: int(value)):
        marker_list = top_markers.get(cluster, [])
        contam = contamination_label(marker_list)
        ref_row = ref_df[ref_df["myeloid_cluster"] == cluster]
        ref_label = "Unknown"
        ref_fraction = 0.0
        if not ref_row.empty:
            ref_label = canonical_myeloid_label(ref_row.iloc[0]["top_reference_label"])
            ref_fraction = float(ref_row.iloc[0]["top_reference_fraction"])

        final_label = canonical_myeloid_label(panel_label.loc[cluster])
        if ref_label in set(MYELOID_PANEL_ORDER) and ref_fraction >= 0.20:
            final_label = ref_label
        if contam is not None:
            final_label = contam
        rows.append(
            {
                "myeloid_cluster": cluster,
                "assigned_label": final_label,
                "top_panel_label": panel_label.loc[cluster],
                "top_reference_label": ref_label,
                "top_reference_fraction": ref_fraction,
                "n_cells": int(cluster_sizes.loc[cluster]),
                "top_markers": ", ".join(marker_list[:8]),
            }
        )

    annotation = pd.DataFrame(rows)
    label_map = annotation.set_index("myeloid_cluster")["assigned_label"]
    myeloid_full.obs["myeloid_label"] = myeloid_full.obs["myeloid_cluster"].map(label_map)

    macrophage = myeloid_full[myeloid_full.obs["myeloid_label"].isin(MACROPHAGE_STATES)].copy()
    macrophage.obs["macrophage_state"] = macrophage.obs["myeloid_label"].astype(str)
    macrophage.obs["macrophage_cluster_label"] = (
        macrophage.obs["macrophage_state"].astype(str) + "_m" + macrophage.obs["myeloid_cluster"].astype(str)
    )

    cluster_order = (
        annotation[annotation["assigned_label"].isin(MACROPHAGE_STATES)]
        .assign(
            state_rank=lambda df: df["assigned_label"].map({name: idx for idx, name in enumerate(MACROPHAGE_ORDER)}),
            cluster_numeric=lambda df: df["myeloid_cluster"].astype(int),
        )
        .sort_values(["state_rank", "cluster_numeric"], kind="stable")["myeloid_cluster"]
        .astype(str)
        .tolist()
    )
    cluster_label_order = [
        f"{row.assigned_label}_m{row.myeloid_cluster}"
        for row in annotation[annotation["assigned_label"].isin(MACROPHAGE_STATES)]
        .assign(
            state_rank=lambda df: df["assigned_label"].map({name: idx for idx, name in enumerate(MACROPHAGE_ORDER)}),
            cluster_numeric=lambda df: df["myeloid_cluster"].astype(int),
        )
        .sort_values(["state_rank", "cluster_numeric"], kind="stable")
        .itertuples()
    ]
    set_categorical_colors(macrophage, "macrophage_state", order=MACROPHAGE_ORDER, palette=MACROPHAGE_PALETTE)
    set_categorical_colors(macrophage, "macrophage_cluster_label", order=cluster_label_order)
    set_categorical_colors(macrophage, "myeloid_cluster", order=cluster_order, palette=macrophage.uns["macrophage_cluster_label_colors"])
    return macrophage, annotation


def save_macrophage_outputs(macrophage: sc.AnnData, annotation: pd.DataFrame) -> None:
    output_dir = ROOT / "results/macrophage"
    figure_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    annotation.to_csv(output_dir / "myeloid_cluster_annotation.csv", index=False)
    (
        macrophage.obs["macrophage_state"]
        .value_counts()
        .rename_axis("macrophage_state")
        .reset_index(name="n_cells")
        .to_csv(output_dir / "macrophage_state_counts.csv", index=False)
    )
    macrophage.write_h5ad(ROOT / "data/processed/macrophage_subclusters.h5ad", compression="gzip")

    save_embedding(
        macrophage,
        ["macrophage_cluster_label", "macrophage_state"],
        figure_dir / "macrophage_umap.png",
        "Macrophage refined subgroups",
        basis="umap",
    )
    save_marker_bubbleplot(
        macrophage,
        groupby="macrophage_state",
        marker_dict={state: MARKERS["macrophage_states"][state] for state in MACROPHAGE_ORDER},
        output_path=figure_dir / "macrophage_marker_bubbleplot.png",
        title="Macrophage subgroup markers",
        groups_order=MACROPHAGE_ORDER,
        marker_group_order=MACROPHAGE_ORDER,
        show_group_labels=True,
    )


def expression_stats(adata: sc.AnnData, groupby: str, genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        raise ValueError("No genes available for interaction analysis.")
    matrix = adata[:, present].X
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    expr = pd.DataFrame(matrix, index=adata.obs_names, columns=present)
    labels = adata.obs[groupby].astype(str)
    mean_expr = expr.groupby(labels, observed=True).mean()
    pct_expr = (expr > 0).groupby(labels, observed=True).mean() * 100.0
    return mean_expr, pct_expr


def interaction_score(mean_lig: float, mean_rec: float, pct_lig: float, pct_rec: float) -> float:
    return float(np.sqrt(max(mean_lig, 0.0) * max(mean_rec, 0.0)) * (pct_lig / 100.0) * (pct_rec / 100.0))


def compute_interactions(
    sender_mean: pd.DataFrame,
    sender_pct: pd.DataFrame,
    receiver_mean: pd.DataFrame,
    receiver_pct: pd.DataFrame,
    direction: str,
) -> pd.DataFrame:
    records = []
    for pair in LIGAND_RECEPTOR_PAIRS:
        ligand = pair["ligand"]
        receptor = pair["receptor"]
        if ligand not in sender_mean.columns or receptor not in receiver_mean.columns:
            continue
        for sender in sender_mean.index:
            for receiver in receiver_mean.index:
                mean_lig = float(sender_mean.loc[sender, ligand])
                mean_rec = float(receiver_mean.loc[receiver, receptor])
                pct_lig = float(sender_pct.loc[sender, ligand])
                pct_rec = float(receiver_pct.loc[receiver, receptor])
                score = interaction_score(mean_lig, mean_rec, pct_lig, pct_rec)
                if score <= 0.01:
                    continue
                records.append(
                    {
                        "direction": direction,
                        "sender": sender,
                        "receiver": receiver,
                        "ligand": ligand,
                        "receptor": receptor,
                        "pathway": pair["pathway"],
                        "mean_ligand_expr": mean_lig,
                        "mean_receptor_expr": mean_rec,
                        "pct_ligand_expr": pct_lig,
                        "pct_receptor_expr": pct_rec,
                        "interaction_score": score,
                    }
                )
    if not records:
        return pd.DataFrame(
            columns=[
                "direction",
                "sender",
                "receiver",
                "ligand",
                "receptor",
                "pathway",
                "mean_ligand_expr",
                "mean_receptor_expr",
                "pct_ligand_expr",
                "pct_receptor_expr",
                "interaction_score",
            ]
        )
    return pd.DataFrame(records).sort_values("interaction_score", ascending=False, kind="stable")


def save_interaction_heatmap(interactions: pd.DataFrame, output_path: Path, title: str) -> None:
    if interactions.empty:
        return
    matrix = (
        interactions.groupby(["sender", "receiver"], observed=True)["interaction_score"]
        .sum()
        .reset_index()
        .pivot(index="sender", columns="receiver", values="interaction_score")
        .fillna(0.0)
    )
    matrix = matrix.loc[[idx for idx in matrix.index if idx in set(matrix.index)], [col for col in matrix.columns if col in set(matrix.columns)]]
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, ax = plt.subplots(figsize=(1.2 * max(4, matrix.shape[1]) + 1.5, 0.75 * max(4, matrix.shape[0]) + 1.5))
    sns.heatmap(matrix, cmap="mako", linewidths=0.5, linecolor="#E5E7EB", ax=ax)
    finalize_axes(ax, title=title, xlabel="", ylabel="", tight=False)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_interaction_dotplot(interactions: pd.DataFrame, output_path: Path, title: str) -> None:
    if interactions.empty:
        return
    top = interactions.head(20).copy()
    top["pair"] = top["ligand"] + "-" + top["receptor"]
    top["receiver"] = pd.Categorical(top["receiver"], categories=list(dict.fromkeys(top["receiver"])), ordered=True)
    top["pair"] = pd.Categorical(top["pair"], categories=list(dict.fromkeys(top["pair"])), ordered=True)

    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig_h = max(4.8, 0.34 * len(top))
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    palette = sns.color_palette("Set2", n_colors=max(top["sender"].nunique(), 3))
    sender_palette = {sender: palette[idx] for idx, sender in enumerate(top["sender"].astype(str).unique())}
    sns.scatterplot(
        data=top,
        x="receiver",
        y="pair",
        size="interaction_score",
        sizes=(40, 360),
        hue="sender",
        palette=sender_palette,
        edgecolor="#374151",
        linewidth=0.4,
        alpha=0.92,
        ax=ax,
    )
    ax.tick_params(axis="x", rotation=35)
    finalize_axes(ax, title=title, xlabel="", ylabel="", tight=False)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1.0), loc="upper left")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_summary_report(cd8_annotation: pd.DataFrame, macrophage_annotation: pd.DataFrame) -> None:
    report_path = ROOT / "results/annotation/subgroup_refinement_summary.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)

    cd8_lines = [
        "- `Tpex_like`: `TCF7/IL7R/SLAMF6` 为主，代表前体样或记忆样 CD8。",
        "- `Early_Tex`: `PDCD1/TOX/TIGIT` 上升，但尚未进入强终末耗竭。",
        "- `Effector_memory`: `GZMK/GZMA/CCL5/CXCR3` 为主，覆盖 Tem 与部分活化效应群。",
        "- `Terminal_Tex`: `HAVCR2/ENTPD1/LAG3/CXCL13/GZMB` 为主，代表终末耗竭。",
        "- `Stress_program`: 主要保留 `IFN_response` 和 `Cycling`，作为叠加程序而不是独立发育终点。",
    ]
    macro_lines = [
        "- `C1QC_TAM`: `C1QA/B/C, APOE, LGMN`，偏吞噬与抗原呈递。",
        "- `FOLR2_TAM`: `FOLR2, LYVE1, PLTP, CD163, F13A1`，偏稳态/驻留样免疫调节。",
        "- `SPP1_TAM`: `SPP1, GPNMB, CTSB, APOC1, CD9`，偏基质重塑与 TAM 程序。",
        "- `Inflammatory_Macro`: `IL1B, NLRP3, NAMPT, PLAUR, VEGFA`，偏炎症与应激。",
    ]
    cd8_merge_csv = "results/cd8/cd8_merged_annotation.csv"
    macro_csv = "results/macrophage/myeloid_cluster_annotation.csv"

    content = "\n".join(
        [
            "# CD8 与 Macrophage 细分亚群标准",
            "",
            "## CD8 细分与合并原则",
            *cd8_lines,
            "",
            "## Macrophage 细分原则",
            *macro_lines,
            "",
            "## CD8 合并结果",
            f"- 详见 `{cd8_merge_csv}`，每个原始 cluster 都映射到一个合并后的 CD8 状态。",
            "",
            "## Macrophage 保留亚群",
            f"- 详见 `{macro_csv}`，保留的 macrophage 群为 `C1QC_TAM / FOLR2_TAM / SPP1_TAM / Inflammatory_Macro`。",
            "",
            "## 通讯分析说明",
            "- 当前环境未安装 R `CellChat`。本次输出采用 CellChat-style 配体-受体打分：`sqrt(mean_ligand * mean_receptor) * pct_ligand * pct_receptor`。",
            "- 该结果适合做 CRC 内部不同亚群的相对比较；若你后续需要严格复现 R CellChat 的数据库和统计检验，可在允许安装依赖后补跑。",
        ]
    )
    report_path.write_text(content + "\n", encoding="utf-8")


def main() -> None:
    logger = setup_logging(ROOT / "logs/annotation/refine_subgroups_and_interactions.log")

    logger.info("Step 1/3: annotate merged CD8 programs")
    cd8_processed = sc.read_h5ad(ROOT / "data/processed/cd8_subclusters.h5ad")
    cd8_plot = annotate_cd8_merged(cd8_processed)
    cd8_annotation = save_cd8_outputs(cd8_plot)
    del cd8_processed, cd8_plot
    gc.collect()

    logger.info("Step 2/3: recluster myeloid compartment and retain macrophage states")
    annotated = sc.read_h5ad(ROOT / "data/processed/crc_integrated_annotated.h5ad")
    myeloid_full, marker_df = run_myeloid_reclustering(annotated, logger)
    macrophage, macrophage_annotation = annotate_myeloid_clusters(myeloid_full, marker_df)
    save_macrophage_outputs(macrophage, macrophage_annotation)
    del annotated, myeloid_full, marker_df, macrophage
    gc.collect()

    logger.info("Step 3/3: compute macrophage-CD8 interaction scores")
    cd8_for_comm = sc.read_h5ad(ROOT / "data/processed/cd8_merged.h5ad")
    macrophage_for_comm = sc.read_h5ad(ROOT / "data/processed/macrophage_subclusters.h5ad")
    interaction_genes = sorted(
        {pair["ligand"] for pair in LIGAND_RECEPTOR_PAIRS}
        | {pair["receptor"] for pair in LIGAND_RECEPTOR_PAIRS}
        | {gene for genes in MARKERS["cd8_merged_states"].values() for gene in genes}
        | {gene for genes in MARKERS["macrophage_states"].values() for gene in genes}
    )

    macro_mean, macro_pct = expression_stats(macrophage_for_comm, "macrophage_state", interaction_genes)
    cd8_mean, cd8_pct = expression_stats(cd8_for_comm, "cd8_merged_subgroup", interaction_genes)

    macro_to_cd8 = compute_interactions(macro_mean, macro_pct, cd8_mean, cd8_pct, "macrophage_to_cd8")
    cd8_to_macro = compute_interactions(cd8_mean, cd8_pct, macro_mean, macro_pct, "cd8_to_macrophage")
    interaction_dir = ROOT / "results/communication"
    interaction_dir.mkdir(parents=True, exist_ok=True)
    macro_to_cd8.to_csv(interaction_dir / "macrophage_to_cd8_interactions.csv", index=False)
    cd8_to_macro.to_csv(interaction_dir / "cd8_to_macrophage_interactions.csv", index=False)

    save_interaction_heatmap(
        macro_to_cd8,
        interaction_dir / "figures/macrophage_to_cd8_heatmap.png",
        "Macrophage to CD8 interaction strength",
    )
    save_interaction_heatmap(
        cd8_to_macro,
        interaction_dir / "figures/cd8_to_macrophage_heatmap.png",
        "CD8 to macrophage interaction strength",
    )
    save_interaction_dotplot(
        macro_to_cd8,
        interaction_dir / "figures/macrophage_to_cd8_top_pairs.png",
        "Top macrophage to CD8 ligand-receptor pairs",
    )
    save_interaction_dotplot(
        cd8_to_macro,
        interaction_dir / "figures/cd8_to_macrophage_top_pairs.png",
        "Top CD8 to macrophage ligand-receptor pairs",
    )

    write_summary_report(cd8_annotation, macrophage_annotation)
    logger.info("Subgroup refinement and interaction analysis finished")


if __name__ == "__main__":
    main()
