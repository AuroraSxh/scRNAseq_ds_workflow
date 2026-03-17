from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

from crc_sc_integration.utils import finalize_axes, set_beautiful_style, setup_logging


ROOT = Path(__file__).resolve().parents[1]
CELLCHAT_DIR = ROOT / "results/communication/cellchat"
OUT_DIR = ROOT / "results/communication/mac_t"
MACROPHAGE_ORDER = ["C1QC_TAM", "FOLR2_TAM", "SPP1_TAM", "Inflammatory_Macro"]
CD8_ORDER = ["Tpex_like", "Early_Tex", "Effector_memory", "Terminal_Tex", "Stress_program"]
LINEAGE_PAIR_LABELS = {
    ("Myeloid", "CD8"): ("Macrophage to CD8", MACROPHAGE_ORDER, CD8_ORDER),
    ("CD8", "Myeloid"): ("CD8 to Macrophage", CD8_ORDER, MACROPHAGE_ORDER),
}
DISPLAY_LABELS = {
    "C1QC_TAM": "C1QC TAM",
    "FOLR2_TAM": "FOLR2 TAM",
    "SPP1_TAM": "SPP1 TAM",
    "Inflammatory_Macro": "Inflammatory\nMacro",
    "Tpex_like": "Tpex-like",
    "Early_Tex": "Early Tex",
    "Effector_memory": "Effector\nmemory",
    "Terminal_Tex": "Terminal Tex",
    "Stress_program": "Stress\nprogram",
}


def pretty_label(value: str) -> str:
    return DISPLAY_LABELS.get(value, value.replace("_", " "))


def load_interactions() -> pd.DataFrame:
    frame = pd.read_csv(CELLCHAT_DIR / "cellchat_myeloid_cd8_interactions.csv")
    frame["source"] = frame["source"].astype(str)
    frame["target"] = frame["target"].astype(str)
    frame["pathway_name"] = frame["pathway_name"].astype(str)
    frame["interaction_name_2"] = frame["interaction_name_2"].astype(str)
    return frame


def save_direction_heatmap(
    interactions: pd.DataFrame,
    *,
    source_lineage: str,
    target_lineage: str,
    output_path: Path,
) -> pd.DataFrame:
    title, row_order, col_order = LINEAGE_PAIR_LABELS[(source_lineage, target_lineage)]
    summary = (
        interactions[
            (interactions["source_lineage"] == source_lineage)
            & (interactions["target_lineage"] == target_lineage)
        ]
        .groupby(["source", "target"], observed=True)
        .agg(
            summed_prob=("prob", "sum"),
            n_interactions=("interaction_name", "nunique"),
            n_pathways=("pathway_name", "nunique"),
        )
        .reset_index()
    )
    matrix = (
        summary.pivot(index="source", columns="target", values="summed_prob")
        .reindex(index=row_order, columns=col_order)
        .fillna(0.0)
    )

    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, ax = plt.subplots(figsize=(8.8, 5.3))
    sns.heatmap(
        matrix,
        cmap="mako",
        linewidths=0.6,
        linecolor="#E5E7EB",
        annot=True,
        fmt=".3g",
        square=True,
        cbar_kws={"label": "Summed communication probability", "shrink": 0.88, "pad": 0.03},
        ax=ax,
    )
    finalize_axes(ax, title=title, xlabel="", ylabel="", tight=False)
    ax.set_xticklabels([pretty_label(label) for label in col_order], rotation=22, ha="right")
    ax.set_yticklabels([pretty_label(label) for label in row_order], rotation=0)
    ax.tick_params(axis="x", pad=6)
    ax.tick_params(axis="y", rotation=0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return summary.sort_values("summed_prob", ascending=False, kind="stable")


def save_pathway_bubble_facets(
    interactions: pd.DataFrame,
    *,
    source_lineage: str,
    target_lineage: str,
    output_path: Path,
    top_n: int = 8,
) -> pd.DataFrame:
    title, source_order, target_order = LINEAGE_PAIR_LABELS[(source_lineage, target_lineage)]
    summary = (
        interactions[
            (interactions["source_lineage"] == source_lineage)
            & (interactions["target_lineage"] == target_lineage)
        ]
        .groupby(["source", "target", "pathway_name"], observed=True)
        .agg(
            summed_prob=("prob", "sum"),
            n_interactions=("interaction_name", "nunique"),
            top_lr=("interaction_name_2", lambda series: series.value_counts().index[0]),
        )
        .reset_index()
    )

    n_panels = len(source_order)
    ncols = 2 if n_panels > 1 else 1
    nrows = int(np.ceil(n_panels / ncols))
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.6 * ncols, 3.9 * nrows))
    axes = np.atleast_1d(axes).ravel()
    fig.subplots_adjust(left=0.10, right=0.88, top=0.90, bottom=0.16, wspace=0.42, hspace=0.46)

    panel_max = max(float(summary["summed_prob"].max()), 1e-9)
    norm = plt.Normalize(vmin=0.0, vmax=panel_max)
    cmap = plt.get_cmap("crest")
    scatter = None

    for ax, source in zip(axes, source_order, strict=False):
        source_frame = summary[summary["source"] == source].copy()
        pathway_order = (
            source_frame.groupby("pathway_name", observed=True)["summed_prob"]
            .sum()
            .sort_values(ascending=False)
            .head(top_n)
            .index.tolist()
        )
        source_frame = source_frame[source_frame["pathway_name"].isin(pathway_order)].copy()
        if source_frame.empty:
            ax.axis("off")
            continue

        source_frame["target"] = pd.Categorical(source_frame["target"], categories=target_order, ordered=True)
        source_frame["pathway_name"] = pd.Categorical(
            source_frame["pathway_name"],
            categories=pathway_order[::-1],
            ordered=True,
        )
        x_vals = source_frame["target"].cat.codes.to_numpy()
        y_vals = source_frame["pathway_name"].cat.codes.to_numpy()
        sizes = np.clip((source_frame["summed_prob"] / panel_max) * 520.0, 34.0, 300.0)

        scatter = ax.scatter(
            x_vals,
            y_vals,
            s=sizes,
            c=source_frame["summed_prob"],
            cmap=cmap,
            norm=norm,
            edgecolors="#374151",
            linewidths=0.45,
            alpha=0.95,
        )
        ax.set_xticks(range(len(target_order)))
        ax.set_xticklabels([pretty_label(label) for label in target_order], rotation=24, ha="right")
        ax.set_yticks(range(len(pathway_order)))
        ax.set_yticklabels(pathway_order[::-1])
        ax.set_xlim(-0.5, len(target_order) - 0.5)
        ax.set_ylim(-0.5, len(pathway_order) - 0.5)
        ax.grid(axis="x", linestyle=":", linewidth=0.6, color="#D1D5DB")
        ax.grid(axis="y", linestyle=":", linewidth=0.6, color="#E5E7EB")
        finalize_axes(ax, title=pretty_label(source), xlabel="", ylabel="", tight=False)
        ax.tick_params(axis="y", pad=4)

    for ax in axes[len(source_order):]:
        ax.axis("off")

    cbar_ax = fig.add_axes([0.905, 0.24, 0.015, 0.50])
    cbar = fig.colorbar(scatter, cax=cbar_ax)
    cbar.set_label("Summed communication probability")

    size_levels = [0.25, 0.5, 0.75, 1.0]
    handles = [
        Line2D(
            [],
            [],
            linestyle="",
            marker="o",
            markersize=np.sqrt(np.clip(level * 520.0, 34.0, 300.0)) / 1.9,
            markerfacecolor="#94A3B8",
            markeredgecolor="#475569",
            alpha=0.9,
            label=f"{panel_max * level:.3g}",
        )
        for level in size_levels
    ]
    fig.legend(
        handles=handles,
        title="Summed probability",
        loc="lower center",
        ncol=len(handles),
        frameon=False,
        bbox_to_anchor=(0.47, 0.03),
        columnspacing=1.6,
        handletextpad=0.6,
    )
    fig.suptitle(f"{title} top pathways", x=0.10, y=0.98, ha="left", fontsize=13, fontweight="semibold")
    fig.text(
        0.10,
        0.94,
        "Bubble size and color reflect pathway-level summed communication probability.",
        ha="left",
        va="center",
        fontsize=9,
        color="#4B5563",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return summary.sort_values(["source", "summed_prob"], ascending=[True, False], kind="stable")


def main() -> None:
    logger = setup_logging(ROOT / "logs/annotation/mac_t_cellchat_summary.log")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Loading CellChat cross-lineage interaction table")
    interactions = load_interactions()

    logger.info("Saving macrophage to CD8 summary heatmap and table")
    mac_to_cd8 = save_direction_heatmap(
        interactions,
        source_lineage="Myeloid",
        target_lineage="CD8",
        output_path=OUT_DIR / "macrophage_to_cd8_heatmap.png",
    )
    mac_to_cd8.to_csv(OUT_DIR / "macrophage_to_cd8_summary.csv", index=False)

    logger.info("Saving CD8 to macrophage summary heatmap and table")
    cd8_to_mac = save_direction_heatmap(
        interactions,
        source_lineage="CD8",
        target_lineage="Myeloid",
        output_path=OUT_DIR / "cd8_to_macrophage_heatmap.png",
    )
    cd8_to_mac.to_csv(OUT_DIR / "cd8_to_macrophage_summary.csv", index=False)

    logger.info("Saving pathway bubble plots")
    mac_pathways = save_pathway_bubble_facets(
        interactions,
        source_lineage="Myeloid",
        target_lineage="CD8",
        output_path=OUT_DIR / "macrophage_to_cd8_top_pathways.png",
    )
    mac_pathways.to_csv(OUT_DIR / "macrophage_to_cd8_pathways.csv", index=False)
    cd8_pathways = save_pathway_bubble_facets(
        interactions,
        source_lineage="CD8",
        target_lineage="Myeloid",
        output_path=OUT_DIR / "cd8_to_macrophage_top_pathways.png",
    )
    cd8_pathways.to_csv(OUT_DIR / "cd8_to_macrophage_pathways.csv", index=False)

    shutil.copyfile(CELLCHAT_DIR / "cellchat.png", OUT_DIR / "cellchat_overview.png")
    interactions.to_csv(OUT_DIR / "cellchat_myeloid_cd8_interactions.csv", index=False)
    logger.info("mac_t CellChat summary finished")


if __name__ == "__main__":
    main()
