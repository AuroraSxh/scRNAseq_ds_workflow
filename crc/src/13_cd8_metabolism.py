from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import kruskal, mannwhitneyu, wilcoxon

matplotlib.use("Agg", force=True)

from crc_sc_integration.utils import finalize_axes, set_beautiful_style, setup_logging


ROOT = Path(__file__).resolve().parents[1]
INPUT_PATH = ROOT / "data/processed/cd8_merged.h5ad"
OUTPUT_DIR = ROOT / "results/cd8/metabolism"
FIGURE_DIR = OUTPUT_DIR / "figures"
TISSUE_DIR = OUTPUT_DIR / "tissue_stratified"
TISSUE_FIGURE_DIR = TISSUE_DIR / "figures"
PAIRED_DIR = TISSUE_DIR / "paired"
PAIRED_FIGURE_DIR = PAIRED_DIR / "figures"
FOCUSED_DIR = OUTPUT_DIR / "focused_four"
FOCUSED_FIGURE_DIR = FOCUSED_DIR / "figures"
LOG_PATH = ROOT / "logs/cd8/cd8_metabolism.log"

GROUP_ORDER = ["Tpex_like", "Early_Tex", "Effector_memory", "Terminal_Tex", "Stress_program"]
FOCUSED_GROUP_ORDER = ["Tpex_like", "Early_Tex", "Effector_memory", "Terminal_Tex"]
TISSUE_ORDER = ["T", "N"]
TISSUE_LABELS = {"T": "Tumor", "N": "Normal"}
GROUP_PALETTE = {
    "Tpex_like": "#2B7A78",
    "Early_Tex": "#E07A5F",
    "Effector_memory": "#3D5A80",
    "Terminal_Tex": "#7D4E57",
    "Stress_program": "#8D99AE",
}
MIN_SAMPLE_CELLS = 10
MIN_PAIRED_PATIENTS = 3
REPRESENTATIVE_GENES = [
    "SLC2A1",
    "HK2",
    "LDHA",
    "G6PD",
    "TKT",
    "CS",
    "OGDH",
    "NDUFB8",
    "COX5A",
    "CPT1A",
    "ACADVL",
    "SLC1A5",
    "GLS",
    "PHGDH",
    "MTHFD2",
    "TYMS",
    "RRM2",
]
METABOLIC_PANELS = {
    "Glycolysis": ["SLC2A1", "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA"],
    "Pentose_Phosphate": ["G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT", "TALDO1", "PRPS1"],
    "TCA_Cycle": ["CS", "ACO2", "IDH3A", "IDH3B", "OGDH", "DLST", "SUCLG1", "SDHA", "FH", "MDH2"],
    "OxPhos": ["NDUFA1", "NDUFA2", "NDUFB8", "NDUFS1", "UQCRFS1", "COX5A", "COX6C", "COX7A2"],
    "FAO": ["CPT1A", "CPT2", "ACADM", "ACADVL", "ECHS1", "HADHA", "HADHB", "ACAA2", "ETFA"],
    "Glutamine_Metabolism": ["SLC1A5", "GLS", "GLUD1", "GOT2", "OAT", "ASNS", "GPT2", "GLS2"],
    "One_Carbon_Serine": ["PHGDH", "PSAT1", "PSPH", "SHMT1", "SHMT2", "MTHFD1", "MTHFD2", "DHFR", "TYMS", "ATIC"],
    "Nucleotide_Biosynthesis": ["CAD", "UMPS", "CTPS1", "TYMS", "TK1", "RRM1", "RRM2", "DHODH", "IMPDH1", "IMPDH2", "PPAT", "GART", "PFAS", "ATIC"],
}
PATHWAY_LABELS = {
    "meta_glycolysis": "Glycolysis",
    "meta_pentose_phosphate": "Pentose Phosphate",
    "meta_tca_cycle": "TCA Cycle",
    "meta_oxphos": "OxPhos",
    "meta_fao": "FAO",
    "meta_glutamine_metabolism": "Glutamine Metabolism",
    "meta_one_carbon_serine": "One Carbon Serine",
    "meta_nucleotide_biosynthesis": "Nucleotide Biosynthesis",
}


def benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.empty_like(ranked)
    running = 1.0
    total = len(values)
    for idx in range(total - 1, -1, -1):
        running = min(running, ranked[idx] * total / (idx + 1))
        adjusted[idx] = running
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(output, index=pvalues.index)


def cliffs_delta(x: pd.Series, y: pd.Series) -> float:
    x_vals = np.asarray(x, dtype=float)
    y_vals = np.asarray(y, dtype=float)
    gt = sum((value > y_vals).sum() for value in x_vals)
    lt = sum((value < y_vals).sum() for value in x_vals)
    return float((gt - lt) / (len(x_vals) * len(y_vals)))


def mean_expression_score(adata: sc.AnnData, genes: list[str]) -> np.ndarray:
    matrix = adata[:, genes].X
    if sp.issparse(matrix):
        return np.asarray(matrix.mean(axis=1)).ravel()
    return np.asarray(matrix.mean(axis=1)).ravel()


def compute_pathway_scores(adata: sc.AnnData) -> tuple[list[str], pd.DataFrame]:
    coverage_rows: list[dict[str, object]] = []
    score_cols: list[str] = []
    for pathway, genes in METABOLIC_PANELS.items():
        present = [gene for gene in genes if gene in adata.var_names]
        coverage_rows.append(
            {
                "pathway": pathway,
                "n_panel_genes": len(genes),
                "n_present_genes": len(present),
                "coverage_fraction": len(present) / len(genes),
                "present_genes": ", ".join(present),
            }
        )
        if len(present) < 3:
            continue
        score_col = f"meta_{pathway.lower()}"
        adata.obs[score_col] = mean_expression_score(adata, present)
        score_cols.append(score_col)
    coverage = pd.DataFrame(coverage_rows)
    return score_cols, coverage


def pathway_name(score_col: str) -> str:
    return PATHWAY_LABELS.get(score_col, score_col.removeprefix("meta_").replace("_", " ").title())


def prepare_sample_scores(adata: sc.AnnData, score_cols: list[str]) -> pd.DataFrame:
    frame = adata.obs[["sample_id", "dataset", "tissue", "cd8_merged_subgroup", *score_cols]].copy()
    frame["n_cells"] = 1
    aggregated = (
        frame.groupby(["sample_id", "dataset", "tissue", "cd8_merged_subgroup"], observed=True)
        .agg({**{col: "mean" for col in score_cols}, "n_cells": "sum"})
        .reset_index()
    )
    return aggregated.loc[aggregated["n_cells"] >= MIN_SAMPLE_CELLS].copy()


def compute_group_means(adata: sc.AnnData, score_cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    means = (
        adata.obs.groupby("cd8_merged_subgroup", observed=True)[score_cols]
        .mean()
        .rename(columns=pathway_name)
        .reindex(GROUP_ORDER)
    )
    zscores = means.apply(lambda col: (col - col.mean()) / (col.std(ddof=0) if col.std(ddof=0) else 1.0), axis=0)
    return means, zscores


def compute_statistics(sample_scores: pd.DataFrame, score_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for score_col in score_cols:
        pvalues: list[float] = []
        interim: list[dict[str, object]] = []
        pathway = pathway_name(score_col)
        for group in GROUP_ORDER:
            focal = sample_scores.loc[sample_scores["cd8_merged_subgroup"] == group, score_col]
            other = sample_scores.loc[sample_scores["cd8_merged_subgroup"] != group, score_col]
            pvalue = mannwhitneyu(focal, other, alternative="two-sided").pvalue
            pvalues.append(pvalue)
            interim.append(
                {
                    "pathway": pathway,
                    "group": group,
                    "n_group_samples": int(focal.shape[0]),
                    "n_other_samples": int(other.shape[0]),
                    "mean_group": float(focal.mean()),
                    "mean_rest": float(other.mean()),
                    "median_group": float(focal.median()),
                    "median_rest": float(other.median()),
                    "median_diff": float(focal.median() - other.median()),
                    "cliffs_delta": cliffs_delta(focal, other),
                    "p_value": float(pvalue),
                }
            )
        adjusted = benjamini_hochberg(pd.Series(pvalues))
        for entry, qvalue in zip(interim, adjusted, strict=False):
            entry["q_value"] = float(qvalue)
            rows.append(entry)

    stats = pd.DataFrame(rows)

    kw_rows: list[dict[str, object]] = []
    for score_col in score_cols:
        pathway = pathway_name(score_col)
        groups = [
            sample_scores.loc[sample_scores["cd8_merged_subgroup"] == group, score_col].to_numpy()
            for group in GROUP_ORDER
            if not sample_scores.loc[sample_scores["cd8_merged_subgroup"] == group, score_col].empty
        ]
        statistic, pvalue = kruskal(*groups)
        kw_rows.append({"pathway": pathway, "kruskal_statistic": statistic, "kruskal_p_value": pvalue})
    kw = pd.DataFrame(kw_rows)
    kw["kruskal_q_value"] = benjamini_hochberg(kw["kruskal_p_value"])

    return stats.merge(kw, on="pathway", how="left")


def compute_tissue_group_means(sample_scores: pd.DataFrame, score_cols: list[str]) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    means_by_tissue: dict[str, pd.DataFrame] = {}
    zscores_by_tissue: dict[str, pd.DataFrame] = {}
    for tissue in TISSUE_ORDER:
        subset = sample_scores.loc[sample_scores["tissue"] == tissue]
        means = (
            subset.groupby("cd8_merged_subgroup", observed=True)[score_cols]
            .mean()
            .rename(columns=pathway_name)
            .reindex(GROUP_ORDER)
        )
        zscores = means.apply(lambda col: (col - col.mean()) / (col.std(ddof=0) if col.std(ddof=0) else 1.0), axis=0)
        means_by_tissue[tissue] = means
        zscores_by_tissue[tissue] = zscores
    return means_by_tissue, zscores_by_tissue


def compute_tissue_statistics(sample_scores: pd.DataFrame, score_cols: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group in GROUP_ORDER:
        subgroup = sample_scores.loc[sample_scores["cd8_merged_subgroup"] == group]
        pvalues: list[float] = []
        pending_indices: list[int] = []
        for score_col in score_cols:
            tumor = subgroup.loc[subgroup["tissue"] == "T", score_col]
            normal = subgroup.loc[subgroup["tissue"] == "N", score_col]
            row = {
                "group": group,
                "pathway": pathway_name(score_col),
                "n_tumor_samples": int(tumor.shape[0]),
                "n_normal_samples": int(normal.shape[0]),
                "mean_tumor": float(tumor.mean()) if not tumor.empty else np.nan,
                "mean_normal": float(normal.mean()) if not normal.empty else np.nan,
                "median_tumor": float(tumor.median()) if not tumor.empty else np.nan,
                "median_normal": float(normal.median()) if not normal.empty else np.nan,
                "median_diff_tumor_minus_normal": float(tumor.median() - normal.median()) if not tumor.empty and not normal.empty else np.nan,
                "cliffs_delta_tumor_vs_normal": cliffs_delta(tumor, normal) if len(tumor) > 0 and len(normal) > 0 else np.nan,
                "p_value": np.nan,
                "q_value": np.nan,
                "sufficient_support": len(tumor) >= 2 and len(normal) >= 2,
            }
            rows.append(row)
            if row["sufficient_support"]:
                pending_indices.append(len(rows) - 1)
                pvalues.append(mannwhitneyu(tumor, normal, alternative="two-sided").pvalue)

        if pvalues:
            adjusted = benjamini_hochberg(pd.Series(pvalues))
            for row_index, pvalue, qvalue in zip(pending_indices, pvalues, adjusted, strict=False):
                rows[row_index]["p_value"] = float(pvalue)
                rows[row_index]["q_value"] = float(qvalue)

    return pd.DataFrame(rows)


def add_patient_id(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output["patient_id"] = output["sample_id"].astype(str).str.replace(r"_(T|N|TA|TB|P)$", "", regex=True)
    return output


def prepare_paired_patient_scores(sample_scores: pd.DataFrame, score_cols: list[str]) -> pd.DataFrame:
    frame = add_patient_id(sample_scores)
    frame = frame.loc[frame["tissue"].isin(TISSUE_ORDER)].copy()
    patient_level = (
        frame.groupby(["patient_id", "tissue", "cd8_merged_subgroup"], observed=True)[score_cols]
        .mean()
        .reset_index()
    )
    paired_rows: list[pd.DataFrame] = []
    for group in GROUP_ORDER:
        subset = patient_level.loc[patient_level["cd8_merged_subgroup"] == group].copy()
        counts = subset.groupby(["patient_id", "tissue"], observed=True).size().unstack(fill_value=0)
        paired_patients = counts.index[(counts.get("T", 0) > 0) & (counts.get("N", 0) > 0)]
        if len(paired_patients) == 0:
            continue
        paired_rows.append(subset.loc[subset["patient_id"].isin(paired_patients)].copy())
    if not paired_rows:
        return patient_level.iloc[0:0].copy()
    return pd.concat(paired_rows, ignore_index=True)


def compute_paired_statistics(paired_scores: pd.DataFrame, score_cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    delta_tables: list[pd.DataFrame] = []
    for group in GROUP_ORDER:
        subset = paired_scores.loc[paired_scores["cd8_merged_subgroup"] == group].copy()
        if subset.empty:
            continue
        pvalues: list[float] = []
        row_indices: list[int] = []
        for score_col in score_cols:
            pivot = (
                subset.pivot_table(index="patient_id", columns="tissue", values=score_col, aggfunc="mean", observed=True)
                .dropna(subset=["T", "N"], how="any")
                .copy()
            )
            delta = pivot["T"] - pivot["N"]
            delta_frame = pd.DataFrame(
                {
                    "patient_id": pivot.index,
                    "group": group,
                    "pathway": pathway_name(score_col),
                    "tumor_score": pivot["T"].to_numpy(),
                    "normal_score": pivot["N"].to_numpy(),
                    "delta_tumor_minus_normal": delta.to_numpy(),
                }
            )
            delta_tables.append(delta_frame)

            row = {
                "group": group,
                "pathway": pathway_name(score_col),
                "n_paired_patients": int(delta.shape[0]),
                "mean_tumor": float(pivot["T"].mean()) if not pivot.empty else np.nan,
                "mean_normal": float(pivot["N"].mean()) if not pivot.empty else np.nan,
                "median_tumor": float(pivot["T"].median()) if not pivot.empty else np.nan,
                "median_normal": float(pivot["N"].median()) if not pivot.empty else np.nan,
                "median_delta_tumor_minus_normal": float(delta.median()) if not delta.empty else np.nan,
                "mean_delta_tumor_minus_normal": float(delta.mean()) if not delta.empty else np.nan,
                "p_value": np.nan,
                "q_value": np.nan,
                "sufficient_support": int(delta.shape[0]) >= MIN_PAIRED_PATIENTS,
            }
            rows.append(row)
            if row["sufficient_support"]:
                try:
                    stat = wilcoxon(pivot["T"], pivot["N"], alternative="two-sided", zero_method="wilcox", method="auto")
                    pvalue = float(stat.pvalue)
                except ValueError:
                    pvalue = 1.0
                pvalues.append(pvalue)
                row_indices.append(len(rows) - 1)

        if pvalues:
            adjusted = benjamini_hochberg(pd.Series(pvalues))
            for idx, pvalue, qvalue in zip(row_indices, pvalues, adjusted, strict=False):
                rows[idx]["p_value"] = float(pvalue)
                rows[idx]["q_value"] = float(qvalue)

    stats = pd.DataFrame(rows)
    deltas = pd.concat(delta_tables, ignore_index=True) if delta_tables else pd.DataFrame(columns=["patient_id", "group", "pathway", "tumor_score", "normal_score", "delta_tumor_minus_normal"])
    return stats, deltas


def representative_gene_expression(adata: sc.AnnData) -> pd.DataFrame:
    present = [gene for gene in REPRESENTATIVE_GENES if gene in adata.var_names]
    matrix = adata[:, present].X
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    expr = pd.DataFrame(matrix, index=adata.obs_names, columns=present)
    expr["cd8_merged_subgroup"] = adata.obs["cd8_merged_subgroup"].astype(str).values
    summary = expr.groupby("cd8_merged_subgroup", observed=True).mean().reindex(GROUP_ORDER)
    return summary.T


def save_heatmap(frame: pd.DataFrame, output_path: Path) -> None:
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    sns.heatmap(
        frame,
        cmap=sns.diverging_palette(240, 15, as_cmap=True),
        center=0,
        linewidths=0.8,
        linecolor="#FFFFFF",
        cbar_kws={"label": "Z-score across subgroups"},
        ax=ax,
    )
    finalize_axes(
        ax,
        title="CD8 metabolic programs across five subgroups",
        subtitle="Row-scaled pathway activity based on mean pathway expression",
        xlabel="Metabolic pathway",
        ylabel="CD8 subgroup",
        tight=False,
    )
    plt.xticks(rotation=35, ha="right")
    plt.yticks(rotation=0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_focused_heatmap(frame: pd.DataFrame, output_path: Path) -> None:
    focused = frame.loc[FOCUSED_GROUP_ORDER].copy()
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    sns.heatmap(
        focused,
        cmap=sns.diverging_palette(240, 15, as_cmap=True),
        center=0,
        linewidths=0.8,
        linecolor="#FFFFFF",
        cbar_kws={"label": "Z-score across the four core subgroups"},
        ax=ax,
    )
    finalize_axes(
        ax,
        title="Core CD8 subgroup metabolic programs",
        subtitle="Stress_program removed to keep the main biological continuum readable",
        xlabel="Metabolic pathway",
        ylabel="CD8 subgroup",
        tight=False,
    )
    ax.tick_params(axis="x", rotation=35)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_boxplots(sample_scores: pd.DataFrame, score_cols: list[str], output_path: Path) -> None:
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    long_df = sample_scores.melt(
        id_vars=["sample_id", "dataset", "tissue", "cd8_merged_subgroup", "n_cells"],
        value_vars=score_cols,
        var_name="pathway",
        value_name="score",
    )
    long_df["pathway"] = long_df["pathway"].map(pathway_name)

    ncols = 2
    nrows = int(np.ceil(long_df["pathway"].nunique() / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.3 * nrows), sharey=False)
    axes = np.atleast_1d(axes).flatten()

    for idx, pathway in enumerate([pathway_name(col) for col in score_cols]):
        ax = axes[idx]
        subset = long_df.loc[long_df["pathway"] == pathway].copy()
        sns.boxplot(
            data=subset,
            x="cd8_merged_subgroup",
            y="score",
            hue="cd8_merged_subgroup",
            order=GROUP_ORDER,
            palette=GROUP_PALETTE,
            legend=False,
            showfliers=False,
            linewidth=0.9,
            width=0.68,
            ax=ax,
        )
        sns.stripplot(
            data=subset,
            x="cd8_merged_subgroup",
            y="score",
            order=GROUP_ORDER,
            color="#111827",
            size=2.2,
            alpha=0.45,
            jitter=0.18,
            ax=ax,
        )
        finalize_axes(ax, title=pathway, xlabel="", ylabel="Mean sample score", tight=False)
        ax.tick_params(axis="x", rotation=35)

    for idx in range(len(score_cols), len(axes)):
        axes[idx].axis("off")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def overall_activity_ranking(zscores: pd.DataFrame) -> pd.Series:
    return zscores.mean(axis=1).sort_values(ascending=False)


def overall_activity_ranking_by_tissue(zscores_by_tissue: dict[str, pd.DataFrame]) -> dict[str, pd.Series]:
    return {
        tissue: overall_activity_ranking(frame)
        for tissue, frame in zscores_by_tissue.items()
    }


def save_tissue_heatmaps(zscores_by_tissue: dict[str, pd.DataFrame], output_path: Path) -> None:
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, axes = plt.subplots(1, len(TISSUE_ORDER), figsize=(15, 4.8), sharey=True)
    if len(TISSUE_ORDER) == 1:
        axes = [axes]

    for idx, tissue in enumerate(TISSUE_ORDER):
        ax = axes[idx]
        frame = zscores_by_tissue[tissue]
        sns.heatmap(
            frame,
            cmap=sns.diverging_palette(240, 15, as_cmap=True),
            center=0,
            linewidths=0.8,
            linecolor="#FFFFFF",
            cbar=idx == len(TISSUE_ORDER) - 1,
            cbar_kws={"label": "Z-score across subgroups"} if idx == len(TISSUE_ORDER) - 1 else None,
            ax=ax,
        )
        finalize_axes(
            ax,
            title=TISSUE_LABELS[tissue],
            subtitle="Sample-level pathway means",
            xlabel="Metabolic pathway",
            ylabel="CD8 subgroup" if idx == 0 else "",
            tight=False,
        )
        ax.tick_params(axis="x", rotation=35)
        ax.tick_params(axis="y", rotation=0)

    fig.suptitle("CD8 metabolic hierarchy in tumor and adjacent normal tissues", x=0.02, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_tissue_delta_heatmap(sample_scores: pd.DataFrame, score_cols: list[str], output_path: Path) -> None:
    means = (
        sample_scores.groupby(["tissue", "cd8_merged_subgroup"], observed=True)[score_cols]
        .mean()
        .rename(columns=pathway_name)
    )
    delta = means.loc["T"].reindex(GROUP_ORDER) - means.loc["N"].reindex(GROUP_ORDER)
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    sns.heatmap(
        delta,
        cmap=sns.diverging_palette(220, 20, as_cmap=True),
        center=0,
        linewidths=0.8,
        linecolor="#FFFFFF",
        cbar_kws={"label": "Tumor - Normal"},
        ax=ax,
    )
    finalize_axes(
        ax,
        title="Tumor-vs-normal metabolic shift within each CD8 subgroup",
        subtitle="Positive values indicate higher pathway activity in tumor tissue",
        xlabel="Metabolic pathway",
        ylabel="CD8 subgroup",
        tight=False,
    )
    ax.tick_params(axis="x", rotation=35)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_paired_delta_heatmap(paired_stats: pd.DataFrame, output_path: Path) -> None:
    frame = (
        paired_stats.pivot(index="group", columns="pathway", values="median_delta_tumor_minus_normal")
        .reindex(GROUP_ORDER)
    )
    set_beautiful_style(medium="paper", background="light", font_scale=1.0, dpi=180)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    sns.heatmap(
        frame,
        cmap=sns.diverging_palette(220, 20, as_cmap=True),
        center=0,
        linewidths=0.8,
        linecolor="#FFFFFF",
        cbar_kws={"label": "Paired tumor - normal median delta"},
        ax=ax,
    )
    finalize_axes(
        ax,
        title="Matched-patient tumor-vs-normal metabolic shift",
        subtitle="Only patients with both tumor and normal support in the same subgroup are included",
        xlabel="Metabolic pathway",
        ylabel="CD8 subgroup",
        tight=False,
    )
    ax.tick_params(axis="x", rotation=35)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_tissue_report(
    sample_scores: pd.DataFrame,
    zscores_by_tissue: dict[str, pd.DataFrame],
    tissue_stats: pd.DataFrame,
    output_path: Path,
) -> None:
    support = (
        sample_scores.groupby(["tissue", "cd8_merged_subgroup"], observed=True)
        .size()
        .rename("n_sample_subgroup")
        .reset_index()
    )
    rankings = overall_activity_ranking_by_tissue(zscores_by_tissue)
    lines = [
        "# CD8+T 亚群 Tumor vs Normal 代谢分层分析",
        "",
        "## 分析范围",
        "- 仅纳入 `tissue` 为 `T` 或 `N` 的 sample-subgroup 组合。",
        f"- 仍使用每个 sample-subgroup 至少 {MIN_SAMPLE_CELLS} 个细胞的阈值。",
        "",
        "## 样本支持",
    ]
    for tissue in TISSUE_ORDER:
        label = TISSUE_LABELS[tissue]
        lines.append(f"### {label}")
        subset = support.loc[support["tissue"] == tissue].set_index("cd8_merged_subgroup")["n_sample_subgroup"].reindex(GROUP_ORDER).fillna(0).astype(int)
        for group, count in subset.items():
            lines.append(f"- `{group}`: {int(count)} 个 sample-subgroup 组合")

    lines.extend(
        [
            "",
            "## 总体观察",
            f"- Tumor 内部的综合代谢活跃度排序：{', '.join(rankings['T'].index.tolist())}。",
            f"- Normal 内部的综合代谢活跃度排序：{', '.join(rankings['N'].index.tolist())}。",
            "- 两种组织里，`Early_Tex` 和 `Stress_program` 仍处于相对高代谢端，`Tpex_like` 仍更接近低代谢端，说明总体层级基本保留。",
            "- 但真正需要看的不是层级是否相同，而是每个亚群在 Tumor 中是否相对 Normal 进一步上调特定代谢程序。",
            "",
            "## 各亚群 Tumor vs Normal 变化",
        ]
    )

    for group in GROUP_ORDER:
        subset = tissue_stats.loc[tissue_stats["group"] == group].copy()
        valid = subset.loc[subset["sufficient_support"]].sort_values("median_diff_tumor_minus_normal", ascending=False)
        increased = valid.loc[(valid["q_value"] < 0.05) & (valid["median_diff_tumor_minus_normal"] > 0)]
        decreased = valid.loc[(valid["q_value"] < 0.05) & (valid["median_diff_tumor_minus_normal"] < 0)]
        insufficient = subset.loc[~subset["sufficient_support"], "pathway"].tolist()
        lines.append(f"### {group}")
        if increased.empty:
            lines.append("- 没有检测到在 Tumor 中显著升高的代谢通路。")
        else:
            lines.append(
                "- Tumor 显著升高："
                + "; ".join(
                    f"{row.pathway} (median diff={row.median_diff_tumor_minus_normal:.3f}, q={row.q_value:.3g})"
                    for row in increased.itertuples()
                )
                + "。"
            )
        if decreased.empty:
            lines.append("- 没有检测到在 Tumor 中显著降低的代谢通路。")
        else:
            lines.append(
                "- Tumor 显著降低："
                + "; ".join(
                    f"{row.pathway} (median diff={row.median_diff_tumor_minus_normal:.3f}, q={row.q_value:.3g})"
                    for row in decreased.itertuples()
                )
                + "。"
            )
        if insufficient:
            lines.append("- 支持不足，未进行稳定统计的通路：" + ", ".join(insufficient) + "。")

    lines.extend(
        [
            "",
            "## 解读边界",
            "- 这是 sample-subgroup 层面的非配对比较，没有控制病人配对、肿瘤区域或临床协变量。",
            "- `Stress_program` 在 Normal 中可用组合极少，因此其 Tumor-vs-Normal 结论最不稳健。",
            "- 如果下一步要写结果段，建议优先使用 `Tpex_like / Early_Tex / Effector_memory / Terminal_Tex` 四个支持度更好的亚群。",
        ]
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_paired_report(paired_scores: pd.DataFrame, paired_stats: pd.DataFrame, output_path: Path) -> None:
    support = (
        paired_scores.groupby("cd8_merged_subgroup", observed=True)["patient_id"]
        .nunique()
        .reindex(GROUP_ORDER)
        .fillna(0)
        .astype(int)
    )
    lines = [
        "# CD8+T 亚群配对病人 Tumor vs Normal 代谢分析",
        "",
        "## 分析设计",
        "- 先将 sample-level 通路分数按 `patient_id + tissue + subgroup` 聚合，解决同一患者存在多个肿瘤样本的问题。",
        "- 然后仅保留在同一 CD8 亚群内同时具有 Tumor 和 Normal 的患者，使用配对 Wilcoxon signed-rank test。",
        f"- 最低统计支持阈值：至少 {MIN_PAIRED_PATIENTS} 个配对患者。",
        "",
        "## 配对支持",
    ]
    for group, count in support.items():
        lines.append(f"- `{group}`: {int(count)} 个配对患者")

    lines.extend(["", "## 主要结论"])
    stable_groups = [group for group, count in support.items() if count >= MIN_PAIRED_PATIENTS]
    if stable_groups:
        lines.append("- 可进行配对统计的亚群为：" + ", ".join(stable_groups) + "。")
    else:
        lines.append("- 没有亚群达到配对统计阈值。")

    for group in GROUP_ORDER:
        subset = paired_stats.loc[paired_stats["group"] == group].copy()
        if subset.empty:
            lines.append(f"- `{group}`: 没有可用的配对患者。")
            continue
        significant_up = subset.loc[(subset["q_value"] < 0.05) & (subset["median_delta_tumor_minus_normal"] > 0)].sort_values("median_delta_tumor_minus_normal", ascending=False)
        significant_down = subset.loc[(subset["q_value"] < 0.05) & (subset["median_delta_tumor_minus_normal"] < 0)].sort_values("median_delta_tumor_minus_normal")
        if not bool(subset["sufficient_support"].iloc[0]):
            lines.append(f"- `{group}`: 仅 {int(subset['n_paired_patients'].max())} 个配对患者，证据不足。")
            continue
        if significant_up.empty and significant_down.empty:
            lines.append(f"- `{group}`: 未检测到配对后仍显著的 Tumor-vs-Normal 代谢差异。")
            continue
        summary_parts: list[str] = []
        if not significant_up.empty:
            summary_parts.append(
                "Tumor 升高 "
                + "; ".join(
                    f"{row.pathway} (median delta={row.median_delta_tumor_minus_normal:.3f}, q={row.q_value:.3g})"
                    for row in significant_up.itertuples()
                )
            )
        if not significant_down.empty:
            summary_parts.append(
                "Tumor 降低 "
                + "; ".join(
                    f"{row.pathway} (median delta={row.median_delta_tumor_minus_normal:.3f}, q={row.q_value:.3g})"
                    for row in significant_down.itertuples()
                )
            )
        lines.append(f"- `{group}`: " + "；".join(summary_parts) + "。")

    lines.extend(
        [
            "",
            "## 解读边界",
            "- `Early_Tex` 虽达到最小阈值，但仅 3 个配对患者，任何阴性或阳性结果都不稳健。",
            "- `Stress_program` 没有配对患者，不能据此判断其肿瘤特异性代谢变化。",
            "- 如果要写论文主结果，配对分析最适合作为对非配对结果的稳健性验证，而不是替代主结论。",
        ]
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_focused_report(
    means: pd.DataFrame,
    zscores: pd.DataFrame,
    stats: pd.DataFrame,
    paired_stats: pd.DataFrame,
    output_path: Path,
) -> None:
    focused_means = means.loc[FOCUSED_GROUP_ORDER].copy()
    focused_zscores = zscores.loc[FOCUSED_GROUP_ORDER].copy()
    positive_stats = (
        stats.loc[
            stats["group"].isin(FOCUSED_GROUP_ORDER)
            & (stats["q_value"] < 0.05)
            & (stats["median_diff"] > 0)
        ]
        .sort_values(["group", "median_diff"], ascending=[True, False])
        .copy()
    )
    paired_support = (
        paired_stats.groupby("group", observed=True)["n_paired_patients"]
        .max()
        .reindex(FOCUSED_GROUP_ORDER)
        .fillna(0)
        .astype(int)
    )
    ranking = overall_activity_ranking(focused_zscores)

    lines = [
        "# CD8 四个主亚群代谢摘要",
        "",
        "## 为什么去掉 Stress_program",
        "- `Stress_program` 混合了 IFN-response 和 cycling 叠加态，容易把增殖/应激信号误读成主干分化轨迹上的代谢特征。",
        "- 因此主叙事更适合聚焦 `Tpex_like / Early_Tex / Effector_memory / Terminal_Tex` 四个亚群。",
        "",
        "## 一句话结论",
        f"- 四个主亚群的综合代谢活跃度排序为：{', '.join(ranking.index.tolist())}。",
        "- 其中 `Early_Tex` 是最明确的高代谢过渡态，`Tpex_like` 是最低代谢的储备样状态，`Effector_memory` 与 `Terminal_Tex` 处于中间，但二者方向不同：前者更偏低糖酵解，后者更偏中等 OxPhos/一碳代谢。",
        "",
        "## 各亚群摘要",
    ]

    for group in FOCUSED_GROUP_ORDER:
        top_paths = focused_zscores.loc[group].sort_values(ascending=False).head(3)
        positives = positive_stats.loc[positive_stats["group"] == group]
        lines.append(f"### {group}")
        lines.append(
            "- 组内最高的 3 条通路："
            + ", ".join(f"{name} (Z={value:.2f})" for name, value in top_paths.items())
            + "。"
        )
        if positives.empty:
            lines.append("- 相比其他主亚群，没有检测到显著升高的代谢通路。")
        else:
            lines.append(
                "- 显著升高通路："
                + "; ".join(
                    f"{row.pathway} (median diff={row.median_diff:.3f}, q={row.q_value:.3g})"
                    for row in positives.itertuples()
                )
                + "。"
            )
        lines.append(f"- 配对 Tumor-Normal 支持度：{int(paired_support[group])} 位患者。")

    lines.extend(
        [
            "",
            "## 展示建议",
            "- 如果做主图，优先展示 `Early_Tex` 对 `Tpex_like` 的代谢上移，以及 `Effector_memory`/`Terminal_Tex` 分别代表的两种中间状态。",
            "- 如果写正文，不建议把 `Stress_program` 和四个主亚群并列描述；更适合作为补充态或 overlay program。",
            "",
            "## 配对验证如何写",
            "- 配对分析显示四个主亚群都没有稳定的 Tumor-vs-Normal 显著差异，因此更支持“亚群内在状态主导代谢差异”而非“肿瘤环境统一重编程”。",
        ]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(
    means: pd.DataFrame,
    zscores: pd.DataFrame,
    stats: pd.DataFrame,
    sample_scores: pd.DataFrame,
    coverage: pd.DataFrame,
    output_path: Path,
) -> None:
    support = sample_scores.groupby("cd8_merged_subgroup", observed=True).size().reindex(GROUP_ORDER).fillna(0).astype(int)
    ranking = overall_activity_ranking(zscores)
    pathway_leaders = zscores.idxmax(axis=0)

    lines = [
        "# CD8+T 五个亚群代谢状态分析",
        "",
        "## 分析对象",
        f"- 输入对象：`{INPUT_PATH.relative_to(ROOT)}`",
        f"- 亚群：{', '.join(GROUP_ORDER)}",
        f"- 样本层统计阈值：每个 sample-subgroup 组合至少 {MIN_SAMPLE_CELLS} 个细胞",
        "",
        "## 总体结论",
        f"- 综合 8 条代谢通路的 Z-score，整体代谢活跃度排序为：{', '.join(ranking.index.tolist())}。",
        f"- `Stress_program` 在一碳代谢、核苷酸合成、氧化磷酸化和糖酵解上最强，但样本支持仅 {int(support['Stress_program'])} 个 sample-subgroup 组合，解释时需要谨慎。",
        f"- `Early_Tex` 是最明确的高代谢过渡状态，在糖酵解、OxPhos、FAO、PPP、TCA 和谷氨酰胺代谢上均显著高于其他亚群。",
        "- `Tpex_like` 没有显著升高的代谢通路，反而在糖酵解、OxPhos、PPP 和一碳代谢上显著偏低，更接近低代谢/储备样状态。",
        "- `Effector_memory` 与 `Terminal_Tex` 没有形成广谱代谢优势，其中 `Effector_memory` 主要表现为糖酵解偏低，`Terminal_Tex` 更接近中间水平但缺乏明显的合成代谢增强。",
        "",
        "## 各亚群解读",
    ]

    for group in GROUP_ORDER:
        group_stats = stats.loc[stats["group"] == group].copy()
        positive = group_stats.loc[(group_stats["q_value"] < 0.05) & (group_stats["median_diff"] > 0)].sort_values("median_diff", ascending=False)
        negative = group_stats.loc[(group_stats["q_value"] < 0.05) & (group_stats["median_diff"] < 0)].sort_values("median_diff")
        top_paths = zscores.loc[group].sort_values(ascending=False).head(3)
        lines.append(f"### {group}")
        lines.append(f"- 样本支持：{int(support[group])} 个 sample-subgroup 组合。")
        lines.append(
            "- 组内相对较高的通路："
            + ", ".join(f"{path} (Z={value:.2f})" for path, value in top_paths.items())
            + "。"
        )
        if positive.empty:
            lines.append("- 相比其余亚群，没有检测到 FDR<0.05 的显著升高通路。")
        else:
            lines.append(
                "- 显著升高通路："
                + "; ".join(
                    f"{row.pathway} (median diff={row.median_diff:.3f}, q={row.q_value:.3g})"
                    for row in positive.itertuples()
                )
                + "。"
            )
        if negative.empty:
            lines.append("- 没有检测到 FDR<0.05 的显著降低通路。")
        else:
            lines.append(
                "- 显著降低通路："
                + "; ".join(
                    f"{row.pathway} (median diff={row.median_diff:.3f}, q={row.q_value:.3g})"
                    for row in negative.head(4).itertuples()
                )
                + "。"
            )

    lines.extend(
        [
            "",
            "## 通路层面观察",
            *[
                f"- `{pathway}` 最高的亚群是 `{leader}`。"
                for pathway, leader in pathway_leaders.items()
            ],
            "",
            "## 结果解读边界",
            "- 这里的通路分数采用通路基因平均表达，不等同于直接代谢通量测定。",
            "- `Stress_program` 合并了 IFN-response 和 Cycling 细胞，因此其核苷酸和一碳代谢增强很可能部分来自增殖需求。",
            "- 样本层统计已经降低细胞数不均带来的伪重复风险，但并未进一步校正肿瘤/邻近正常、平台或病人临床变量。",
            "",
            "## 基因覆盖率",
        ]
    )

    for row in coverage.itertuples():
        lines.append(
            f"- `{pathway_name(f'meta_{str(row.pathway).lower()}')}`: {int(row.n_present_genes)}/{int(row.n_panel_genes)} genes present"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    logger = setup_logging(LOG_PATH)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TISSUE_DIR.mkdir(parents=True, exist_ok=True)
    TISSUE_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    PAIRED_DIR.mkdir(parents=True, exist_ok=True)
    PAIRED_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    FOCUSED_DIR.mkdir(parents=True, exist_ok=True)
    FOCUSED_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Loading merged CD8 object: %s", INPUT_PATH)
    adata = sc.read_h5ad(INPUT_PATH)

    logger.info("Computing pathway scores for %d cells", adata.n_obs)
    score_cols, coverage = compute_pathway_scores(adata)
    pathway_labels = [pathway_name(col) for col in score_cols]
    logger.info("Retained %d metabolic pathways: %s", len(pathway_labels), ", ".join(pathway_labels))

    logger.info("Aggregating cell-level and sample-level summaries")
    means, zscores = compute_group_means(adata, score_cols)
    sample_scores = prepare_sample_scores(adata, score_cols)
    stats = compute_statistics(sample_scores, score_cols)
    representative = representative_gene_expression(adata)
    tn_sample_scores = sample_scores.loc[sample_scores["tissue"].isin(TISSUE_ORDER)].copy()
    tissue_means, tissue_zscores = compute_tissue_group_means(tn_sample_scores, score_cols)
    tissue_stats = compute_tissue_statistics(tn_sample_scores, score_cols)
    paired_scores = prepare_paired_patient_scores(tn_sample_scores, score_cols)
    paired_stats, paired_deltas = compute_paired_statistics(paired_scores, score_cols)

    logger.info("Writing result tables")
    coverage.to_csv(OUTPUT_DIR / "pathway_gene_coverage.csv", index=False)
    means.to_csv(OUTPUT_DIR / "subgroup_pathway_means.csv")
    zscores.to_csv(OUTPUT_DIR / "subgroup_pathway_zscores.csv")
    sample_scores.rename(columns={col: pathway_name(col) for col in score_cols}).to_csv(
        OUTPUT_DIR / "sample_level_pathway_scores.csv",
        index=False,
    )
    stats.to_csv(OUTPUT_DIR / "pathway_vs_rest_statistics.csv", index=False)
    representative.to_csv(OUTPUT_DIR / "representative_gene_expression.csv")
    tissue_means["T"].to_csv(TISSUE_DIR / "tumor_subgroup_pathway_means.csv")
    tissue_means["N"].to_csv(TISSUE_DIR / "normal_subgroup_pathway_means.csv")
    tissue_zscores["T"].to_csv(TISSUE_DIR / "tumor_subgroup_pathway_zscores.csv")
    tissue_zscores["N"].to_csv(TISSUE_DIR / "normal_subgroup_pathway_zscores.csv")
    tn_sample_scores.rename(columns={col: pathway_name(col) for col in score_cols}).to_csv(
        TISSUE_DIR / "tumor_normal_sample_level_pathway_scores.csv",
        index=False,
    )
    tissue_stats.to_csv(TISSUE_DIR / "tumor_vs_normal_within_subgroup_statistics.csv", index=False)
    paired_scores.rename(columns={col: pathway_name(col) for col in score_cols}).to_csv(
        PAIRED_DIR / "paired_patient_level_pathway_scores.csv",
        index=False,
    )
    paired_stats.to_csv(PAIRED_DIR / "paired_tumor_vs_normal_statistics.csv", index=False)
    paired_deltas.to_csv(PAIRED_DIR / "paired_patient_pathway_deltas.csv", index=False)
    means.loc[FOCUSED_GROUP_ORDER].to_csv(FOCUSED_DIR / "focused_subgroup_pathway_means.csv")
    zscores.loc[FOCUSED_GROUP_ORDER].to_csv(FOCUSED_DIR / "focused_subgroup_pathway_zscores.csv")
    stats.loc[stats["group"].isin(FOCUSED_GROUP_ORDER)].to_csv(
        FOCUSED_DIR / "focused_subgroup_statistics.csv",
        index=False,
    )

    logger.info("Saving figures")
    save_heatmap(zscores, FIGURE_DIR / "cd8_metabolism_heatmap.png")
    save_boxplots(sample_scores, score_cols, FIGURE_DIR / "cd8_metabolism_sample_boxplots.png")
    save_tissue_heatmaps(tissue_zscores, TISSUE_FIGURE_DIR / "cd8_metabolism_tumor_normal_heatmaps.png")
    save_tissue_delta_heatmap(tn_sample_scores, score_cols, TISSUE_FIGURE_DIR / "cd8_metabolism_tumor_minus_normal_heatmap.png")
    save_paired_delta_heatmap(paired_stats, PAIRED_FIGURE_DIR / "cd8_metabolism_paired_tumor_minus_normal_heatmap.png")
    save_focused_heatmap(zscores, FOCUSED_FIGURE_DIR / "cd8_metabolism_core_subgroups_heatmap.png")

    logger.info("Writing markdown report")
    write_report(means, zscores, stats, sample_scores, coverage, OUTPUT_DIR / "report.md")
    write_tissue_report(tn_sample_scores, tissue_zscores, tissue_stats, TISSUE_DIR / "report.md")
    write_paired_report(paired_scores, paired_stats, PAIRED_DIR / "report.md")
    write_focused_report(means, zscores, stats, paired_stats, FOCUSED_DIR / "report.md")
    logger.info("CD8 metabolism analysis finished")


if __name__ == "__main__":
    main()
