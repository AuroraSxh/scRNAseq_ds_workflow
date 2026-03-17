from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc

from crc_sc_integration.utils import (
    cluster_mean_scores,
    load_yaml,
    reference_major_label,
    save_heatmap,
    save_umap,
    score_gene_panels,
    setup_logging,
)


ROOT = Path(__file__).resolve().parents[1]
MARKERS = load_yaml(ROOT / "config/marker_panels.yaml")


def main() -> None:
    logger = setup_logging(ROOT / "logs/annotate/annotate_major.log")
    results_dir = ROOT / "results/annotation"
    adata = sc.read_h5ad(ROOT / "data/processed/crc_integrated_clusters.h5ad")

    logger.info("Scoring major cell-type marker panels")
    score_cols = score_gene_panels(adata, MARKERS["major_cell_types"], prefix="major")
    cluster_scores = cluster_mean_scores(adata, "cluster", score_cols)
    score_to_label = {}
    for panel_name in MARKERS["major_cell_types"]:
        panel_col = f"major_{panel_name.lower().replace(' ', '_').replace('/', '_')}"
        score_to_label[panel_col] = panel_name
    adata.obs["major_cell_type_score_only"] = adata.obs[score_cols].idxmax(axis=1).map(score_to_label)

    adata.obs["reference_major"] = adata.obs["source_global_cluster"].map(reference_major_label)
    adata.obs["major_cell_type"] = adata.obs["major_cell_type_score_only"]
    ref_mask = adata.obs["reference_major"].notna() & (adata.obs["reference_major"] != "Unknown")
    adata.obs.loc[ref_mask, "major_cell_type"] = adata.obs.loc[ref_mask, "reference_major"]

    reference = (
        adata.obs.loc[ref_mask]
        .groupby(["cluster", "reference_major"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    ref_totals = reference.groupby("cluster", observed=True)["n_cells"].sum().rename("ref_total")
    reference = reference.merge(ref_totals, on="cluster", how="left")
    reference["ref_fraction"] = reference["n_cells"] / reference["ref_total"]
    ref_best = reference.sort_values(["cluster", "ref_fraction"], ascending=[True, False]).drop_duplicates("cluster")

    cluster_annotations = []
    for cluster, scores in cluster_scores.iterrows():
        top_score_col = scores.idxmax()
        marker_label = score_to_label[top_score_col]
        ref_row = ref_best[ref_best["cluster"] == cluster]
        final_label = marker_label
        ref_label = "Unknown"
        ref_fraction = 0.0
        if not ref_row.empty:
            ref_label = ref_row.iloc[0]["reference_major"]
            ref_fraction = float(ref_row.iloc[0]["ref_fraction"])
            if ref_label != "Unknown" and ref_fraction >= 0.60:
                final_label = ref_label
        cluster_annotations.append(
            {
                "cluster": cluster,
                "major_cell_type": final_label,
                "top_marker_score_panel": marker_label,
                "reference_major": ref_label,
                "reference_fraction": ref_fraction,
            }
        )

    cluster_annotation_df = pd.DataFrame(cluster_annotations)

    logger.info("Writing annotation outputs")
    cluster_annotation_df.to_csv(results_dir / "cluster_major_annotation.csv", index=False)
    adata.obs["major_cell_type"].value_counts().rename_axis("major_cell_type").reset_index(name="n_cells").to_csv(
        results_dir / "major_cell_type_counts.csv", index=False
    )
    adata.write_h5ad(ROOT / "data/processed/crc_integrated_annotated.h5ad", compression="gzip")
    cluster_heatmap = cluster_scores.copy()
    cluster_heatmap.columns = [col.removeprefix("major_") for col in cluster_heatmap.columns]
    save_heatmap(
        cluster_heatmap,
        results_dir / "figures/major_marker_score_heatmap.png",
        "Cluster-level major cell-type scores",
    )
    save_umap(
        adata,
        ["cluster", "major_cell_type", "dataset"],
        results_dir / "figures/major_annotation_umap.png",
        "Major cell-type annotation",
    )
    logger.info("Major annotation step finished")


if __name__ == "__main__":
    main()
