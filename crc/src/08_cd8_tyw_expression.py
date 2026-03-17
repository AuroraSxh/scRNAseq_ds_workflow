from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc

from crc_sc_integration.utils import (
    get_present_genes,
    save_gene_bubbleplot,
    save_gene_expression_embedding,
    save_gene_violinplot,
    setup_logging,
)


ROOT = Path(__file__).resolve().parents[1]
GENES = ["TYW1", "TYW2", "TYW3", "TYW4"]


def main() -> None:
    logger = setup_logging(ROOT / "logs/cd8/cd8_tyw_expression.log")
    output_dir = ROOT / "results/cd8/tyw_genes"
    figure_dir = output_dir / "figures"

    logger.info("Loading CD8 object for TYW gene visualization")
    adata = sc.read_h5ad(ROOT / "data/processed/cd8_subclusters.h5ad")

    gene_table = pd.DataFrame(
        {
            "gene": GENES,
            "present_in_cd8_object": [gene in adata.var_names for gene in GENES],
        }
    )
    present_genes = get_present_genes(adata, GENES)
    gene_table["used_for_plotting"] = gene_table["gene"].isin(present_genes)

    output_dir.mkdir(parents=True, exist_ok=True)
    gene_table.to_csv(output_dir / "tyw_gene_availability.csv", index=False)
    logger.info("Requested genes present in CD8 object: %s", ", ".join(present_genes) if present_genes else "none")

    if not present_genes:
        raise RuntimeError("None of TYW1/TYW2/TYW3/TYW4 are present in the CD8 object.")

    cluster_order = (
        adata.obs["cd8_cluster_label"].cat.categories.tolist()
        if hasattr(adata.obs["cd8_cluster_label"], "cat")
        else sorted(adata.obs["cd8_cluster_label"].astype(str).unique().tolist())
    )

    logger.info("Saving TYW bubble plot")
    save_gene_bubbleplot(
        adata=adata,
        groupby="cd8_cluster_label",
        genes=GENES,
        output_path=figure_dir / "tyw_genes_bubbleplot.png",
        title="TYW gene expression across CD8 subgroups",
        groups_order=cluster_order,
    )

    basis = "umap"
    logger.info("Saving TYW UMAP expression panels using basis=%s", basis)
    save_gene_expression_embedding(
        adata=adata,
        genes=GENES,
        output_path=figure_dir / "tyw_genes_umap_expression.png",
        title="TYW gene expression on CD8 UMAP",
        basis=basis,
        ncols=2,
    )

    logger.info("Saving TYW violin plots")
    save_gene_violinplot(
        adata=adata,
        groupby="cd8_cluster_label",
        genes=GENES,
        output_path=figure_dir / "tyw_genes_violin.png",
        title="TYW gene expression across CD8 subgroups",
        groups_order=cluster_order,
    )
    logger.info("TYW gene visualization finished")


if __name__ == "__main__":
    main()
