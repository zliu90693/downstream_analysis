# %%
import scanpy as sc
import anndata
# %%
Hsal = sc.read_h5ad("../Sheng_SA_2020_Hsal/data/7_cluster-output/concat.h5ad")
Acer = sc.read_h5ad("../Acer/data/7_cluster-output/concat.h5ad")
Amel = sc.read_h5ad("../Zhang_iScience_2022_Amel/data/7_cluster-output/concat.h5ad")
# %%
# def sc.pl.embedding(
#     adata: anndata.AnnData,
#     basis: str,
#     color,
# ) -> None:
#     sc.pl.embedding(
#         adata,
#         basis=basis,
#         color=color,
#     )
# %%
# Hsal
sc.pl.embedding(Hsal, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Hsal, basis="X_umap_harmony_scran", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Hsal, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
# Acer
sc.pl.embedding(Acer, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Acer, basis="X_umap_harmony_scran", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Acer, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
# Amel
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
