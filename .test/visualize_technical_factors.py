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
# def visualize(
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
# sc.pl.embedding(Hsal, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.25", "leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00"])
# sc.pl.embedding(Acer, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.25", "leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00"])
# sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.25", "leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00"])
# %%
sc.pl.embedding(Hsal, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50"])
sc.pl.embedding(Acer, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50"])
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50"])
# %%
sc.pl.embedding(Acer, basis="X_umap_harmony_log1p", color=["leiden_harmony_scran_res0.50"])
sc.pl.embedding(Acer, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50"])
sc.pl.embedding(Acer, basis="X_umap_harmony_pearson", color=["leiden_harmony_scran_res0.50"])

# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["leiden_harmony_log1p_res0.50"])
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50"])
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["leiden_harmony_pearson_res0.50"])
# %%
# Amel
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])

# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
import scrublet as scr

Amel.X = Amel.layers["counts"].copy()
sc.pp.scrublet(Amel)
# %%
