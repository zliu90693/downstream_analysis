# 解决哑铃形问题: 
# %%
import scanpy as sc
# %%
project_name = "Zhang_iScience_2022_Amel"
Amel_h5ad = sc.read_h5ad(f"./{project_name}/data/7_cluster-output/concat.h5ad")
# %%
sc.pl.embedding(Amel_h5ad, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.25", "leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00"], legend_loc="on data")
sc.pl.embedding(Amel_h5ad, basis="X_umap_harmony_pearson", color=["leiden_harmony_pearson_res0.25", "leiden_harmony_pearson_res0.50", "leiden_harmony_pearson_res1.00"], legend_loc="on data")
sc.pl.embedding(Amel_h5ad, basis="X_umap_harmony_log1p", color=["leiden_harmony_log1p_res0.25", "leiden_harmony_log1p_res0.50", "leiden_harmony_log1p_res1.00"], legend_loc="on data")
# %%
sc.pl.embedding(Amel_h5ad, basis="X_umap_harmony_scran", color=["log1p_total_counts", "log1p_n_genes_by_counts", "n_genes_by_counts", "pct_counts_mt", "batch", "size_factors"])
sc.pl.embedding(Amel_h5ad, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "n_genes_by_counts", "pct_counts_mt", "batch", "size_factors"])
sc.pl.embedding(Amel_h5ad, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "n_genes_by_counts", "pct_counts_mt", "batch", "size_factors"])
# %%
