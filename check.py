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
# res: 0.5
# 19
# 22
# 17

# res: 1.0
# 17 31
# 28 40
# 26 37
# %%
import numpy as np
from scipy.stats import spearmanr, pearsonr
import anndata
# %%
target_cluster = '22'  # 替换为您的哑铃簇ID
mask = Amel_h5ad.obs['leiden_harmony_log1p_res0.50'] == target_cluster
# 计算n_genes与UMAP坐标的相关性
umap_coords = Amel_h5ad.obsm['X_umap_harmony_log1p'][mask]
n_genes = Amel_h5ad.obs.loc[mask, 'n_genes_by_counts']
# Spearman相关（对非线性关系更鲁棒）
rho_umap1, p1 = spearmanr(umap_coords[:, 0], n_genes)
rho_umap2, p2 = spearmanr(umap_coords[:, 1], n_genes)
# %%
print(f"哑铃簇内 n_genes vs UMAP1: ρ={rho_umap1:+.3f}, p={p1:.2e}")
print(f"哑铃簇内 n_genes vs UMAP2: ρ={rho_umap2:+.3f}, p={p2:.2e}")
# %%
def check_correlation(
    adata: anndata.AnnData,
    target_cluster: str,
    obs_cluster_key: str,
    obsm_cluster_key: str,
    param2be_checked: str,
) -> None:
    mask = adata.obs[obs_cluster_key] == target_cluster
    umap_coords = adata.obsm[obsm_cluster_key][mask]
    param = adata.obs.loc[mask, param2be_checked]

    rho_umap1, p1 = spearmanr(umap_coords[:, 0], param)
    rho_umap2, p2 = spearmanr(umap_coords[:, 1], param)

    print(f"{param2be_checked} vs UMAP1: ρ={rho_umap1:+.3f}, p={p1:.2e}")
    print(f"{param2be_checked} vs UMAP2: ρ={rho_umap2:+.3f}, p={p2:.2e}")
# %%
# ------------------------------------- 22 ------------------------------------
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='22',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='log1p_n_genes_by_counts'
)
# log1p_n_genes_by_counts vs UMAP1: ρ=-0.919, p=0.00e+00
# log1p_n_genes_by_counts vs UMAP2: ρ=-0.780, p=0.00e+00
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='22',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='log1p_total_counts'
)
# log1p_total_counts vs UMAP1: ρ=-0.884, p=0.00e+00
# log1p_total_counts vs UMAP2: ρ=-0.760, p=0.00e+00
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='22',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='pct_counts_mt'
)
# pct_counts_mt vs UMAP1: ρ=+0.111, p=7.01e-09
# pct_counts_mt vs UMAP2: ρ=+0.135, p=2.14e-12
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='22',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='size_factors'
)
# size_factors vs UMAP1: ρ=-0.923, p=0.00e+00
# size_factors vs UMAP2: ρ=-0.774, p=0.00e+00
# %%
# ------------------------------------- 28 ------------------------------------
check_correlation(
    adata=Amel_h5ad,
    target_cluster='28',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='log1p_n_genes_by_counts'
)
# n_genes_by_counts vs UMAP1: ρ=+0.776, p=0.00e+00
# n_genes_by_counts vs UMAP2: ρ=+0.122, p=1.18e-24
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='28',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='log1p_total_counts'
)
# log1p_total_counts vs UMAP1: ρ=+0.775, p=0.00e+00
# log1p_total_counts vs UMAP2: ρ=+0.058, p=1.51e-06
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='28',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='pct_counts_mt'
)
# pct_counts_mt vs UMAP1: ρ=-0.104, p=4.27e-18
# pct_counts_mt vs UMAP2: ρ=-0.431, p=0.00e+00
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='28',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='size_factors'
)
# size_factors vs UMAP1: ρ=+0.731, p=0.00e+00
# size_factors vs UMAP2: ρ=+0.136, p=3.87e-30
# %%
# ------------------------------------- 12 ------------------------------------
check_correlation(
    adata=Amel_h5ad,
    target_cluster='12',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='log1p_n_genes_by_counts'
)
# log1p_n_genes_by_counts vs UMAP1: ρ=+0.829, p=0.00e+00
# log1p_n_genes_by_counts vs UMAP2: ρ=-0.895, p=0.00e+00
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='12',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='log1p_total_counts'
)
# log1p_total_counts vs UMAP1: ρ=+0.787, p=0.00e+00
# log1p_total_counts vs UMAP2: ρ=-0.856, p=0.00e+00
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='12',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='pct_counts_mt'
)
# pct_counts_mt vs UMAP1: ρ=-0.407, p=3.71e-87
# pct_counts_mt vs UMAP2: ρ=+0.352, p=3.95e-64
# %%
check_correlation(
    adata=Amel_h5ad,
    target_cluster='12',
    obs_cluster_key='leiden_harmony_log1p_res0.50',
    obsm_cluster_key='X_umap_harmony_log1p',
    param2be_checked='size_factors'
)
# size_factors vs UMAP1: ρ=+0.877, p=0.00e+00
# size_factors vs UMAP2: ρ=-0.932, p=0.00e+00
# %%
