# %%
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.sparse as sp
# %%
Amel = sc.read_h5ad("./data/Amel_1.25.h5ad")
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00", "leiden_harmony_scran_res1.25"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["leiden_harmony_pearson_res0.50", "leiden_harmony_pearson_res1.00", "leiden_harmony_pearson_res1.25"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["leiden_harmony_log1p_res0.50", "leiden_harmony_log1p_res1.00", "leiden_harmony_log1p_res1.25"],legend_loc="on data")
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=[ "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
# %%
def check_sparsity_pair(adata, cluster_key, cluster_a, cluster_b):
    # mask_a = adata.obs[cluster_key] == cluster_a
    # mask_b = adata.obs[cluster_key] == cluster_b
    mask_a = (adata.obs[cluster_key] == cluster_a).values
    mask_b = (adata.obs[cluster_key] == cluster_b).values
    #! 加上 .values 将 Pandas Series 转换为 NumPy 数组，下文的 counts 是一个 SciPy 稀疏矩阵，无法通过 bool 索引取子集
    
    counts = adata.layers["counts"]  # 原始计数
    
    X_a = counts[mask_a]
    X_b = counts[mask_b]
    
    if sp.issparse(X_a):
        X_a = X_a.toarray()
        X_b = X_b.toarray()
    
    # 每个细胞的基因检出数（非零基因数）
    genes_per_cell_a = (X_a > 0).sum(axis=1)
    genes_per_cell_b = (X_b > 0).sum(axis=1)
    
    # 每个细胞的dropout率（零值比例）
    dropout_a = (X_a == 0).mean(axis=1)
    dropout_b = (X_b == 0).mean(axis=1)
    
    print(f"Average number of genes detected in cluster {cluster_a}: {genes_per_cell_a.mean():.1f} ± {genes_per_cell_a.std():.1f}")
    print(f"Average number of genes detected in cluster {cluster_b}: {genes_per_cell_b.mean():.1f} ± {genes_per_cell_b.std():.1f}")
    print(f"Average dropout rate for cluster {cluster_a}: {dropout_a.mean():.3f}")
    print(f"Average dropout rate for cluster {cluster_b}: {dropout_b.mean():.3f}")
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(genes_per_cell_a, bins=50, alpha=0.6, label=f"cluster {cluster_a}")
    axes[0].hist(genes_per_cell_b, bins=50, alpha=0.6, label=f"cluster {cluster_b}")
    axes[0].set_xlabel("Number of genes detected per cell")
    axes[0].legend()
    
    axes[1].hist(dropout_a, bins=50, alpha=0.6, label=f"cluster {cluster_a}")
    axes[1].hist(dropout_b, bins=50, alpha=0.6, label=f"cluster {cluster_b}")
    axes[1].set_xlabel("dropout rate")
    axes[1].legend()
    
    plt.tight_layout()
    plt.show()
# %%
#! leiden_harmony_log1p_res1.25
#! taillight & head
#! 27 & 46
#! 28 & 33
#! 44 & 51
#! 24 & 41
#! 36 & 47
#! 34 & 35
#! 22 & 19
#! 20 & 37
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="27",
    cluster_b="46"
)
# Average number of genes detected in cluster 27: 3051.7 ± 594.9
# Average number of genes detected in cluster 46: 1067.8 ± 480.3
# Average dropout rate for cluster 27: 0.736
# Average dropout rate for cluster 46: 0.908
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="28",
    cluster_b="33"
)
# Average number of genes detected in cluster 28: 3378.6 ± 565.4
# Average number of genes detected in cluster 33: 1369.0 ± 674.8
# Average dropout rate for cluster 28: 0.707
# Average dropout rate for cluster 33: 0.881
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="44",
    cluster_b="51"
)
# Average number of genes detected in cluster 44: 2606.0 ± 604.0
# Average number of genes detected in cluster 51: 877.0 ± 384.7
# Average dropout rate for cluster 44: 0.774
# Average dropout rate for cluster 51: 0.924
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="24",
    cluster_b="41"
)
# Average number of genes detected in cluster 24: 3217.0 ± 554.4
# Average number of genes detected in cluster 41: 1284.5 ± 611.4
# Average dropout rate for cluster 24: 0.721
# Average dropout rate for cluster 41: 0.889
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="36",
    cluster_b="47"
)
# Average number of genes detected in cluster 36: 2722.0 ± 608.2
# Average number of genes detected in cluster 47: 1035.9 ± 498.2
# Average dropout rate for cluster 36: 0.764
# Average dropout rate for cluster 47: 0.910
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="34",
    cluster_b="35"
)
# Average number of genes detected in cluster 34: 2552.4 ± 655.5
# Average number of genes detected in cluster 35: 957.2 ± 420.9
# Average dropout rate for cluster 34: 0.779
# Average dropout rate for cluster 35: 0.917
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="22",
    cluster_b="19"
)
# Average number of genes detected in cluster 22: 3316.8 ± 650.4
# Average number of genes detected in cluster 19: 1141.8 ± 609.4
# Average dropout rate for cluster 22: 0.713
# Average dropout rate for cluster 19: 0.901
# %%
check_sparsity_pair(
    Amel,
    cluster_key="leiden_harmony_log1p_res1.25",
    cluster_a="20",
    cluster_b="37"
)
# Average number of genes detected in cluster 22: 3316.8 ± 650.4
# Average number of genes detected in cluster 19: 1141.8 ± 609.4
# Average dropout rate for cluster 22: 0.713
# Average dropout rate for cluster 19: 0.901
# %%
