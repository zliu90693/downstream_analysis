# %%
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
# %%
Amel = sc.read_h5ad("./data/Amel_dumbbell_marked.h5ad")
# %%
def run_leiden(
    adata: anndata.AnnData,
    neo_key: str,
    neighbors_key: str,
    resolution: float = 1.0,
    flavor: str = "igraph",
    n_iterations: int = 2,
) -> None:
    sc.tl.leiden(
        adata, 
        key_added=f"leiden_{neo_key}_res{resolution:.2f}",  # 防止浮点数转字符串时可能出现 res0.6000000000000001
        neighbors_key=neighbors_key,
        resolution=resolution, 
        flavor=flavor, 
        n_iterations=n_iterations
    )
# %%
run_leiden(Amel, neo_key="harmony_scran", neighbors_key="neighbors_harmony_scran", resolution=1.25, flavor="igraph", n_iterations=2)
run_leiden(Amel, neo_key="harmony_pearson", neighbors_key="neighbors_harmony_pearson", resolution=1.25, flavor="igraph", n_iterations=2)
run_leiden(Amel, neo_key="harmony_log1p", neighbors_key="neighbors_harmony_log1p", resolution=1.25, flavor="igraph", n_iterations=2)

# %%
# sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00", "leiden_harmony_scran_res1.25"],legend_loc="on data")
# sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["leiden_harmony_pearson_res0.50", "leiden_harmony_pearson_res1.00", "leiden_harmony_pearson_res1.25"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["leiden_harmony_log1p_res0.50", "leiden_harmony_log1p_res1.00", "leiden_harmony_log1p_res1.25"],legend_loc="on data")
# %%
cluster_col = 'leiden_harmony_log1p_res1.00'
batch_col = 'batch'

plt.figure(figsize=(12, 6))
# multiple='fill' 是关键，它会将每个柱子的总高度归一化为 1 (即 100%)
# discrete=True 表示 X 轴是离散类别
sns.histplot(
    data=Amel.obs, 
    x=cluster_col, 
    hue=batch_col, 
    multiple="fill", 
    discrete=True, 
    shrink=0.8,      # 让柱子之间有点缝隙，更好看
    palette="tab20"  # 如果 batch 很多，用 tab20 颜色区分度高
)

plt.ylabel("Proportion")
plt.title("Batch Composition per Cluster (100% Stacked)")
# plt.xticks(rotation=45, ha='right')
# plt.legend(title='Batch', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
# %%
