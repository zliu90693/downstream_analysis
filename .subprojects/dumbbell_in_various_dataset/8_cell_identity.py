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
Amel.write_h5ad("./data/Amel_1.25.h5ad")
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["leiden_harmony_pearson_res0.50", "leiden_harmony_pearson_res1.00"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["leiden_harmony_log1p_res0.50", "leiden_harmony_log1p_res1.00"],legend_loc="on data")
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["dumbbells"], legend_loc="on data")
# %%
#! 寻找细胞簇marker，以确定哑铃中的细胞是单一细胞类型被拆分，还是两个细胞类型被合并
#! -------------------------------------------- 使用log1p_0.5和scran_0.5 -------------------------------------------- 
#! 对应关系：
#! log1p_0.5, scran_0.5
#! 14, 35 & 23
#! 18, 4 & 31
#! 25, 39 & 27
#! 15, 32 & 18
#! 16, 20 & 22
#! 22, 29 & 17
#! 21, 36 & 28
#! 24, 33
#! 19, 37 & 24
#! 13, 12 & 34
#! 26, 38
# %%
#? log1p 0.50
Amel.X = Amel.layers["log1p_norm"]
sc.tl.rank_genes_groups(
    Amel, groupby="leiden_harmony_log1p_res0.50", method="wilcoxon", key_added="dea_log1p_0.50"
)
# %%
sc.tl.dendrogram(
    Amel,
    groupby="leiden_harmony_log1p_res0.50",
)
sc.pl.rank_genes_groups_dotplot(
    Amel, groupby="leiden_harmony_log1p_res0.50", standard_scale="var", n_genes=5, key="dea_log1p_0.50"
)
# %%
sc.tl.filter_rank_genes_groups(
    Amel,
    min_in_group_fraction=0.2,
    max_out_group_fraction=0.2,
    key="dea_log1p_0.50",
    key_added="dea_log1p_0.50_filtered",
)
# %%
