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
dumbell_log1p_cluster = ["18","14","25","16","15","22","21","19","13","26","24"]
no_dumbbell_log1p_cluster = list(set(Amel.obs["leiden_harmony_log1p_res0.50"].unique()) - set(dumbell_log1p_cluster))
# %%
Amel.obs["total_counts_scran"] = Amel.layers["scran_normalization"].sum(1)
Amel.obs["total_counts_log1p"] = Amel.layers["log1p_norm"].sum(1)
Amel.obs["total_counts_pearson"] = Amel.layers["analytic_pearson_residuals"].sum(1)
# %%
# %%
def visualize_cluster_counts_better(
    adata: anndata.AnnData,
    normalize_total_counts: str,
    dumbell_list: list,
    cluster_clip: str = "leiden_harmony_log1p_res0.50",
    col_wrap: int = 4,
) -> None:
    # 1. 高效构建 tidy DataFrame（避免逐行 append）
    frames = []
    for single_cluster in dumbell_list:
        cell_sums = np.array(adata[adata.obs[cluster_clip] == single_cluster].obs[normalize_total_counts]).flatten()
        frames.append(pd.DataFrame({
            "cluster": single_cluster,
            "cell_sum": cell_sums,
        }))
    df = pd.concat(frames, ignore_index=True)

    # 2. 一行 displot 搞定分面
    g = sns.displot(
        data=df,
        x="cell_sum",
        col="cluster",
        col_wrap=col_wrap,
        bins=100,
        kde=True,
        height=3,
        aspect=1.3,
    )
    g.set_titles("Cluster {col_name}")
    g.set_axis_labels("UMI counts per cell", "Frequency")
    plt.tight_layout()
    plt.show()
# %%
# %%
#! ---------------------------------- 四种总深度中的哑铃簇 ----------------------------------
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts_log1p", 
    dumbell_list=dumbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts_scran",
    dumbell_list=dumbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts",
    dumbell_list=dumbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts_pearson",
    dumbell_list=dumbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
#! ---------------------------------- 四种总深度中的非哑铃簇 ----------------------------------
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts_log1p", 
    dumbell_list=no_dumbbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts_scran",
    dumbell_list=no_dumbbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts",
    dumbell_list=no_dumbbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
visualize_cluster_counts_better(
    Amel, 
    normalize_total_counts = "total_counts_pearson",
    dumbell_list=no_dumbbell_log1p_cluster, 
    cluster_clip = "leiden_harmony_log1p_res0.50"
)
# %%
# %%
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_n_genes_by_counts', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_total_counts', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
# %%
#! ---------------------------------------------- All in Violin ---------------------------------------------- 
#? log1p
sc.pl.embedding(Amel, 
                basis="X_umap_harmony_log1p", 
                color=["leiden_harmony_log1p_res0.50"],
                legend_loc="on data")
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_n_genes_by_counts', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_total_counts', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='total_counts_log1p', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()

# %%
#? scran
sc.pl.embedding(Amel, 
                basis="X_umap_harmony_scran", 
                color=["leiden_harmony_scran_res0.50"],
                legend_loc="on data")
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_n_genes_by_counts', 
    groupby='leiden_harmony_scran_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_total_counts', 
    groupby='leiden_harmony_scran_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='total_counts_scran', 
    groupby='leiden_harmony_scran_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
#? pearson
sc.pl.embedding(Amel, 
                basis="X_umap_harmony_pearson", 
                color=["leiden_harmony_pearson_res0.50"],
                legend_loc="on data")
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_n_genes_by_counts', 
    groupby='leiden_harmony_pearson_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_total_counts', 
    groupby='leiden_harmony_pearson_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='total_counts_pearson', 
    groupby='leiden_harmony_pearson_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
