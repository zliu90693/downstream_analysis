# %%
#! 检测哑铃形簇在pca中的出现情况

import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# %%
Amel = sc.read_h5ad("./data/Amel_dumbbell_marked.h5ad")
# %%
sc.pl.embedding(Amel, basis="X_pca_harmony_scran", color=["dumbbells"])
sc.pl.embedding(Amel, basis="X_pca_harmony_log1p", color=["dumbbells"])
sc.pl.embedding(Amel, basis="X_pca_harmony_pearson", color=["dumbbells"])
# %%
#! 小提琴图是以聚类簇为分组的，存在循环论证的嫌疑（聚类本身可能受深度影响）。建议直接对全体西方蜜蜂细胞画log1p_n_genes_by_counts的全局直方图
# sns.histplot(
#     Amel.obs["log1p_n_genes_by_counts"], 
#     # x="log1p_n_genes_by_counts", 
#     bins=100,
#     kde=True,)
# %%
# # 绘制直方图并叠加核密度估计曲线 (kde=True)
# sns.histplot(
#     data=Amel.obs["log1p_n_genes_by_counts"], 
#     bins=100,       # 设置分箱数，单细胞数据通常较多，100左右比较合适
#     kde=True,       # 叠加密度曲线，更容易看清分布形态（如单峰/双峰）
#     color="skyblue",
#     edgecolor="black"
# )
# %%
sns.histplot(
    Amel.obs["log1p_n_genes_by_counts"], 
    bins=100,
    kde=True,)
# %%
sns.displot(
    Amel.obs, 
    x="log1p_n_genes_by_counts",
    col="batch",
    bins=100,
    kde=True,)
# %%
Acer = sc.read_h5ad("../../Acer/data/7_cluster-output/concat.h5ad")
# %%
sns.histplot(
    Acer.obs["log1p_n_genes_by_counts"], 
    bins=100,
    kde=True,)
sns.displot(
    Acer.obs, 
    x="log1p_n_genes_by_counts",
    col="batch",
    bins=100,
    kde=True,)
# %%
#! 三种归一化方法得到的数据集在全部细胞中的直方图
Amel.obs["total_counts_scran"] = Amel.layers["scran_normalization"].sum(1)
Amel.obs["total_counts_log1p"] = Amel.layers["log1p_norm"].sum(1)
Amel.obs["total_counts_pearson"] = Amel.layers["analytic_pearson_residuals"].sum(1)
# %%
sns.histplot(
    Amel.obs["total_counts_scran"], 
    bins=100,
    kde=True,)
# %%
sns.displot(
    Amel.obs, 
    x="total_counts_scran",
    col="batch",
    bins=100,
    kde=True,)
# %%
sns.histplot(
    Amel.obs["total_counts_pearson"], 
    bins=100,
    kde=True,)
# %%
sns.displot(
    Amel.obs, 
    x="total_counts_pearson",
    col="batch",
    bins=100,
    kde=True,)
# %%
sns.histplot(
    Amel.obs["total_counts_log1p"], 
    bins=100,
    kde=True,)
# %%
sns.displot(
    Amel.obs, 
    x="total_counts_log1p",
    col="batch",
    bins=100,
    kde=True,)
# %%
