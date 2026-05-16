# %%

import logging
import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.sparse import issparse

# %%

F1_adata = sc.read_h5ad("./Zhang_iScience_2022_Amel/h5_QC/ft_F1.h5ad")

# %% 

# filt genes appear in < 20 Cells

print(f"Total number of genes: {F1_adata.n_vars}")

# Min 20 cells - filters out 0 count genes
sc.pp.filter_genes(F1_adata, min_cells=20)
print(f"Number of genes after cell filter: {F1_adata.n_vars}")

# %%

sns.histplot(F1_adata.obs["total_counts"], bins=100, kde=False)

# %%
# ---------------------------------------- Shifted logarithm ----------------------------------------

scales_counts = sc.pp.normalize_total(F1_adata, target_sum=None, inplace=False) # inplace=False, 防止Normalize覆盖X
F1_adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

# %%

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(F1_adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(F1_adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
# 注意, Normalize的对象不是total_counts, 而是X(每个细胞每个基因的UMI计数), 所以dim是1
axes[1].set_title("Shifted logarithm")
plt.show()

# %%

# ---------------------------------------- Scran ----------------------------------------

from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import subprocess

F1_adata_pp = F1_adata.copy()
sc.pp.normalize_total(F1_adata_pp)
sc.pp.log1p(F1_adata_pp)
sc.pp.pca(F1_adata_pp, n_comps=15)
sc.pp.neighbors(F1_adata_pp)
sc.tl.leiden(
    F1_adata_pp, key_added="groups", flavor="igraph", n_iterations=2, directed=False
)

# %%

F1_data_mat = F1_adata_pp.X.T
if issparse(F1_data_mat):
    if F1_data_mat.nnz > 2**31 - 1:
        F1_data_mat = F1_data_mat.tocoo()
    else:
        F1_data_mat = F1_data_mat.tocsc()
    # Convert sparse matrix to dense numpy array
    # F1_data_mat = F1_data_mat.toarray()
mmwrite("./Zhang_iScience_2022_Amel/metadata/F1_data_mat.mtx", F1_data_mat)

# %%

F1_adata_pp.obs[["groups"]].to_csv("./Zhang_iScience_2022_Amel/metadata/input_groups.csv")

# %%

del F1_adata_pp

# %%

### 2.1_run_scran.R

# %%

import pandas as pd

size_factors = pd.read_csv("./Zhang_iScience_2022_Amel/metadata/size_factors.csv", index_col=0).squeeze()

# %%

F1_adata.obs["size_factors"] = size_factors.values
scran = F1_adata.X / F1_adata.obs["size_factors"].values[:, None]
scran_logged = np.log1p(scran)
F1_adata.layers["scran_normalization"] = csr_matrix(scran_logged)

# ---------------------------------------- Pearson ---------------------------------------- 

# %%

analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(F1_adata, inplace=False)  # inplace=False, 防止Normalize覆盖X
F1_adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])

# %%

# check
mat = F1_adata.layers["analytic_pearson_residuals"]
# 检查 inf 和 nan 的数量
print("inf 数量:", np.isinf(mat.data).sum())
print("nan 数量:", np.isnan(mat.data).sum())

# 检查全0基因数量
zero_genes = np.asarray(F1_adata.X.sum(0)).flatten() == 0
print("全零基因数量:", zero_genes.sum())
print("全零基因名称:", F1_adata.var_names[zero_genes].tolist())

# %%

fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(F1_adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(
    F1_adata.layers["analytic_pearson_residuals"].sum(1), bins=100, kde=False, ax=axes[1]
)
axes[1].set_title("Analytic Pearson residuals")
plt.show()

# %%

F1_adata.write("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm_F1.h5ad")

# %%
