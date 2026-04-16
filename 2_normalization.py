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

sns.histplot(F1_adata.obs["total_counts"], bins=100, kde=False)

# %%
# ---------------------------------------- Shifted logarithm ----------------------------------------

scales_counts = sc.pp.normalize_total(F1_adata, target_sum=None, inplace=False)
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

analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(F1_adata, inplace=False)

