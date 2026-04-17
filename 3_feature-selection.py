# %%

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.sparse import issparse
from scipy.io import mmwrite

# %%

F1_adata = sc.read_h5ad("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm_F1.h5ad")

# %%

X_sparse = F1_adata.X.T

# %%

if issparse(X_sparse):
    if X_sparse.nnz > 2**31 - 1:
        X_sparse = X_sparse.tocoo()
    else:
        X_sparse = X_sparse.tocsc()

# %%

mmwrite("./Zhang_iScience_2022_Amel/metadata/X_sparse.mtx", X_sparse)

# %%

F1_adata.obs.to_csv("./Zhang_iScience_2022_Amel/metadata/obs.csv")
F1_adata.var.to_csv("./Zhang_iScience_2022_Amel/metadata/var.csv")

# %%

### run 3.1_run_devianceFeatureSelection.R

# %%

import pandas as pd
deviance_results = pd.read_csv("./Zhang_iScience_2022_Amel/metadata/deviance_results.csv", index_col=0)
# F1_adata.var["deviance"] = deviance_results["deviance"].values

# %%

# 将偏差提取为数组
binomial_deviance = deviance_results["deviance"].values
idx = binomial_deviance.argsort()[-4000:] # 偏差最大的4000个基因
mask = np.zeros(F1_adata.var_names.shape, dtype=bool)
mask[idx] = True

# %%

F1_adata.var["highly_deviant"] = mask
F1_adata.var["binomial_deviance"] = binomial_deviance

# %%

sc.pp.highly_variable_genes(F1_adata, layer="scran_normalization")
# 计算每个基因在所有细胞中的均值和离散度, 计算完成后, F1_adata.var增加四列

# %%

ax = sns.scatterplot(
    data=F1_adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
)
ax.set_xlim(None, 1.5)
ax.set_ylim(None, 3)
plt.show()

# %%

F1_adata.write("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm-Fs_F1.h5ad")

# %%


