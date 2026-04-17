# %%

import logging
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
F1_adata.var["deviance"] = deviance_results["deviance"].values

# %%

idx = deviance_results.argsort()[-4000:]
