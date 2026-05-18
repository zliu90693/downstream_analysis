# %%

import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.sparse import issparse
import anndata

from scipy.sparse import csr_matrix
from scipy.io import mmwrite



# %%

# combined_h5ad = sc.read_h5ad("./Sheng_SA_2020_Hsal/3_concated-output/concat.h5ad")

# %%

# ----------------------------------- Shifted logarithm -----------------------------------

def do_shift_log(
    adata: anndata.AnnData
) -> None:
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

def visualize_shift_log(
    adata: anndata.AnnData
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
    axes[1].set_title("Shifted logarithm")
    plt.show()

# %%

# ----------------------------------- Scran -----------------------------------

def 

# %%

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/3_concated-output/concat.h5ad")

# %%

do_shift_log(combined_h5ad)
visualize_shift_log(combined_h5ad)

# %%
