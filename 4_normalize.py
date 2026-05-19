# %%

import numpy as np
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.sparse import issparse
import anndata

from scipy.sparse import csr_matrix
from scipy.io import mmwrite

import pandas as pd
import subprocess



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

def pre_scran(
    proj_name: str,
    adata: anndata.AnnData
) -> None:
    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(
        adata_pp, key_added="groups", flavor="igraph", n_iterations=2, directed=False
    )

    data_mat = adata_pp.X.T
    if issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
    
    mmwrite(f"./{proj_name}/metadata/4_scran-data-mat.mtx", data_mat)
    cells = pd.DataFrame({"cell_id": adata_pp.obs_names.tolist()})
    cells.to_csv(f"./{proj_name}/metadata/4_cellnames.csv", index=False)
    adata_pp.obs[["groups"]].to_csv(f"./{proj_name}/metadata/4_scran-groups.csv", index_label="cell_id")

def do_scran(
    proj_name: str
) -> None:
    subprocess.run(["./4._run_scran.R", proj_name], check=True)

def combine_scran_output(
    proj_name: str,
    adata: anndata.AnnData
) -> None:
    size_factors_df = pd.read_csv(f"./{proj_name}/metadata/4_scran-size-factors.csv")
    size_factors = size_factors_df.set_index("cell_id")["size_factors"]

    adata.obs["size_factors"] = size_factors.loc[adata.obs_names]

    assert adata.obs["size_factors"].isna().sum() == 0, \
        "NaN detected! cell_id is not fully aligned; please check the indices."
    assert (adata.obs["size_factors"] <= 0).sum() == 0, \
        "Non-positive size factors detected; scran calculation may be erroneous."

    scran = adata.X / adata.obs["size_factors"].values[:, None]
    scran_logged = np.log1p(scran)
    adata.layers["scran_normalization"] = csr_matrix(scran_logged)

def visualize_scran(
    adata: anndata.AnnData
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(
        adata.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1]
    )
    axes[1].set_title("log1p with Scran estimated size factors")
    plt.show()

# ---------------------------------------- Pearson correlation ----------------------------------------


    
# %%

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/3_concated-output/concat.h5ad")

# %%

do_shift_log(combined_h5ad)
visualize_shift_log(combined_h5ad)

# %%

pre_scran(project_name, combined_h5ad)

# %%

do_scran(project_name) # Use rm_ambient_doublet environment!!!

# %%

combine_scran_output(project_name, combined_h5ad)

# %%

visualize_scran(combined_h5ad)

# %%
