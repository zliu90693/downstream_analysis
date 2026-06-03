# %%

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.sparse import issparse
from scipy.io import mmwrite
import anndata
import subprocess
import pandas as pd

# %%

def feature_selection_export(
    project_name: str,
    adata: anndata.AnnData
) -> None:
    X_sparse = adata.X.T.tocoo()
    mmwrite(f"./{project_name}/metadata/5_fs-data-mat.mtx", X_sparse)
    adata.obs.to_csv(
        f"./{project_name}/metadata/5_fs-cellnames.csv", 
        index_label="cell_id"
    )
    adata.var.to_csv(
        f"./{project_name}/metadata/5_fs-genenames.csv", 
        index_label="gene_id"
    )

def do_feature_selection(
    project_name: str,
) -> None:
    subprocess.run(["./5._run_devianceFeatureSelection.R", project_name], check=True)

def combine_feature_selection_output(
    project_name: str,
    adata: anndata.AnnData
) -> None:
    deviance_results = pd.read_csv(f"{project_name}/metadata/5_fs_deviance_results.csv", index_col=0)
    aligned_df = deviance_results.reindex(adata.var_names)
    binomial_deviance = aligned_df["deviance"].values
    adata.var["binomial_deviance"] = binomial_deviance

    if np.any(np.isnan(binomial_deviance)):
        missing_num = np.isnan(binomial_deviance).sum()
        raise ValueError(f"Found {missing_num} genes with missing deviance values!")

    idx = binomial_deviance.argsort()[-4000:] # 偏差最大的4000个基因
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    adata.var["highly_deviant"] = mask

def do_hvg(
    adata: anndata.AnnData
) -> None:
    sc.pp.highly_variable_genes(adata, layer="scran_normalization") # 计算每个基因在所有细胞中的均值和离散度, 计算完成后, adata.var增加四列

def visualize_feature_selection(
    adata: anndata.AnnData
) -> None:
    ax = sns.scatterplot(
        data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
    ) # 根据 "均值 x 离散度" 进行可视化, 同时标注高偏差基因
    ax.set_xlim(None, 1.5)
    ax.set_ylim(None, 3)
    plt.show()

# %%
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# ---------------------------------------- Harpegnathos venator ----------------------------------------

# %%
project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/4_normalize-output/concat.h5ad")
# %%
feature_selection_export(project_name, combined_h5ad)
# %%
do_feature_selection(project_name)
# %%
combine_feature_selection_output(project_name, combined_h5ad)
do_hvg(combined_h5ad)
visualize_feature_selection(combined_h5ad)
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/5_feature-selection-output/concat.h5ad")
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "Acer"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/4_normalize-output/concat.h5ad")
# %%
feature_selection_export(project_name, combined_h5ad)
# %%
do_feature_selection(project_name)
# %%
combine_feature_selection_output(project_name, combined_h5ad)
do_hvg(combined_h5ad)
visualize_feature_selection(combined_h5ad)
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/5_feature-selection-output/concat.h5ad")
# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
# %%
project_name = "Zhang_iScience_2022_Amel"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/4_normalize-output/concat.h5ad")
# %%
feature_selection_export(project_name, combined_h5ad)
# %%
do_feature_selection(project_name)
# %%
combine_feature_selection_output(project_name, combined_h5ad)
do_hvg(combined_h5ad)
visualize_feature_selection(combined_h5ad)
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/5_feature-selection-output/concat.h5ad")
# %%
