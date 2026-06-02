# %%
import scanpy as sc
import numpy as np
from scipy.stats import spearmanr, pearsonr
import anndata
import scvi
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "0"
# %%
project_name = "Zhang_iScience_2022_Amel"
Amel_h5ad = sc.read_h5ad(f"../{project_name}/data/3_concated-output/concat.h5ad")
# %%
Amel_h5ad.layers['counts'] = Amel_h5ad.X.copy()
sc.pp.normalize_total(Amel_h5ad)
sc.pp.log1p(Amel_h5ad)
Amel_h5ad.layers['logcounts'] = Amel_h5ad.X.copy()

sc.pp.highly_variable_genes(Amel_h5ad, n_top_genes=2000, flavor="cell_ranger", batch_key='batch')

Amel_hvg = Amel_h5ad[:, Amel_h5ad.var["highly_variable"]].copy()
Amel_hvg
# %%
scvi.model.SCVI.setup_anndata(Amel_hvg, layer="counts", batch_key='batch')
# %%
Amel_scvi = scvi.model.SCVI(Amel_hvg)
Amel_scvi.view_anndata_setup()
# %%
Amel_scvi.train(early_stopping=True)
# %%
Amel_hvg.obsm["X_scVI"] = Amel_scvi.get_latent_representation()
# %%
# -------------------------------作校正后的批次效应---------------------------------
sc.pp.neighbors(Amel_hvg, use_rep="X_scVI")
sc.tl.umap(Amel_hvg)
fig = sc.pl.umap(Amel_hvg, color=["batch"], wspace=1, return_fig=True)
fig.savefig(f"../{project_name}/figures/.test/Amel_integration.pdf", bbox_inches='tight')
# %%
# -------------------------------未作校正的批次效应---------------------------------
sc.tl.pca(Amel_h5ad)
sc.pp.neighbors(Amel_h5ad, use_rep='X_pca')
sc.tl.umap(Amel_h5ad)
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(Amel_h5ad, color=["batch"], wspace=1, return_fig=True)
# %%
Amel_hvg.write_h5ad('../Zhang_iScience_2022_Amel/data/.test/Amel_hvg_integration.h5ad', compression='gzip')
Amel_h5ad.write_h5ad('../Zhang_iScience_2022_Amel/data/.test/Amel_all_concat.h5ad', compression='gzip')
# %%
Amel_h5ad.obsm["X_scVI"] = Amel_hvg.obsm["X_scVI"] # 本来是要这么做吗???
sc.pp.neighbors(Amel_h5ad, use_rep='X_scVI')
sc.tl.umap(Amel_h5ad)
sc.tl.leiden(Amel_h5ad, resolution=0.25)
fig = sc.pl.umap(Amel_h5ad, legend_loc="on data", color='leiden', return_fig=True)
# %%
