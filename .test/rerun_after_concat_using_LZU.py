# %%
import scanpy as sc
import numpy as np
from scipy.stats import spearmanr, pearsonr
import anndata
import scvi
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
scvi.model.SCVI.setup_anndata(Amel_hvg, layer="counts", batch_key='batch', TF_CPP_MIN_LOG_LEVEL=0)