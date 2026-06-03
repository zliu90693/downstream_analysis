# %%
import scanpy as sc
# import anndata
# %%
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "species/Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/3_concated-output/concat.h5ad")
# %%
combined_h5ad.X = combined_h5ad.layers['counts'].copy()
sc.pp.normalize_total(combined_h5ad)
sc.pp.log1p(combined_h5ad)
combined_h5ad.layers['logcounts'] = combined_h5ad.X.copy()
# %%
if combined_h5ad.X.max() == combined_h5ad.layers["logcounts"].max(): # flavor='cell_ranger' 需要 log1p 之后的表达矩阵, 如果使用seurat_v3 则需要原始表达矩阵
    sc.pp.highly_variable_genes(combined_h5ad, n_top_genes=4000, flavor='cell_ranger', batch_key='batch')
else:
    print("Warning: The maximum value of the log-transformed data is different from the original data, which may lead to incorrect results when using flavor='cell_ranger'. Please check the data and consider using flavor='seurat_v3' instead.")
# %%
# hvg = combined_h5ad[:, combined_h5ad.var['highly_variable']].copy()
# hvg
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/4_feature-selection-output/concat.h5ad")
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "species/Acer"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/3_concated-output/concat.h5ad")
# %%
combined_h5ad.X = combined_h5ad.layers['counts'].copy()
sc.pp.normalize_total(combined_h5ad)
sc.pp.log1p(combined_h5ad)
combined_h5ad.layers['logcounts'] = combined_h5ad.X.copy()
# %%
if combined_h5ad.X.max() == combined_h5ad.layers["logcounts"].max(): # flavor='cell_ranger' 需要 log1p 之后的表达矩阵, 如果使用seurat_v3 则需要原始表达矩阵
    sc.pp.highly_variable_genes(combined_h5ad, n_top_genes=4000, flavor='cell_ranger', batch_key='batch')
else:
    print("Warning: The maximum value of the log-transformed data is different from the original data, which may lead to incorrect results when using flavor='cell_ranger'. Please check the data and consider using flavor='seurat_v3' instead.")
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/4_feature-selection-output/concat.h5ad")

# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
project_name = "species/Zhang_iScience_2022_Amel"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/3_concated-output/concat.h5ad")
# %%
combined_h5ad.X = combined_h5ad.layers['counts'].copy()
sc.pp.normalize_total(combined_h5ad)
sc.pp.log1p(combined_h5ad)
combined_h5ad.layers['logcounts'] = combined_h5ad.X.copy()
# %%
if combined_h5ad.X.max() == combined_h5ad.layers["logcounts"].max(): # flavor='cell_ranger' 需要 log1p 之后的表达矩阵, 如果使用seurat_v3 则需要原始表达矩阵
    sc.pp.highly_variable_genes(combined_h5ad, n_top_genes=4000, flavor='cell_ranger', batch_key='batch')
else:
    print("Warning: The maximum value of the log-transformed data is different from the original data, which may lead to incorrect results when using flavor='cell_ranger'. Please check the data and consider using flavor='seurat_v3' instead.")
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/4_feature-selection-output/concat.h5ad")