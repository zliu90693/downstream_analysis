# %%
# use envrionment scvi!!! see .env/scvi.yml
import scvi
import scanpy as sc
import time
# %%
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------
# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "species/Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/4_feature-selection-output/concat.h5ad")
# %%
hvg = combined_h5ad[:, combined_h5ad.var['highly_variable']].copy()
scvi.model.SCVI.setup_anndata(hvg, layer='counts', batch_key='batch') #! 指定使用原始计数矩阵, 是在这里指定的!!!

model = scvi.model.SCVI(hvg)
model.view_anndata_setup()

model.train(early_stopping=True, accelerator='cpu')
# %%
combined_h5ad.X = combined_h5ad.layers['logcounts'].copy()
sc.tl.pca(combined_h5ad)
sc.pp.neighbors(combined_h5ad, use_rep='X_pca')
sc.tl.umap(combined_h5ad)
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(combined_h5ad, color=["batch"], wspace=1, return_fig=True)
fig.savefig(f"./{project_name}/figures/5_scVI/umap_without_integrate_batch.pdf", bbox_inches='tight')
# %%
combined_h5ad.obsm['X_scVI'] = model.get_latent_representation()
sc.pp.neighbors(combined_h5ad, use_rep="X_scVI")
sc.tl.umap(combined_h5ad)
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(combined_h5ad, color=["batch"], wspace=1, return_fig=True)
fig.savefig(f"./{project_name}/figures/5_scVI/umap_integrated_batch.pdf", bbox_inches='tight')
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/5_scVI-output/concat.h5ad")
hvg.write_h5ad(f"./{project_name}/data/5_scVI-output/hvg.h5ad")
model.save(f"./{project_name}/model/scVI_model_{time.strftime('%Y%m%d')}", save_anndata=False)

# %%
# ---------------------------------------- Apis cerana ----------------------------------------
project_name = "species/Acer"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/4_feature-selection-output/concat.h5ad")
# %%
hvg = combined_h5ad[:, combined_h5ad.var['highly_variable']].copy()
scvi.model.SCVI.setup_anndata(hvg, layer='counts', batch_key='batch')

model = scvi.model.SCVI(hvg)
model.view_anndata_setup()

model.train(early_stopping=True, accelerator='cpu')
# %%
combined_h5ad.X = combined_h5ad.layers['logcounts'].copy()
sc.tl.pca(combined_h5ad)
sc.pp.neighbors(combined_h5ad, use_rep='X_pca')
sc.tl.umap(combined_h5ad)
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(combined_h5ad, color=["batch"], wspace=1, return_fig=True)
fig.savefig(f"./{project_name}/figures/5_scVI/umap_without_integrate_batch.pdf", bbox_inches='tight')
# %%
combined_h5ad.obsm['X_scVI'] = model.get_latent_representation()
sc.pp.neighbors(combined_h5ad, use_rep="X_scVI")
sc.tl.umap(combined_h5ad)
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(combined_h5ad, color=["batch"], wspace=1, return_fig=True)
fig.savefig(f"./{project_name}/figures/5_scVI/umap_integrated_batch.pdf", bbox_inches='tight')
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/5_scVI-output/concat.h5ad")
hvg.write_h5ad(f"./{project_name}/data/5_scVI-output/hvg.h5ad")
model.save(f"./{project_name}/model/scVI_model_{time.strftime('%Y%m%d')}", save_anndata=False)

# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
# %%
project_name = "species/Zhang_iScience_2022_Amel"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/4_feature-selection-output/concat.h5ad")
# %%
hvg = combined_h5ad[:, combined_h5ad.var['highly_variable']].copy()
scvi.model.SCVI.setup_anndata(hvg, layer='counts', batch_key='batch')

model = scvi.model.SCVI(hvg)
model.view_anndata_setup()

model.train(early_stopping=True, accelerator='cpu')
# %%
combined_h5ad.X = combined_h5ad.layers['logcounts'].copy()
sc.tl.pca(combined_h5ad)
sc.pp.neighbors(combined_h5ad, use_rep='X_pca')
sc.tl.umap(combined_h5ad)
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(combined_h5ad, color=["batch"], wspace=1, return_fig=True)
fig.savefig(f"./{project_name}/figures/5_scVI/umap_without_integrate_batch.pdf", bbox_inches='tight')
# %%
combined_h5ad.obsm['X_scVI'] = model.get_latent_representation()
sc.pp.neighbors(combined_h5ad, use_rep="X_scVI")
sc.tl.umap(combined_h5ad)
sc.set_figure_params(figsize=[6, 6], dpi_save=300)
fig = sc.pl.umap(combined_h5ad, color=["batch"], wspace=1, return_fig=True)
fig.savefig(f"./{project_name}/figures/5_scVI/umap_integrated_batch.pdf", bbox_inches='tight')
# %%
combined_h5ad.write_h5ad(f"./{project_name}/data/5_scVI-output/concat.h5ad")
hvg.write_h5ad(f"./{project_name}/data/5_scVI-output/hvg.h5ad")
model.save(f"./{project_name}/model/scVI_model_{time.strftime('%Y%m%d')}", save_anndata=False)
