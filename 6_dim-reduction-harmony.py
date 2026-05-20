# %%

import scanpy as sc
import anndata

# %%

def run_pca(
    adata: anndata.AnnData,
    layer: str = "scran_normalization",
) -> None:
    adata.X = adata.layers[layer]
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable")

def visualize_pca(
    adata: anndata.AnnData,
    color
) -> None:
    sc.pl.pca_scatter(adata, color=color)

# %%

def run_tsne(
    adata: anndata.AnnData,
    use_rep: str = "X_pca",
) -> None:
    sc.tl.tsne(adata, use_rep=use_rep)

def visualize_tsne(
    adata: anndata.AnnData,
    color
) -> None:
    sc.pl.tsne(adata, color=color)

# %%

def run_umap(
    adata: anndata.AnnData,
) -> None:
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

def visualize_umap(
    adata: anndata.AnnData,
    color,
) -> None:
    sc.pl.umap(adata, color=color)

# %%
def run_harmony(
    adata: anndata.AnnData,
    batch_key: str = "batch",
    use_rep: str = "X_pca",
    adjusted_basis: str = 'X_pca_harmony'
) -> None:
    sc.external.pp.harmony_integrate(adata, key=batch_key, basis=use_rep, adjusted_basis=adjusted_basis)


# %%

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/5_feature-selection-output/concat.h5ad")

# %%

run_pca(combined_h5ad, layer="scran_normalization")
# %%

visualize_pca(combined_h5ad, color="total_counts") # ????
visualize_pca(combined_h5ad, color="batch")
sc.pl.pca_loadings(combined_h5ad, components=[1, 2])
visualize_pca(combined_h5ad, color="pct_counts_mt")

# %%

run_tsne(combined_h5ad, use_rep="X_pca")

# %%
visualize_tsne(combined_h5ad, color="total_counts")
visualize_tsne(combined_h5ad, color="batch")
visualize_tsne(combined_h5ad, color="pct_counts_mt")

# %%
run_umap(combined_h5ad)
# %%
visualize_umap(combined_h5ad, color=["total_counts", "pct_counts_mt", "batch"])

# %%
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca", adjusted_basis='X_pca_harmony')
# %%
