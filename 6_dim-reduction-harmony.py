# %%

import scanpy as sc
import anndata
import harmonypy as hm

# %%

def run_pca(
    adata: anndata.AnnData,
    layer: str = "scran_normalization",
) -> None:
    adata.X = adata.layers[layer]
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", mask_var="highly_variable")

# %%

def run_tsne(
    adata: anndata.AnnData,
    use_rep: str = "X_pca",
) -> None:
    sc.tl.tsne(adata, use_rep=use_rep)


# %%

def run_umap(
    adata: anndata.AnnData,
) -> None:
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)


# %%
def run_harmony(
    adata: anndata.AnnData,
    batch_key: str = "batch",
    use_rep: str = "X_pca",
    adjusted_basis: str = 'X_pca_harmony'
) -> None:
    # sc.external.pp.harmony_integrate(adata, key=batch_key, basis=use_rep, adjusted_basis=adjusted_basis)

    pcs = adata.obsm[use_rep] # (n_cells, n_pcs)
    # Run Harmony
    harmony_out = hm.run_harmony(pcs, adata.obs, batch_key)
    # Z_corr is already a numpy array with shape (n_cells, n_pcs)
    adata.obsm[adjusted_basis] = harmony_out.Z_corr
    # Use harmonized PCs for downstream analysis
    # sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    # sc.tl.umap(adata)

# %%

def visualize(
    adata: anndata.AnnData,
    basis: str,
    color,
) -> None:
    sc.pl.embedding(
        adata,
        basis=basis,
        color=color,
    )

# %%

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/5_feature-selection-output/concat.h5ad")

# %%

run_pca(combined_h5ad, layer="scran_normalization")
# %%

visualize(combined_h5ad, basis="X_pca", color=["total_counts", "pct_counts_mt", "batch"])
sc.pl.pca_loadings(combined_h5ad, components=[1, 2])

# %%

run_tsne(combined_h5ad, use_rep="X_pca")

# %%
visualize(combined_h5ad, basis="X_tsne", color=["total_counts", "pct_counts_mt", "batch"])

# %%
run_umap(combined_h5ad)
# %%
visualize(combined_h5ad, basis="X_umap", color=["total_counts", "pct_counts_mt", "batch"])

# %%
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca", adjusted_basis='X_pca_harmony')
# %%
