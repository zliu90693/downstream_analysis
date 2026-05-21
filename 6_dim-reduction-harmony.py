# %%

import scanpy as sc
import anndata
import harmonypy as hm

# %%

def run_pca(
    adata: anndata.AnnData,
    neo_key: str,
    layer: str = "scran_normalization"
) -> None:
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(
        adata, 
        svd_solver="arpack", 
        mask_var="highly_variable", 
        layer=layer)
    adata.obsm[f"X_pca_{neo_key}"] = adata.obsm.pop("X_pca")
    adata.varm[f"PCs_{neo_key}"] = adata.varm.pop("PCs")
    adata.uns[f"pca_{neo_key}"] = adata.uns.pop("pca")


# %%

def run_tsne(
    adata: anndata.AnnData,
    neo_key: str,
    use_rep: str = "X_pca",
) -> None:
    sc.tl.tsne(adata, use_rep=use_rep)
    adata.obsm[f"X_tsne_{neo_key}"] = adata.obsm.pop("X_tsne")
    adata.uns[f"tsne_{neo_key}"] = adata.uns.pop("tsne")


# %%

def run_umap(
    adata: anndata.AnnData,
    neo_key: str,
    use_rep: str = "X_pca",
    n_pcs: int = None,
) -> None:
    sc.pp.neighbors(adata, use_rep=use_rep, key_added=f"neighbors_{neo_key}", n_pcs=n_pcs)
    sc.tl.umap(adata, neighbors_key=f"neighbors_{neo_key}")
    adata.obsm[f"X_umap_{neo_key}"] = adata.obsm.pop("X_umap")
    adata.uns[f"umap_{neo_key}"] = adata.uns.pop("umap")

# %%

# uns: 'neighbors*', 'umap'
# obsm: 'X_umap'
# obsp: 'distances*', 'connectivities*'

"""
AnnData object with n_obs × n_vars = 18932 × 9814
    obs: 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_20_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 'outlier', 'mt_outlier', 'batch', 'size_factors'
    var: 'gene_ids', 'feature_types', 'genome', 'mt', 'n_cells', 'binomial_deviance', 'highly_deviant', 'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'hvg', 'pca_scran', 'pca_pearson', 'pca_log1p', 'batch_colors', 'tsne_scran', 'tsne_pearson', 'tsne_log1p', 'neighbors', 'umap'
    obsm: 'X_pca_scran', 'X_pca_pearson', 'X_pca_log1p', 'X_tsne_scran', 'X_tsne_pearson', 'X_tsne_log1p', 'X_umap'
    varm: 'PCs_scran', 'PCs_pearson', 'PCs_log1p'
    layers: 'analytic_pearson_residuals', 'counts', 'log1p_norm', 'scran_normalization'
    obsp: 'distances', 'connectivities'
"""

"""
sc.pp.neighbors
key_added: 
If not specified, the neighbors data is stored in .uns['neighbors'], distances and connectivities are stored in .obsp['distances'] and .obsp['connectivities'] 
respectively. If specified, the neighbors data is added to .uns[key_added], distances are stored in .obsp[key_added+'_distances'] and connectivities in .obsp[key_added+'_connectivities'].
"""

"""
sc.tl.umap
neighbors_key
Umap looks in ~anndata.AnnData.uns\ [neighbors_key] for neighbors settings and ~anndata.AnnData.obsp\ [.uns[neighbors_key]['connectivities_key']] for connectivities.
"""


# %%
def run_harmony( # From https://github.com/slowkow/harmonypy/issues/49
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

# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# %%

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/5_feature-selection-output/concat.h5ad")
# %%
run_pca(combined_h5ad, layer="scran_normalization", neo_key="scran")
run_pca(combined_h5ad, layer="analytic_pearson_residuals", neo_key="pearson")
run_pca(combined_h5ad, layer="log1p_norm", neo_key="log1p")
# %%
visualize(combined_h5ad, basis="X_pca_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# sc.pl.pca_loadings(combined_h5ad, components=[1, 2])
# %%
run_tsne(combined_h5ad, use_rep="X_pca_scran", neo_key="scran")
run_tsne(combined_h5ad, use_rep="X_pca_pearson", neo_key="pearson")
run_tsne(combined_h5ad, use_rep="X_pca_log1p", neo_key="log1p") 
# %%
visualize(combined_h5ad, basis="X_tsne_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_umap(combined_h5ad, use_rep="X_pca_scran", neo_key="scran")
run_umap(combined_h5ad, use_rep="X_pca_pearson", neo_key="pearson")
run_umap(combined_h5ad, use_rep="X_pca_log1p", neo_key="log1p")
# %%
visualize(combined_h5ad, basis="X_umap_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca_scran", adjusted_basis='X_pca_harmony_scran')
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca_pearson", adjusted_basis='X_pca_harmony_pearson')
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca_log1p", adjusted_basis='X_pca_harmony_log1p')
# %%
visualize(combined_h5ad, basis="X_pca_harmony_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_harmony_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_harmony_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_tsne(combined_h5ad, use_rep="X_pca_harmony_scran", neo_key="harmony_scran")
run_tsne(combined_h5ad, use_rep="X_pca_harmony_pearson", neo_key="harmony_pearson")
run_tsne(combined_h5ad, use_rep="X_pca_harmony_log1p", neo_key="harmony_log1p") 
# %%
visualize(combined_h5ad, basis="X_tsne_harmony_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_harmony_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_harmony_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_umap(combined_h5ad, use_rep="X_pca_harmony_scran", neo_key="harmony_scran")
run_umap(combined_h5ad, use_rep="X_pca_harmony_pearson", neo_key="harmony_pearson")
run_umap(combined_h5ad, use_rep="X_pca_harmony_log1p", neo_key="harmony_log1p")
# %%
visualize(combined_h5ad, basis="X_umap_harmony_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_harmony_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_harmony_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
combined_h5ad.write(f"./{project_name}/6_dim-reduction-output/concat.h5ad")
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
# %%

project_name = "Acer"
combined_h5ad = sc.read_h5ad(f"./{project_name}/5_feature-selection-output/concat.h5ad")
# %%
run_pca(combined_h5ad, layer="scran_normalization", neo_key="scran")
run_pca(combined_h5ad, layer="analytic_pearson_residuals", neo_key="pearson")
run_pca(combined_h5ad, layer="log1p_norm", neo_key="log1p")
# %%
visualize(combined_h5ad, basis="X_pca_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# sc.pl.pca_loadings(combined_h5ad, components=[1, 2])
# %%
run_tsne(combined_h5ad, use_rep="X_pca_scran", neo_key="scran")
run_tsne(combined_h5ad, use_rep="X_pca_pearson", neo_key="pearson")
run_tsne(combined_h5ad, use_rep="X_pca_log1p", neo_key="log1p") 
# %%
visualize(combined_h5ad, basis="X_tsne_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_umap(combined_h5ad, use_rep="X_pca_scran", neo_key="scran")
run_umap(combined_h5ad, use_rep="X_pca_pearson", neo_key="pearson")
run_umap(combined_h5ad, use_rep="X_pca_log1p", neo_key="log1p")
# %%
visualize(combined_h5ad, basis="X_umap_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca_scran", adjusted_basis='X_pca_harmony_scran')
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca_pearson", adjusted_basis='X_pca_harmony_pearson')
run_harmony(combined_h5ad, batch_key="batch", use_rep="X_pca_log1p", adjusted_basis='X_pca_harmony_log1p')
# %%
visualize(combined_h5ad, basis="X_pca_harmony_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_harmony_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_pca_harmony_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_tsne(combined_h5ad, use_rep="X_pca_harmony_scran", neo_key="harmony_scran")
run_tsne(combined_h5ad, use_rep="X_pca_harmony_pearson", neo_key="harmony_pearson")
run_tsne(combined_h5ad, use_rep="X_pca_harmony_log1p", neo_key="harmony_log1p") 
# %%
visualize(combined_h5ad, basis="X_tsne_harmony_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_harmony_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_tsne_harmony_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
run_umap(combined_h5ad, use_rep="X_pca_harmony_scran", neo_key="harmony_scran")
run_umap(combined_h5ad, use_rep="X_pca_harmony_pearson", neo_key="harmony_pearson")
run_umap(combined_h5ad, use_rep="X_pca_harmony_log1p", neo_key="harmony_log1p")
# %%
visualize(combined_h5ad, basis="X_umap_harmony_scran", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_harmony_pearson", color=["total_counts", "pct_counts_mt", "batch"])
visualize(combined_h5ad, basis="X_umap_harmony_log1p", color=["total_counts", "pct_counts_mt", "batch"])
# %%
combined_h5ad.write(f"./{project_name}/6_dim-reduction-output/concat.h5ad")