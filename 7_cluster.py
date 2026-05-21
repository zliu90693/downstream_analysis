# %%
import scanpy as sc
import anndata
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

def run_leiden(
    adata: anndata.AnnData,
    neo_key: str,
    neighbors_key: str,
    resolution: float = 1.0,
    flavor: str = "igraph",
    n_iterations: int = 2,
) -> None:
    sc.tl.leiden(
        adata, 
        key_added=f"leiden_{neo_key}_res{resolution:.2f}",  # 防止浮点数转字符串时可能出现 res0.6000000000000001
        neighbors_key=neighbors_key,
        resolution=resolution, 
        flavor=flavor, 
        n_iterations=n_iterations
    )
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
        legend_loc="on data"
    )
# %%
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/6_dim-reduction-output/concat.h5ad")
# %%
run_umap(combined_h5ad, use_rep="X_pca_harmony_scran", neo_key="harmony_scran", n_pcs=30)
run_umap(combined_h5ad, use_rep="X_pca_harmony_pearson", neo_key="harmony_pearson", n_pcs=30)
run_umap(combined_h5ad, use_rep="X_pca_harmony_log1p", neo_key="harmony_log1p", n_pcs=30)
# %%
for reso in [0.25, 0.5, 1.0]:
    run_leiden(combined_h5ad, neo_key="harmony_scran", neighbors_key="neighbors_harmony_scran", resolution=reso, flavor="igraph", n_iterations=2)
    run_leiden(combined_h5ad, neo_key="harmony_pearson", neighbors_key="neighbors_harmony_pearson", resolution=reso, flavor="igraph", n_iterations=2)
    run_leiden(combined_h5ad, neo_key="harmony_log1p", neighbors_key="neighbors_harmony_log1p", resolution=reso, flavor="igraph", n_iterations=2)
# %%
visualize(combined_h5ad, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.25", "leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00"])
visualize(combined_h5ad, basis="X_umap_harmony_pearson", color=["leiden_harmony_pearson_res0.25", "leiden_harmony_pearson_res0.50", "leiden_harmony_pearson_res1.00"])
visualize(combined_h5ad, basis="X_umap_harmony_log1p", color=["leiden_harmony_log1p_res0.25", "leiden_harmony_log1p_res0.50", "leiden_harmony_log1p_res1.00"])
# %%
combined_h5ad.write(f"./{project_name}/7_cluster-output/concat.h5ad")
# %%
