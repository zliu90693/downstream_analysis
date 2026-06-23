# %%
import scanpy as sc
from scipy.stats import spearmanr
import anndata
import pandas as pd
# %%
def check_cluster_spearman(
    adata: anndata.AnnData,
    cluster_id: str,
    normalize_reso_cluster_key: str,
    normalize_umap_key: str,
    technical_param: str
) -> tuple:
    mask = adata.obs[ normalize_reso_cluster_key] == cluster_id
    umap_coords = adata.obsm[normalize_umap_key][mask]
    param = adata.obs.loc[mask, technical_param]

    rho_umap1, p1 = spearmanr(umap_coords[:, 0], param)
    rho_umap2, p2 = spearmanr(umap_coords[:, 1], param)

    return rho_umap1, p1, rho_umap2, p2

def mk_dataframe(
    species: str,
    adata: anndata.AnnData,
    # cluster_id: str,
    normalize_reso_cluster_key: str,
    normalize_umap_key: str,
    technical_param_list: list,
) -> tuple:
    clusters_list = list(adata.obs[normalize_reso_cluster_key].unique())
    cluster = []
    rho_umap1 = []
    p1 = []
    rho_umap2 = []
    p2 = []
    param = []
    normalize_func = normalize_umap_key.split("_")[-1]
    leiden_reso = normalize_reso_cluster_key[-4:]
    for cluster_id in clusters_list:
        for single_param in technical_param_list:
            spearmans = check_cluster_spearman(
                adata=adata,
                cluster_id=cluster_id,
                normalize_reso_cluster_key=normalize_reso_cluster_key,
                normalize_umap_key=normalize_umap_key,
                technical_param=single_param
            )
            cluster.append(cluster_id)
            param.append(single_param)
            rho_umap1.append(spearmans[0])
            p1.append(spearmans[1])
            rho_umap2.append(spearmans[2])
            p2.append(spearmans[3])
    pearson_corr_df = pd.DataFrame({
        "species": species,
        "normalize_func": normalize_func,
        "leiden_reso": leiden_reso,
        "cluster": cluster,
        "param": param,
        "rho_umap1": rho_umap1,
        "p1": p1,
        "rho_umap2": rho_umap2,
        "p2": p2
    })
    file_name = f"{species}_{normalize_func}_reso{leiden_reso}.csv"
    return pearson_corr_df, file_name

# species, normalize_func, leiden_reso, cluster, param, rho_umap1, p1, rho_umap2, p2
# %%
# %%
Amel_h5ad = sc.read_h5ad(f"../Zhang_iScience_2022_Amel/data/7_cluster-output/concat.h5ad")
# %%
# output = check_cluster_spearman(
#     adata=Amel_h5ad,
#     cluster_id='22',
#     normalize_reso_cluster_key='leiden_harmony_log1p_res0.50',
#     normalize_umap_key='X_umap_harmony_log1p',
#     technical_param='log1p_n_genes_by_counts'
# )
# # %%
# output[1]
# %%
Amel_df = mk_dataframe(
    species="Amel",
    adata=Amel_h5ad,
    normalize_reso_cluster_key="leiden_harmony_log1p_res0.50",
    normalize_umap_key="X_umap_harmony_log1p",
    technical_param_list=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"]
)[0]
Amel_name = mk_dataframe(
    species="Amel",
    adata=Amel_h5ad,
    normalize_reso_cluster_key="leiden_harmony_log1p_res0.50",
    normalize_umap_key="X_umap_harmony_log1p",
    technical_param_list=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"]
)[1]
Amel_df.to_csv(f"./metadata/{Amel_name}")
# %%
