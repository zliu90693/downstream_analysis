# %%
import scanpy as sc
from scipy.stats import spearmanr
import anndata
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
import numpy as np
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
technical_param_list = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"]
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
# %%
Amel_log1p_050 = mk_dataframe(
    species="Amel",
    adata=Amel_h5ad,
    normalize_reso_cluster_key="leiden_harmony_log1p_res0.50",
    normalize_umap_key="X_umap_harmony_log1p",
    technical_param_list=technical_param_list
)
Amel_df = Amel_log1p_050[0]
Amel_name = Amel_log1p_050[1]
Amel_df.to_csv(f"./metadata/{Amel_name}")
# %%
Amel_scran_050 = mk_dataframe(
    species="Amel",
    adata=Amel_h5ad,
    normalize_reso_cluster_key="leiden_harmony_scran_res0.50",
    normalize_umap_key="X_umap_harmony_scran",
    technical_param_list=technical_param_list
)
Amel_df = Amel_scran_050[0]
Amel_name = Amel_scran_050[1]
Amel_df.to_csv(f"./metadata/{Amel_name}")
# %%
Amel_pearson_050 = mk_dataframe(
    species="Amel",
    adata=Amel_h5ad,
    normalize_reso_cluster_key="leiden_harmony_pearson_res0.50",
    normalize_umap_key="X_umap_harmony_pearson",
    technical_param_list=technical_param_list
)
Amel_df = Amel_pearson_050[0]
Amel_name = Amel_pearson_050[1]
Amel_df.to_csv(f"./metadata/{Amel_name}")
# %%
Hsal_h5ad = sc.read_h5ad(f"../Sheng_SA_2020_Hsal/data/7_cluster-output/concat.h5ad")
# %%
Hsal_log1p_050 = mk_dataframe(
    species="Hsal",
    adata=Hsal_h5ad,
    normalize_reso_cluster_key="leiden_harmony_log1p_res0.50",
    normalize_umap_key="X_umap_harmony_log1p",
    technical_param_list=technical_param_list
)
Hsal_df = Hsal_log1p_050[0]
Hsal_name = Hsal_log1p_050[1]
Hsal_df.to_csv(f"./metadata/{Hsal_name}")
# %%
Hsal_scran_050 = mk_dataframe(
    species="Hsal",
    adata=Hsal_h5ad,
    normalize_reso_cluster_key="leiden_harmony_scran_res0.50",
    normalize_umap_key="X_umap_harmony_scran",
    technical_param_list=technical_param_list
)
Hsal_df = Hsal_scran_050[0]
Hsal_name = Hsal_scran_050[1]
Hsal_df.to_csv(f"./metadata/{Hsal_name}")
# %%
Hsal_pearson_050 = mk_dataframe(
    species="Hsal",
    adata=Hsal_h5ad,
    normalize_reso_cluster_key="leiden_harmony_pearson_res0.50",
    normalize_umap_key="X_umap_harmony_pearson",
    technical_param_list=technical_param_list
)
Hsal_df = Hsal_pearson_050[0]
Hsal_name = Hsal_pearson_050[1]
Hsal_df.to_csv(f"./metadata/{Hsal_name}")
# %%
Acer_h5ad = sc.read_h5ad(f"../Acer/data/7_cluster-output/concat.h5ad")
# %%
Acer_log1p_050 = mk_dataframe(
    species="Acer",
    adata=Acer_h5ad,
    normalize_reso_cluster_key="leiden_harmony_log1p_res0.50",
    normalize_umap_key="X_umap_harmony_log1p",
    technical_param_list=technical_param_list
)
Acer_df = Acer_log1p_050[0]
Acer_name = Acer_log1p_050[1]
Acer_df.to_csv(f"./metadata/{Acer_name}")
# %%
Acer_scran_050 = mk_dataframe(
    species="Acer",
    adata=Acer_h5ad,
    normalize_reso_cluster_key="leiden_harmony_scran_res0.50",
    normalize_umap_key="X_umap_harmony_scran",
    technical_param_list=technical_param_list
)
Acer_df = Acer_scran_050[0]
Acer_name = Acer_scran_050[1]
Acer_df.to_csv(f"./metadata/{Acer_name}")
# %%
Acer_pearson_050 = mk_dataframe(
    species="Acer",
    adata=Acer_h5ad,
    normalize_reso_cluster_key="leiden_harmony_pearson_res0.50",
    normalize_umap_key="X_umap_harmony_pearson",
    technical_param_list=technical_param_list
)
Acer_df = Acer_pearson_050[0]
Acer_name = Acer_pearson_050[1]
Acer_df.to_csv(f"./metadata/{Acer_name}")

#! --------------------------------------------------------------------------------------------------
#! --------------------------------------- Analyse output csv ---------------------------------------
#! --------------------------------------------------------------------------------------------------
# %%
# def rho_heatmap(
#     df: pd.DataFrame,
#     filename: str,
# ) -> None:
#     pivot_umap1 = df.pivot_table(
#         index='cluster', 
#         columns='param', 
#         values='rho_umap1'
#     )
#     fig, ax = plt.subplots(figsize=(8, max(6, len(pivot_umap1) * 0.5)))
#     sns.heatmap(
#         pivot_umap1,
#         cmap="RdBu_r",        # 红蓝双色：红=正相关，蓝=负相关
#         center=0,             # 以0为中心对称着色
#         vmin=-0.8, vmax=0.8,  # 固定色标范围，确保跨图可比
#         annot=True,           # 在格子内显示数值
#         fmt=".2f",            # 保留两位小数
#         linewidths=0.5,       # 格子边框线
#         linecolor='white',
#         cbar_kws={'label': 'Spearman ρ (UMAP1)', 'shrink': 0.8},
#         ax=ax
#     )
#     ax.set_title("UMAP1 vs Technical Covariates", fontsize=14)
#     ax.set_ylabel("Leiden Cluster")
#     ax.set_xlabel("Technical Covariate")
#     plt.tight_layout()
#     plt.savefig(f"./figures/spearman_heatmap_{filename}_UMAP1_rho.png", dpi=150, bbox_inches='tight')
#     plt.show()

#     pivot_umap2 = df.pivot_table(
#         index='cluster', 
#         columns='param', 
#         values='rho_umap2'
#     )
#     fig, ax = plt.subplots(figsize=(8, max(6, len(pivot_umap2) * 0.5)))
#     sns.heatmap(
#         pivot_umap2,
#         cmap="RdBu_r",        
#         center=0,             
#         vmin=-0.8, vmax=0.8,  
#         annot=True,           
#         fmt=".2f",            
#         linewidths=0.5,       
#         linecolor='white',
#         cbar_kws={'label': 'Spearman ρ (UMAP2)', 'shrink': 0.8},
#         ax=ax
#     )

#     ax.set_title("UMAP2 vs Technical Covariates", fontsize=14)
#     ax.set_ylabel("Leiden Cluster")
#     ax.set_xlabel("Technical Covariate")
#     plt.tight_layout()
#     plt.savefig(f"./figures/spearman_heatmap_{filename}_UMAP2_rho.png", dpi=150, bbox_inches='tight')
#     plt.show()
# # %%
# Amel_log1p_50_df = pd.read_csv("./metadata/Amel_log1p_reso0.50.csv")
# rho_heatmap(Amel_log1p_50_df, "Amel_log1p_50")
# # %%
# Amel_scran_50_df = pd.read_csv("./metadata/Amel_scran_reso0.50.csv")
# rho_heatmap(Amel_scran_50_df, "Amel_scran_50")
# # %%
# Amel_pearson_50_df = pd.read_csv("./metadata/Amel_pearson_reso0.50.csv")
# rho_heatmap(Amel_pearson_50_df, "Amel_pearson_50")
# # %%
# Acer_log1p_50_df = pd.read_csv("./metadata/Acer_log1p_reso0.50.csv")
# rho_heatmap(Acer_log1p_50_df, "Acer_log1p_50")
# # %%
# Acer_scran_50_df = pd.read_csv("./metadata/Acer_scran_reso0.50.csv")
# rho_heatmap(Acer_scran_50_df, "Acer_scran_50")
# # %%
# Acer_pearson_50_df = pd.read_csv("./metadata/Acer_pearson_reso0.50.csv")
# rho_heatmap(Acer_pearson_50_df, "Acer_pearson_50")
# # %%
# Hsal_log1p_50_df = pd.read_csv("./metadata/Hsal_log1p_reso0.50.csv")
# rho_heatmap(Hsal_log1p_50_df, "Hsal_log1p_50")
# # %%
# Hsal_scran_50_df = pd.read_csv("./metadata/Hsal_scran_reso0.50.csv")
# rho_heatmap(Hsal_scran_50_df, "Hsal_scran_50")
# # %%
# Hsal_pearson_50_df = pd.read_csv("./metadata/Hsal_pearson_reso0.50.csv")
# rho_heatmap(Hsal_pearson_50_df, "Hsal_pearson_50")
# # %%
# Amel_log1p_50_df
# %%
#! 例： 使用Amel_log1p_50_df重写 rho_heatmap 函数整套流程
# %%
def do_fdr_bh(
    df: pd.DataFrame,
) -> pd.DataFrame:
    _, p1_adj, _, _ = multipletests(df["p1"], method="fdr_bh")
    _, p2_adj, _, _ = multipletests(df["p2"], method="fdr_bh")
    df["p1_fdr"] = p1_adj
    df["p2_fdr"] = p2_adj
    return df

def make_mask(
    df: pd.DataFrame,
    threshold: float = 0.05,
) -> tuple:
    df = df.copy()
    rho1 = df.pivot(index="cluster", columns="param", values="rho_umap1")
    rho2 = df.pivot(index="cluster", columns="param", values="rho_umap2")
    p1_mat = df.pivot(index="cluster", columns="param", values="p1_fdr")
    p2_mat = df.pivot(index="cluster", columns="param", values="p2_fdr")
    rho1_masked = rho1.where(p1_mat < threshold, other=np.nan)
    rho2_masked = rho2.where(p2_mat < threshold, other=np.nan)
    return rho1_masked, rho2_masked

def visualize_heat(
    rho_masked: tuple,
    filename: str
) -> None:
    rho1_masked = rho_masked[0]
    rho2_masked = rho_masked[1]
    vmax = np.nanmax([
        np.nanmax(np.abs(rho1_masked.values)),
        np.nanmax(np.abs(rho2_masked.values))
    ])
    fig, axes = plt.subplots(1, 2, figsize=(14, max(4, len(rho1_masked) * 0.5 + 2)))
    common_kwargs = dict(
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="grey",
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Spearman ρ"},
    )
    for ax, rho_mat, title in zip(
        axes,
        [rho1_masked, rho2_masked],
        ["UMAP1", "UMAP2"]
    ):
        sns.heatmap(rho_mat, ax=ax, **common_kwargs)
        ax.set_xlabel("Technical factor")
        ax.set_ylabel("Leiden cluster")
    plt.tight_layout()
    plt.savefig(f"./figures/{filename}_UMAP2_rho.png", dpi=150, bbox_inches='tight')
    plt.show()
# %%
Amel_log1p_50_df = pd.read_csv("./metadata/Amel_log1p_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Amel_log1p_50_df)), "Amel_log1p_50")
# %%
Amel_scran_50_df = pd.read_csv("./metadata/Amel_scran_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Amel_scran_50_df)), "Amel_scran_50")
# %%
Amel_pearson_50_df = pd.read_csv("./metadata/Amel_pearson_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Amel_pearson_50_df)), "Amel_pearson_50")
# %%
Acer_log1p_50_df = pd.read_csv("./metadata/Acer_log1p_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Acer_log1p_50_df)), "Acer_log1p_50_df")
# %%
Acer_scran_50_df = pd.read_csv("./metadata/Acer_scran_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Acer_scran_50_df)), "Acer_scran_50_df")
# %%
Acer_pearson_50_df = pd.read_csv("./metadata/Acer_pearson_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Acer_pearson_50_df)), "Acer_pearson_50_df")
# %%
Hsal_log1p_50_df = pd.read_csv("./metadata/Hsal_log1p_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Hsal_log1p_50_df)), "Hsal_log1p_50_df")
# %%
Hsal_scran_50_df = pd.read_csv("./metadata/Hsal_scran_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Hsal_scran_50_df)), "Hsal_scran_50_df")
# %%
Hsal_pearson_50_df = pd.read_csv("./metadata/Hsal_pearson_reso0.50.csv")
visualize_heat(make_mask(do_fdr_bh(Hsal_pearson_50_df)), "Hsal_pearson_50_df")