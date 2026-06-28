# %%
import scanpy as sc
import numpy as np
import pandas as pd
# %%
Amel = sc.read_h5ad("../Zhang_iScience_2022_Amel/data/7_cluster-output/concat.h5ad")
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["leiden_harmony_scran_res0.25", "leiden_harmony_scran_res0.50", "leiden_harmony_scran_res1.00"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["leiden_harmony_pearson_res0.25", "leiden_harmony_pearson_res0.50", "leiden_harmony_pearson_res1.00"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["leiden_harmony_log1p_res0.25", "leiden_harmony_log1p_res0.50", "leiden_harmony_log1p_res1.00"],legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
# %%
Amel.obs["log1p_Taillight"] = np.nan
Amel.obs.loc[(Amel.obs["leiden_harmony_log1p_res1.00"] =="22") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "24") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "27") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "21")  | (Amel.obs["leiden_harmony_log1p_res1.00"] == "17"), "log1p_Taillight"] = "Taillight"
Amel.obs.loc[(Amel.obs["leiden_harmony_log1p_res1.00"] =="34") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "41") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "28") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "39")  | (Amel.obs["leiden_harmony_log1p_res1.00"] == "18"), "log1p_Taillight"] = "Head"

# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["log1p_Taillight"])
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["log1p_Taillight"])
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["log1p_Taillight"])
# %%
cross = pd.crosstab(
    Amel.obs["batch"],
    Amel.obs["log1p_Taillight"],
    margins=True
)
print(cross)
# %%
for layer in ['analytic_pearson_residuals', 'log1p_norm', 'scran_normalization']:
    print(f"--------------------------------------{layer}--------------------------------------")
    Amel_subset = Amel[Amel.obs["log1p_Taillight"].notna()].copy()
    sc.tl.rank_genes_groups(
        Amel_subset,
        groupby="log1p_Taillight",
        groups=["Taillight"],
        reference="Head",
        method="wilcoxon",
        layer=layer,
        use_raw=False
    )
    sc.pl.rank_genes_groups(Amel_subset, n_genes=20, sharey=False)
# %%
log1p_Taillight = Amel.obs["log1p_Taillight"]
log1p_Taillight
# %%
Amel_scVI = sc.read_h5ad("../../downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/6_cluster-output/concat.h5ad")
Amel_scVI
# %%
Amel_scVI.obs
# %%
Amel_scVI.obs["log1p_Taillight"] = log1p_Taillight
# %%
sc.pl.embedding(Amel_scVI, basis="X_umap", color=["log1p_Taillight"]) #! 可视化Taillight在scVI中的分布
# %%
#! --------------------------------------- Head 和 Taillight 都是什么？ ---------------------------------------
Amel_scVI.X = Amel_scVI.layers["logcounts"].copy()
sc.pl.embedding(Amel_scVI, basis="X_umap", color=["LOC408372", "LOC408804"], cmap="Reds") # KC
# %%
sc.pl.embedding(Amel_scVI, basis="X_umap", color=["LOC410151", "LOC408480"], cmap="Reds") # glia
# %%
sc.pl.embedding(Amel_scVI, basis="X_umap", color=["LOC724282", "LOC724148", "LOC410657", "LOC413466"], cmap="Reds") # OPN
# %%
sc.pl.embedding(Amel_scVI, basis="X_umap", color=["LOC551684", "LOC411597"], cmap="Reds") # Hemo
# %%
sc.pl.embedding(Amel_scVI, basis="X_umap", color=["LOC409937", "LOC724608"], cmap="Reds")
# %%
#! --------------------------------------- 在harmony前，log1p_Taillight的簇是什么样子的？ ---------------------------------------