# %%
import scanpy as sc
import numpy as np
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
# %%
Amel.obs["log1p_Taillight"] = np.nan
Amel.obs.loc[(Amel.obs["leiden_harmony_log1p_res1.00"] =="22") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "24") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "27") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "21")  | (Amel.obs["leiden_harmony_log1p_res1.00"] == "17"), "log1p_Taillight"] = "Taillight"
Amel.obs.loc[(Amel.obs["leiden_harmony_log1p_res1.00"] =="34") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "41") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "28") | (Amel.obs["leiden_harmony_log1p_res1.00"] == "39")  | (Amel.obs["leiden_harmony_log1p_res1.00"] == "18"), "log1p_Taillight"] = "Head"

# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["log1p_Taillight"])
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["log1p_Taillight"])
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["log1p_Taillight"])
# %%
