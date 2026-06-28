# %%
import scanpy as sc
import numpy as np
import pandas as pd
# %%
Amel = sc.read_h5ad("../../Zhang_iScience_2022_Amel/data/7_cluster-output/concat.h5ad")
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
# sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["size_factors", "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
Amel.obs["dumbbells"] = np.nan
for dumbbell_cluster in [18,14,25,16,15,22,21,19,13,26,24]:
    Amel.obs.loc[(Amel.obs["leiden_harmony_log1p_res0.50"] == str(dumbbell_cluster)), "dumbbells"] = str(dumbbell_cluster)
# 18,14,25,16,15,22,21,19,13,26,24,20
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["dumbbells"], legend_loc="on data")
# %%
Amel.write_h5ad("./data/Amel_dumbbell_marked.h5ad")
# %%
dumbbells = Amel.obs["dumbbells"]
dumbbells.to_csv("./metadata/log1p_dumbbells.csv")
# %%
# -----------------------------------------------------------------------------------------
Amel = sc.read_h5ad("./data/Amel_dumbbell_marked.h5ad")
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["dumbbells"], legend_loc="on data")
# %%
