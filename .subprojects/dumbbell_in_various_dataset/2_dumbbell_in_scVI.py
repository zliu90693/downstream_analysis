# %%
import scanpy as sc
import numpy as np
import pandas as pd
# %%
Amel_scVI = sc.read_h5ad("../../../downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/6_cluster-output/concat.h5ad")
Amel_scVI
# %%
dumbbells = pd.read_csv("./metadata/log1p_dumbbells.csv", index_col=0)["dumbbells"]
# %%
dumbbells = pd.to_numeric(dumbbells, errors='coerce')
dumbbells = dumbbells.astype('Int64')
dumbbells = dumbbells.reindex(Amel_scVI.obs_names)
# %%
Amel_scVI.obs["dumbbells"] = dumbbells.astype("category") # 防止作为连续变量处理
# %%
sc.pl.embedding(Amel_scVI, basis="X_umap", color=["dumbbells"], legend_loc="on data")
# %%
