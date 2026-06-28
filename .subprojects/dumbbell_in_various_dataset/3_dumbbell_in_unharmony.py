# %%
import scanpy as sc
import numpy as np
import pandas as pd
# %%
Amel = sc.read_h5ad("./data/Amel_dumbbell_marked.h5ad")
# %%
sc.pl.embedding(Amel, basis="X_umap", color=["dumbbells"], legend_loc="on data")
# %%
sc.pl.embedding(Amel, basis="X_umap_scran", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_log1p", color=["dumbbells"], legend_loc="on data")
sc.pl.embedding(Amel, basis="X_umap_pearson", color=["dumbbells"], legend_loc="on data")
# %%
