# %%

import scanpy as sc

# %%

F1_adata = sc.read_h5ad("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm-Fs-Dr_F1.h5ad")
F1_adata

# %%

sc.pp.neighbors(F1_adata, n_pcs=30)
sc.tl.umap(F1_adata)

# %%

# sc.tl.leiden(F1_adata, flavor="igraph", n_iterations=2) 

# n_iterations=2, 细化过程执行次数
# 默认resolution参数是1

# %%

sc.tl.leiden(
    F1_adata, key_added="leiden_res0_25", resolution=0.25, flavor="igraph", n_iterations=2
)
sc.tl.leiden(
    F1_adata, key_added="leiden_res0_5", resolution=0.5, flavor="igraph", n_iterations=2
)
sc.tl.leiden(
    F1_adata, key_added="leiden_res1", resolution=1.0, flavor="igraph", n_iterations=2
)

# %%

sc.pl.umap(
    F1_adata, 
    color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"],
    legend_loc="on data",
)

# %%

F1_adata.write("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm-Fs-Dr-Cl_F1.h5ad")

# %%
