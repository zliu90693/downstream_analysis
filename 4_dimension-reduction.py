# %%

import scanpy as sc

# %%

F1_adata = sc.read_h5ad("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm-Fs_F1.h5ad") # 在进入降维等步骤之前，需要将目标归一化层的内容赋值给 adata.X
F1_adata

# %%

F1_adata.X = F1_adata.layers["scran_normalization"] # 替换(覆盖)原有的F1_adata.X, 因为降维需要Normalize后的数据

# %%

### ------------------ PCA, 线性降维, 选前 10-50 个主成分 (PC)，用于后续的分析任务 ------------------

# 将highly_variable赋值为highly_deviant, 以便在sc.pp.pca 中使用 scanpy 的 'use_highly_variable' 参数

F1_adata.var["highly_variable"] = F1_adata.var["highly_deviant"]
sc.pp.pca(F1_adata, svd_solver="arpack", mask_var="highly_variable")

# %%

sc.pl.pca_scatter(F1_adata, color="total_counts")

# %%

# ------------------ t-SNE ------------------

sc.tl.tsne(F1_adata, use_rep="X_pca")
sc.pl.tsne(F1_adata, color="total_counts")

# %%

# ------------------ UMAP ------------------

sc.pp.neighbors(F1_adata) # 先做neighbors再做UMAP
sc.tl.umap(F1_adata)

# %%

sc.pl.umap(F1_adata, color="total_counts")

# %%

F1_adata

# %%

F1_adata.write("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm-Fs-Dr_F1.h5ad")

# %%
