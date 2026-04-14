# %%

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import median_abs_deviation

import scanpy as sc
import scvi

# sc.settings.verbosity = 0
# sc.settings.set_figure_params(dpi=150, facecolor="white", frameon=False)

# %%

mito_genes = pd.read_csv("./Zhang_iScience_2022_Amel/metadata/mito.txt", header=None, names=["gene_id"])

# %%

F1_adata = sc.read_10x_h5("./Zhang_iScience_2022_Amel/h5_from_fastq2matrix/F1.h5")
F1_adata.var_names_make_unique()
F1_adata

# %%

# F2_adata = sc.read_10x_h5("./Zhang_iScience_2022_Amel/h5_from_fastq2matrix/F2.h5")
# F2_adata

# %%

F1_adata.var["mt"] = F1_adata.var["gene_ids"].isin(mito_genes["gene_id"])


# %%

sc.pp.calculate_qc_metrics(
        F1_adata, 
        qc_vars=['mt'], 
        inplace=True, # 直接修改原对象, 不返回新对象
        percent_top=[20], # 计算表达量排名前20基因占总表达量百分比, 判断表达量是否集中于少数基因
        log1p=True # 对部分指标进行log1p变换
    )

# %%

print(F1_adata.obs['total_counts'].describe()) # 检查其中是否包含空液滴

# %%
# 绘制_genes_by_counts、total_counts 和 pct_counts_mt, 三个QC相关协变量

import numpy as np
import seaborn as sns
from scipy.stats import median_abs_deviation

def check_3_QC_covariates(adata):
    sns.displot(adata.obs["total_counts"], bins=100, kde=False) 
    sc.pl.violin(adata, 'total_counts')
    sc.pl.violin(adata, "pct_counts_mt")
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# %%

check_3_QC_covariates(F1_adata)

# %%

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

# %%

def add_outlier_column(adata, nmad, nmad_mt, pct_counts_mt):
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", nmad) # 注意, 这里的标准是log!
        | is_outlier(adata, "log1p_n_genes_by_counts", nmad)
        | is_outlier(adata, "pct_counts_in_top_20_genes", nmad)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", nmad_mt) | (
        adata.obs["pct_counts_mt"] > pct_counts_mt 
    )
    print(adata.obs["outlier"].value_counts())
    print(adata.obs["mt_outlier"].value_counts())


# %%

add_outlier_column(F1_adata, nmad=5, nmad_mt=3, pct_counts_mt=15) # Outliner: 3693; 6035(mt)

# %%

sc.pl.scatter(F1_adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# %%

def filt_with_outliner_mtoutliner(adata):
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy() # 注意, 这里不是原地修改, 这是创建了一个新对象, 如果不返回, 这个对象在函数内会被丢弃
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    return adata

# %%

F1_adata = filt_with_outliner_mtoutliner(F1_adata)

# %%

F1_adata.obs["outlier"].value_counts() # 确认完成过滤

# %%

check_3_QC_covariates(F1_adata)

# %%

F1_adata.write("./Zhang_iScience_2022_Amel/h5_QC/ft_F1.h5ad")

# %%

