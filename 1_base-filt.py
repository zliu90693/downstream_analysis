# %%
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import median_abs_deviation
import anndata

import scanpy as sc
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor # 用多进程!! hdf5库和多线程配合极易引发死锁!!!
import matplotlib.pyplot as plt

# %%
def _read_single(file_path: str) -> tuple[str, anndata.AnnData]: #_worker 函数, 返回 (文件名, AnnData)
    name = Path(file_path).stem
    adata = sc.read_10x_h5(file_path)
    adata.var_names_make_unique()
    return name, adata

def load_h5_parallel(
    project_name: str,
    dir_name: str,
    suffix: str = ".h5ad",
    max_workers: int = 8  # 进程数不要超过CPU核心数
):
    directory = f"{project_name}/{dir_name}"
    path = Path(directory)
    files = [str(f) for f in path.glob(f"*{suffix}") if f.is_file()]
    
    print(f"Found {len(files)} files: {[Path(f).name for f in files]}")  
    
    # 文件小时，串行反而最快
    if len(files) <= 4:
        return {f: sc.read_h5ad(f) for f in files}
    
    # 文件多且大时，用进程池绕过GIL
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(_read_single, files)
        return dict(results)

def add_mito(
    project_name: str, 
    adata: anndata.AnnData,
    mito_prefix: str = "MT-"
) -> None:
    if adata.var_names.str.startswith(mito_prefix).sum() == 13:
        adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    elif adata.var_names.str.startswith(mito_prefix).sum() == 0:
        mito_genes = pd.read_csv(f"./{project_name}/metadata/mito.txt", header=None, names=["gene_id"])
        adata.var["mt"] = adata.var["gene_ids"].isin(mito_genes["gene_id"])
    else:
        raise ValueError("Unable to determine the identification format for mitochondrial genes. Please check the data or provide a correct list of mitochondrial genes.")

def cal_metrics(adata: anndata.AnnData) -> None:
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt'], 
        inplace=True, # 直接修改原对象, 不返回新对象
        percent_top=[20], # 计算表达量排名前20基因占总表达量百分比, 判断表达量是否集中于少数基因
        log1p=True # 对部分指标进行log1p变换
    )

def check_3_QC_covariates(
    project_name: str, 
    file_name: str,
    dir_name: str,
    adata: anndata.AnnData
) -> None:
    out_dir = Path(f"./{project_name}/figures/{dir_name}")
    out_dir.mkdir(parents=True, exist_ok=True)

    g = sns.displot(adata.obs["total_counts"], bins=100, kde=False) 
    g.figure.savefig(out_dir / f"{file_name}_total_counts_hist.png", dpi=300, bbox_inches='tight')
    plt.close(g.figure)

    sc.pl.violin(adata, 'total_counts', show=False)
    plt.gcf().savefig(out_dir / f"{file_name}_total_counts_violin.png", dpi=300, bbox_inches='tight')
    plt.close()

    sc.pl.violin(adata, "pct_counts_mt", show=False)
    plt.gcf().savefig(out_dir / f"{file_name}_pct_counts_mt_violin.png", dpi=300, bbox_inches='tight')
    plt.close()

    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False)
    plt.gcf().savefig(out_dir / f"{file_name}_scatter.png", dpi=300, bbox_inches='tight')
    plt.close()


def _is_outlier(
    adata: anndata.AnnData,
    metric: str, 
    nmads: int
) -> pd.Series:
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def add_outlier_column(
    adata: anndata.AnnData, 
    nmad: int = 5, 
    nmad_mt: int = 3, 
    pct_counts_mt: int = 20
) -> None:
    adata.obs["outlier"] = (
        _is_outlier(adata, "log1p_total_counts", nmad) # 注意, 这里的标准是log!
        | _is_outlier(adata, "log1p_n_genes_by_counts", nmad)
        | _is_outlier(adata, "pct_counts_in_top_20_genes", nmad)
    )
    adata.obs["mt_outlier"] = _is_outlier(adata, "pct_counts_mt", nmad_mt) | (
        adata.obs["pct_counts_mt"] > pct_counts_mt
    )
    print(adata.obs["outlier"].value_counts())
    print(adata.obs["mt_outlier"].value_counts())

def filter_outliers(adata: anndata.AnnData) -> anndata.AnnData:
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[~adata.obs["outlier"] & ~adata.obs["mt_outlier"]].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    return adata

# def check_3_QC_covariates(dir_name = "2_after_filt")

# %%

# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "Sheng_SA_2020_Hsal"
ann_dic = load_h5_parallel(project_name, dir_name="./data/0_h5_from_fastq2matrix", suffix=".h5", max_workers=8)
# %%
for key, adata in ann_dic.items():
    add_mito(project_name, adata)
    cal_metrics(adata)
    check_3_QC_covariates(project_name, key, "1-1_before-filt", adata)

# %%
for key, adata in ann_dic.items():
    add_outlier_column(adata, nmad=5, nmad_mt=3, pct_counts_mt=20)
    adata = filter_outliers(adata)
    ann_dic[key] = adata # 更新字典中的对象!!! 非常重要!!!
    # 否则, 局部变量 adata 确实指向了新对象，但字典 ann_dic[key] 的引用从未改变
    # adata = xxx 只是把标签 adata 贴到新对象上，不会反向修改原来持有该对象的地方（字典、列表、全局变量等）
# %%
for key, adata in ann_dic.items():
    check_3_QC_covariates(project_name, key, "1-2_after-filt", adata)
# %%
out_dir = Path(f"./{project_name}/data/1_base-filt-output")
out_dir.mkdir(parents=True, exist_ok=True)
for key, adata in ann_dic.items():
    adata.write_h5ad(f"./{project_name}/data/1_base-filt-output/{key}.h5ad")

# %%

# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "Acer"
ann_dic = load_h5_parallel(project_name, dir_name="./data/0_h5_from_fastq2matrix", suffix=".h5", max_workers=8)
# %%
for key, adata in ann_dic.items():
    add_mito(project_name, adata)
    cal_metrics(adata)
    check_3_QC_covariates(project_name, key, "1-1_before-filt", adata)

# %%
for key, adata in ann_dic.items():
    add_outlier_column(adata, nmad=5, nmad_mt=3, pct_counts_mt=10)
    adata = filter_outliers(adata)
    ann_dic[key] = adata # 更新字典中的对象!!! 非常重要!!!
    # 否则, 局部变量 adata 确实指向了新对象，但字典 ann_dic[key] 的引用从未改变
    # adata = xxx 只是把标签 adata 贴到新对象上，不会反向修改原来持有该对象的地方（字典、列表、全局变量等）
# %%
for key, adata in ann_dic.items():
    check_3_QC_covariates(project_name, key, "1-2_after-filt", adata)
# %%
out_dir = Path(f"./{project_name}/data/1_base-filt-output")
out_dir.mkdir(parents=True, exist_ok=True)
for key, adata in ann_dic.items():
    adata.write_h5ad(f"./{project_name}/data/1_base-filt-output/{key}.h5ad")

# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
# %%
project_name = "Zhang_iScience_2022_Amel"
ann_dic = load_h5_parallel(project_name, dir_name="./data/0_h5_from_fastq2matrix", suffix=".h5", max_workers=8)
# %%
for key, adata in ann_dic.items():
    add_mito(project_name, adata)
    cal_metrics(adata)
    check_3_QC_covariates(project_name, key, "1-1_before-filt", adata)

# %%
for key, adata in ann_dic.items():
    add_outlier_column(adata, nmad=5, nmad_mt=3, pct_counts_mt=15)
    adata = filter_outliers(adata)
    ann_dic[key] = adata # 更新字典中的对象!!! 非常重要!!!
    # 否则, 局部变量 adata 确实指向了新对象，但字典 ann_dic[key] 的引用从未改变
    # adata = xxx 只是把标签 adata 贴到新对象上，不会反向修改原来持有该对象的地方（字典、列表、全局变量等）
# %%
for key, adata in ann_dic.items():
    check_3_QC_covariates(project_name, key, "1-2_after-filt", adata)
# %%
out_dir = Path(f"./{project_name}/data/1_base-filt-output")
out_dir.mkdir(parents=True, exist_ok=True)
for key, adata in ann_dic.items():
    adata.write_h5ad(f"./{project_name}/data/1_base-filt-output/{key}.h5ad")
# %%
