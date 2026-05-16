# %%
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import median_abs_deviation
import anndata

import scanpy as sc
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor  # ← 用线程池
import matplotlib.pyplot as plt

# %%

def _read_single(file_path: str) -> tuple[str, anndata.AnnData]: #_worker 函数, 返回 (文件名, AnnData)
    name = Path(file_path).stem
    adata = sc.read_10x_h5(file_path)
    adata.var_names_make_unique()
    return name, adata

def load_h5_parallel(
    project_name: str, 
    suffix: str = ".h5", 
    max_workers: int = 8  # 线程可以开更多，开销小
):
    directory = f"{project_name}/h5_from_fastq2matrix"
    path = Path(directory)
    files = [str(f) for f in path.glob(f"*{suffix}") if f.is_file()]
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(_read_single, files)
        # 确保完全消费迭代器
        return dict(results)

# %%

def add_mito(
    project_name: str, 
    adata: anndata.AnnData,
    mito_prefix: str = "MT"
) -> anndata.AnnData:
    if adata.var_names.str.startswith(mito_prefix).sum() == 13:
        adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    elif adata.var_names.str.startswith(mito_prefix).sum() == 0:
        mito_genes = pd.read_csv(f"./{project_name}/metadata/mito.txt", header=None, names=["gene_id"])
        adata.var["mt"] = adata.var_names.isin(mito_genes["gene_id"])
    else:
        raise ValueError("Unable to determine the identification format for mitochondrial genes. Please check the data or provide a correct list of mitochondrial genes.")
    return adata

def cal_metrics(adata: anndata.AnnData) -> anndata.AnnData:
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt'], 
        inplace=True, # 直接修改原对象, 不返回新对象
        percent_top=[20], # 计算表达量排名前20基因占总表达量百分比, 判断表达量是否集中于少数基因
        log1p=True # 对部分指标进行log1p变换
    )
    return adata

def check_3_QC_covariates(
    project_name: str, 
    file_name: str,
    adata: anndata.AnnData
) -> None:
    out_dir = Path(f"./{project_name}/figures/1_before-filt")
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

# %%

Hsal_ann_dic = load_h5_parallel("Sheng_SA_2020_Hsal")
for key, adata in Hsal_ann_dic.items():
    adata = add_mito("Sheng_SA_2020_Hsal", adata)
    adata = cal_metrics(adata)
    check_3_QC_covariates("Sheng_SA_2020_Hsal", key, adata)


# %%
