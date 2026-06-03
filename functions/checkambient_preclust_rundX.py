import anndata
import scanpy as sc
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor # 用多进程!! hdf5库和多线程配合极易引发死锁!!!
import subprocess

def preprocess_adata(adata: anndata.AnnData) -> None:
    """仅执行一次的数据预处理与降维"""
    if "X_pca" in adata.obsm:
        return  # 已预处理则跳过
        
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
        
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(adata, n_pcs=20)

def run_leiden(adata: anndata.AnnData, reso: float = 0.5) -> None:
    """仅运行聚类，可反复调用不同分辨率"""
    sc.tl.leiden(adata, key_added=f"leiden_{reso:.2f}", resolution=reso)
    n_clusters = adata.obs[f"leiden_{reso:.2f}"].nunique()
    print(f"Clustering Complete: {n_clusters} clusters, resolution: {reso} (Recommended range: 10~30)")

def run_decontX(
    project_name: str,
    cluster_col: str = "leiden_1.0"
) -> None:
    subprocess.run([
        "./2._checkambient-decontX.R", project_name, cluster_col
    ], check=True)