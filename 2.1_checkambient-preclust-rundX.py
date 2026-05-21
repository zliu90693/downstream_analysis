# %%
import anndata
import scanpy as sc
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor # 用多进程!! hdf5库和多线程配合极易引发死锁!!!
import subprocess

# %%

def _read_single_h5ad(file_path: str) -> tuple[str, anndata.AnnData]: #_worker 函数, 返回 (文件名, AnnData)
    name = Path(file_path).stem
    adata = sc.read_h5ad(file_path)
    # adata.var_names_make_unique()
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
        results = executor.map(_read_single_h5ad, files)
        return dict(results)

# %%

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

# %%
def run_decontX(
    project_name: str,
    cluster_col: str = "leiden_1.0"
) -> None:
    subprocess.run([
        "./2._checkambient-decontX.R", project_name, cluster_col
    ], check=True)

# %%

# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
h5ad_dic = load_h5_parallel(project_name, dir_name="1_base-filt-output", suffix=".h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    print(key)
    preprocess_adata(adata)
    for reso in [0.25, 0.5, 0.8, 1.0]:
        run_leiden(adata, reso)
# %%
for key, adata in h5ad_dic.items():
    adata.write_h5ad(f"./{project_name}/2_checkambient-output/{key}.h5ad")
# %%

# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "Acer"
h5ad_dic = load_h5_parallel(project_name, dir_name="1_base-filt-output", suffix=".h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    print(key)
    preprocess_adata(adata)
    for reso in [0.2, 0.25, 0.5, 0.8, 1.0]:
        run_leiden(adata, reso)
# %%
run_decontX(project_name, cluster_col="leiden_0.20")
# %%
for key, adata in h5ad_dic.items():
    adata.write_h5ad(f"./{project_name}/2_checkambient-output/{key}.h5ad")

# %%
