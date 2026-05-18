# %%
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, Any

# %%
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE" 

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
    max_workers: int = 8  # 线程可以开更多，开销小
):
    directory = f"{project_name}/{dir_name}"
    path = Path(directory)
    files = [str(f) for f in path.glob(f"*{suffix}") if f.is_file()]
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(_read_single_h5ad, files)
        return dict(results)

def load_h5ad_parallel(
    project_name: str, 
    dir_name: str,
    suffix: str = ".h5ad", 
    max_workers: int = 8
) -> Dict[str, Any]:
    directory = Path(project_name) / dir_name
    if not directory.is_dir():
        raise FileNotFoundError(f"directory not exist: {directory}")
        
    files = [f for f in directory.glob(f"*{suffix}") if f.is_file()]
    if not files:
        raise ValueError(f"No files matching *{suffix} found in {directory}")

    results = {}
    # 使用 as_completed 便于逐文件捕获异常，避免单点失败丢弃全部结果
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 提交任务，返回 Future 对象
        future_to_path = {executor.submit(_read_single_h5ad, str(f)): f for f in files}
        
        for future in as_completed(future_to_path):
            file_path = future_to_path[future]
            try:
                key, adata = future.result()  # 必须返回 (key, value)
                results[key] = adata
            except Exception as e:
                print(f"Read failed {file_path.name}: {e}")
                # 可选：results[file_path.stem] = None 或记录日志
                
    return results

def visualize(
    project_name: str, 
    file_name: str,
    adata: anndata.AnnData
) -> None:
    out_dir = Path(f"./{project_name}/figures/2.3-1_checkambient")

    sc.pl.embedding(adata, basis="decontX_UMAP", color="decontX_contamination", cmap="Reds", show=False)
    plt.gcf().savefig(out_dir / f"{file_name}_contamination_umap.png", dpi=300, bbox_inches='tight')
    plt.close()

    sc.pl.embedding(adata, basis="decontX_UMAP", color="decontX_clusters", show=False)
    plt.gcf().savefig(out_dir / f"{file_name}_clusters_umap.png", dpi=300, bbox_inches='tight')
    plt.close()

    plt.hist(adata.obs["decontX_contamination"], bins=50)
    plt.xlabel("Contamination Score")
    plt.ylabel("Cell Count")
    plt.savefig(out_dir / f"{file_name}_contamination_hist.png", dpi=300, bbox_inches='tight')
    plt.close()

    sc.pl.violin(adata, keys="decontX_contamination", groupby="decontX_clusters")
    plt.gcf().savefig(out_dir / f"{file_name}_contamination_violin.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(adata.obs["decontX_contamination"].describe())

# %%

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
Hsal_h5ad_dic = load_h5ad_parallel(project_name, dir_name="2_checkambient-output", suffix="_decontX.h5ad", max_workers=8)

# %%


# %%

for key, adata in Hsal_h5ad_dic.items():
    visualize(project_name, key, adata)
