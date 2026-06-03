from concurrent.futures import ProcessPoolExecutor # 用多进程!! hdf5库和多线程配合极易引发死锁!!!
from pathlib import Path
import scanpy as sc
import anndata
from functions.concat_libs import *

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
    if len(files) <= 2:
        return {f: sc.read_h5ad(f) for f in files}
    
    # 文件多且大时，用进程池绕过GIL
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(_read_single_h5ad, files)
        return dict(results)