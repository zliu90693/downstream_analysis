# %%
from concurrent.futures import ProcessPoolExecutor # 用多进程!! hdf5库和多线程配合极易引发死锁!!!
from pathlib import Path
import scanpy as sc
import anndata
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
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------
# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "Sheng_SA_2020_Hsal"
h5ad_dic = load_h5_parallel(project_name, dir_name="1_base-filt-output", suffix=".h5ad", max_workers=8)
# %%
sample_keys = list(h5ad_dic.keys())
adata_list = list(h5ad_dic.values())
# %%
h5ad_concat = sc.concat(
    adata_list,
    keys=sample_keys,
    index_unique='_',
    label='batch',
    merge='same'
)
# %%
h5ad_concat.layers["counts"] = h5ad_concat.X.copy()
# %%
sc.pp.filter_genes(h5ad_concat, min_cells=20) # 过滤掉至少 20 个细胞中未检测到的基因，因为这些基因不具有参考价值, 在极少数细胞中检测到的基因通常是技术噪声、环境 RNA 污染或随机低水平转录的结果
# %%
h5ad_concat.write(f"./{project_name}/3_concated-output/concat.h5ad")
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
project_name = "Acer"
h5ad_dic = load_h5_parallel(project_name, dir_name="1_base-filt-output", suffix=".h5ad", max_workers=8)
# %%
sample_keys = list(h5ad_dic.keys())
adata_list = list(h5ad_dic.values())
# %%
h5ad_concat = sc.concat(
    adata_list,
    keys=sample_keys,
    index_unique='_',
    label='batch',
    merge='same'
)
# %%
h5ad_concat.layers["counts"] = h5ad_concat.X.copy()
# %%
sc.pp.filter_genes(h5ad_concat, min_cells=20) # 过滤掉至少 20 个细胞中未检测到的基因，因为这些基因不具有参考价值, 在极少数细胞中检测到的基因通常是技术噪声、环境 RNA 污染或随机低水平转录的结果
# %%
h5ad_concat.write(f"./{project_name}/3_concated-output/concat.h5ad")
# %%
