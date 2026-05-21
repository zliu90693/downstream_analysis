# %%
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor # 用多进程!! hdf5库和多线程配合极易引发死锁!!!


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

    sc.pl.violin(adata, keys="decontX_contamination", groupby="decontX_clusters", show=False)
    plt.gcf().savefig(out_dir / f"{file_name}_contamination_violin.png", dpi=300, bbox_inches='tight')
    plt.close()

    print(adata.obs["decontX_contamination"].describe())

# %%

# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
Hsal_h5ad_dic = load_h5_parallel(project_name, dir_name="2_checkambient-output", suffix="_decontX.h5ad", max_workers=8)

# %%

for key, adata in Hsal_h5ad_dic.items():
    visualize(project_name, key, adata)

# %%
