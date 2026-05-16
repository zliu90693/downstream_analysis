# %%
import anndata
import scanpy as sc
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor 

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

# %%

def pre_clust(
    adata: anndata.AnnData,
    reso: float = 0.5,
) -> None:
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_pcs=20)
    
    sc.tl.leiden(adata, key_added=f"leiden_{reso}", resolution=reso)

    n_clusters = adata.obs[f"leiden_{reso}"].nunique()
    print(f"Clustering Complete: {n_clusters} clusters (Recommended range: 10~30)")
    return adata

# %%

Hsal_h5ad_dic = load_h5_parallel("Sheng_SA_2020_Hsal", dir_name="1_base-filt-output", suffix=".h5ad", max_workers=8)

# %%

for key, adata in Hsal_h5ad_dic.items():
    pre_clust(adata, reso=1.0)

# %%

for key, adata in Hsal_h5ad_dic.items():
    adata.write_h5ad(f"./Sheng_SA_2020_Hsal/2_checkambient-output/{key}.h5ad")

# %%
