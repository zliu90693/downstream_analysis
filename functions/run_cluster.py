import scanpy as sc
import anndata

def run_leiden(
    adata: anndata.AnnData,
    resolution: float = 0.5,
    flavor: str = "leidenalg",
    n_iterations: int = 2,
) -> None:
    sc.tl.leiden(
        adata, 
        key_added=f"leiden_res{resolution:.2f}",  # 防止浮点数转字符串时可能出现 res0.6000000000000001
        resolution=resolution, 
        flavor=flavor, 
        n_iterations=n_iterations
    )