import matplotlib.pyplot as plt
import scanpy as sc
import anndata
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor # 用多进程!! hdf5库和多线程配合极易引发死锁!!!

def visualize_dX(
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