# %%
import scanpy as sc
from functions.run_cluster import run_leiden
# %%
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "species/Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/5_scVI-output/concat.h5ad")
# %%
# sc.tl.leiden(combined_h5ad, resolution=0.5)
# fig = sc.pl.umap(combined_h5ad, legend_loc="on data", color='leiden', return_fig=True)
# %%
for reso in [0.25, 0.50, 1.00]:
    run_leiden(
        adata=combined_h5ad,
        resolution=reso,
        flavor="leidenalg",
        n_iterations=2
    )
    fig = sc.pl.umap(combined_h5ad, legend_loc="on data", color=f'leiden_res{reso:.2f}', return_fig=True)
    fig.savefig(f"./{project_name}/figures/6_cluster/umap_integrated_res{reso}.pdf", bbox_inches='tight')

combined_h5ad.write_h5ad(f"./{project_name}/data/6_cluster-output/concat.h5ad")

# sc.tl.leiden使用的是 adata.obsp['connectivities'], 是 sc.pp.neighbors() 产生的
# %%
# fig = sc.pl.umap(combined_h5ad, legend_loc="on data", color='leiden_res0.25', return_fig=True)
# fig.savefig(f"./{project_name}/figures/5_scVI/umap_integrated_batch.pdf", bbox_inches='tight')
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
project_name = "species/Acer"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/5_scVI-output/concat.h5ad")
# %%
for reso in [0.25, 0.50, 1.00]:
    run_leiden(
        adata=combined_h5ad,
        resolution=reso,
        flavor="leidenalg",
        n_iterations=2
    )
    fig = sc.pl.umap(combined_h5ad, legend_loc="on data", color=f'leiden_res{reso:.2f}', return_fig=True)
    fig.savefig(f"./{project_name}/figures/6_cluster/umap_integrated_res{reso}.pdf", bbox_inches='tight')

combined_h5ad.write_h5ad(f"./{project_name}/data/6_cluster-output/concat.h5ad")

# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
project_name = "species/Zhang_iScience_2022_Amel"
combined_h5ad = sc.read_h5ad(f"./{project_name}/data/5_scVI-output/concat.h5ad")
# %%
for reso in [0.25, 0.50, 1.00]:
    run_leiden(
        adata=combined_h5ad,
        resolution=reso,
        flavor="leidenalg",
        n_iterations=2
    )
    fig = sc.pl.umap(combined_h5ad, legend_loc="on data", color=f'leiden_res{reso:.2f}', return_fig=True)
    fig.savefig(f"./{project_name}/figures/6_cluster/umap_integrated_res{reso}.pdf", bbox_inches='tight')

combined_h5ad.write_h5ad(f"./{project_name}/data/6_cluster-output/concat.h5ad")

# %%
