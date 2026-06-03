# %%
from functions.concat_libs import *
from functions.checkambient_preclust_rundX import *

# %%

# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------
# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------

project_name = "Sheng_SA_2020_Hsal"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/1_base-filt-output", suffix=".h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    print(key)
    preprocess_adata(adata)
    for reso in [0.25, 0.5, 0.8, 1.0]:
        run_leiden(adata, reso)
# %%
for key, adata in h5ad_dic.items():
    adata.write_h5ad(f"./{project_name}/data/2_checkambient-output/{key}.h5ad")
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "Acer"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/1_base-filt-output", suffix=".h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    print(key)
    preprocess_adata(adata)
    for reso in [0.2, 0.25, 0.5, 0.8, 1.0]:
        run_leiden(adata, reso)
# %%
for key, adata in h5ad_dic.items():
    adata.write_h5ad(f"./{project_name}/data/2_checkambient-output/{key}.h5ad")
# %%
run_decontX(project_name, cluster_col="leiden_0.20")
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "Zhang_iScience_2022_Amel"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/1_base-filt-output", suffix=".h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    print(key)
    preprocess_adata(adata)
    for reso in [0.2, 0.25, 0.5, 0.8, 1.0]:
        run_leiden(adata, reso)
# %%
for key, adata in h5ad_dic.items():
    adata.write_h5ad(f"./{project_name}/data/2_checkambient-output/{key}.h5ad")
# %%
run_decontX(project_name, cluster_col="leiden_0.50")
# %%
