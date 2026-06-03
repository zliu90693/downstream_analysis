# %%
from functions.concat_libs import *
from functions.checkambient_visualize import *

# %%

# %%
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------
# %%
project_name = "Sheng_SA_2020_Hsal"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/2_checkambient-output", suffix="_decontX.h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    visualize_dX(project_name, key, adata)
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "Acer"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/2_checkambient-output", suffix="_decontX.h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    visualize_dX(project_name, key, adata)
# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
# %%
project_name = "Zhang_iScience_2022_Amel"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/2_checkambient-output", suffix="_decontX.h5ad", max_workers=8)
# %%
for key, adata in h5ad_dic.items():
    visualize_dX(project_name, key, adata)