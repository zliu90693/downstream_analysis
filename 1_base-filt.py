# %%
from pathlib import Path
from functions.base_filt import *
# def check_3_QC_covariates(dir_name = "2_after_filt")

# %%

# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------

# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "species/Sheng_SA_2020_Hsal"
ann_dic = load_h5_parallel(project_name, dir_name="./data/0_h5_from_fastq2matrix", suffix=".h5", max_workers=8)
# %%
for key, adata in ann_dic.items():
    add_mito(project_name, adata)
    cal_metrics(adata)
    check_3_QC_covariates(project_name, key, "1-1_before-filt", adata)

# %%
for key, adata in ann_dic.items():
    add_outlier_column(adata, nmad=5, nmad_mt=3, pct_counts_mt=20)
    adata = filter_outliers(adata)
    ann_dic[key] = adata # 更新字典中的对象!!! 非常重要!!!
    # 否则, 局部变量 adata 确实指向了新对象，但字典 ann_dic[key] 的引用从未改变
    # adata = xxx 只是把标签 adata 贴到新对象上，不会反向修改原来持有该对象的地方（字典、列表、全局变量等）
# %%
for key, adata in ann_dic.items():
    check_3_QC_covariates(project_name, key, "1-2_after-filt", adata)
# %%
out_dir = Path(f"./{project_name}/data/1_base-filt-output")
out_dir.mkdir(parents=True, exist_ok=True)
for key, adata in ann_dic.items():
    adata.write_h5ad(f"./{project_name}/data/1_base-filt-output/{key}.h5ad")

# %%

# ---------------------------------------- Apis cerana ----------------------------------------
# %%
project_name = "species/Acer"
ann_dic = load_h5_parallel(project_name, dir_name="./data/0_h5_from_fastq2matrix", suffix=".h5", max_workers=8)
# %%
for key, adata in ann_dic.items():
    add_mito(project_name, adata)
    cal_metrics(adata)
    check_3_QC_covariates(project_name, key, "1-1_before-filt", adata)

# %%
for key, adata in ann_dic.items():
    add_outlier_column(adata, nmad=5, nmad_mt=3, pct_counts_mt=10)
    adata = filter_outliers(adata)
    ann_dic[key] = adata # 更新字典中的对象!!! 非常重要!!!
    # 否则, 局部变量 adata 确实指向了新对象，但字典 ann_dic[key] 的引用从未改变
    # adata = xxx 只是把标签 adata 贴到新对象上，不会反向修改原来持有该对象的地方（字典、列表、全局变量等）
# %%
for key, adata in ann_dic.items():
    check_3_QC_covariates(project_name, key, "1-2_after-filt", adata)
# %%
out_dir = Path(f"./{project_name}/data/1_base-filt-output")
out_dir.mkdir(parents=True, exist_ok=True)
for key, adata in ann_dic.items():
    adata.write_h5ad(f"./{project_name}/data/1_base-filt-output/{key}.h5ad")

# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
# %%
project_name = "species/Zhang_iScience_2022_Amel"
ann_dic = load_h5_parallel(project_name, dir_name="./data/0_h5_from_fastq2matrix", suffix=".h5", max_workers=8)
# %%
for key, adata in ann_dic.items():
    add_mito(project_name, adata)
    cal_metrics(adata)
    check_3_QC_covariates(project_name, key, "1-1_before-filt", adata)

# %%
for key, adata in ann_dic.items():
    add_outlier_column(adata, nmad=5, nmad_mt=3, pct_counts_mt=15)
    adata = filter_outliers(adata)
    ann_dic[key] = adata # 更新字典中的对象!!! 非常重要!!!
    # 否则, 局部变量 adata 确实指向了新对象，但字典 ann_dic[key] 的引用从未改变
    # adata = xxx 只是把标签 adata 贴到新对象上，不会反向修改原来持有该对象的地方（字典、列表、全局变量等）
# %%
for key, adata in ann_dic.items():
    check_3_QC_covariates(project_name, key, "1-2_after-filt", adata)
# %%
out_dir = Path(f"./{project_name}/data/1_base-filt-output")
out_dir.mkdir(parents=True, exist_ok=True)
for key, adata in ann_dic.items():
    adata.write_h5ad(f"./{project_name}/data/1_base-filt-output/{key}.h5ad")
# %%
