# %%
import scanpy as sc
from functions.concat_libs import *
# %%
# ------------------------------------------------------------------------------------------------------
# --------------------------------------------- Pipeline -----------------------------------------------
# ------------------------------------------------------------------------------------------------------
# %%
# ---------------------------------------- Harpegnathos venator ----------------------------------------
project_name = "species/Sheng_SA_2020_Hsal"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/1_base-filt-output", suffix=".h5ad", max_workers=8)
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
h5ad_concat.write(f"./{project_name}/data/3_concated-output/concat.h5ad")
# %%
# ---------------------------------------- Apis cerana ----------------------------------------
project_name = "species/Acer"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/1_base-filt-output", suffix=".h5ad", max_workers=8)
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
h5ad_concat.write(f"./{project_name}/data/3_concated-output/concat.h5ad")
# %%
# ---------------------------------------- Apis mellifera ----------------------------------------
# %%
project_name = "species/Zhang_iScience_2022_Amel"
h5ad_dic = load_h5_parallel(project_name, dir_name="data/1_base-filt-output", suffix=".h5ad", max_workers=8)
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
h5ad_concat.write(f"./{project_name}/data/3_concated-output/concat.h5ad")
# %%
