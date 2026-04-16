# %%

library(scran)
library(BiocParallel)
library(Matrix)

# %%

setwd("/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Zhang_iScience_2022_Amel")
# ── 读取 Python 传来的文件 ──────────────────────────────────
data_mat <- readMM("./metadata/F1_data_mat.mtx")
input_groups_df <- read.csv("./metadata/input_groups.csv", row.names = 1)
input_groups <- input_groups_df$groups

# %%

# ── 计算 size factors ───────────────────────────────────────
sce <- SingleCellExperiment(list(counts = data_mat))
size_factors <- sizeFactors(
    computeSumFactors(
        sce,
        clusters = input_groups,
        min.mean = 0.1,
        BPPARAM = MulticoreParam()
    )
)

# 对于每个细胞, size_factors的含义是该细胞相对于数据集"平均水平"的测序深度;
# 大于1, 说明深度高于平均水平, 需向下缩放;
# 小于1, 深度低于平均水平, 向上缩放;

# %%

length(size_factors)
summary(size_factors) # 确保无NaN

# %%

# ── 输出结果 ────────────────────────────────────────────────
write.csv(
    data.frame(size_factors = size_factors),
    "./metadata/size_factors.csv"
)
