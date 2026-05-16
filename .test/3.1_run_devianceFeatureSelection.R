# %%

R

# %%

library(Matrix)
library(SingleCellExperiment)
library(scry)

# %%

setwd("/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Zhang_iScience_2022_Amel/metadata")

# %%

# ── 读取文件 ────────────────────────────────────────────────
X <- as(readMM("X_sparse.mtx"), "CsparseMatrix")
obs <- read.csv("obs.csv", row.names = 1)
var <- read.csv("var.csv", row.names = 1)

# %%

# ── 构建 SCE 并运行 ─────────────────────────────────────────
sce <- SingleCellExperiment(
    assays = list(X = X),
    colData = obs,
    rowData = var
)
sce <- devianceFeatureSelection(sce, assay = "X")

# %%

# ── 输出结果 ────────────────────────────────────────────────
deviance_df <- data.frame(
    deviance = rowData(sce)$binomial_deviance,
    row.names = rownames(sce)
)
write.csv(deviance_df, "deviance_results.csv")
