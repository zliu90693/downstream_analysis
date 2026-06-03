#!/home/liuzhiyu/Software/miniconda3/envs/rm_ambient_doublet/bin/Rscript

.libPaths(c("/home/liuzhiyu/Software/miniconda3/envs/rm_ambient_doublet/lib/R/library", .libPaths()))

library(Matrix)
library(SingleCellExperiment)
library(scry)

args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]

X_sparse <- readMM(paste0("./", project_name, "/metadata/5_fs-data-mat.mtx"))
obs_df <- read.csv(
    paste0("./", project_name, "/metadata/5_fs-cellnames.csv"),
    row.names = 1,          # 将第一列作为行名
    stringsAsFactors = FALSE
)
var_df <- read.csv(
    paste0("./", project_name, "/metadata/5_fs-genenames.csv"),
    row.names = 1,          # 将第一列作为行名
    stringsAsFactors = FALSE
)

rownames(X_sparse) <- rownames(var_df)
colnames(X_sparse) <- rownames(obs_df)

stopifnot(identical(rownames(X_sparse), rownames(var_df)))
stopifnot(identical(colnames(X_sparse), rownames(obs_df)))
# 赋名之后立即验证赋名是否成功且一致

sce <- SingleCellExperiment(
    assays = list(X = X_sparse),
    colData = obs_df,
    rowData = var_df
)

sce <- devianceFeatureSelection(sce, assay = "X")

deviance_df <- data.frame(
    deviance = rowData(sce)$binomial_deviance,
    row.names = rownames(sce)
)
write.csv(deviance_df, paste0("./", project_name, "/metadata/5_fs_deviance_results.csv"))