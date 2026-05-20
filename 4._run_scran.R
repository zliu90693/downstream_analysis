#!/home/liuzhiyu/Software/miniconda3/envs/rm_ambient_doublet/bin/Rscript

.libPaths(c("/home/liuzhiyu/Software/miniconda3/envs/rm_ambient_doublet/lib/R/library", .libPaths()))

library(scran)
library(BiocParallel)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
project_name <- args[1]

data_mat <- readMM(paste0("./", project_name, "/metadata/4_scran-data-mat.mtx"))
cells_df <- read.csv(paste0("./", project_name, "/metadata/4_cellnames.csv"), stringsAsFactors = FALSE)
colnames(data_mat) <- cells_df$cell_id
input_groups_df <- read.csv(paste0("./", project_name, "/metadata/4_scran-groups.csv"), row.names = 1)

if (!all(colnames(data_mat) %in% rownames(input_groups_df))) {
    stop("Groups file missing cells present in matrix!")
}

input_groups <- input_groups_df[colnames(data_mat), "groups"]

sce <- SingleCellExperiment(list(counts = data_mat))
size_factors <- sizeFactors(
    computeSumFactors(
        sce,
        clusters = input_groups,
        min.mean = 0.1,
        BPPARAM = MulticoreParam()
    )
)

if (any(is.na(size_factors))) {
    stop("Detected NA size factors - checking clusters...")
}
if (any(size_factors <= 0)) {
    stop("Detected non-positive size factors - checking clusters...")
}

write.csv(
    data.frame(
        cell_id = colnames(data_mat), 
        size_factors = size_factors,
        stringsAsFactors = FALSE
    ),
    paste0("./", project_name, "/metadata/4_scran-size-factors.csv"),
    row.names = FALSE
)

