# %%

R

# %%

library(celda)
library(SingleCellExperiment)
library(Seurat)
library(anndataR)

library(zellkonverter)
library(future.apply)
library(tools)

# %%

load_h5ad_parallel <- function(project_name, suffix = ".h5ad", max_workers = 8) {
    # 1. 构建目录路径 & 获取文件列表
    directory <- file.path(project_name, "2_checkambient-output")
    files <- list.files(directory, pattern = paste0("\\", suffix, "$"), full.names = TRUE)
    
    if (length(files) == 0) stop("No matching files found. Please check the path or file extension.")
    
    # 2. 配置并行后端（多进程，类似 Python 的 ProcessPool）
    plan(multisession, workers = max_workers)
    on.exit(plan(sequential)) # 函数退出时自动恢复串行模式
    
    # 3. 并行读取 .h5ad → 自动转为 Anndata 对象
    adata_list <- future_lapply(files, function(f) {
        read_h5ad(f)
    })
    
    # 4. 提取文件名（去后缀）作为列表的“键”
    names(adata_list) <- file_path_sans_ext(basename(files))
    
    return(adata_list)
}

# %%

Anndata_2_SCE <- function(anndata) {
    counts_raw <- anndata$layers[["counts"]]
    if (is.null(counts_raw)) {
        stop("'counts' not found in `anndata$layers`. Please ensure that `adata.layers['counts'] = adata.X.copy()` has been executed on the Python side.")
    }

    counts_mtx <- Matrix::t(counts_raw) # 根据pbmc的例子, decontX需要输入原始计数矩阵(未经过normalize)
    if (!inherits(counts_mtx, "dgCMatrix")) {
        counts_mtx <- as(counts_mtx, "CsparseMatrix")
    }
    
    sce <- SingleCellExperiment(
        assays = list(counts = counts_mtx),
        colData = as.data.frame(anndata$obs),
        rowData = as.data.frame(anndata$var)
    )
    # 显式同步行列名（防 anndataR 版本差异导致名称错位）
    rownames(sce) <- rownames(anndata$var)
    colnames(sce) <- rownames(anndata$obs)
    return(sce)
}

# %%

materialize <- function(sce, assay_name) {
    # Matrix::t(as(assay(sce, assay_name), "CsparseMatrix"))
    mtx <- assay(sce, assay_name)
    if (!inherits(mtx, "dgCMatrix")) {
        mtx <- as(mtx, "CsparseMatrix")
    }
    return(Matrix::t(mtx))
}

SCE_2_Anndata <- function(sce) {
    counts_t <- materialize(sce, "counts")
    decontXcounts_t <- materialize(sce, "decontXcounts")
    decontX_meta <- metadata(sce)$decontX
    decontX_meta$runParams$logfile <- NULL
    # decontX_meta$runParams$z <- NULL
    anndata <- AnnData(
        X      = counts_t,
        obs    = as.data.frame(colData(sce)),
        var    = as.data.frame(rowData(sce)),
        layers = list(decontXcounts = decontXcounts_t),
        obsm   = list(X_decontX_UMAP = reducedDim(sce, "decontX_UMAP")),
        uns    = list(decontX = decontX_meta)
    )
    return(anndata)
}

# %%

setwd("/home/liuzhiyu/Projects/neo_caste/downstream_analysis")

# %%

# -----------------------------------------------------------------------------
# ---------------------------    Workflow    ----------------------------------
# -----------------------------------------------------------------------------


Hsal_h5ad_dic <- load_h5ad_parallel("Sheng_SA_2020_Hsal")

# %%

Hsal_sce_dic <- lapply(Hsal_h5ad_dic, Anndata_2_SCE)  # 直接覆盖原对象
gc() # 手动触发垃圾回收，释放未使用的内存

# %%

cluster_col <- "leiden_0.8"
Hsal_decont_dic <- lapply(Hsal_sce_dic, function(sce) {
    # 1. 安全检查：确认聚类列存在于 colData 中
    if (!cluster_col %in% colnames(colData(sce))) {
        stop(paste("Sample", names(Hsal_sce_dic)[which(Hsal_sce_dic == sce)], 
                "without ", cluster_col, "column in colData"))
    }
    
    # 2. 提取聚类标签（转为字符向量，避免因子水平干扰 decontX 内部逻辑）
    z <- as.character(colData(sce)[[cluster_col]])
    
    # 3. 运行 decontX（自动读取 "counts" assay，使用自定义 z 标签）
    sce <- decontX(x = sce, z = z)
    
    return(sce)
})

# %%

lapply(Hsal_decont_dic, function(sce) {
    summary(sce$decontX_contamination)
})

# %%

Hsal_deconth5ad_dic <- lapply(Hsal_decont_dic, SCE_2_Anndata)

# %%

output_dir <- file.path("Sheng_SA_2020_Hsal", "2_checkambient-output")
for (key in names(Hsal_deconth5ad_dic)) {
    adata <- Hsal_deconth5ad_dic[[key]]
    filename <- file.path(output_dir, paste0(key, "_decontX.h5ad"))
    
    if (file.exists(filename)) {
        message(sprintf("Skip, : %s exists", filename))
        next
    }
    
    tryCatch({
        write_h5ad(adata, filename, mode = "w")
        message(sprintf("%s Saved", filename))
    }, error = function(e) {
        warning(sprintf("Error while saving [%s]: %s", key, e$message))
    })
    }
