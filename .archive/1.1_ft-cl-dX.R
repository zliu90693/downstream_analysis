# %%

R
library(celda)
library(SingleCellExperiment)
library(Seurat)
library(anndataR)

# %%

setwd("/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Zhang_iScience_2022_Amel/h5_QC")

# %%

F1_adata <- read_h5ad("ft-cl_F1.h5ad")

# %%

Anndata_2_SCE <- function(anndata) {
    counts_mtx <- Matrix::t(anndata$layers[["counts"]]) # 根据pbmc的例子, decontX需要输入原始计数矩阵(未经过normalize)
    
    if (!inherits(counts_mtx, "dgCMatrix")) {
        counts_mtx <- as(counts_mtx, "CsparseMatrix")
    }
    
    sce <- SingleCellExperiment(
        assays = list(counts = counts_mtx),
        colData = as.data.frame(anndata$obs),
        rowData = as.data.frame(anndata$var)
    )

    rownames(sce) <- rownames(anndata$var)
    colnames(sce) <- rownames(anndata$obs)
    
    return(sce)
}

# %%

F1_sce <- Anndata_2_SCE(F1_adata)

# %%

F1_sce_dX <- decontX(x = F1_sce, z = F1_sce$leiden_0.5)

# %%

summary(F1_sce_dX$decontX_contamination)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 8.061e-05 6.794e-03 4.889e-02 1.135e-01 1.521e-01 9.855e-01 

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

F1_sce_dX_adata <- SCE_2_Anndata(F1_sce_dX)

# %%

write_h5ad(F1_sce_dX_adata, "ft-cl-dX_F1_0.5.h5ad", mode = "w")
