# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html

library(Seurat)
library(SeuratData)
library(SeuratDisk)

# TS_Blood.h5ad filtered to ~20k genes with raw counts in AnnData.X
Convert("tabula_filtered.h5ad", dest = "h5seurat", overwrite = TRUE)
tabula <- LoadH5Seurat("tabula_filtered.h5seurat", misc=FALSE, meta.data=FALSE, tools=FALSE, neighbors=FALSE, graphs=FALSE)
tabula
saveRDS(tabula, file = "TS_Blood_filtered.rds")

