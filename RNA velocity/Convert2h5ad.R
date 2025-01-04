library(Seurat)
library(SeuratDisk)
setwd('/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/RNA velocity')
# Read processed Seurat objects
d6.object.cdcs <- readRDS('../scRNA-seq analysis/Outputs/rds/Dnr6_cDCs.rds')
d10.object.cdcs <- readRDS('../scRNA-seq analysis/Outputs/rds/Dnr10_cDCs.rds')
d12.object.cdcs <- readRDS('../scRNA-seq analysis/Outputs/rds/Dnr12_cDCs.rds')
d6.object.cdc2s <- readRDS('../scRNA-seq analysis/Outputs/rds/Dnr6_cDC2s.rds')

# Donor 6 all cDCs
SaveH5Seurat(d6.object.cdcs, filename = "Outputs/h5ad/d6.object.cdcs.h5Seurat")
Convert("Outputs/h5ad/d6.object.cdcs.h5Seurat", dest = "h5ad")
# Donor 6 cDC2s
SaveH5Seurat(d6.object.cdc2s, filename = "Outputs/h5ad/d6.object.cdc2s.h5Seurat")
Convert("Outputs/h5ad/d6.object.cdc2s.h5Seurat", dest = "h5ad")

# Donor 10 all cDCs
SaveH5Seurat(d10.object.cdcs, filename = "Outputs/h5ad/d10.object.cdcs.h5Seurat")
Convert("Outputs/h5ad/d10.object.cdcs.h5Seurat", dest = "h5ad")

# Donor 12 all cDCs
SaveH5Seurat(d12.object.cdcs, filename = "Outputs/h5ad/d12.object.cdcs.h5Seurat")
Convert("Outputs/h5ad/d12.object.cdcs.h5Seurat", dest = "h5ad")