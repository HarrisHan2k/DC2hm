setwd("/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/CITE-seq")
library(pacman)
p_load(Seurat, ggplot2, ggplotify, ggpubr, harmony, 
       DoubletFinder, future, dplyr, tidyverse,ggrastr, openxlsx,
       RColorBrewer, ggrastr,CytoTRACE, monocle3,ClusterGVis,
       clusterProfiler, ggsci, ComplexHeatmap, circlize,scDblFinder, dsb,
       ggridges, monocle3, SeuratWrappers)
set.seed(0822)
output_dir <- ''
# QC for Cellranger Outputs and create filtered RDS files -----------------------------------------------
work.list <- c("D6","D23")

# Define the function to compute and store initial criteria for QC
PlotQC <- function(project){
  matrix <- Read10X(paste0('01_CellRanger/',project,'/filtered_feature_bc_matrix'))
  object <- CreateSeuratObject(matrix$`Gene Expression`, project = project, min.cells = 3, 
                               min.features = 200)
  adt.assay <- matrix$`Antibody Capture`
  adt.assay <- adt.assay[,colnames(object)]
  object[['ADT']] <- CreateAssayObject(adt.assay)
  object[['percent.mt']] <- PercentageFeatureSet(object ,pattern = '^MT-')
  vlnplots <- VlnPlot(object, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
                      combine = F)
  figure <- ggarrange(plotlist = vlnplots, ncol = 3, common.legend = T, legend = 'none')
  qc.vlnplot <- annotate_figure(figure, top = paste0('QC_Violin_Plot_sample_', project))
  scatter.1 <- FeatureScatter(object, feature1 = 'nCount_RNA', 
                              feature2 = 'percent.mt')
  scatter.2 <- FeatureScatter(object, feature1 = 'nCount_RNA', 
                              feature2 = 'nFeature_RNA')
  pdf(paste0(output_dir,project,'_QC_Plots.pdf'))
  plot(qc.vlnplot)
  plot(scatter.1)
  plot(scatter.2)
  dev.off()
}
for (i in seq(1,length(work.list))) {
  PlotQC(project = work.list[[i]])
}
# Define thresholds for QC
th.nfeature.mins <- rep(500, 2)
th.nfeature.maxs <- rep(6000,2)
th.mt.pcts <- rep(10,2)
th.count.mins <- rep(2000,2)
threshold.df <- data.frame(th.nfeature.mins=th.nfeature.mins,
                           th.nfeature.maxs=th.nfeature.maxs,
                           th.mt.pct=th.mt.pcts,
                           th.count.mins=th.count.mins,
                           row.names = work.list)
write.csv(threshold.df, paste0(output_dir,'QC_Thresholds.csv'))

# Actually perform QC, store RDS files, and thresholds used for QC
cell.numbers.before <- c()
cell.numbers.after <- c()
# Define the function
DoQCAndStore <- function(project,th.nfeature.min,
                         th.nfeature.max, th.mt.pct,th.count.min){
  matrix <- Read10X(paste0('01_CellRanger/',project,'/filtered_feature_bc_matrix'))
  object <- CreateSeuratObject(matrix$`Gene Expression`, project = project, min.cells = 3, 
                               min.features = 200)
  cell.number.before <- ncol(object)
  object[['percent.mt']] <- PercentageFeatureSet(object ,pattern = '^MT-')
  object <- object[,object$nFeature_RNA > th.nfeature.min &
                     object$nFeature_RNA < th.nfeature.max & 
                     object$percent.mt < th.mt.pct &
                     object$nCount_RNA > th.count.min]
  cell.number.after <- ncol(object)
  saveRDS(object, paste0(output_dir, project, '.object.rds'))
  return(c(cell.number.before, cell.number.after))
}
# Operation
for (i in seq(1, length(work.list))) {
  cell.numbers.before.and.after <- DoQCAndStore(work.list[[i]],
                                                th.nfeature.max = th.nfeature.maxs[[i]],
                                                th.nfeature.min = th.nfeature.mins[[i]],
                                                th.mt.pct = th.mt.pcts[[i]],
                                                th.count.min=th.count.mins[[i]])
  cell.numbers.before <- c(cell.numbers.before, cell.numbers.before.and.after[[1]])
  cell.numbers.after <- c(cell.numbers.after, cell.numbers.before.and.after[[2]])
}

write.csv(data.frame('Donor'=c('Donor6', 'Donor 23'),
                     'Before'=cell.numbers.before,
                     'After'=cell.numbers.after),
          paste0(output_dir,'Cell_numbers.csv'))



# Clustering analysis for Donor 6 -----------------------------------------
d6.object <- readRDS('02_Clustering/D6.object.rds')
d6.object <- NormalizeData(d6.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
sce <- as.SingleCellExperiment(d6.object)
sce <- scDblFinder(sce)
d6.object$doublet_class <- sce$scDblFinder.class
d6.object$doublet_score <- sce$scDblFinder.score
# Filter out doublets
d6.object <- d6.object[,d6.object$doublet_class!='doublet']
DefaultAssay(d6.object) <- 'RNA'
d6.object <- NormalizeData(d6.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d6.object <- RunUMAP(d6.object , verbose = F,dims = 1:10)
d6.object <- FindNeighbors(d6.object)
d6.object <- FindClusters(d6.object ,resolution = 0.5)
rasterize(DimPlot(d6.object, label = T, pt.size = .1), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D6_RNA.pdf'), width = 6, height = 4)
DefaultAssay(d6.object) <- 'RNA'
markers.d6.rna <- FindAllMarkers(d6.object)

# Clustering analysis for Donor 23 -----------------------------------------
d23.object <- readRDS('02_Clustering/D23.object.rds')
d23.object <- NormalizeData(d23.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
sce <- as.SingleCellExperiment(d23.object)
sce <- scDblFinder(sce)
d23.object$doublet_class <- sce$scDblFinder.class
d23.object$doublet_score <- sce$scDblFinder.score
# Filter out doublets
d23.object <- d23.object[,d23.object$doublet_class!='doublet']
d6.object <- NormalizeData(d6.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d23.object <- RunUMAP(d23.object , verbose = F,dims = 1:10)
d23.object <- FindNeighbors(d23.object)
d23.object <- FindClusters(d23.object ,resolution = 0.5)
rasterize(DimPlot(d23.object, label = T, pt.size = .1), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D23_RNA.pdf'), width = 6, height = 4)
DefaultAssay(d23.object) <- 'RNA'
markers.d23.rna <- FindAllMarkers(d23.object)

write.xlsx(list('D6_RNA'=markers.d6.rna,
                'D23_RNA'=markers.d23.rna),
           paste0(output_dir,'Clustering_markers_RNA.xlsx'))


# D6 DC clustering --------------------------------------------------------
d6.object.dc <- d6.object[,Idents(d6.object) %in% c(0,3,6,5,7,10,11,13)]
d6.object.dc <- NormalizeData(d6.object.dc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d6.object.dc <- RunUMAP(d6.object.dc, verbose = F,dims = 1:10)
d6.object.dc <- FindNeighbors(d6.object.dc)
d6.object.dc <- FindClusters(d6.object.dc)
rasterize(DimPlot(d6.object.dc, pt.size = .1, label = T), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D6_DC_RNA.pdf'), width = 6, height = 4)
DefaultAssay(d6.object.dc) <- 'RNA'
# Display markers for a general DC clustering
DotPlot(d6.object.dc, features = c('ITGAX','CLEC9A','XCR1','CADM1',
                                  'CD1C','CLEC10A','CLEC4A',
                                  'AXL','SIGLEC6','CD5',
                                  'CCR7','IL7R','CRLF2',
                                  'RORC','PRDM16',
                                  'IL3RA','CLEC4C',
                                  'CD14','CD163'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave(paste0(output_dir,'Dotplot_DC_D6_unannotated.pdf'),
       width = 10 ,height = 4)
markers.d6.dc.rna <- FindAllMarkers(d6.object.dc)
Idents(d6.object.dc) <- 'seurat_clusters'
d6.object.dc <- RenameIdents(d6.object.dc,
                             '0'='cDC2','1'='DC3','2'='pDC',
                             '3'='cDC2','4'='cDC2','5'='DC2pre-hm',
                             '6'='cDC2','7'='CD11c-hi AS DC','8'='CD11c-lo AS DC',
                             '9'='pDC','10'='cDC1',
                             '11'='DC2hm')
rasterize(DimPlot(d6.object.dc, pt.size = .1, label = T), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D6_DC_RNA_Annotated.pdf'), width = 6, height = 4)
DefaultAssay(d6.object.dc) <- 'RNA'
DotPlot(d6.object.dc, features = c('ITGAX','FLT3',
                                   'CLEC9A','XCR1','CADM1',
                                   'CD1C','CLEC10A','CLEC4A',
                                   'AXL','SIGLEC6','CD5',
                                   'CCR7','IL7R','CRLF2',
                                   'RORC','PRDM16','AIRE',
                                   'IL3RA','CLEC4C',
                                   'CD14','CD163'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave(paste0(output_dir,'Dotplot_DC_D6_annotated.pdf'),
       width = 10 ,height = 4)
# D23 DC clustering -------------------------------------------------------
d23.object.dc <- d23.object[,Idents(d23.object) %in% c(0,1,3,4,5,9,2,11,12)]
d23.object.dc <- NormalizeData(d23.object.dc) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d23.object.dc <- RunUMAP(d23.object.dc, verbose = F,dims = 1:10)
d23.object.dc <- FindNeighbors(d23.object.dc)
d23.object.dc <- FindClusters(d23.object.dc)
rasterize(DimPlot(d23.object.dc, pt.size = .1, label = T), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D23_DC_RNA.pdf'), width = 6, height = 4)
DefaultAssay(d23.object.dc) <- 'RNA'
DotPlot(d23.object.dc, features = c('ITGAX','FLT3',
                                    'CLEC9A','XCR1','CADM1',
                                   'CD1C','CLEC10A','CLEC4A',
                                   'AXL','SIGLEC6','CD5',
                                   'CCR7','IL7R','CRLF2',
                                   'RORC','PRDM16',
                                   'IL3RA','CLEC4C',
                                   'CD14','CD163','MKI67'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave(paste0(output_dir,'Dotplot_DC_D23_unannotated.pdf'),
       width = 10 ,height = 4)
#markers.d23.dc.rna <- FindAllMarkers(d23.object.dc)
Idents(d23.object.dc) <- 'seurat_clusters'
d23.object.dc <- RenameIdents(d23.object.dc,
                             '0'='cDC2','1'='cDC2','2'='cDC2',
                             '3'='pDC','4'='CD11c-hi AS DC','5'='DC3',
                             '6'='DC3','7'='cDC2','8'='CD11c-lo AS DC',
                             '9'='DC2hm','10'='pDC','11'='cDC1',
                             '12'='RORgt+ APC','13'='CD11c-hi AS DC',
                             '14'='Proliferating')
rasterize(DimPlot(d23.object.dc, pt.size = .1, label = T), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D23_DC_RNA_Annotated.pdf'), width = 6, height = 4)
DefaultAssay(d23.object.dc) <- 'RNA'
DotPlot(d23.object.dc, features = c('ITGAX','FLT3',
                                    'CLEC9A','XCR1','CADM1',
                                   'CD1C','CLEC10A','CLEC4A',
                                   'AXL','SIGLEC6','CD5',
                                   'CCR7','IL7R','CRLF2',
                                   'RORC','PRDM16','AIRE',
                                   'IL3RA','CLEC4C',
                                   'CD14','CD163','MKI67'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave(paste0(output_dir,'Dotplot_DC_D23_annotated.pdf'),
       width = 10 ,height = 4)
saveRDS(d23.object, paste0(output_dir,'d23.object.all.rds'))
saveRDS(d6.object, paste0(output_dir,'d6.object.all.rds'))

# Normalize ADT assays for Donor 6 ----------------------------------------------------
# Standard dsb workgflow  
# Define drop outs
raw = Seurat::Read10X("01_CellRanger/D6_protein/raw_feature_bc_matrix/")
cells = Seurat::Read10X("01_CellRanger/D6_protein/sample_filtered_feature_bc_matrix/")
prot = raw$`Antibody Capture`
rna = raw$`Gene Expression`
# define cell-containing barcodes and separate cells and empty drops
stained_cells = colnames(cells$`Gene Expression`)
background = setdiff(colnames(raw$`Gene Expression`), stained_cells)
mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below
md = data.frame(
  rna.size = log10(Matrix::colSums(rna)), 
  prot.size = log10(Matrix::colSums(prot)), 
  n.gene = Matrix::colSums(rna > 0), 
  mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
)
# add indicator for barcodes Cell Ranger called as cells
md$drop.class = ifelse(rownames(md) %in% stained_cells, 'cell', 'background')

# remove barcodes with no evidence of capture in the experiment
md = md[md$rna.size > 0 & md$prot.size > 0, ]
background_drops <- rownames(
  md[ md$prot.size > 1.5 & 
        md$prot.size < 3 & 
        md$rna.size < 2.5, ]
) 
background.adt.mtx = as.matrix(prot[ , background_drops])

cellmd = md[md$drop.class == 'cell', ]

# filter drops with + / - 3 median absolute deviations from the median library size
rna.mult = (3*mad(cellmd$rna.size))
prot.mult = (3*mad(cellmd$prot.size))
rna.lower = median(cellmd$rna.size) - rna.mult
rna.upper = median(cellmd$rna.size) + rna.mult
prot.lower = median(cellmd$prot.size) - prot.mult
prot.upper = median(cellmd$prot.size) + prot.mult

# filter rows based on droplet qualty control metrics
qc_cells = rownames(
  cellmd[cellmd$prot.size > prot.lower & 
           cellmd$prot.size < prot.upper & 
           cellmd$rna.size > rna.lower & 
           cellmd$rna.size < rna.upper & 
           cellmd$mt.prop < 0.14, ]
)
cell.adt.raw = as.matrix(prot[ , qc_cells])
cell.rna.raw = rna[ ,qc_cells]
isotype.ab <- c("Isotype_MOPC.21",  "Isotype_MOPC.173",
                "Isotype_MPC.11",   "Isotype_RTK4530",
                "Isotype_RTK2071", "Isotype_RTK2758")
cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = cell.adt.raw, 
  empty_drop_matrix = background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.ab,
  scale.factor = 'standardize',
  quantile.clipping = TRUE,
  quantile.clip = c(0.01, 0.99)
)
# split the data into separate matrices for RNA and ADT
d6.object.dc.adt.filtered <- d6.object.dc[,intersect(colnames(d6.object.dc),
                                                     colnames(cells.dsb.norm))]
d6.object.dc.adt.filtered[['ADTNormalized']] <- CreateAssayObject(cells.dsb.norm[,
                                                                          intersect(colnames(d6.object.dc),
                                                                                    colnames(cells.dsb.norm))])
DefaultAssay(d6.object.dc.adt.filtered) <- 'ADTNormalized'

pdf(paste0(output_dir, 'All_feature_ADT_D6_DC_new.pdf'), width = 8, height = 6)
DefaultAssay(d6.object.dc.adt.filtered) <- 'ADTNormalized'
for (i in rownames(d6.object.dc.adt.filtered)) {
  plot(rasterize(FeaturePlot(d6.object.dc.adt.filtered, pt.size = .1, features = i), 
                 dpi = 300)+
         scale_color_viridis_c())
}
dev.off()
rasterize(DimPlot(d6.object.dc.adt.filtered, pt.size = .1, label = T), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D6_DC_RNA_Annotated_ADT_filtered.pdf'),
       width = 6, height = 4)
DefaultAssay(d6.object.dc.adt.filtered) <- 'RNA'
DotPlot(d6.object.dc.adt.filtered, features = c('ITGAX','FLT3','HLA-DRA',
                                                 'CLEC9A','XCR1','CADM1','BATF3',
                                                 
                                                'IRF8','CD1C','CLEC10A','CLEC4A',
                                                 'IRF4','AXL','SIGLEC6','CD5',
                                                 'CCR7','IL7R','CRLF2','CD83',
                                                'LAMP3','CD40','CD274','CD200',
                                                'ICOSLG','CCL22',
                                                 'IL3RA','CLEC4C','NRP1','TCF4',
                                                 'CD14','CD163','PLAUR'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  scale_color_viridis_c()
ggsave(paste0(output_dir,'Dotplot_DC_D6_annotated_ADT_filtered.pdf'),
       width = 10 ,height = 4)
DefaultAssay(d6.object.dc.adt.filtered) <- 'ADTNormalized'
markers.d6.dc.adtnorm <- FindAllMarkers(d6.object.dc.adt.filtered)
DefaultAssay(d6.object.dc.adt.filtered) <- 'RNA'
markers.d6.dc.rna.annotated <- FindAllMarkers(d6.object.dc.adt.filtered)

# Normalize ADT assays for Donor 23 ----------------------------------------------------------------
raw = Seurat::Read10X("01_CellRanger/D23_protein/raw_feature_bc_matrix/")
cells = Seurat::Read10X("01_CellRanger/D23_protein/sample_filtered_feature_bc_matrix/")
prot = raw$`Antibody Capture`
rna = raw$`Gene Expression`
# define cell-containing barcodes and separate cells and empty drops
stained_cells = colnames(cells$`Gene Expression`)
background = setdiff(colnames(raw$`Gene Expression`), stained_cells)
mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below

md = data.frame(
  rna.size = log10(Matrix::colSums(rna)), 
  prot.size = log10(Matrix::colSums(prot)), 
  n.gene = Matrix::colSums(rna > 0), 
  mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
)
# add indicator for barcodes Cell Ranger called as cells
md$drop.class = ifelse(rownames(md) %in% stained_cells, 'cell', 'background')

# remove barcodes with no evidence of capture in the experiment
md = md[md$rna.size > 0 & md$prot.size > 0, ]
background_drops <- rownames(
  md[ md$prot.size > 1.5 & 
        md$prot.size < 3 & 
        md$rna.size < 2.5, ]
) 
background.adt.mtx = as.matrix(prot[ , background_drops])

cellmd = md[md$drop.class == 'cell', ]

# filter drops with + / - 3 median absolute deviations from the median library size
rna.mult = (3*mad(cellmd$rna.size))
prot.mult = (3*mad(cellmd$prot.size))
rna.lower = median(cellmd$rna.size) - rna.mult
rna.upper = median(cellmd$rna.size) + rna.mult
prot.lower = median(cellmd$prot.size) - prot.mult
prot.upper = median(cellmd$prot.size) + prot.mult

# filter rows based on droplet qualty control metrics
qc_cells = rownames(
  cellmd[cellmd$prot.size > prot.lower & 
           cellmd$prot.size < prot.upper & 
           cellmd$rna.size > rna.lower & 
           cellmd$rna.size < rna.upper & 
           cellmd$mt.prop < 0.14, ]
)
cell.adt.raw = as.matrix(prot[ , qc_cells])
cell.rna.raw = rna[ ,qc_cells]
isotype.ab <- c("Isotype_MOPC.21",  "Isotype_MOPC.173",
                "Isotype_MPC.11",   "Isotype_RTK4530",
                "Isotype_RTK2071", "Isotype_RTK2758")
cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = cell.adt.raw, 
  empty_drop_matrix = background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.ab,
  scale.factor = 'standardize',
  quantile.clipping = TRUE,
  quantile.clip = c(0.01, 0.99)
)
# split the data into separate matrices for RNA and ADT
d23.object.dc.adt.filtered <- d23.object.dc[,intersect(colnames(d23.object.dc),
                                                     colnames(cells.dsb.norm))]
d23.object.dc.adt.filtered[['ADTNormalized']] <- CreateAssayObject(cells.dsb.norm[,
                                                                          intersect(colnames(d23.object.dc),
                                                                                    colnames(cells.dsb.norm))])
DefaultAssay(d23.object.dc.adt.filtered) <- 'ADTNormalized'

pdf(paste0(output_dir, 'All_feature_ADT_D23_DC_new.pdf'), width = 8, height = 6)
DefaultAssay(d23.object.dc.adt.filtered) <- 'ADTNormalized'
for (i in rownames(d23.object.dc.adt.filtered)) {
  plot(rasterize(FeaturePlot(d23.object.dc.adt.filtered, pt.size = .1, features = i), 
                 dpi = 300)+
         scale_color_viridis_c())
}
dev.off()
rasterize(DimPlot(d23.object.dc.adt.filtered, pt.size = .1, label = T), dpi = 300)
ggsave(paste0(output_dir,'UMAP_D23_DC_RNA_Annotated_ADT_filtered.pdf'),
       width = 6, height = 4)
DefaultAssay(d23.object.dc.adt.filtered) <- 'RNA'
DotPlot(d23.object.dc.adt.filtered, features = c('ITGAX','FLT3',
                                    'CLEC9A','XCR1','CADM1',
                                    'CD1C','CLEC10A','CLEC4A',
                                    'AXL','SIGLEC6','CD5',
                                    'CCR7','IL7R','CRLF2',
                                    'RORC','PRDM16','AIRE',
                                    'IL3RA','CLEC4C',
                                    'CD14','CD163','MKI67'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  scale_color_viridis_c()
ggsave(paste0(output_dir,'Dotplot_DC_D23_annotated_ADT_filtered.pdf'),
       width = 10 ,height = 4)
DefaultAssay(d23.object.dc.adt.filtered) <- 'ADTNormalized'
markers.d23.dc.adtnorm <- FindAllMarkers(d23.object.dc.adt.filtered)
DefaultAssay(d23.object.dc.adt.filtered) <- 'RNA'
markers.d23.dc.rna.annotated <- FindAllMarkers(d23.object.dc.adt.filtered)
DefaultAssay(d23.object.dc.adt.filtered) <- 'RNA'
FeaturePlot(d23.object.dc.adt.filtered,
            features = c('CD14','CD163'), pt.size = .1)
ggsave(paste0(output_dir,'Featureplot_CD14_CD163.pdf'),
       width = 8, height = 4)
FeaturePlot(d6.object.dc.adt.filtered,
            features = c('CD14','CD163'), pt.size = .1)
ggsave(paste0(output_dir,'D6_Featureplot_CD14_CD163.pdf'),
       width = 8, height = 4)
write.xlsx(list('Donor 23 RNA'=markers.d23.dc.rna.annotated,
                'Donor 23 ADTnorm'=markers.d23.dc.adtnorm,
                'Donor 6 RNA'=markers.d6.dc.rna.annotated,
                'Donor 6 ADTnorm'=markers.d6.dc.adtnorm),
           paste0(output_dir,'DC_RNA&ADT_markers_annotated_clusters.xlsx'),
           rowNames=T)


# Display All DC1/2 marker in all Donor 6 DCs
dc.marker <- c("FITC","Hu.CD370","Hu.CD141",
               "Hu.CD1c","Hu.CD172a", "Hu.CD301")
DefaultAssay(d6.object.dc.adt.filtered) <- 'ADTNormalized'
expression_data_compare <- as.data.frame(t(as.data.frame(GetAssayData(d6.object.dc.adt.filtered[dc.marker, ]))))
# Add the cluster information
expression_data_compare$Cluster <- Idents(d6.object.dc.adt.filtered)
expression_data_compare$Cluster <- factor(expression_data_compare$Cluster)
plot_data_compare <- expression_data_compare %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Marker",
    values_to = "Expression"
  )
plot_data_compare$Marker <- factor(plot_data_compare$Marker,
                                   levels = c("FITC","Hu.CD370","Hu.CD141",
                                              "Hu.CD1c","Hu.CD172a", "Hu.CD301"))
ggplot(plot_data_compare, aes(x = Expression, y = Cluster, 
                              fill = Cluster, color = Cluster)) +
  geom_density_ridges2(alpha = 0.7) +
  facet_wrap(~ Marker, ncol = 3, scales = "free_x") +

  # scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  
  labs(
    title = "Marker Expression Across Clusters in Donor 6",
    x = "Expression Level",
    y = "Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" # Hide legend as y-axis is labeled
  ) 

ggsave(paste0(output_dir,'DC1_2_marker_in_all_DC_expression_ridges_D6.pdf'),
       width = 8, height = 4)
# DC1 and DC2hm only
plot_data_compare <- plot_data_compare[plot_data_compare$Cluster %in% c('cDC1',
                                                                        'DC2hm'),]
ggplot(plot_data_compare, aes(x = Expression, y = Cluster, 
                              fill = Cluster, color = Cluster)) +
  geom_density_ridges2(alpha = 0.7) +
  
  facet_wrap(~ Marker, ncol = 3, scales = "free_x") +
  # scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  
  labs(
    title = "Marker Expression Across Clusters in Donor 6",
    x = "Expression Level",
    y = "Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" # Hide legend as y-axis is labeled
  ) 
ggsave(paste0(output_dir,'DC1_2_marker_in_selected_DC_expression_ridges_D6.pdf'),
       width = 8, height = 4)

# All DC1/2 marker in all Donor 23 DCs
dc.marker <- c("FITC","Hu.CD370","Hu.CD141",
               "Hu.CD1c","Hu.CD172a", "Hu.CD301")
DefaultAssay(d23.object.dc.adt.filtered) <- 'ADTNormalized'
expression_data_compare <- as.data.frame(t(as.data.frame(GetAssayData(d23.object.dc.adt.filtered[dc.marker, ]))))
# Add the cluster information
expression_data_compare$Cluster <- Idents(d23.object.dc.adt.filtered)
# Convert Cluster to a factor with desired order for plotting, if necessary
expression_data_compare$Cluster <- factor(expression_data_compare$Cluster)

plot_data_compare <- expression_data_compare %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Marker",
    values_to = "Expression"
  )
plot_data_compare$Marker <- factor(plot_data_compare$Marker,
                                   levels = c("FITC","Hu.CD370","Hu.CD141",
                                              "Hu.CD1c","Hu.CD172a", "Hu.CD301"))
ggplot(plot_data_compare, aes(x = Expression, y = Cluster, 
                              fill = Cluster, color = Cluster)) +
  geom_density_ridges2(alpha = 0.7) +

  facet_wrap(~ Marker, ncol = 3, scales = "free_x") +
  
  
  # scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  
  labs(
    title = "Marker Expression Across Clusters in Donor 23",
    x = "Expression Level",
    y = "Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" # Hide legend as y-axis is labeled
  ) 
ggsave(paste0(output_dir,'DC1_2_marker_in_all_DC_expression_ridges_D23.pdf'),
       width = 8, height = 4)
  # DC1 and DC2hm only
plot_data_compare <- plot_data_compare[plot_data_compare$Cluster %in% c('cDC1',
                                                                        'DC2hm'),]
ggplot(plot_data_compare, aes(x = Expression, y = Cluster, 
                              fill = Cluster, color = Cluster)) +
  geom_density_ridges2(alpha = 0.7) +
  
  # --- THIS IS THE MISSING LINE ---
  # Add this to create a separate panel for each marker
  facet_wrap(~ Marker, ncol = 3, scales = "free_x") +
  
  # scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  
  labs(
    title = "Marker Expression Across Clusters in Donor 23",
    x = "Expression Level",
    y = "Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" # Hide legend as y-axis is labeled
  ) 
ggsave(paste0(output_dir,'DC1_2_marker_in_selected_DC_expression_ridges_D23.pdf'),
       width = 8, height = 4)


# Save cell-type markers --------------------------------------------------
DefaultAssay(d6.object.dc.adt.filtered) <- 'ADTNormalized'
d6.adt.markers <- FindAllMarkers(d6.object.dc.adt.filtered)
DefaultAssay(d6.object.dc.adt.filtered) <- 'RNA'
d6.rna.markers <- FindAllMarkers(d6.object.dc.adt.filtered)
DefaultAssay(d23.object.dc.adt.filtered) <- 'ADTNormalized'
d23.adt.markers <- FindAllMarkers(d23.object.dc.adt.filtered)
DefaultAssay(d23.object.dc.adt.filtered) <- 'RNA'
d23.rna.markers <- FindAllMarkers(d23.object.dc.adt.filtered)
write.xlsx(list('Donor 6 RNA'=d6.rna.markers,
                'Donor 6 ADT'=d6.adt.markers,
                'Donor 23 RNA'=d23.rna.markers,
                'Donor 23 ADT'=d23.adt.markers),
           paste0(output_dir,
                  'Supplemental table CITE-seq cluster markers.xlsx'))
# Save RDS ----------------------------------------------------------------
saveRDS(d23.object.dc.adt.filtered, paste0(output_dir,
                                           'd23.object.dc.adt.filtered_new.rds'))
saveRDS(d6.object.dc.adt.filtered, paste0(output_dir,
                                           'd6.object.dc.adt.filtered_new.rds'))

# Monocle 3 Donor 23 ---------------------------------------------------------------
  #  Donor 23 all cells
counts_matrix <- GetAssayData(
  d23.object.dc.adt.filtered,
  assay = "RNA",
  slot = "counts"
)
d23.object.dc.adt.filtered$celltype <- Idents(d23.object.dc.adt.filtered)
cell_meta <- d23.object.dc.adt.filtered@meta.data
gene_meta <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix)
)
message("Creating new cell_data_set (CDS)...")
cds <- new_cell_data_set(
  expression_data = counts_matrix,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
cds <- preprocess_cds(cds, 
                      num_dim = 100)
cds <- cluster_cells(cds, reduction_method = "PCA")
# Use Seurat UMAP coordinates
seurat_umap_coords <- Embeddings(d23.object.dc.adt.filtered, reduction = "umap")

reducedDims(cds)[["UMAP"]] <- seurat_umap_coords

cds <- cluster_cells(cds, 
                     reduction_method = "UMAP")
cds <- learn_graph(cds,use_partition = T,close_loop = F)
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "celltype", 
                     label_cell_groups = TRUE,label_roots = F,label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10)+
            , dpi = 300)
ggsave('02_Clustering/Monocle3_inferred trajectory_D23.pdf',
       width = 4, height = 3)

cyto_results <- CytoTRACE(as.matrix(counts_matrix), ncores = 6,
                          enableFast = F)
d23.object.dc.adt.filtered$CytoTRACE_Score <- cyto_results$CytoTRACE
d23.object.dc.adt.filtered$CytoTRACE_Rank <- cyto_results$CytoTRACErank
FeaturePlot(d23.object.dc.adt.filtered, 
            features = "CytoTRACE_Score", 
            reduction = "umap",pt.size = .1) +
  ggtitle("CytoTRACE Score (Differentiation Potential)") +
  scale_color_viridis_c()
ggsave(paste0(output_dir,'CytoTRACE_Donor23_all_cell_feature.pdf'),
       width = 4.5, height = 3)
cds@colData$celltype <- factor(cds@colData$celltype,
                               levels = levels(d23.object.dc.adt.filtered))
VlnPlot(d23.object.dc.adt.filtered, 
        features = "CytoTRACE_Score", 
        pt.size = 0) +
  ggtitle("CytoTRACE Score by Cluster") +
  NoLegend()
ggsave(paste0(output_dir,'CytoTRACE_Donor23_all_cell_vlnplot.pdf'),
       width = 4.5, height = 3)
cds <- order_cells(cds)
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "pseudotime",
                     label_cell_groups = TRUE,label_roots = F,label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10), dpi = 300)
ggsave('02_Clustering/Monocle3_inferred pseudotime_D23.pdf',
       width = 4, height = 3)
saveRDS(cds, paste0(output_dir, 'Donor23_all_cell_monocle3_cds.rds'))
  # cDC2 trajectory
DefaultAssay(d23.object.dc.adt.filtered) <- 'RNA'
d23.object.dc.adt.filtered.cdc2 <- d23.object.dc.adt.filtered[,d23.object.dc.adt.filtered$celltype %in% c('cDC2',
                                                                                                      'DC2hm')]
d23.object.dc.adt.filtered.cdc2 <- NormalizeData(d23.object.dc.adt.filtered.cdc2)   %>% 
  FindVariableFeatures() %>%
   ScaleData() %>%
   RunPCA()
d23.object.dc.adt.filtered.cdc2 <- RunUMAP(d23.object.dc.adt.filtered.cdc2, 
                                          dims = 1:10)
# Donor 23 cDC2 monocle3 trajectory, similar workflow

counts_matrix <- GetAssayData(
  d23.object.dc.adt.filtered.cdc2,
  assay = "RNA",
  slot = "counts"
)
cell_meta <- d23.object.dc.adt.filtered.cdc2@meta.data

gene_meta <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix)
)
cds <- new_cell_data_set(
  expression_data = counts_matrix,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
cds <- preprocess_cds(cds, num_dim = 100)
message("Running Monocle 3 cluster_cells (using PCA)...")
cds <- cluster_cells(cds, 
                     reduction_method = "PCA")
seurat_umap_coords <- Embeddings(d23.object.dc.adt.filtered.cdc2,
                                 reduction = "umap")

reducedDims(cds)[["UMAP"]] <- seurat_umap_coords

cds <- cluster_cells(cds, 
                     reduction_method = "UMAP")
cds <- learn_graph(cds,use_partition = F,
                   close_loop = F,
                   learn_graph_control = list(minimal_branch_len=5))
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "celltype", 
                     label_cell_groups = TRUE,
                     label_roots = F,
                     label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10), dpi = 300)
ggsave('02_Clustering/Monocle3_inferred_trajectory_cDC2_D23.pdf',
       width = 4, height = 3)
FeaturePlot(d23.object.dc.adt.filtered.cdc2, 
            features = "CytoTRACE_Score", 
            reduction = "umap",pt.size = .1) +
  ggtitle("CytoTRACE Score (Differentiation Potential)") +
  scale_color_viridis_c()
ggsave(paste0(output_dir,'CytoTRACE_Donor23_cDC2_cell_feature.pdf'),
       width = 4.5, height = 3)
VlnPlot(d23.object.dc.adt.filtered.cdc2, 
        features = "CytoTRACE_Score", 
        pt.size = 0) +
  ggtitle("CytoTRACE Score by Cluster") +
  NoLegend()
ggsave(paste0(output_dir,'CytoTRACE_Donor23_cDC2_vlnplot.pdf'),
       width = 4.5, height = 3)
cds <- order_cells(cds)
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "pseudotime", 
                     label_cell_groups = TRUE,label_roots = F,label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10), dpi = 300)
ggsave('02_Clustering/Monocle3_inferred pseudotime_cDC2_D23.pdf',
       width = 4, height = 3)
saveRDS(cds, paste0(output_dir, 'Donor23_cDC2_monocle3_cds.rds'))
mono_pseudotime <- pseudotime(cds)

pseudotime_df <- data.frame(
  cell_name = names(mono_pseudotime),
  pseudotime = mono_pseudotime
)

d23.object.dc.adt.filtered.cdc2$pseudotime <- pseudotime_df[,2]

metadata_df <- d23.object.dc.adt.filtered.cdc2@meta.data
d23.object.dc.adt.filtered.cdc2 <- ScaleData(d23.object.dc.adt.filtered.cdc2,
                                            features = rownames(d23.object.dc.adt.filtered.cdc2))

genes_to_plot <- read.csv('../Code deposit/scRNA-seq analysis/Highlighted_gene_selected.csv')$Gene[4:31]
full_expression_matrix <- GetAssayData(d23.object.dc.adt.filtered.cdc2, 
                                       assay = "RNA", 
                                       slot = "data")

pdf(paste0(output_dir,
           "Gene_Expression_vs_Pseudotime_Donor23_RNA.pdf"), 
    width = 7, height = 5)


# Showing per-gene expression in pseudotime
for (gene in genes_to_plot) {
  
  
  if (!gene %in% rownames(full_expression_matrix)) {
    message(paste("Warning: Gene", gene, "not found. Skipping."))
    next # Skip to the next gene
  }
  
  
  plot_df <- data.frame(
    expression = full_expression_matrix[gene,],
    pseudotime = metadata_df$pseudotime,
    cluster = metadata_df$celltype
  )
  
  
  p <- ggplot(plot_df, aes(x = pseudotime, y = expression)) +
    geom_point(aes(color = cluster), size = 0.5, alpha = 0.4) + 
    geom_smooth(method = "gam", color = "black", se = FALSE) + 
    ggtitle(paste(gene, "Expression Over Pseudotime")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Pseudotime", y = "Log-Normalized Expression")
  
  
  # 'print()' is necessary inside a loop to render the plot
  print(p)
  
  message(paste("Plotted", gene))
}
dev.off()

# ADT plot

genes_to_plot <- rownames(d23.object.dc.adt.filtered.cdc2@assays$ADTNormalized)
full_expression_matrix <- GetAssayData(d23.object.dc.adt.filtered.cdc2, 
                                       assay = "ADTNormalized", 
                                       slot = "data")

pdf(paste0(output_dir,
           "Gene_Expression_vs_Pseudotime_Donor23_ADT.pdf"), 
    width = 7, height = 5)

message("Starting PDF generation...")


for (gene in genes_to_plot) {
  
  # Check if the gene is in our expression matrix
  if (!gene %in% rownames(full_expression_matrix)) {
    message(paste("Warning: Gene", gene, "not found. Skipping."))
    next # Skip to the next gene
  }
  
  
  plot_df <- data.frame(
    expression = full_expression_matrix[gene,],
    pseudotime = metadata_df$pseudotime,
    cluster = metadata_df$celltype
  )
  
  
  p <- ggplot(plot_df, aes(x = pseudotime, y = expression)) +
    geom_point(aes(color = cluster), size = 0.5, alpha = 0.4) + 
    geom_smooth(method = "gam", color = "black", se = FALSE) + 
    ggtitle(paste(gene, "Expression Over Pseudotime")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Pseudotime", y = "Log-Normalized Expression")
  
  
  # 'print()' is necessary inside a loop to render the plot
  print(p)
  
  message(paste("Plotted", gene))
}
dev.off()



# Monocle 3 Donor 6 ---------------------------------------------------------------
#  Donor 6 all cells, similar workflow
counts_matrix <- GetAssayData(
  d6.object.dc.adt.filtered,
  assay = "RNA",
  slot = "counts"
)
d6.object.dc.adt.filtered$celltype <- Idents(d6.object.dc.adt.filtered)

cell_meta <- d6.object.dc.adt.filtered@meta.data

gene_meta <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix)
)

cds <- new_cell_data_set(
  expression_data = counts_matrix,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)

cds <- preprocess_cds(cds, num_dim = 100)


cds <- cluster_cells(cds, reduction_method = "PCA")

seurat_umap_coords <- Embeddings(d6.object.dc.adt.filtered, reduction = "umap")

reducedDims(cds)[["UMAP"]] <- seurat_umap_coords

cds <- cluster_cells(cds, 
                     reduction_method = "UMAP")
cds <- learn_graph(cds,use_partition = T,close_loop = F)
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "celltype", 
                     label_cell_groups = TRUE,label_roots = F,label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10, cell_size = 0.5,
                      alpha = 1)+
            scale_color_manual(values = c('cDC1'='#C5E524',
                                          'pDC'='#8965A9',
                                          'cDC2'='#38C2E5',
                                          'DC2pre-hm'='#384C94',
                                          'DC2hm'='#E679C5',
                                          'CD11c-hi AS DC'="#E62A8A",
                                          'CD11c-lo AS DC'='#EC6968',
                                          'DC3'='#F0A42F',
                                          'pDC'='#9966CC')), dpi = 300)
ggsave('02_Clustering/Monocle3_inferred trajectory.pdf',
       width = 4, height = 3)

rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "pseudotime", 
                     label_cell_groups = TRUE,label_roots = F,label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10,
                     cell_size = 0.5,
                     alpha = 1), dpi = 300)
ggsave('02_Clustering/Monocle3_inferred pseudotime.pdf', width = 4.2, height = 3)
cyto_results <- CytoTRACE(as.matrix(counts_matrix), ncores = 6,
                          enableFast = F)
d6.object.dc.adt.filtered$CytoTRACE_Score <- cyto_results$CytoTRACE
d6.object.dc.adt.filtered$CytoTRACE_Rank <- cyto_results$CytoTRACErank
FeaturePlot(d6.object.dc.adt.filtered, 
            features = "CytoTRACE_Score", 
            reduction = "umap",pt.size = .1) +
  ggtitle("CytoTRACE Score (Differentiation Potential)") +
  scale_color_viridis_c()
ggsave(paste0(output_dir,'CytoTRACE_Donor6_all_cell_feature.pdf'),
       width = 4.5, height = 3)
cds@colData$celltype <- factor(cds@colData$celltype,
                               levels = levels(d6.object.dc.adt.filtered))
VlnPlot(d6.object.dc.adt.filtered, 
        features = "CytoTRACE_Score", 
        pt.size = 0) +
  ggtitle("CytoTRACE Score by Cluster") +
  NoLegend()
ggsave(paste0(output_dir,'CytoTRACE_Donor6_all_cell_vlnplot.pdf'),
       width = 4.5, height = 3)
cds <- order_cells(cds)
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "pseudotime", 
                     label_cell_groups = TRUE,label_roots = F,label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10),dpi = 300)
ggsave('02_Clustering/Monocle3_inferred pseudotime.pdf',
       width = 4, height = 3)
saveRDS(cds, paste0(output_dir, 'Donor6_all_cell_monocle3_cds.rds'))
# cDC2 trajectory
d6.object.dc.adt.filtered.cdc2 <- d6.object.dc.adt.filtered[,d6.object.dc.adt.filtered$celltype %in% c('cDC2',
                                                                                                       'DC2pre-hm',
                                                                                                       'DC2hm')]
d6.object.dc.adt.filtered.cdc2 <- NormalizeData(d6.object.dc.adt.filtered.cdc2)   %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
d6.object.dc.adt.filtered.cdc2 <- RunUMAP(d6.object.dc.adt.filtered.cdc2, 
                                          dims = 1:5)
# Donor 6 cDC2 monocle3 trajectory

counts_matrix <- GetAssayData(
  d6.object.dc.adt.filtered.cdc2,
  assay = "RNA",
  slot = "counts"
)

cell_meta <- d6.object.dc.adt.filtered.cdc2@meta.data

gene_meta <- data.frame(
  gene_short_name = rownames(counts_matrix),
  row.names = rownames(counts_matrix)
)

cds <- new_cell_data_set(
  expression_data = counts_matrix,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
cds <- preprocess_cds(cds, num_dim = 100)

cds <- cluster_cells(cds, 
                     reduction_method = "PCA")
seurat_umap_coords <- Embeddings(d6.object.dc.adt.filtered.cdc2,
                                 reduction = "umap")

reducedDims(cds)[["UMAP"]] <- seurat_umap_coords

cds <- cluster_cells(cds, 
                     reduction_method = "UMAP")
cds <- learn_graph(cds,use_partition = F,
                   close_loop = F)
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "celltype", 
                     label_cell_groups = TRUE,
                     label_roots = F,
                     label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10,
                     cell_size = 0.5,
                     alpha = 1)+
            scale_color_manual(values = c('cDC1'='#C5E524',
                                          'pDC'='#8965A9',
                                          'cDC2'='#38C2E5',
                                          'DC2pre-hm'='#384C94',
                                          'DC2hm'='#E679C5',
                                          'CD11c-hi AS DC'="#E62A8A",
                                          'CD11c-lo AS DC'='#EC6968',
                                          'DC3'='#F0A42F',
                                          'pDC'='#9966CC')), dpi = 300)
ggsave('02_Clustering/Monocle3_inferred_trajectory_cDC2.pdf',
       width = 4, height = 3)
FeaturePlot(d6.object.dc.adt.filtered.cdc2, 
            features = "CytoTRACE_Score", 
            reduction = "umap",pt.size = .1) +
  ggtitle("CytoTRACE Score (Differentiation Potential)") +
  scale_color_viridis_c()
ggsave(paste0(output_dir,'CytoTRACE_Donor6_cDC2_cell_feature.pdf'),
       width = 4.5, height = 3)
VlnPlot(d6.object.dc.adt.filtered.cdc2, 
        features = "CytoTRACE_Score", 
        pt.size = 0) +
  ggtitle("CytoTRACE Score by Cluster") +
  NoLegend()
ggsave(paste0(output_dir,'CytoTRACE_Donor6_cDC2_vlnplot.pdf'),
       width = 4.5, height = 3)
cds <- order_cells(cds)
rasterize(plot_cells(cds,
                     reduction_method = "UMAP",
                     color_cells_by = "pseudotime", 
                     label_cell_groups = TRUE,label_roots = F,label_branch_points = F,
                     label_leaves = F,
                     #label_branch_points = TRUE,
                     graph_label_size = 10,
                     cell_size = 0.5,
                     alpha = 1), dpi = 300)
ggsave('02_Clustering/Monocle3_inferred pseudotime_cDC2.pdf',
       width = 4.2, height = 3)
saveRDS(cds, paste0(output_dir, 'Donor6_cDC2_monocle3_cds.rds'))

mono_pseudotime <- pseudotime(cds)

pseudotime_df <- data.frame(
  cell_name = names(mono_pseudotime),
  pseudotime = mono_pseudotime
)

d6.object.dc.adt.filtered.cdc2$pseudotime <- pseudotime_df[,2]

metadata_df <- d6.object.dc.adt.filtered.cdc2@meta.data
d6.object.dc.adt.filtered.cdc2 <- ScaleData(d6.object.dc.adt.filtered.cdc2,
                                            features = rownames(d6.object.dc.adt.filtered.cdc2))

genes_to_plot <- read.csv('Highlighted_gene_selected.csv')$Gene[4:31]
full_expression_matrix <- GetAssayData(d6.object.dc.adt.filtered.cdc2, 
                                       assay = "RNA", 
                                       slot = "data")

pdf(paste0(output_dir,
           "Gene_Expression_vs_Pseudotime_Donor6_RNA.pdf"), 
    width = 6, height = 4)

message("Starting PDF generation...")

# 3.2 --- Loop through each gene ---
for (gene in c(genes_to_plot,'CD83')) {
  
  # Check if the gene is in our expression matrix
  if (!gene %in% rownames(full_expression_matrix)) {
    message(paste("Warning: Gene", gene, "not found. Skipping."))
    next # Skip to the next gene
  }
  
  # 3.3. Create the data frame for this one gene
  plot_df <- data.frame(
    expression = full_expression_matrix[gene,],
    pseudotime = metadata_df$pseudotime,
    cluster = metadata_df$celltype
  )
  
  # 3.4. Create the ggplot object
  p <- rasterize(ggplot(plot_df, aes(x = pseudotime, y = expression)) +
                   geom_point(aes(color = cluster), size = 1.5, alpha = 0.8) + 
                   geom_smooth(method = "gam", color = "grey60", se = FALSE) + 
                   ggtitle(paste(gene, "Expression Over Pseudotime")) +
                   theme_classic() +
                   theme(plot.title = element_text(hjust = 0.5)) +
                   labs(x = "Pseudotime", y = "Log-Normalized Expression")+
                   scale_color_manual(values = c('cDC1'='#C5E524',
                                                 'pDC'='#8965A9',
                                                 'cDC2'='#38C2E5',
                                                 'DC2pre-hm'='#384C94',
                                                 'DC2hm'='#E679C5',
                                                 'CD11c-hi AS DC'="#E62A8A",
                                                 'CD11c-lo AS DC'='#EC6968',
                                                 'DC3'='#F0A42F',
                                                 'pDC'='#9966CC')),dpi=300)
  
  # 3.5. Print the plot to the PDF
  # 'print()' is necessary inside a loop to render the plot
  print(p)
  
  message(paste("Plotted", gene))
}
dev.off()

# ADT plot

genes_to_plot <- rownames(d6.object.dc.adt.filtered.cdc2@assays$ADTNormalized)
full_expression_matrix <- GetAssayData(d6.object.dc.adt.filtered.cdc2, 
                                       assay = "ADTNormalized", 
                                       slot = "data")

pdf(paste0(output_dir,
           "Gene_Expression_vs_Pseudotime_Donor6_ADT.pdf"), 
    width = 7, height = 5)


# 3.2 --- Loop through each gene ---
for (gene in genes_to_plot) {
  
  if (!gene %in% rownames(full_expression_matrix)) {
    message(paste("Warning: Gene", gene, "not found. Skipping."))
    next # Skip to the next gene
  }
  
  plot_df <- data.frame(
    expression = full_expression_matrix[gene,],
    pseudotime = metadata_df$pseudotime,
    cluster = metadata_df$celltype
  )
  
  p <- ggplot(plot_df, aes(x = pseudotime, y = expression)) +
    geom_point(aes(color = cluster), size = 0.5, alpha = 0.4) + 
    geom_smooth(method = "gam", color = "black", se = FALSE) + 
    ggtitle(paste(gene, "Expression Over Pseudotime")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Pseudotime", y = "Log-Normalized Expression")
  # 'print()' is necessary inside a loop to render the plot
  print(p)
  
  message(paste("Plotted", gene))
}
dev.off()
