library(pacman)
p_load(Seurat, ggplot2, ggplotify, ggpubr, harmony, 
       DoubletFinder, future, dplyr, tidyverse,ggrastr, openxlsx,
       RColorBrewer, ggrastr,CytoTRACE, monocle3,ClusterGVis,
       clusterProfiler, ggsci, ComplexHeatmap, circlize)
set.seed(0822)
setwd("/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/scRNA-seq analysis")
# QC for Cellranger Outputs and create filtered RDS files -----------------------------------------------
work.list <- c("Dnr6","Dnr10","Dnr12")

# Define the function to compute and store initial criteria for QC
PlotQC <- function(project){
  matrix <- Read10X(paste0('CellrangerOutputs/',project,'/filtered_feature_bc_matrix'))
  object <- CreateSeuratObject(matrix, project = project, min.cells = 3, 
                               min.features = 200)
  object[['percent.mt']] <- PercentageFeatureSet(object ,pattern = '^MT-')
  vlnplots <- VlnPlot(object, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
                      combine = F)
  figure <- ggarrange(plotlist = vlnplots, ncol = 3, common.legend = T, legend = 'none')
  qc.vlnplot <- annotate_figure(figure, top = paste0('QC_Violin_Plot_sample_', project))
  scatter.1 <- FeatureScatter(object, feature1 = 'nCount_RNA', 
                              feature2 = 'percent.mt')
  scatter.2 <- FeatureScatter(object, feature1 = 'nCount_RNA', 
                              feature2 = 'nFeature_RNA')
  pdf(paste0('Outputs/figures/',project,'_QC_Plots.pdf'))
  plot(qc.vlnplot)
  plot(scatter.1)
  plot(scatter.2)
  dev.off()
}
for (i in seq(1,length(work.list))) {
  PlotQC(project = work.list[[i]])
}
# Define thresholds for QC
th.nfeature.mins <- rep(500, 3)
th.nfeature.maxs <- rep(7500,3)
th.mt.pcts <- rep(10,3)
th.count.mins <- rep(2000,3)
threshold.df <- data.frame(th.nfeature.mins=th.nfeature.mins,
                           th.nfeature.maxs=th.nfeature.maxs,
                           th.mt.pct=th.mt.pcts,
                           th.count.mins=th.count.mins,
                           row.names = work.list)
write.csv(threshold.df, 'Outputs/tables/thresholds.csv')

# Actually perform QC, store RDS files, and thresholds used for QC
cell.numbers.before <- c()
cell.numbers.after <- c()
  # Define the function
DoQCAndStore <- function(project,th.nfeature.min,
                         th.nfeature.max, th.mt.pct,th.count.min){
  matrix <- Read10X(paste0('CellrangerOutputs/',project,'/filtered_feature_bc_matrix'))
  object <- CreateSeuratObject(matrix, project = project, min.cells = 3, 
                               min.features = 200)
  cell.number.before <- ncol(object)
  object[['percent.mt']] <- PercentageFeatureSet(object ,pattern = '^MT-')
  object <- subset(object, subset = nFeature_RNA > th.nfeature.min &
                     nFeature_RNA <th.nfeature.max & percent.mt < th.mt.pct &
                     nCount_RNA > th.count.min)
  cell.number.after <- ncol(object)
  saveRDS(object, paste0('Outputs/rds/', project, '.object.rds'))
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

write.csv(data.frame('Donor'=c('Donor6', 'Donor 10', 'Donor 12'),
                     'Before'=cell.numbers.before,
                     'After'=cell.numbers.after),
          'Outputs/tables/cell_numbers.csv')


# Clustering analysis on the Dnr6 object ----------------------------------

d6.object <- readRDS('Outputs/rds/Dnr6.object.rds')
d6.object  <- NormalizeData(d6.object ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d6.object  <- RunUMAP(d6.object , verbose = F,dims = 1:10)
d6.object  <- FindNeighbors(d6.object )
d6.object  <- FindClusters(d6.object ,resolution = 1.1)
rasterize(DimPlot(d6.object)+
            ggtitle('Donor 6'),
          dpi = 300)
ggsave('Outputs/figures/UMAP_dnr6_all_cells.pdf',
       width = 4,height = 3.5)
# General annotation
DotPlot(d6.object , features = c('CLEC9A','XCR1','BATF3',
                             'CD1C','ITGAX','CLEC10A',
                             'CCR7','IL7R','LAMP3','IDO1',
                             'CLEC4C','IL3RA','TCF4',
                             'AXL','SIGLEC6',
                             'HSPA6','G0S2','NKG7',
                             'CD14','LYZ','CD68','FCGR3A',
                             'CD19','MS4A1','CD79A',
                             'CD38','CD27',
                             'MKI67','CD34'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Outputs/figures/Donor6_General_markers.pdf', width = 10, height = 5)
d6.markers <- FindAllMarkers(d6.object, only.pos = T)
# Save markers
marker.cluster.list <- list()
for (i in names(table(markers$cluster))) {
  marker.cluster <- markers[markers$cluster==i & markers$p_val_adj<0.05,]
  marker.cluster.list[[i]] <- marker.cluster[order(marker.cluster$avg_log2FC,
                                                   decreasing = T),]
}
write.xlsx(marker.cluster.list, 
           'Outputs/tables/D6_Markers_by_clusters_all_cells.xlsx',
           rowNames=T)
  # Annotate object
d6.object <- RenameIdents(d6.object,
                       c('0'='cDC2','1'='DC2hm','2'='cDC2',
                         '3'='DC2pre-hm','4'='cDC2','5'='cDC2_UPR',
                         '6'='cDC2','7'='AXL DC','8'='DC2hm','9'='Low reads',
                         '10'='cDC1','11'='Mono/Macro',
                         '12'='B','13'='NK', '14'='DC2pre-hm',
                         '15'='Low reads','16'='Mono/Macro',
                         '17'='B proliferating'))
saveRDS(d6.object, 'Outputs/rds/Dnr6_all_cells.rds')
d6.object$celltype_hh_1 <- Idents(d6.object )
d6.object.filtered <- d6.object [,!d6.object $celltype_hh_1 %in% c('Low reads',
                                                        'Doublets')]
  # Figure S2B, Donor 6 UMAP of all cells
rasterize(DimPlot(d6.object.filtered, label = T)+NoLegend()+
            ggtitle('Donor 6'),dpi = 300)
ggsave('Outputs/figures/UMAP_dnr6_annotated.filtered.pdf',
       width = 4,height = 3.5)
# Fine-tune clustering of cDCs

d6.object.cdcs <- d6.object[,Idents(d6.object) %in% c('cDC1','cDC2','DC2pre-hm',
                                             'DC2hm')]
d6.object.cdcs <- NormalizeData(d6.object.cdcs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d6.object.cdcs <- RunUMAP(d6.object.cdcs, dims=1:8)
d6.object.cdcs <- FindNeighbors(d6.object.cdcs) %>% FindClusters(resolution = 0.5)
d6.object.cdcs <- RenameIdents(d6.object.cdcs, c('0'='DC2hm','1'='cDC2','2'='cDC2',
                                           '3'='cDC2','4'='DC2pre-hm',
                                           '5'='cDC1','6'='cDC2',
                                           '7'='DC2pre-hm'))
levels(d6.object.cdcs) <- c('cDC1','cDC2','DC2pre-hm','DC2hm')
d6.object.cdcs$celltype_hh_1 <- Idents(d6.object.cdcs)

  # Figure 1D UMAP of Donor 6 cDCs
rasterize(DimPlot(d6.object.cdcs, pt.size = .1,
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5')),
          dpi = 300)
ggsave('Outputs/figures/Donor_6_UMAP_cDCs.pdf', width = 4.5, height = 3)
  # Figure 1F Dot plot showing marker expression in Donor 6 cDC subpopulations
DotPlot(d6.object.cdcs, 
        features = read.csv('Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Outputs/figures/Dotplot_of_cDC_markers_Donor6.pdf',
       width = 10, height = 2.7)
saveRDS(d6.object.cdcs,
        'Outputs/rds/Dnr6_cDCs.rds')

  # Figure 5A, CRLF2 feature plot
rasterize(FeaturePlot(d6.object.cdcs, features = 'CRLF2', pt.size = .1)+
            scale_color_gradient(low = 'lightgrey', high='#FF2B72'),
          dpi = 300)
ggsave('Outputs/figures/CRLF2_featureplot.pdf', width = 4.5, height = 3)

# Clustering analysis on the Dnr10 object ---------------------------------
d10.object <- readRDS('Outputs/rds/Dnr10.object.rds')
d10.object  <- NormalizeData(d10.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d10.object  <- RunUMAP(d10.object , verbose = F,dims = 1:10)
d10.object  <- FindNeighbors(d10.object )
d10.object  <- FindClusters(d10.object)
rasterize(DimPlot(d10.object)+
            ggtitle('Donor 10'),
          dpi = 300)
ggsave('Outputs/figures/UMAP_dnr10_all_cells.pdf',
       width = 4,height = 3.5)
# General annotation
DotPlot(d10.object, features = c('CLEC9A','XCR1','BATF3',
                             'CD1C','ITGAX','CLEC10A',
                             'CCR7','IL7R','LAMP3','IDO1',
                             'CLEC4C','IL3RA','TCF4',
                             'AXL','SIGLEC6',
                             'HSPA6','G0S2','DNJAB1',
                             'CD14','LYZ','CD68','FCGR3A',
                             'CD19','MS4A1','CD79A',
                             'CD38','CD27',
                             'MKI67','CD34'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Outputs/figures/Donor_10_General_markers_.pdf', width = 10, height = 5)
markers <- FindAllMarkers(d10.object, only.pos = T)
# Save markers
marker.cluster.list <- list()
for (i in names(table(markers$cluster))) {
  marker.cluster <- markers[markers$cluster==i & markers$p_val_adj<0.05,]
  marker.cluster.list[[i]] <- marker.cluster[order(marker.cluster$avg_log2FC,
                                                   decreasing = T),]
}
write.xlsx(marker.cluster.list, 
           'Outputs/tables/D10_Markers_by_clusters_all_cells.xlsx',
           rowNames=T)
# Annotate object
d10.object <- RenameIdents(d10.object,
                       c('0'='cDC2','1'='cDC2_UPR','2'='cDC2',
                         '3'='cDC2_UPR','4'='cDC2','5'='cDC2',
                         '6'='DC2pre-hm','7'='cDC2','8'='cDC2','9'='cDC1',
                         '10'='Proliferating','11'='Low reads',
                         '12'='DC2hm','13'='Mono/Macro', '14'='pDC'))
saveRDS(d10.object, 'Outputs/rds/Dnr10_all_cells.rds')
d10.object$celltype_hh_1 <- Idents(d10.object )
d10.object.filtered <- d10.object [,!d10.object $celltype_hh_1 %in% c('Low reads',
                                                                   'Doublets')]
# Figure S2B, Donor 10 UMAP of all cells
rasterize(DimPlot(d10.object.filtered, label = T)+NoLegend()+
            ggtitle('Donor 10'),dpi = 300)
ggsave('Outputs/figures/UMAP_dnr10_annotated.filtered.pdf',
       width = 4,height = 3.5)
# Fine-tune clustering of cDCs
d10.object.cdcs <- d10.object[,Idents(d10.object) %in% c('cDC1','cDC2',
                                                         'DC2pre-hm',
                                                         'DC2hm')]
d10.object.cdcs <- NormalizeData(d10.object.cdcs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d10.object.cdcs <- RunUMAP(d10.object.cdcs, dims=1:7)
d10.object.cdcs <- FindNeighbors(d10.object.cdcs) %>% FindClusters(resolution = 0.5)
d10.object.cdcs <- RenameIdents(d10.object.cdcs, c('0'='cDC2','1'='cDC2',
                                           '2'='cDC2','3'='cDC2',
                                           '4'='DC2pre-hm',
                                           '5'='cDC1', '6'='DC2hm'))

d10.object.cdcs$celltype_hh_1 <- Idents(d10.object.cdcs)
rasterize(DimPlot(d10.object.cdcs, label = T)+NoLegend()+
            ggtitle('Donor 10'),dpi = 300)
levels(d10.object.cdcs) <- c('cDC1','cDC2','DC2pre-hm', 'DC2hm')
# UMAP of Donor 10 cDCs
rasterize(DimPlot(d10.object.cdcs, pt.size = .1,
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5')),
          dpi = 300)
ggsave('Outputs/figures/Donor_10_UMAP_cDCs.pdf', width = 4.5, height = 3)
# Figure 1F Dot plot showing marker expression in Donor 6 cDC subpopulations
DotPlot(d10.object.cdcs, 
        features = read.csv('Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Outputs/figures/Dotplot_of_cDC_markers_Donor10.pdf',
       width = 10, height = 2.7)
saveRDS(d10.object.cdcs,
        'Outputs/rds/Dnr10_cDCs.rds')
# Clustering analysis on the Dnr12 object ---------------------------------
d12.object <- readRDS('Outputs/rds/Dnr12.object.rds')
d12.object  <- NormalizeData(d12.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d12.object  <- RunUMAP(d12.object , verbose = F,dims = 1:10)
d12.object  <- FindNeighbors(d12.object )
d12.object  <- FindClusters(d12.object)
rasterize(DimPlot(d12.object)+
            ggtitle('Donor 12'),
          dpi = 300)
ggsave('Outputs/figures/UMAP_dnr12_all_cells.pdf',
       width = 4,height = 3.5)
# General annotation
# General annotation
DotPlot(d12.object, features = c('CLEC9A','XCR1','BATF3',
                                 'CD1C','ITGAX','CLEC10A',
                                 'CCR7','IL7R','LAMP3','IDO1',
                                 'CLEC4C','IL3RA','TCF4',
                                 'AXL','SIGLEC6',
                                 'HSPA6','G0S2','DNJAB1',
                                 'CD14','LYZ','CD68','FCGR3A',
                                 'CD19','MS4A1','CD79A',
                                 'CD38','CD27',
                                 'MKI67','CD34'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Outputs/figures/Donor_12_General_markers_.pdf', width = 10, height = 5)
markers <- FindAllMarkers(d12.object, only.pos = T)
# Save markers
marker.cluster.list <- list()
for (i in names(table(markers$cluster))) {
  marker.cluster <- markers[markers$cluster==i & markers$p_val_adj<0.05,]
  marker.cluster.list[[i]] <- marker.cluster[order(marker.cluster$avg_log2FC,
                                                   decreasing = T),]
}
write.xlsx(marker.cluster.list, 
           'Outputs/tables/D12_Markers_by_clusters_all_cells.xlsx',
           rowNames=T)
# Annotate object
d12.object <- RenameIdents(d12.object,
                       c('0'='cDC2','1'='cDC2','2'='cDC2',
                         '3'='cDC2','4'='cDC2','5'='cDC1',
                         '6'='cDC2','7'='cDC2','8'='DC2pre-hm','9'='pDC',
                         '10'='cDC2','11'='DC2hm','12'='cDC2','13'='ILC',
                         '14'='B','15'='Low reads','16'='Doublets','17'='pDC'))
saveRDS(d12.object, 'Outputs/rds/Dnr12_all_cells.rds')
d12.object$celltype_hh_1 <- Idents(d12.object )
d12.object.filtered <- d12.object [,!d12.object $celltype_hh_1 %in% c('Low reads',
                                                                      'ILC')]
# Figure S2B, Donor 12 UMAP of all cells
rasterize(DimPlot(d12.object.filtered, label = T)+NoLegend()+
            ggtitle('Donor 12'),dpi = 300)
ggsave('Outputs/figures/UMAP_dnr12_annotated.filtered.pdf',
       width = 4,height = 3.5)
d12.object.cdcs <- d12.object[,Idents(d12.object) %in% c('cDC1','cDC2',
                                                         'DC2pre-hm',
                                                         'DC2hm')]
d12.object.cdcs <- NormalizeData(d12.object.cdcs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
d12.object.cdcs <- RunUMAP(d12.object.cdcs, dims=1:10)
rasterize(DimPlot(d12.object.cdcs, label = T)+NoLegend()+
            ggtitle('Donor 12'),dpi = 300)
levels(d12.object.cdcs) <- c('cDC1','cDC2','DC2pre-hm', 'DC2hm')
# UMAP of Donor 12 cDCs
rasterize(DimPlot(d12.object.cdcs, pt.size = .1,
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5')),
          dpi = 300)
ggsave('Outputs/figures/Donor_12_UMAP_cDCs.pdf', width = 4.5, height = 3)
# Figure 1F Dot plot showing marker expression in Donor 6 cDC subpopulations
DotPlot(d12.object.cdcs, 
        features = read.csv('Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Outputs/figures/Dotplot_of_cDC_markers_Donor12.pdf',
       width = 10, height = 2.7)
saveRDS(d12.object.cdcs,
        'Outputs/rds/Dnr12_cDCs.rds')


# CytoTRACE analysis of Dnr6 cDCs -----------------------------------------
  # Subset Donor 6 cDC2s
d6.object.cdc2s <- d6.object.cdcs[,d6.object.cdcs$celltype_hh_1!='cDC1']
d6.object.cdc2s <- NormalizeData(d6.object.cdc2s) %>% FindVariableFeatures() %>% ScaleData(features = rownames(d6.object.cdc2s)) %>% RunPCA(verbose=FALSE) 
d6.object.cdc2s <- RunUMAP(d6.object.cdc2s, dims=1:10)
DimPlot(d6.object.cdc2s)
d6.object.cdc2s <- FindNeighbors(d6.object.cdc2s) %>% FindClusters(resolution = 0.2)
# Filter out cluster 3 as it highly express C1QC and APOE, may be Macro
d6.object.cdc2s <- d6.object.cdc2s[,
                                 d6.object.cdc2s$RNA_snn_res.0.2!='3']
d6.object.cdc2s <- NormalizeData(d6.object.cdc2s) %>% FindVariableFeatures() %>% ScaleData(features = rownames(d6.object.cdc2s)) %>% RunPCA(verbose=FALSE) 
d6.object.cdc2s <- RunUMAP(d6.object.cdc2s, dims=1:10)
d6.object.cdc2s <- RenameIdents(d6.object.cdc2s,c('0'='cDC2','1'='DC2hm',
                                                '2'='DC2pre-hm',
                                                '4'='DC2pre-hm'))
d6.object.cdc2s$celltype_hh_1 <- Idents(d6.object.cdc2s)
levels(d6.object.cdc2s) <- c('cDC2', 'DC2pre-hm','DC2hm')
# This cDC2-to-DC2hm UMAP is used for Figure 1H-1K
rasterize(DimPlot(d6.object.cdc2s, pt.size = .1,
                  cols = c('#38C2E5','#384C94','#E679C5')),
          dpi = 300)
ggsave('Outputs/figures/Donor_6_UMAP_cDC2s.pdf', width = 4.5, height = 3)
# Save the dnr6 cDC2s
saveRDS(d6.object.cdc2s, 'Outputs/rds/Dnr6_cDC2s.rds')
  # CytoTRACE analysis
exp.mat <- as.matrix(d6.object.cdc2s@assays$RNA@counts)
phe <- as.character(d6.object.cdc2s$celltype_hh_1)
names(phe) <- rownames(d6.object.cdc2s@meta.data)
cytotrace.result <- CytoTRACE(exp.mat,ncores = 5)
# Figure 1I and 1J
plotCytoTRACE(cytotrace.result,phenotype = phe,
              emb = d6.object.cdc2s@reductions$umap@cell.embeddings,outputDir = 'Outputs/figures/')


# Monocle3 analysis -------------------------------------------------------
count.matrix <- GetAssayData(d6.object.cdc2s, assay = 'RNA',
                             slot = 'counts')
cell_metadata <- d6.object.cdc2s@meta.data
gene.names <- data.frame(gene_short_name=rownames(count.matrix),
                         row.names = rownames(count.matrix))
# Build CDS object
cds <- new_cell_data_set(count.matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene.names)
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds, preprocess_method = 'PCA',
                        reduction_method = 'UMAP')
cds <- cluster_cells(cds)
# Use Seurat UMAP
cds@int_colData$reducedDims$UMAP <- d6.object.cdc2s@reductions$umap@cell.embeddings
plot_cells(cds, color_cells_by = 'celltype_hh_1')
# Monocle 3 analysis
cds <- learn_graph(cds,use_partition = F, 
                   learn_graph_control=list('minimal_branch_len'=22,
                                            'orthogonal_proj_tip'=F))
plot_cells(cds)
cds <- order_cells(cds)
# Plot pseudotime in UMAP space, Figure 1K
rasterize(plot_cells(cds, color_cells_by = 'pseudotime',label_cell_groups = F,
                     label_roots = F,label_branch_points = F,label_leaves = F,
                     alpha = 1,cell_stroke = 1,cell_size = 0.1),
          dpi = 300)
ggsave('Outputs/figures/Pseudotime_cDC2_traject.pdf',
       width = 4.5, height = 3)

# Modulated gene expression in pseudotime ---------------------------------

# Modulating genes  
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", 
                              cores = 8)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = unique(c(genes, 
                                                  VariableFeatures(d6.object.cdc2s))))
ck <- clusterData(exp = mat,
                  cluster.method = "kmeans",
                  cluster.num = 3)
saveRDS(ck,'Outputs/rds/ck.rds')
ck <- readRDS('ck.rds')
# Enrichment
  # Define enrichment function
IDTransandenrichWP <- function(gene.list){
  result <- enrichWP(bitr(gene.list,
                          fromType = 'SYMBOL', toType = 'ENTREZID',
                          OrgDb = 'org.Hs.eg.db')$ENTREZID,
                     organism = 'Homo sapiens')
  return(result)
}
IDTransandenrichKEGG <- function(gene.list){
  result <- enrichKEGG(bitr(gene.list,
                            fromType = 'SYMBOL', toType = 'ENTREZID',
                            OrgDb = 'org.Hs.eg.db')$ENTREZID,
                       organism = 'hsa')
  return(result)
}
enrich.result.list <- list()
enrich_go <- data.frame()
for (i in seq(1,3)) {
  go.enrich <- enrichGO(unique(ck$wide.res[ck$wide.res$cluster==i,
                                           'gene']),
                        OrgDb = 'org.Hs.eg.db',keyType = 'SYMBOL',
                        ont = 'ALL')@result
  wiki.enrich <- IDTransandenrichWP(unique(ck$wide.res[ck$wide.res$cluster==i,
                                                       'gene']))@result
  wiki.enrich$ONTOLOGY <- rep('Wikipathway',nrow(wiki.enrich))
  enrich.result.list[[i]] <- rbind(go.enrich, wiki.enrich,kegg.enrichment)
  enrich.result.list.filtered <- enrich.result.list[[i]]
  enrich.result.list.filtered <- enrich.result.list.filtered[enrich.result.list.filtered$Count > 5 & enrich.result.list.filtered$p.adjust<0.05,]
  pp <- data.frame(id=paste0('C',i),
                   term=enrich.result.list.filtered$Description,
                   pval=enrich.result.list.filtered$p.adjust)
  enrich_go <- rbind(enrich_go,
                     pp)
}
enrich_go.clean <- enrich_go[!duplicated(enrich_go$term),]
write.csv(enrich_go.clean,'Outputs/tables/Enrich_cluster_cleaned.csv')
# Clean enrichment
# Save Enrich wikipathway genes
write.xlsx(enrich.result.list, 'Outputs/tables/Enrich_cluster.xlsx')
write.csv(enrich_go,'Outputs/tables/Enrich_cluster_filtered_combined.csv')
write.xlsx(list('C1'=unique(ck$wide.res[ck$wide.res$cluster==1,
                                        'gene']),
                'C2'=unique(ck$wide.res[ck$wide.res$cluster==2,
                                        'gene']),
                'C3'=unique(ck$wide.res[ck$wide.res$cluster==3,
                                        'gene'])),
           'Outputs/tables/Gene_clusters.xlsx')
dc3.markers <- read.csv('Highlighted_gene_selected_with_isgs.csv')
# Figure 1L pseudotime heatmap, note that final version abbreviate some pathways for layout's sake
pdf('Outputs/figures/monocle3_pseudotime_heatmap.pdf',height = 10,width = 9,
    onefile = F)
visCluster(object = ck,
           plot.type = "both",
           add.sampleanno = F,
           markGenes = c('OASL','ITGAX','TLR2','NLRP3','IL1B','CD1C','ITGAX','FCGR1A',
                         'IL18','IFI6','IFITM2','IFITM1','PYCARD','IFIT1',
                         'MX1','MX2','ZEB2','IFI44','RPL26','RPL24','RPS3A',
                         intersect(dc3.markers$Gene, 
                                   unique(ck$wide.res[ck$wide.res$cluster==2,'gene'])),
                         'HLA-DRA','HLA-DRB1','CD86','CXCL10','CTLA4',
                         'CD14','CD36'),
           annoTerm.data = enrich_go.clean[c(1225,1139,1007,1000,983,
                                             1296,1275,
                                             484,272,340,333, 291, 846),],
           go.size = 10, markGenes.side = 'left',
           line.side = 'left', show_row_dend = F, go.col = 'black',
           ht.col.list = list(col_range = c(-2, 0, 2),
                              col_color = c('#2A9FF3','white' ,'#FF2B72')),
           ctAnno.col = c('#38C2E5','#E679C5','#384C94'),cluster.order = c(2,3,1))
dev.off()

# Integrate cDCs from all 3 donors ------------------------------------
combined.object <- merge(d6.object.cdcs,c(d10.object.cdcs,d12.object.cdcs))
combined.object <- NormalizeData(combined.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunHarmony('orig.ident', 
                                                                                                                                      plot_convergence = T,
                                                                                                                                      max.iter.harmony = 50)
combined.object <- RunUMAP(combined.object, verbose = F, 
                           reduction = 'harmony', dims = 1:6)
combined.object <- FindNeighbors(combined.object, reduction = 'harmony')
combined.object <- FindClusters(combined.object, resolution = 0.5)
combined.object <- combined.object[,!combined.object$seurat_clusters%in% c('4',
                                                                           '8',
                                                                           '9',
                                                                           '10')]
combined.object <- NormalizeData(combined.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunHarmony('orig.ident', 
                                                                                                                                      plot_convergence = T,
                                                                                                                                      max.iter.harmony = 50)
combined.object <- RunUMAP(combined.object, verbose = F, 
                           reduction = 'harmony', dims = 1:5)
combined.object <- RenameIdents(combined.object,c('0'='cDC2','1'='cDC2',
                                                  '2'='cDC2',
                                                  '3'='DC2hm', 
                                                  '5'='DC2pre-hm',
                                                  '6'='cDC2',
                                                  '7'='cDC1'))
levels(combined.object) <- c('cDC1','cDC2','DC2pre-hm','DC2hm')
combined.object$orig.ident <- factor(combined.object$orig.ident,
                                     levels =  c('Dnr6', 'Dnr10', 'Dnr12'))
  # Figure S2C merged cDC UMAP
rasterize(DimPlot(combined.object,
                  pt.size = .1,
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5')),
          dpi = 300)+
  rasterize(DimPlot(combined.object,
                    pt.size = .1,
                    group.by = 'orig.ident')+
              scale_color_npg(),
            dpi = 300)
ggsave('Outputs/figures/UMAP_merged_cDCs.pdf', width = 9, height = 3)
combined.object$celltype_hh_1 <- Idents(combined.object)
DotPlot(combined.object, 
        features = read.csv('Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  # Figure S2D, Dot plot showing the expression of merged cDCs
ggsave('Outputs/figures/Dotplot_of_cDC_markers_merged.pdf',
       width = 10, height = 2.7)
  # Figure S2E, heatmap showing the expression of merfed cDCs
marker.genes <- FindAllMarkers(combined.object,only.pos = T)
marker.genes.ranked <- marker.genes[order(marker.genes$avg_log2FC, decreasing = T),]
# Scale all the genes for illustration
combined.object <- ScaleData(combined.object, features = rownames(combined.object))
# Get the top 100 markers for each cluster
genes.to.show <- c(marker.genes.ranked[marker.genes.ranked$cluster=='cDC1','gene'],
                   marker.genes.ranked[marker.genes.ranked$cluster=='cDC2','gene'],
                   marker.genes.ranked[marker.genes.ranked$cluster=='DC2hm','gene'])
  # Get matrix for complex heatmap
heatmap.matrix <- GetAssayData(combined.object[genes.to.show,],
                               slot = 'scale.data')
heatmap.matrix <- as.matrix(heatmap.matrix[genes.to.show,names(sort(Idents(combined.object)))])
  # Make annotation of highlighted marker genes
genes.to.anno <- intersect(read.csv('Highlighted_gene_selected_extended.csv')$Gene,
                           genes.to.show)
gene.pos <- which(rownames(heatmap.matrix) %in% genes.to.anno)
row.anno <- rowAnnotation(gene=anno_mark(at=gene.pos,
                                         labels = genes.to.anno))
top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=c('#C5E524',
                                                                '#38C2E5',
                                                                '#384C94',
                                                                '#E679C5')),
                                                 labels = levels(sort(Idents(combined.object))),
                                                 labels_gp = gpar(cex=0.5, col='white')))
col_fun = colorRamp2(c(-2, 1, 2), c('#2A9FF3','white' ,'#FF2B72'))
pdf('Outputs/figures/Heatmap_of_cDC_markers.pdf', width = 6, height = 10)
Heatmap(heatmap.matrix, cluster_rows = F, cluster_columns = F,
        show_column_names = F, show_row_names = F, right_annotation = row.anno,
        use_raster = T, raster_quality =10, raster_device = 'png',
        col = col_fun, raster_by_magick = F, column_split = sort(Idents(combined.object)),
        top_annotation = top_anno)
dev.off()
# Digital CCR7 sorting ----------------------------------------------------
combined.object$ccr7_exp <- ifelse(combined.object@assays$RNA@counts['CCR7',]==0,
                                   'CCR7-', 'CCR7+')
  # Figure S4A, UMAP colored by CCR7 RNA detection
rasterize(DimPlot(combined.object,group.by = 'ccr7_exp',pt.size = 0.1,)+
            scale_color_aaas(),
          dpi = 300)
ggsave('Outputs/figures/UMAP_CCR7_sorting_splited.pdf', width = 4.5,
       height = 3)
  # Calculate proportion and do the sorting
proportion.table <- as.data.frame(prop.table(table(combined.object$ccr7_exp,
                                                   Idents(combined.object)),
                                             margin = 2))
colnames(proportion.table) <- c('Sorted Cell Type',
                                'Annotated Cell Type',
                                'Proportion')
  # Figure S4B, roportion bar plot
ggplot(proportion.table, aes(x=`Annotated Cell Type`,
                             y=Proportion, fill=`Sorted Cell Type`))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_aaas()
ggsave('Outputs/figures/CCR7_sorting_sort_by_cluster.pdf', 
       width = 4, height = 3)
  # Figure S4C, DC2hm marker expression 
VlnPlot(combined.object[, Idents(combined.object)=='cDC2'], 
        features = c('CCR7', 'CD83', 'CD274', 'FSCN1',
                     'FAS', 'ICOSLG', 'CCL22'), group.by = 'ccr7_exp',
        pt.size = 0, ncol = 7, cols = c('#3e498d','#da2f20'))
ggsave('Outputs/figures/Violin_plot_DC2hm_exp.pdf', width = 10, height = 3)
# Hierachical clustering with fetus and hu-mice ---------------------------
# Fetal spleen object were obtained from xxx
fetal.spleen <- readRDS('object.cDCs.rds') #cDCs from all fetal tissues
fetal.spleen <- fetal.spleen[,fetal.spleen$tissue=='Spleen']
fetal.spleen <- RenameIdents(fetal.spleen,
                             c('cDC1'='cDC1', 'cDC2'='cDC2',
                               'DC3'='DC2hm'))
  # Figure S2H, UMAP of fetal spleens
rasterize(DimPlot(fetal.spleen, pt.size = .1,
                  cols = c('#C5E524', '#38C2E5','#E679C5')),
          dpi = 300)
ggsave('Outputs/figures/Fetal_spleen_UMAP_cDCs.pdf', width = 4.5, height = 3)
# Figure 1I Dot plot showing marker expression in fetal spleen cDC subpopulations
DotPlot(fetal.spleen, 
        features = read.csv('Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Outputs/figures/Dotplot_of_cDC_markers_fetal_spleen.pdf',
       width = 10, height = 2.7)
# Hu-mice spleen object were obtained from xxx
humice.spleen <- readRDS('Mock_humice.rds') 
# Figure S2F, UMAP of fetal spleens
rasterize(DimPlot(humice.spleen, pt.size = .1,
                  cols = c('#C5E524', '#38C2E5','#E679C5')),
          dpi = 300)
ggsave('Outputs/figures/Humice_spleen_UMAP_cDCs.pdf', width = 4.5, height = 3)
# Figure S2G Dot plot showing marker expression in humice cDC subpopulations
DotPlot(humice.spleen, 
        features = read.csv('Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Outputs/figures/Dotplot_of_cDC_markers_humice.pdf',
       width = 10, height = 2.7)
# Hierachical clustering
  # Prepare pseudobulk profiles
# Prepare pseudobulk
markers <- FindAllMarkers(combined.object, only.pos = T)
marker.cluster.list <- list()
for (i in names(table(markers$cluster))) {
  marker.cluster <- markers[markers$cluster==i & markers$p_val_adj<0.05,]
  marker.cluster.list[[i]] <- marker.cluster[order(marker.cluster$avg_log2FC,
                                                   decreasing = T),]
}
write.xlsx(marker.cluster.list, 
           'Outputs/tables/cDC_markers_per_cluster_merged.xlsx',
           rowNames=T)
# Prepare pseudobulk
  # Prepare gene lists for analysis
genes <- unique(c(marker.cluster.list[[1]][1:100,'gene'],
                  marker.cluster.list[[2]][1:100,'gene'],
                  marker.cluster.list[[4]][1:100,'gene']))
genes <- intersect(genes, rownames(fetal.spleen))
genes <- intersect(genes, rownames(humice.spleen))
# To scale all features before generating pseudobulks
combined.object <- ScaleData(combined.object, features = rownames(combined.object))
fetal.spleen <- ScaleData(fetal.spleen , features = rownames(fetal.spleen))
humice.spleen <- ScaleData(humice.spleen, features = rownames(humice.spleen))
pseudo.bulk <- data.frame(list('human spleen cDC1'=AverageExpression(combined.object[genes,Idents(combined.object)=='cDC1']),
                               'human fetal spleen cDC1'=AverageExpression(fetal.spleen[genes,Idents(fetal.spleen)=='cDC1']),
                               'humanzied mice cDC1'=AverageExpression(humice.spleen[genes,Idents(humice.spleen)=='cDC1']),
                               'human spleen cDC2'=AverageExpression(combined.object[genes,Idents(combined.object)=='cDC2']),
                               'human fetal spleen cDC2'=AverageExpression(fetal.spleen[genes,Idents(fetal.spleen)=='cDC2']),
                               'humanzied mice cDC2'=AverageExpression(humice.spleen[genes,Idents(humice.spleen)=='cDC2']),
                               'human spleen DC2hm'=AverageExpression(combined.object[genes,Idents(combined.object)=='DC2hm']),
                               'human fetal spleen DC2hm'=AverageExpression(fetal.spleen[genes,Idents(fetal.spleen)=='DC2hm']),
                               'humanzied mice DC2hm'=AverageExpression(humice.spleen[genes,Idents(humice.spleen)=='DC2hm'])),
                          row.names = genes)
colnames(pseudo.bulk) <- names(list('human spleen cDC1'=AverageExpression(combined.object[genes,Idents(combined.object)=='cDC1']),
                                    'human fetal spleen cDC1'=AverageExpression(fetal.spleen[genes,Idents(fetal.spleen)=='cDC1']),
                                    'humanzied mice cDC1'=AverageExpression(humice.spleen[genes,Idents(humice.spleen)=='cDC1']),
                                    'human spleen cDC2'=AverageExpression(combined.object[genes,Idents(combined.object)=='cDC2']),
                                    'human fetal spleen cDC2'=AverageExpression(fetal.spleen[genes,Idents(fetal.spleen)=='cDC2']),
                                    'humanzied mice cDC2'=AverageExpression(humice.spleen[genes,Idents(humice.spleen)=='cDC2']),
                                    'human spleen DC2hm'=AverageExpression(combined.object[genes,Idents(combined.object)=='DC2hm']),
                                    'human fetal spleen DC2hm'=AverageExpression(fetal.spleen[genes,Idents(fetal.spleen)=='DC2hm']),
                                    'humanzied mice DC2hm'=AverageExpression(humice.spleen[genes,Idents(humice.spleen)=='DC2hm'])))
col.anno <- data.frame('Cell_type'=c(rep('cDC1',3),
                                     rep('cDC2',3),
                                     rep('DC2hm',3)),
                       row.names = colnames(pseudo.bulk))
row.anno <- data.frame('Gene_group'=c(rep('cDC1 markers',90),
                                      rep('cDC2 markers',88),
                                      rep('DC2hm markers',88)),
                       row.names = rownames(pseudo.bulk))
as.ggplot(pheatmap::pheatmap(pseudo.bulk, cluster_rows = F,cluster_cols = T,
                             scale = 'row', show_rownames = F, cellwidth = 25,
                             annotation_col = col.anno, 
                             annotation_colors = list('Cell_type'=c('cDC1'='#C5E524',
                                                                    'cDC2'='#38C2E5',
                                                                    'DC2hm'='#E679C5'),
                                                      'Gene_group'=c('cDC1 markers'='#C5E524',
                                                                     'cDC2 markers'='#38C2E5',
                                                                     'DC2hm markers'='#E679C5')),
                             annotation_row = row.anno,
                             clustering_distance_cols = "canberra"))
ggsave('Outputs/figures/Hierarchical_clustering_heatmap.pdf')


# Save the workspace
save.image(file = 'scRNA-seq_analysis.RData')
