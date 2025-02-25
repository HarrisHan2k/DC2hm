library(Seurat)
library(ggplot2)
library(ggpubr)
library(harmony)
library(ggrastr)
library(viridis)
setwd("/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/Hu-mice")

# QC ----------------------------------------------------------------------
work.list <- c("MOCK_1","MOCK_3","TSLP_1", "TSLP_2_3")
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
  PlotQC(work.list[[i]])
}
th.nfeature.mins <- rep(500,4)
th.nfeature.maxs <- rep(7500,4)
th.mt.pcts <- rep(10,4)
th.count.mins <- rep(2000,4)
cell.numbers.before <- c()
cell.numbers.after <- c()
DoQCAndStore <- function(project,th.nfeature.min,
                         th.nfeature.max, th.mt.pct,th.count.min){
  matrix <- Read10X(paste0('CellrangerOutputs/',project,'/filtered_feature_bc_matrix'))
  object <- CreateSeuratObject(matrix, project = project, min.cells = 3, 
                               min.features = 200)
  cell.number.before <- ncol(object)
  object[['percent.mt']] <- PercentageFeatureSet(object ,pattern = '^mt-')
  object <- subset(object, subset = nFeature_RNA > th.nfeature.min &
                     nFeature_RNA <th.nfeature.max & percent.mt < th.mt.pct &
                     nCount_RNA > th.count.min)
  cell.number.after <- ncol(object)
  saveRDS(object, paste0('Outputs/rds/', project, '.object.rds'))
  return(c(cell.number.before, cell.number.after))
}
for (i in seq(1, length(work.list))) {
  cell.numbers.before.and.after <- DoQCAndStore(work.list[[i]],
                                                th.nfeature.max = th.nfeature.maxs[[i]],
                                                th.nfeature.min = th.nfeature.mins[[i]],
                                                th.mt.pct = th.mt.pcts[[i]],
                                                th.count.min=th.count.mins[[i]])
  cell.numbers.before <- c(cell.numbers.before, cell.numbers.before.and.after[[1]])
  cell.numbers.after <- c(cell.numbers.after, cell.numbers.before.and.after[[2]])
}
# Clustering all cells--------------------------------------------------------------

object.mock1 <- readRDS('Outputs/rds/MOCK_1.object.rds')
object.mock3 <- readRDS('Outputs/rds/MOCK_3.object.rds')
object.tslp.1 <- readRDS('Outputs/rds/TSLP_1.object.rds')
object.tslp.2.3 <- readRDS('Outputs/rds/TSLP_2_3.object.rds')

combined.object <- merge(object.mock1, 
                         c(object.mock3,object.tslp.1,object.tslp.2.3))
combined.object <- NormalizeData(combined.object) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunHarmony('orig.ident', 
                                                                                                                                      plot_convergence = T)
combined.object <- RunUMAP(combined.object, verbose = F, 
                           reduction = 'harmony', dims = 1:15)
combined.object <- FindNeighbors(combined.object, reduction = 'harmony')
combined.object <- FindClusters(combined.object, resolution = 0.5)
rasterize(DimPlot(combined.object, label = T,
                  pt.size = .1),dpi=300)
ggsave('Outputs/figures/UMAP_of_all_cells.pdf', width = 4, height = 3)
DotPlot(combined.object, features = c('CLEC9A','XCR1','CLEC10A',
                                      'CD1C','CCR7','LAMP3'))
ggsave('Outputs/figures/Dotplot_DC_markers.pdf', width = 5, height = 5)
saveRDS(combined.object,
        'Outputs/rds/Combined_clustered_object_all.rds') # This object is available @ GEO

# Clustering DCs and compare proportions ----------------------------------
  # cluster 12 and 18 contains DCs
combined.object.dcs <- combined.object[,Idents(combined.object) %in% c('12',
                                                                       '18')]
combined.object.dcs <- NormalizeData(combined.object.dcs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunHarmony('orig.ident', 
                                                                                                                                      plot_convergence = T)
combined.object.dcs <- RunUMAP(combined.object.dcs, verbose = F, 
                           reduction = 'harmony', dims = 1:10)
combined.object.dcs <- FindNeighbors(combined.object.dcs, reduction = 'harmony')
combined.object.dcs <- FindClusters(combined.object.dcs, resolution = 0.5)
DotPlot(combined.object.dcs, features = c('CLEC9A','XCR1','CLEC10A',
                                          'CD1C','CCR7','LAMP3'))
  # Cluster 2,3,4 represent DC2hm, cDC2, and cDC1s respectively
combined.object.dcs <- combined.object.dcs[,Idents(combined.object.dcs) %in% c('2',
                                                                               '3',
                                                                               '4')]
combined.object.dcs <- NormalizeData(combined.object.dcs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE) %>% RunHarmony('orig.ident', 
                                                                                                                                              plot_convergence = T)
combined.object.dcs <- RunUMAP(combined.object.dcs, verbose = F, 
                               reduction = 'harmony', dims = 1:10)
combined.object.dcs <- FindNeighbors(combined.object.dcs)
combined.object.dcs <- FindClusters(combined.object.dcs, resolution = .5)
combined.object.dcs <- RenameIdents(combined.object.dcs, c('3'='cDC1',
                                                           '0'='DC2hm',
                                                           '1'='cDC2',
                                                           '2'='cDC1'))
levels(combined.object.dcs) <- c('cDC1','cDC2','DC2hm')
rasterize(DimPlot(combined.object.dcs), dpi = 300)
ggsave('Outputs/figures/UMAP_cDCs.pdf', width = 4, height = 3)
DotPlot(combined.object.dcs[,combined.object.dcs$orig.ident %in% c('TSLP_1',
                                                                   'TSLP_2_3')], 
        features = read.csv('Highlighted_gene_selected_with_isgs.csv')[,1])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  scale_color_viridis_c()
ggsave('Outputs/figures/Dotplot_markers_cDCs.pdf', width = 10, height = 2.5)  

combined.object.dcs$group <- ifelse(combined.object.dcs$orig.ident %in% c('MOCK_1',
                                                                          'MOCK_3'),
                                    'MOCK','TSLP')
proportions <- data.frame(prop.table(table(Idents(combined.object.dcs),
                 combined.object.dcs$group), margin = 2))
colnames(proportions) <- c('Celltype', 'Group', 'Freq')
ggplot(proportions, aes(x=Group, y=Freq, fill=Celltype))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_fill_manual(values = c('#C5E524', '#38C2E5','#E679C5'))+
  theme_classic()
ggsave('Outputs/figures/Proportions_barplot.pdf', width = 3, height = 3)
saveRDS(combined.object.dcs, 'Outputs/rds/Combined_clustered_object_cdcs.rds')
