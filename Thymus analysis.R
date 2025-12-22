library(ggrastr)
library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(scCustomize)
library(openxlsx)
library(ggplotify)
library(pheatmap)
library(lsa)
set.seed(0822)
setwd("Analysis on E-MTAB-8581(thymus atlas)")
  # The RDS file could be obtained @https://developmental.cellatlas.io/thymus-development and converveted from h5ad to rds
thymus.object <- readRDS('../../../thymus data/hu_thymus.rds')
# mTEC clustering analysis ------------------------------------------------
# Subset mTECs
thymus.object.mTEC <- thymus.object[, thymus.object$Anno_level_fig1 %in% c('mTEC(I)',
                                                                   'mTEC(II)',
                                                                   'mTEC(III)',
                                                                   'mTEC(IV)')]
thymus.object.mTEC <- NormalizeData(thymus.object.mTEC)  %>% FindVariableFeatures() %>% ScaleData()%>% RunPCA()
thymus.object.mTEC$Sample <- droplevels(factor(thymus.object.mTEC$Sample))
thymus.object.mTEC <- RunHarmony(thymus.object.mTEC,
                                 group.by.vars = 'Sample')
thymus.object.mTEC <- RunUMAP(thymus.object.mTEC, dims = 1:5,
                              reduction = 'harmony')
rasterize(DimPlot(thymus.object.mTEC,
                  group.by = 'Anno_level_fig1',
                  pt.size = .1)+
            scale_color_simpsons(),dpi = 300)
ggsave('Outputs/figures/_mTEC.pdf', width = 4, height = 3)
  # mTEC subtype marker
DotPlot(thymus.object.mTEC, features = c('CD80','CD40','AIRE',
                                         'IVL','BIK','GADD45G'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  scale_color_viridis_c()
ggsave('Outputs/figures/Dotplot_mTECs.pdf', width = 8, height = 4)
 # DC2hm signature scores
thymus.object.mTEC <- AddModuleScore(thymus.object.mTEC,
                                            list(read.xlsx('../Analysis on sorted cDCs/Outputs/tables/Signature_genes.xlsx',
                                                          sheet = 'IL7RpCCR7p')$genes),
                                            name = 'IL7RpCCR7p.score') # Run the signature score in Analysis on sorted cDCs first
Idents(thymus.object.mTEC) <- 'Anno_level_fig1'
VlnPlot(thymus.object.mTEC, features = 'IL7RpCCR7p.score1', pt.size = 0)+
  scale_fill_simpsons()
ggsave('Outputs/figures/IL7R+CCR7+_score_mTECs.pdf', width = 4, height = 3)
rasterize(FeaturePlot(thymus.object.mTEC, features = 'IL7RpCCR7p.score1',
                      pt.size = .1)+
            scale_color_viridis_c(option = 4), dpi = 300)
ggsave('Outputs/figures/Featureplot_IL7R+CCR7+_score_mTECs.pdf',
       width = 4, height = 3)
# mTEC signature expression in cDCs -----------------------------------------------
  # canonical mTEC marker expression in cDCs
# Show mTEC marker dotplot
combined.object <- readRDS('../scRNA-seq analysis/Outputs/rds/combined.object.rds')
gene_names <- c(
  "FAS", "CCL19", "RELB", "NFKBIA", "SOCS3", "CD40", "TNFAIP2", "HLA-DQA1", "ICAM1",
  "LGALS3", "ALDH2", "SYNGR2", "PRRG4", "TRIP10", "VASN", "ZFP36L1", "CHCHD10", "JAG1",
  "BLVRA", "BIRC3", "TMEM150C", "TYMP", "SCARA3", "JUP", "DSC2", "ANXA2", "GBP1", 
  "KREMEN2", "SQSTM1", "RBP1", "STAP2", "CDC42EP4", "ATF3", "PNRC1", "DDR1", "PTPRF", 
  "PTGDS", "LINCO1551", "SYNPO2", "BHLHE40", "TP63", "LAMC2", "OSMR", "KLF5", "TRIB1", 
  "FXYD3", "DMKN", "ERRFI1", "PERP", "KRT5")
DotPlot(combined.object, features = gene_names)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  scale_color_viridis_c()
ggsave('Outputs/figures/mTEC_marker_expression_in_cDCs.pdf', width = 12, height = 3.5)
  # mTEC subtype markers
markers.list <- list()
thymus.object.mTEC$Anno_level_fig1 <- droplevels(thymus.object.mTEC$Anno_level_fig1)
for (i in levels(thymus.object.mTEC$Anno_level_fig1)) {
  cluster.marker <- FindMarkers(thymus.object.mTEC, ident.1 = i,
                                group.by = 'Anno_level_fig1')
  markers.list[[i]] <- cluster.marker[order(cluster.marker$avg_log2FC,
                                            decreasing = T),]
}
openxlsx::write.xlsx(markers.list,
                     'Outputs/tables/mTECs_cluster_markers.xlsx',rowNames=T)
  # mTEC II signature score
mtec2.markers <- read.xlsx('Outputs/tables/mTECs_cluster_markers.xlsx',
                           sheet = "mTEC(II)", rowNames = T)
d6.object <- AddModuleScore(d6.object,
                            features = list(rownames(mtec2.markers)[1:200]),
                            name = 'mTEC(II)score', seed = 0822)
VlnPlot(d6.object, features = 'mTEC.II.score1', cols = c('cDC1'='#C5E524',
                                                         'cDC2'='#38C2E5',
                                                         'DC2pre-hm'='#384C94',
                                                         'DC2hm'='#E679C5'),
        pt.size = 0)
ggsave('Outputs/figures/mTEC2.score_in_cDCs.pdf', 
       width = 4, height = 3)
rasterize(FeaturePlot(d6.object, features = 'mTEC.II.score1',
                      pt.size = .1)+
            scale_color_viridis_c(option = 4), dpi = 300)
ggsave('Outputs/figures/Featureplot_mTECII_score_cDCs.pdf',
       width = 4, height = 3)
rasterize(FeaturePlot(d6.object, features = 'mTEC.II.score1',
                      pt.size = .1)+
            scale_color_viridis_c(option = 4), dpi = 300)

# Densities of mTECs positve for DC2hm markers
plt.list <- list()
for (i in c('AIRE','CCR7','IL7R', 'CRLF2','STAT5A')) {
  plt.list[[i]] <- rasterize(as.ggplot(Plot_Density_Custom(thymus.object.mTEC, features = i,
                                                           pt.size = 1)),
                             dpi=300)
}
ggarrange(plotlist = plt.list, ncol = 5, nrow = 1)
ggsave('Outputs/figures/Marker_density_mTEC.pdf', width = 24, height = 4)
Plot_Density_Joint_Only(seurat_object = thymus.object.mTEC, features = c('AIRE','CCR7','IL7R', 'CRLF2','STAT5A'))


# Cosine similarity -------------------------------------------------------
cdc.pseudubulk <- AverageExpression(combined.object, assays = "RNA", slot = "data")$RNA
mtec.pseudobulk <- AverageExpression(thymus.object.mTEC, assays = "RNA", slot = "data")$RNA
common_genes <- intersect(VariableFeatures(combined.object), rownames(thymus.object.dc.and.mtec))
cdc.pseudubulk <- cdc.pseudubulk[common_genes, ]
mtec.pseudobulk <- mtec.pseudobulk[common_genes, ]
similarity_matrix <- matrix(0, ncol = ncol(cdc.pseudubulk), nrow = ncol(mtec.pseudobulk),
                            dimnames = list(colnames(mtec.pseudobulk), colnames(cdc.pseudubulk)))

# Calculate cosine similarity for each pair of clusters
for (i in colnames(cdc.pseudubulk)) {
  for (j in colnames(mtec.pseudobulk)) {
    similarity_matrix[j, i] <- cosine(cdc.pseudubulk[, i], mtec.pseudobulk[, j])
  }
}
palette.heatmap <- colorRampPalette(c('#2A9FF3','white' ,'#FF2B72'))
as.ggplot(pheatmap(similarity_matrix, cluster_rows = F,
                   cluster_cols = F, scale = 'row', 
                   display_numbers = T, cellwidth = 50, cellheight = 50,
                   color =  palette.heatmap(200)))
ggsave('Outputs/figures/Cosine_similarity_heatmap.pdf')

