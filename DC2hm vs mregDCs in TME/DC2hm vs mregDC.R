library(Seurat)
library(ggrastr)
library(ggsci)
library(ggpubr)
library(clusterProfiler)
library(openxlsx)
library(enrichplot)
set.seed(0822)
setwd('/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/DC2hm vs mregDCs in TME')
# Visualization of scRNA-seq profiles -------------------------------------
  # Read Seurat object
combined.object <- readRDS('combined.object.dc2hm.mregdc.rds')

  # Figure S8B, t-SNE plot of binary-grouped profiles
rasterize(DimPlot(combined.object, group.by = 'binary.group',
                  pt.size = 0.1,reduction = 'tsne',
                  cols = rev(c('#70d6ff','#ff70a6'))),
          dpi = 300)
ggsave('Outputs/figures/TNSE_group.pdf', 
       width = 4.5, height = 3)
  # Figure S8B, t-SNE plot of sample-grouped profiles
rasterize(DimPlot(combined.object, group.by = 'orig.ident',
                  pt.size = 0.1,reduction = 'tsne')+
            scale_color_npg(), dpi = 300)
ggsave('Outputs/figures/TNSE_samplep.pdf', 
       width = 4.5, height = 3)
  # Figure S8C, showing CCR7 and LAMP3 expression to ensure they are mature DCs
plt.list <- list()
for (i in c('FLT3','CCR7','LAMP3','CRLF2')) {
  p <- rasterize(FeaturePlot(combined.object, features = i,
                             pt.size = 0.1,
                             reduction = 'tsne')+
                   scale_color_gradient(low = 'lightgrey', high='#FF2B72'),
                 dpi=300)
  plt.list[[i]] <- p
}
ggarrange(plotlist = plt.list, ncol = 4)
ggsave('Outputs/figures/Feature_plots_DC_maturation.pdf',
       width = 16, height = 4)
  # Figure S8D, dot plot showing selected gene expression in the 2 groups
DotPlot(combined.object, features = c('IDO2','AHR','CRLF2','IL7R',
                                      'ITGAX','HLA-DPA1','HLA-DRA','HLA-DPB1',
                                      'ZEB2','ETS2','HIVEP1','HIVEP2','ETV6','TET2',
                                      'BCL2',
                                      'CCL5','CCL19','CXCL9','CXCL10',
                                      'IL12B','IL18','IL32',
                                      'IFI6','IFI27','IFITM1',
                                      
                                      ,'OXA1L',
                                      ),
        group.by = 'binary.group')+
  coord_flip()+
  scale_color_gradient(low = 'lightgrey', high='#FF2B72')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Outputs/figures/Tumor_vs_spleen_genes_Dotplot.pdf',
       width = 4.5, height = 6)
  # Figure S8F, violin plot showing the CRLF2 expression inthe 2 groups
VlnPlot(combined.object, features = c('CRLF2'), group.by = 'binary.group',
        pt.size = 0,cols = rev(c('#70d6ff','#ff70a6')))
ggsave('Outputs/figures/Vlnplot_CRLF2.pdf', width = 4, height = 3)

# GSEA analysis -----------------------------------------------------------
  # Spleen vs Tumor DE analysis
tumor.vs.spleen.markers <- FindMarkers(combined.object,
                                       ident.1 = 'Tumor', ident.2 = 'HC_spleen',
                                       group.by = 'binary.group')
  # Remove false pleen.vs.tumorse-positive DEGs which may due to gene name discrepancy
tumor.vs.spleen.markers <- tumor.vs.spleen.markers[!(tumor.vs.spleen.markers$pct.1 > 0.5 & tumor.vs.spleen.markers$pct.2==0) |
                                                     !(tumor.vs.spleen.markers$pct.2 > 0.5 & tumor.vs.spleen.markers$pct.1==0),]
tumor.vs.spleen.markers <- tumor.vs.spleen.markers[order(tumor.vs.spleen.markers$avg_log2FC,
                                                         decreasing = T),]
write.xlsx(list('DE resuslts'=tumor.vs.spleen.markers,
                'mregDC Up'=tumor.vs.spleen.markers[tumor.vs.spleen.markers$avg_log2FC>1&
                                                     tumor.vs.spleen.markers$p_val_adj<0.05,],
                'mregDC Down'=tumor.vs.spleen.markers[tumor.vs.spleen.markers$avg_log2FC <(-1)&
                                                       tumor.vs.spleen.markers$p_val_adj<0.05,]),
           'Outputs/tables/Tumor_vs_spleen_markers.xlsx', rowNames=T)

  # Rank DEGs for further GSEA
tumor.vs.spleen.markers.ranked <- tumor.vs.spleen.markers$avg_log2FC
names(tumor.vs.spleen.markers.ranked) <- rownames(tumor.vs.spleen.markers)
  # Hallmark gene sets, obtained from MSigDB
tumor.vs.spleen.dc3.hallmark <- GSEA(geneList = tumor.vs.spleen.markers.ranked, 
                                     TERM2GENE = read.gmt('h.all.v2023.2.Hs.symbols.gmt'))
gsea.hallmark.show <- tumor.vs.spleen.dc3.hallmark@result
gsea.hallmark.show$Description <- factor(gsea.hallmark.show$Description,
                                         levels = gsea.hallmark.show[order(gsea.hallmark.show$NES),
                                                                     'Description'])
  # Figure S8E, bar plot showing the Hallmark GSEA results
ggplot(gsea.hallmark.show, aes(x=NES, y=Description, fill=-log10(p.adjust)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  scale_fill_gradient(high = '#fff7b7', low ='#e15383')
ggsave('Outputs/figures/GSEA_Tumor_vs_spleen_barplot.pdf',
       width = 6, height = 4)
write.csv(tumor.vs.spleen.dc3.hallmark@result,
          'Outputs/tables/Tumor_vs_spleen_Hallmark.csv')
