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
  # Read Seurat object which is available @ FigShare
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
DotPlot(combined.object, features = c('CRLF2','IL7R',
                                          'HLA-DPA1','HLA-DRA','HLA-DRB5','HLA-DPB1',
                                          'ISG15','IFI27','IFI44L','MX1','IRF7',
                                          'UQCR10','COX5B','NDUFS3','CYB5R3','OXA1L'),
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
spleen.vs.tumor.markers <- FindMarkers(combined.object,
                                       ident.2 = 'Tumor', ident.1 = 'HC_spleen',
                                       group.by = 'binary.group')
  # Remove false-positive DEGs which may due to gene name discrepancy
spleen.vs.tumor.markers <- spleen.vs.tumor.markers[!(spleen.vs.tumor.markers$pct.1 > 0.5 & spleen.vs.tumor.markers$pct.2==0) |
                                                     !(spleen.vs.tumor.markers$pct.2 > 0.5 & spleen.vs.tumor.markers$pct.1==0),]
spleen.vs.tumor.markers <- spleen.vs.tumor.markers[order(spleen.vs.tumor.markers$avg_log2FC,
                                                         decreasing = T),]
write.xlsx(list('DE resuslts'=spleen.vs.tumor.markers,
                'Spleen Up'=spleen.vs.tumor.markers[spleen.vs.tumor.markers$avg_log2FC>1&
                                                     spleen.vs.tumor.markers$p_val_adj<0.05,],
                'Spleen Down'=spleen.vs.tumor.markers[spleen.vs.tumor.markers$avg_log2FC <(-1)&
                                                       spleen.vs.tumor.markers$p_val_adj<0.05,]),
           'Outputs/tables/Spleen_vs_tumor_markers.xlsx', rowNames=T)

  # Rank DEGs for further GSEA
spleen.vs.tumor.markers.ranked <- spleen.vs.tumor.markers$avg_log2FC
names(spleen.vs.tumor.markers.ranked) <- rownames(spleen.vs.tumor.markers)
  # Hallmark gene sets, obtained from MSigDB
spleen.vs.tumor.dc3.hallmark <- GSEA(geneList = spleen.vs.tumor.markers.ranked, 
                                     TERM2GENE = read.gmt('h.all.v2023.2.Hs.symbols.gmt'))
gsea.hallmark.show <- spleen.vs.tumor.dc3.hallmark@result
gsea.hallmark.show$Description <- factor(gsea.hallmark.show$Description,
                                         levels = gsea.hallmark.show[order(gsea.hallmark.show$NES),
                                                                     'Description'])
  # Figure S8E, bar plot showing the Hallmark GSEA results
ggplot(gsea.hallmark.show, aes(x=NES, y=Description, fill=-log10(p.adjust)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  scale_fill_gradient(high = '#fff7b7', low ='#e15383')
ggsave('Outputs/figures/GSEA_spleen_vs_tumor_dc3_barplot.pdf',
       width = 6, height = 4)
  # Figure S8G, GSEA of selected pathways
gsea.selected.tslp.signaling <- GSEA(spleen.vs.tumor.markers.ranked, 
                               TERM2GENE = rbind(read.gmt('WP_THYMIC_STROMAL_LYMPHOPOIETIN_TSLP_SIGNALING_PATHWAY.v2023.2.Hs.gmt'),
                                                 read.gmt('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING.v2023.2.Hs.gmt'),
                                                 read.gmt('LINDSTEDT_DENDRITIC_CELL_MATURATION_B.v2023.2.Hs.gmt')),pvalueCutoff = 1,minGSSize = 1)
gsea.selected.tslp.signaling@result$Description <- factor(gsea.selected.tslp.signaling@result$Description,
                                                    levels = c("LINDSTEDT_DENDRITIC_CELL_MATURATION_B",
                                                               "WP_THYMIC_STROMAL_LYMPHOPOIETIN_TSLP_SIGNALING_PATHWAY",
                                                               "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"))
gseaplot2(gsea.selected.tslp.signaling, geneSetID = c(1,2,3),
          color = c('#FFBE0B','#ff5e5b','#00BBF9'))
ggsave('Outputs/figures/GSEA_spleen_vs_tumor.pdf',
       width = 4, height = 3.5)
write.xlsx(list('Hallmark'=spleen.vs.tumor.dc3.hallmark@result,
                'Mature_IFN_TSLP'=gsea.selected.tslp.signaling@result),
           'Outputs/tables/Spleen_vs_Tumor_mature_DC_vs_DC2hm_GSEA.xlsx',
           rowNames=T)
