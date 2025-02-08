library(CATALYST)
library(ggpubr)
library(ggrastr)
library(ggplotify)
setwd('/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/CyTOF analysis')
  # The objects are processed by Qiu-Chen and are avivaible @ FigShare
sce_all_cells <- readRDS('sce_all.rds')
sce_mdc <- readRDS('sce_cdc_annotated.rds')
  # Figure S1B, all cells t-SNE
tsne_all <- plotDR(sce_all_cells, dr = 'TSNE', 'celltype2')
tsne_all <- rasterize(tsne_all, dpi = 300)
tsne_all
ggsave('Outputs/figures/CyTOF_tsne_alls.pdf', width = 4, height = 3)
  # Figure S1C, All cell heatmap
proteins <- c("CD63_LAMP_3","CD11c","CD141", "CD1c","CD172a_SIRPa","FceRIa",
              "CD197_CCR7","HLA_DR","CD14","CD16", "HLA_A_B_C","CD127_IL_7Ra",
              "CD3","CD4", "CD11b","CD303_BDCA_2","CD304_BDCA4","CD66b",
              "CD117_c_kit","CD86","CD19","CD335_NKp46","CD123_IL_3Ra","CD34")
as.ggplot(plotExprHeatmap(sce_all_cells, features = proteins, 
                  by = "cluster_id", k = "celltype2"))
ggsave('Outputs/figures/All_cell_exp_heatmap.pdf')
  # Figure 1B t-SNE plot of cyTOF cDCs
tsne <- plotDR(sce_mdc, dr = 'TSNE', 'meta15')+
  scale_color_manual(values = c('#C5E524','#38C2E5','#E679C5'))
tsne <- rasterize(tsne, layer='Point', dpi=300)
tsne
ggsave('Outputs/figures/CyTOF_tsne_cDCs.pdf', width = 4, height = 3)
  # Figure 1C, cDC marker heatmap
as.ggplot(plotExprHeatmap(sce_mdc,  features = c("CD11c","HLA_DR","CD1c","CD172a_SIRPa",
                                                 "FceRIa","CD103","CD141","CD127_IL_7Ra",
                                                 'CD95_Fas',
                                                 "CD197_CCR7","CD83","CD40_TNFRSF5",
                                                 "CD273_PD_L2","CD86","CD274_PD_L1"),
                          by = 'cluster_id', k = 'meta15', row_clust = F,
                          col_clust = F))
ggsave('Outputs/figures/Heatmap_cDCs_markers.pdf', width=8, height=3)

# Figure S1,feature plots
proteins_all <- rownames(sce_mdc@assays@data$counts)[!rownames(sce_mdc@assays@data$counts) %in% c("89Y","102Pd", "103Rh", "104_barcode" ,  "105_barcode",
                                                                                                  "106_barcode" ,  "107Ag", "108_barcode",   "109Ag",
                                                                                                  "110_barcode","111Cd", "112Cd", "113In", "114Cd",
                                                                                                  "116Cd", "139La", "140Ce","DNA1",  "192Pt", "DNA2",  "cisplatin",
                                                                                                  "195Pt", "196Pt")]
for (i in seq(1,length(proteins_all))) {
  p <- plotDR(sce_mdc, 'TSNE', proteins_all[[i]])+
    scale_color_gradientn(colors = c('#2A9FF3','white' ,'#FF2B72'))
  p <- rasterize(p, layers = 'Point', dpi = 300)
  pdf(paste0('Outputs/figures/',proteins_all[[i]],'_featureplot.pdf'),
      width = 4,height = 4)
  plot(p)
  dev.off()
}

