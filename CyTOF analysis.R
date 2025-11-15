library(CATALYST)
library(ggpubr)
library(ggrastr)
library(ggplotify)
library(dplyr)
library(ggsci)
library(colorspace)
library(tidyr)
library(paletteer)
library(SingleCellExperiment) 
library(scran)                
library(scater)               
library(uwot)                 
library(Rtsne)
library(viridis)
library(ggridges)
setwd('CyTOF analysis')
  # The objects are processed by Qiu-Chen and are avivaible @ FigShare
sce_all_cells <- readRDS('sce_all.rds')
sce_mdc <- readRDS('sce_cdc_annotated.rds')
  # Foundamental visualizations below
tsne_all <- plotDR(sce_all_cells, dr = 'TSNE', 'celltype2')
tsne_all <- rasterize(tsne_all, dpi = 300)
tsne_all
ggsave('Outputs/figures/CyTOF_tsne_alls.pdf', width = 8, height = 6)
tsne_all_ccr7 <- plotDR(sce_all_cells, dr = 'TSNE', "CD197_CCR7")+
  scale_color_viridis_c()
tsne_all_ccr7 <- rasterize(tsne_all_ccr7, dpi = 300)
tsne_all_ccr7 
ggsave('Outputs/figures/CyTOF_tsne_alls_CCR7.pdf', width = 8, height = 6)

tsne_all_il7r <- plotDR(sce_all_cells, dr = 'TSNE', "CD127_IL_7Ra")+
  scale_color_viridis_c()
tsne_all_il7r <- rasterize(tsne_all_il7r, dpi = 300)
tsne_all_il7r
ggsave('Outputs/figures/CyTOF_tsne_alls_IL7R.pdf', width = 8, height = 6)

tsne_all_cd1c <- plotDR(sce_all_cells, dr = 'TSNE', "CD1c")+
  scale_color_viridis_c()
tsne_all_cd1c <- rasterize(tsne_all_cd1c, dpi = 300)
tsne_all_cd1c
ggsave('Outputs/figures/CyTOF_tsne_alls_CD1c.pdf', width = 8, height = 6)

tsne_all_hladr <- plotDR(sce_all_cells, dr = 'TSNE', "HLA_DR")+
  scale_color_viridis_c()
tsne_all_hladr <- rasterize(tsne_all_hladr, dpi = 300)
tsne_all_hladr
ggsave('Outputs/figures/CyTOF_tsne_alls_HLA-DR.pdf', width = 8, height = 6)

tsne_all_cd141 <- plotDR(sce_all_cells, dr = 'TSNE', "CD141")+
  scale_color_viridis_c()
tsne_all_cd141 <- rasterize(tsne_all_cd141, dpi = 300)
tsne_all_cd141
ggsave('Outputs/figures/CyTOF_tsne_alls_CD141.pdf', width = 8, height = 6)

tsne_all_cd3 <- plotDR(sce_all_cells, dr = 'TSNE', "CD3")+
  scale_color_viridis_c()
tsne_all_cd3 <- rasterize(tsne_all_cd3, dpi = 300)
tsne_all_cd3
ggsave('Outputs/figures/CyTOF_tsne_alls_CD3.pdf', width = 8, height = 6)

tsne_all_cd11c <- plotDR(sce_all_cells, dr = 'TSNE', "CD11c")+
  scale_color_viridis_c()
tsne_all_cd11c <- rasterize(tsne_all_cd11c, dpi = 300)
tsne_all_cd11c
ggsave('Outputs/figures/CyTOF_tsne_alls_CD11c.pdf', width = 8, height = 6)

tsne_all_cd103 <- plotDR(sce_all_cells, dr = 'TSNE', "CD103")+
  scale_color_viridis_c()
tsne_all_cd103 <- rasterize(tsne_all_cd103, dpi = 300)
tsne_all_cd103
ggsave('Outputs/figures/CyTOF_tsne_alls_CD103.pdf', width = 8, height = 6)
tsne_all_cd83 <- plotDR(sce_all_cells, dr = 'TSNE', "CD83")+
  scale_color_viridis_c()
tsne_all_cd83 <- rasterize(tsne_all_cd83, dpi = 300)
tsne_all_cd83
ggsave('Outputs/figures/CyTOF_tsne_alls_CD83.pdf', width = 8, height = 6)
proteins <- c("CD63_LAMP_3","CD11c","CD141", "CD1c","CD172a_SIRPa","FceRIa",
              "CD197_CCR7","HLA_DR","CD14","CD16", "HLA_A_B_C","CD127_IL_7Ra",
              "CD3","CD4", "CD11b","CD303_BDCA_2","CD304_BDCA4","CD66b",
              "CD117_c_kit","CD86","CD19","CD335_NKp46","CD123_IL_3Ra","CD34")
as.ggplot(plotExprHeatmap(sce_all_cells, features = proteins, 
                  by = "cluster_id", k = "celltype2"))
ggsave('Outputs/figures/All_cell_exp_heatmap.pdf')
tsne <- plotDR(sce_mdc, dr = 'TSNE', 'meta15')+
  scale_color_manual(values = c('#C5E524','#38C2E5','#E679C5'))
tsne <- rasterize(tsne, layer='Point', dpi=300)
tsne
ggsave('Outputs/figures/CyTOF_tsne_cDCs.pdf', width = 4, height = 3)

as.ggplot(plotExprHeatmap(sce_mdc,  features = c("CD11c","HLA_DR","CD1c","CD172a_SIRPa",
                                                 "FceRIa","CD103","CD141","CD127_IL_7Ra",
                                                 'CD95_Fas',
                                                 "CD197_CCR7","CD83","CD40_TNFRSF5",
                                                 "CD273_PD_L2","CD86","CD274_PD_L1"),
                          by = 'cluster_id', k = 'meta15', row_clust = F,
                          col_clust = F))
ggsave('Outputs/figures/Heatmap_cDCs_markers.pdf', width=8, height=3)

proteins_all <- rownames(sce_mdc@assays@data$counts)[!rownames(sce_mdc@assays@data$counts) %in% c("89Y","102Pd", "103Rh", "104_barcode" ,  "105_barcode",
                                                                                                  "106_barcode" ,  "107Ag", "108_barcode",   "109Ag",
                                                                                                  "110_barcode","111Cd", "112Cd", "113In", "114Cd",
                                                                                                  "116Cd", "139La", "140Ce","DNA1",  "192Pt", "DNA2",  "cisplatin",
                                                                                                  "195Pt", "196Pt")]


for (i in seq(1,length(proteins_all))) {
  p <- plotDR(sce_mdc, 'TSNE', proteins_all[[i]])+
    scale_color_viridis_c()
  p <- rasterize(p, layers = 'Point', dpi = 300)
  pdf(paste0('Outputs/figures/',proteins_all[[i]],'_featureplot.pdf'),
      width = 4,height = 4)
  plot(p)
  dev.off()
}


# Compare cDC1/cDC2 marker in DC2hm ---------------------------------------

expression_data <- assay(sce_mdc)[genes_of_interest, ]

meta_mapping <- metadata(sce_mdc)$cluster_codes
meta15_lookup <- as.character(meta_mapping[, 'meta15']) 
names(meta15_lookup) <- as.character(meta_mapping[, 'som100']) 
sce_mdc$meta15_cluster <- meta15_lookup[as.character(sce_mdc$cluster_id)]


# Define the markers you want to plot
markers_to_plot <- c("CD1c", "CD301_CLEC10A",'CD172a_SIRPa',
                            "CD141",'CD103',
                            'CD197_CCR7',
                            'HLA_DR',
                            "CD127_IL_7Ra",
                     "CD205_DEC_205")

# Extract expression data for selected markers
expression_data <- as.data.frame(t(assay(sce_mdc, "exprs")[markers_to_plot, ]))
expression_data$Cluster <- sce_mdc$meta15_cluster

plot_data_multi <- expression_data %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Marker",
    values_to = "Expression"
  )
plot_data_multi$Marker <- factor(plot_data_multi$Marker,
                                 levels = markers_to_plot)
ggplot(plot_data_multi, aes(x = Cluster, y = Expression, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.7,scale = "width") +
  geom_boxplot(width = 0.03, outlier.shape = NA, fill = "white", color = "black") +
  facet_wrap(~ Marker, scales = "free_y", ncol = 4) + # Separate plots for each marker, free y-scales
  labs(
    x = "meta15 Cluster",
    y = "Expression (Normalized)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold") # Bold facet titles (marker names)
  ) +
  scale_fill_manual(values = c('#C5E524','#38C2E5','#E679C5'))
ggsave('Outputs/figures/Marker_expression_violins.pdf',
       width = 10, height = 4)

# DC2hm backgating-like analysis ------------------------------------------

# Identify cells belonging to the "DC2hm" cluster
cells_in_DC2hm_cluster_idx <- which(sce_mdc$meta15_cluster =='DC2hm')

# Create a new SCE object containing only cells from the "DC2hm" cluster
sce_dc2hm <- sce_mdc[, cells_in_DC2hm_cluster_idx]
markers_to_plot_subset <- c("CD1c", "CD301_CLEC10A",'CD172a_SIRPa',
                            "CD141",'CD103')
expression_data_subset <- as.data.frame(t(assay(sce_mdc, "exprs")[markers_to_plot_subset, ]))
expression_data_subset$Clusters <- sce_mdc$meta15_cluster
plot_data_subset <- expression_data_subset %>%
  pivot_longer(
    cols = Clusters, # Take all columns since we don't have a 'Cluster' column for faceting now
    names_to = "Marker",
    values_to = "Expression"
  )
plot_data_subset$Marker <- factor(plot_data_subset$Marker,
                                  levels =  c("CD1c", "CD301_CLEC10A",'CD172a_SIRPa',
                                              "CD141",'CD103'))

ggplot(plot_data_subset, aes(x = Expression)) +
  geom_density(fill = Clusters, color = Clusters, alpha = 0.7) +
  facet_wrap(~ Marker, scales = "free_y", ncol = 1) + # <--- Changed to scales = "free"
  labs(
    title = paste0("Expression of cDC1/2 marker within the DC2hm"),
    x = "Expression (Normalized)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )
ggsave('Outputs/figures/Density_cDC1_and_cDC2_marker_in_DC2hm.pdf', 
       width = 6, height = 6)


# DC subpopulation marker expression  -------------------------------------------


markers_to_compare <- c('CD197_CCR7',
                        "CD127_IL_7Ra",
                        'CD172a_SIRPa',
                        "CD1c",
                        'CD103',
                        'CD141')
expression_data_compare <- as.data.frame(t(assay(sce_mdc, "exprs")[markers_to_compare, ]))
expression_data_compare$Cluster <- sce_mdc$meta15_cluster
expression_data_compare$Cluster <- factor(expression_data_compare$Cluster)
plot_data_compare <- expression_data_compare %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Marker",
    values_to = "Expression"
  )

ggplot(plot_data_compare, aes(x = Expression, y = Cluster, 
                              fill = Cluster, color = Cluster)) +
  geom_density_ridges2(alpha = 0.7) +
  
  facet_wrap(~ Marker, ncol = 3, scales = "free_x") +
  # scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  
  labs(
    title = "Marker Expression Across Clusters",
    x = "Expression Level",
    y = "Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none" # Hide legend as y-axis is labeled
  ) +
  scale_fill_manual(values = c('#C5E524', '#38C2E5', '#E679C5')) +
  scale_color_manual(values = c('#C5E524', '#38C2E5', '#E679C5'))

ggsave('Outputs/figures/Marker_expression_ridges.pdf',
       width = 10, height = 4)


