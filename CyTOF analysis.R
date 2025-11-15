library(CATALYST)
library(ggpubr)
library(ggrastr)
library(ggplotify)
library(dplyr)
library(ggsci)
library(colorspace)
library(tidyr)
library(paletteer)
library(SingleCellExperiment) # Core object for single-cell data
library(scran)                # For highly variable gene detection
library(scater)               # For plotting and some DR wrappers
library(uwot)                 # For UMAP algorithm
library(Rtsne)
library(viridis)
library(ggridges)
setwd('/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/CyTOF analysis')
  # The objects are processed by Qiu-Chen and are avivaible @ FigShare
sce_all_cells <- readRDS('sce_all.rds')
sce_mdc <- readRDS('sce_cdc_annotated.rds')
  # Figure S1B, all cells t-SNE
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
meta15_lookup <- as.character(meta_mapping[, 'meta15']) # Assuming meta15 is the second column
names(meta15_lookup) <- as.character(meta_mapping[, 'som100']) # Assuming original cluster ID is the first column
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

# Reshape the data from wide to long format for ggplot2's facetting
plot_data_multi <- expression_data %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Marker",
    values_to = "Expression"
  )
plot_data_multi$Marker <- factor(plot_data_multi$Marker,
                                 levels = markers_to_plot)
# Generate the plot for multiple markers
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

# Density Plot (often preferred for continuous data distributions)
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





# DC2hm gating in cDC1 and cDC2 -------------------------------------------


markers_to_compare <- c('CD197_CCR7',
                        "CD127_IL_7Ra",
                        'CD172a_SIRPa',
                        "CD1c",
                        'CD103',
                        'CD141')
# Extract expression data for selected markers
expression_data_compare <- as.data.frame(t(assay(sce_mdc, "exprs")[markers_to_compare, ]))
# Add the cluster information
expression_data_compare$Cluster <- sce_mdc$meta15_cluster
# Convert Cluster to a factor with desired order for plotting, if necessary
expression_data_compare$Cluster <- factor(expression_data_compare$Cluster)

# Reshape the data from wide to long format for ggplot2's facetting
plot_data_compare <- expression_data_compare %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Marker",
    values_to = "Expression"
  )

ggplot(plot_data_compare, aes(x = Expression, y = Cluster, 
                              fill = Cluster, color = Cluster)) +
  geom_density_ridges2(alpha = 0.7) +
  
  # --- THIS IS THE MISSING LINE ---
  # Add this to create a separate panel for each marker
  facet_wrap(~ Marker, ncol = 3, scales = "free_x") +
  
  # Optional: Use a pseudo-log scale if your data has negative values
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


# Digital sorting --------------------------------------------------
# CCR7
ccr7_normalized_expr <- assay(sce_mdc, "exprs")['CD197_CCR7', ]
sce_mdc$CCR7_status <- ifelse(
  ccr7_normalized_expr > 0,
  "CCR7+",
  "CCR7-"
)
sce_mdc$CCR7_status <- factor(sce_mdc$CCR7_status, levels = c("CCR7+", "CCR7-"))
plot_data_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, CCR7_status) %>%
  # Remove any rows where meta15_cluster might be NA (though unlikely with CATALYST)
  filter(!is.na(meta15_cluster))

# Calculate counts for each combination of meta15_cluster and CCR7_status
proportion_counts <- plot_data_proportions %>%
  group_by(CCR7_status, meta15_cluster) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each CCR7_status group
proportion_data <- proportion_counts %>%
  group_by(CCR7_status) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportion_data$CCR7_status <- factor(proportion_data$CCR7_status,
                                      levels = c("CCR7-",
                                                 "CCR7+"
                                                 ))
  # Sorted by cell type
ggplot(proportion_data, aes(x = meta15_cluster, y = Proportion, fill = CCR7_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.5) + # position="fill" makes bars sum to 100%
  labs(
    title = "Proportion CCR7- vs CCR7+ Cells in DC clusters",
    x = "CCR7 Status",
    y = "Proportion of Cells",
    fill = "meta15 Cluster"
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_aaas()
ggsave('Outputs/figures/Proportions_CCR7_cells_in_clusters.pdf',
       width = 5, height = 3)



# IL7R
il7r_normalized_expr <- assay(sce_mdc, "exprs")['CD127_IL_7Ra', ]
sce_mdc$il7r_status <- ifelse(
  il7r_normalized_expr > 0,
  "IL7R+",
  "IL7R-"
)
sce_mdc$il7r_status <- factor(sce_mdc$il7r_status, levels = c("IL7R+", "IL7R-"))
plot_data_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, il7r_status) %>%
  # Remove any rows where meta15_cluster might be NA (though unlikely with CATALYST)
  filter(!is.na(meta15_cluster))

# Calculate counts for each combination of meta15_cluster and il7r_status
proportion_counts <- plot_data_proportions %>%
  group_by(il7r_status, meta15_cluster) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each il7r_status group
proportion_data <- proportion_counts %>%
  group_by(il7r_status) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportion_data$il7r_status <- factor(proportion_data$il7r_status,
                                      levels = c("IL7R-",
                                                 "IL7R+"
                                                 ))

ggplot(proportion_data, aes(x = meta15_cluster, y = Proportion, fill = il7r_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.5) + # position="fill" makes bars sum to 100%
  labs(
    title = "Proportion IL7R- vs IL7R+ Cells in DC clusters",
    x = "IL7R Status",
    y = "Proportion of Cells",
    fill = "meta15 Cluster"
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_aaas()
ggsave('Outputs/figures/Proportions_IL7R_cells_in_clusters.pdf',
       width = 5, height = 3)

# CD1c
cd1c_normalized_expr <- assay(sce_mdc, "exprs")['CD1c', ]
sce_mdc$cd1c_status <- ifelse(
  cd1c_normalized_expr > 0,
  "CD1c+",
  "CD1c-"
)
sce_mdc$cd1c_status <- factor(sce_mdc$cd1c_status, levels = c("CD1c+", "CD1c-"))
plot_data_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, cd1c_status) %>%
  # Remove any rows where meta15_cluster might be NA (though unlikely with CATALYST)
  filter(!is.na(meta15_cluster))

# Calculate counts for each combination of meta15_cluster and cd1c_status
proportion_counts <- plot_data_proportions %>%
  group_by(cd1c_status, meta15_cluster) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each cd1c_status group
proportion_data <- proportion_counts %>%
  group_by(cd1c_status) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportion_data$cd1c_status <- factor(proportion_data$cd1c_status,
                                      levels = c("CD1c-",
                                                 "CD1c+"
                                      ))
ggplot(proportion_data, aes(x = meta15_cluster, y = Proportion, fill = cd1c_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.5) + # position="fill" makes bars sum to 100%
  labs(
    title = "Proportion CD1c- vs CD1c+ Cells in DC clusters",
    x = "CD1c Status",
    y = "Proportion of Cells",
    fill = "meta15 Cluster"
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_aaas()
ggsave('Outputs/figures/Proportions_CD1c_cells_in_clusters.pdf',
       width = 5, height = 3)

# CD103
cd103_normalized_expr <- assay(sce_mdc, "exprs")['CD103', ]
sce_mdc$cd103_status <- ifelse(
  cd103_normalized_expr > 0,
  "CD103+",
  "CD103-"
)
sce_mdc$cd103_status <- factor(sce_mdc$cd103_status, levels = c("CD103+", "CD103-"))
plot_data_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, cd103_status) %>%
  # Remove any rows where meta15_cluster might be NA (though unlikely with CATALYST)
  filter(!is.na(meta15_cluster))

# Calculate counts for each combination of meta15_cluster and cd103_status
proportion_counts <- plot_data_proportions %>%
  group_by(cd103_status, meta15_cluster) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each cd103_status group
proportion_data <- proportion_counts %>%
  group_by(cd103_status) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportion_data$cd103_status <- factor(proportion_data$cd103_status,
                                      levels = c("CD103-",
                                                 "CD103+"
                                      ))
ggplot(proportion_data, aes(x = meta15_cluster, y = Proportion, fill = cd103_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.5) + # position="fill" makes bars sum to 100%
  labs(
    title = "Proportion CD103- vs CD103+ Cells in DC clusters",
    x = "CD103 Status",
    y = "Proportion of Cells",
    fill = "meta15 Cluster"
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_aaas()
ggsave('Outputs/figures/Proportions_CD103_cells_in_clusters.pdf',
       width = 5, height = 3)
# CD141
cd141_normalized_expr <- assay(sce_mdc, "exprs")['CD141', ]
sce_mdc$cd141_status <- ifelse(
  cd141_normalized_expr > 0,
  "CD141+",
  "CD141-"
)
sce_mdc$cd141_status <- factor(sce_mdc$cd141_status, levels = c("CD141+", "CD141-"))
plot_data_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, cd141_status) %>%
  # Remove any rows where meta15_cluster might be NA (though unlikely with CATALYST)
  filter(!is.na(meta15_cluster))

# Calculate counts for each combination of meta15_cluster and cd141_status
proportion_counts <- plot_data_proportions %>%
  group_by(cd141_status, meta15_cluster) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each cd141_status group
proportion_data <- proportion_counts %>%
  group_by(cd141_status) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportion_data$cd141_status <- factor(proportion_data$cd141_status,
                                       levels = c("CD141-",
                                                  "CD141+"
                                       ))
ggplot(proportion_data, aes(x = meta15_cluster, y = Proportion, fill = cd141_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.5) + # position="fill" makes bars sum to 100%
  labs(
    title = "Proportion CD141- vs CD141+ Cells in DC clusters",
    x = "CD141 Status",
    y = "Proportion of Cells",
    fill = "meta15 Cluster"
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_aaas()
ggsave('Outputs/figures/Proportions_CD141_cells_in_clusters.pdf',
       width = 5, height = 3)

# Sirpa
sirpa_normalized_expr <- assay(sce_mdc, "exprs")['CD172a_SIRPa', ]
sce_mdc$sirpa_status <- ifelse(
  sirpa_normalized_expr > 0,
  "SIRPA+",
  "SIRPA-"
)
sce_mdc$sirpa_status <- factor(sce_mdc$sirpa_status, levels = c("SIRPA+", "SIRPA-"))
plot_data_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, sirpa_status) %>%
  # Remove any rows where meta15_cluster might be NA (though unlikely with CATALYST)
  filter(!is.na(meta15_cluster))

# Calculate counts for each combination of meta15_cluster and sirpa_status
proportion_counts <- plot_data_proportions %>%
  group_by(sirpa_status, meta15_cluster) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each sirpa_status group
proportion_data <- proportion_counts %>%
  group_by(sirpa_status) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportion_data$sirpa_status <- factor(proportion_data$sirpa_status,
                                       levels = c("SIRPA-",
                                                  "SIRPA+"
                                       ))
ggplot(proportion_data, aes(x = meta15_cluster, y = Proportion, fill = sirpa_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.5) + # position="fill" makes bars sum to 100%
  labs(
    title = "Proportion SIRPA- vs SIRPA+ Cells in DC clusters",
    x = "SIRPA Status",
    y = "Proportion of Cells",
    fill = "meta15 Cluster"
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_aaas()
ggsave('Outputs/figures/Proportions_SIRPA_cells_in_clusters.pdf',
       width = 5, height = 3)

# CLEC10A
clec10a_normalized_expr <- assay(sce_mdc, "exprs")['CD301_CLEC10A', ]
sce_mdc$clec10a_status <- ifelse(
  clec10a_normalized_expr > 0,
  "CLEC10A+",
  "CLEC10A-"
)
sce_mdc$clec10a_status <- factor(sce_mdc$clec10a_status, levels = c("CLEC10A+", "CLEC10A-"))
plot_data_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, clec10a_status) %>%
  # Remove any rows where meta15_cluster might be NA (though unlikely with CATALYST)
  filter(!is.na(meta15_cluster))

# Calculate counts for each combination of meta15_cluster and clec10a_status
proportion_counts <- plot_data_proportions %>%
  group_by(clec10a_status, meta15_cluster) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each clec10a_status group
proportion_data <- proportion_counts %>%
  group_by(clec10a_status) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
proportion_data$clec10a_status <- factor(proportion_data$clec10a_status,
                                       levels = c("CLEC10A-",
                                                  "CLEC10A+"
                                       ))
ggplot(proportion_data, aes(x = meta15_cluster, y = Proportion, fill = clec10a_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.5) + # position="fill" makes bars sum to 100%
  labs(
    title = "Proportion CLEC10A- vs CLEC10A+ Cells in DC clusters",
    x = "CLEC10A Status",
    y = "Proportion of Cells",
    fill = "meta15 Cluster"
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  scale_fill_aaas()
ggsave('Outputs/figures/Proportions_CLEC10A_cells_in_clusters.pdf',
       width = 5, height = 3)

# Double marker sorting ---------------------------------------------------
  # CCR7 and IL7R
sce_mdc$ccr7_il7r_status <- factor(paste0(sce_mdc$CCR7_status, sce_mdc$il7r_status),
                                levels = c("CCR7-IL7R-","CCR7+IL7R-","CCR7-IL7R+", "CCR7+IL7R+"))
message("Distribution of the new ccr7_il7r_status attribute:")
print(table(sce_mdc$ccr7_il7r_status, useNA = "always"))

# Create a data frame with the relevant columns
plot_data_double_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, ccr7_il7r_status) %>%
  # Remove any rows where meta15_cluster or ccr7_il7r_status might be NA
  filter(!is.na(meta15_cluster) & !is.na(ccr7_il7r_status))

# Calculate counts for each combination of meta15_cluster and ccr7_il7r_status
double_proportion_counts <- plot_data_double_proportions %>%
  group_by(meta15_cluster, ccr7_il7r_status) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each meta15_cluster group
# This ensures each meta15_cluster bar sums to 100%
double_proportion_data <- double_proportion_counts %>%
  group_by(meta15_cluster) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

double_proportion_data$meta15_cluster <- factor(double_proportion_data$meta15_cluster, levels = order_by_main_subset)

ggplot(double_proportion_data, aes(x = meta15_cluster, y = Proportion, fill = ccr7_il7r_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3) + # position="fill" makes bars sum to 100%
  labs(
    title = "Composition of DC clusters by CCR7 & IL7R Expression",
    x = "meta15 Cluster",
    y = "Proportion of Cells",
    fill = "CCR7/IL7R Status" # Legend title
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), # Rotate & bold cluster names
    legend.title = element_text(face = "bold"),
    legend.position = "right" # Place legend to the right
  ) +
  scale_fill_npg()
ggsave('Outputs/figures/Proportions_IL7R_CCR7cells_in_clusters.pdf')

# CD141 and CD1c
sce_mdc$ccr7_il7r_status <- factor(paste0(sce_mdc$CCR7_status, sce_mdc$il7r_status),
                                   levels = c("CCR7-IL7R-","CCR7+IL7R-","CCR7-IL7R+", "CCR7+IL7R+"))
message("Distribution of the new ccr7_il7r_status attribute:")
print(table(sce_mdc$ccr7_il7r_status, useNA = "always"))

# Create a data frame with the relevant columns
plot_data_double_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, ccr7_il7r_status) %>%
  # Remove any rows where meta15_cluster or ccr7_il7r_status might be NA
  filter(!is.na(meta15_cluster) & !is.na(ccr7_il7r_status))

# Calculate counts for each combination of meta15_cluster and ccr7_il7r_status
double_proportion_counts <- plot_data_double_proportions %>%
  group_by(meta15_cluster, ccr7_il7r_status) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each meta15_cluster group
# This ensures each meta15_cluster bar sums to 100%
double_proportion_data <- double_proportion_counts %>%
  group_by(meta15_cluster) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

double_proportion_data$meta15_cluster <- factor(double_proportion_data$meta15_cluster, levels = order_by_main_subset)

ggplot(double_proportion_data, aes(x = meta15_cluster, y = Proportion, fill = ccr7_il7r_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3) + # position="fill" makes bars sum to 100%
  labs(
    title = "Composition of DC clusters by CCR7 & IL7R Expression",
    x = "meta15 Cluster",
    y = "Proportion of Cells",
    fill = "CCR7/IL7R Status" # Legend title
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), # Rotate & bold cluster names
    legend.title = element_text(face = "bold"),
    legend.position = "right" # Place legend to the right
  ) +
  scale_fill_npg()
ggsave('Outputs/figures/Proportions_IL7R_CCR7cells_in_clusters.pdf')

# CD1c and CD103
sce_mdc$ccr7_il7r_status <- factor(paste0(sce_mdc$CCR7_status, sce_mdc$il7r_status),
                                   levels = c("CCR7-IL7R-","CCR7+IL7R-","CCR7-IL7R+", "CCR7+IL7R+"))
message("Distribution of the new ccr7_il7r_status attribute:")
print(table(sce_mdc$ccr7_il7r_status, useNA = "always"))

# Create a data frame with the relevant columns
plot_data_double_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, ccr7_il7r_status) %>%
  # Remove any rows where meta15_cluster or ccr7_il7r_status might be NA
  filter(!is.na(meta15_cluster) & !is.na(ccr7_il7r_status))

# Calculate counts for each combination of meta15_cluster and ccr7_il7r_status
double_proportion_counts <- plot_data_double_proportions %>%
  group_by(meta15_cluster, ccr7_il7r_status) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each meta15_cluster group
# This ensures each meta15_cluster bar sums to 100%
double_proportion_data <- double_proportion_counts %>%
  group_by(meta15_cluster) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

double_proportion_data$meta15_cluster <- factor(double_proportion_data$meta15_cluster, levels = order_by_main_subset)

ggplot(double_proportion_data, aes(x = meta15_cluster, y = Proportion, fill = ccr7_il7r_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3) + # position="fill" makes bars sum to 100%
  labs(
    title = "Composition of DC clusters by CCR7 & IL7R Expression",
    x = "meta15 Cluster",
    y = "Proportion of Cells",
    fill = "CCR7/IL7R Status" # Legend title
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), # Rotate & bold cluster names
    legend.title = element_text(face = "bold"),
    legend.position = "right" # Place legend to the right
  ) +
  scale_fill_npg()
ggsave('Outputs/figures/Proportions_IL7R_CCR7cells_in_clusters.pdf')

# CD141 and CD1c
sce_mdc$ccr7_il7r_status <- factor(paste0(sce_mdc$CCR7_status, sce_mdc$il7r_status),
                                   levels = c("CCR7-IL7R-","CCR7+IL7R-","CCR7-IL7R+", "CCR7+IL7R+"))
message("Distribution of the new ccr7_il7r_status attribute:")
print(table(sce_mdc$ccr7_il7r_status, useNA = "always"))

# Create a data frame with the relevant columns
plot_data_double_proportions <- as.data.frame(colData(sce_mdc)) %>%
  dplyr::select(meta15_cluster, ccr7_il7r_status) %>%
  # Remove any rows where meta15_cluster or ccr7_il7r_status might be NA
  filter(!is.na(meta15_cluster) & !is.na(ccr7_il7r_status))

# Calculate counts for each combination of meta15_cluster and ccr7_il7r_status
double_proportion_counts <- plot_data_double_proportions %>%
  group_by(meta15_cluster, ccr7_il7r_status) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate proportions within each meta15_cluster group
# This ensures each meta15_cluster bar sums to 100%
double_proportion_data <- double_proportion_counts %>%
  group_by(meta15_cluster) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

double_proportion_data$meta15_cluster <- factor(double_proportion_data$meta15_cluster, levels = order_by_main_subset)

ggplot(double_proportion_data, aes(x = meta15_cluster, y = Proportion, fill = ccr7_il7r_status)) +
  geom_col(position = "fill", color = "black", linewidth = 0.3) + # position="fill" makes bars sum to 100%
  labs(
    title = "Composition of DC clusters by CCR7 & IL7R Expression",
    x = "meta15 Cluster",
    y = "Proportion of Cells",
    fill = "CCR7/IL7R Status" # Legend title
  ) +
  scale_y_continuous(labels = scales::percent) + # Format y-axis as percentage
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), # Rotate & bold cluster names
    legend.title = element_text(face = "bold"),
    legend.position = "right" # Place legend to the right
  ) +
  scale_fill_npg()
ggsave('Outputs/figures/Proportions_IL7R_CCR7cells_in_clusters.pdf')


