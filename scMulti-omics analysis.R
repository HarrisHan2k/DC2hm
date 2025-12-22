library(Seurat)
library(viridis)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(Signac)
library(ComplexHeatmap)
library(ggVennDiagram)
library(circlize)
library(ggplotify)
library(clusterProfiler)
library(openxlsx)
library(ggpubr)
library(future)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(ggrastr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
library(ggsci)
set.seed(0822)
setwd('scMulti-omics analysis')
# The preprocessing and clustering analysis was performed by Qiu-Chen
# Start with this processed object which could be obtained @ GEO: GSE267255
# Clustering and marker analysis on cDCs ----------------------------------
multiomic.object <- readRDS('hu_spleen_multiomics.rds')
# Update nomenclature
multiomic.object <- RenameIdents(multiomic.object, c('mDC1'='cDC1',
                                                     'mDC2'='cDC2',
                                                     'mDC2_stressed'='cDC2_UPR',
                                                     'mDC2_activated'='DC2pre-hm',
                                                     'DC3'='DC2hm', # Note that DC3 is an in-house code for DC2hm, never confuse with CD14+ DCs
                                                     'cycling_B'='Proliferating',
                                                     'NKs'='NK',
                                                     'Macrophages'='Mono/Macro',
                                                     'Monocytes'='Mono/Macro'))
multiomic.object$celltype_hh <- Idents(multiomic.object)
levels(multiomic.object) <- c("cDC1","cDC2",
                              "DC2pre-hm",
                              "DC2hm",
                              "cDC2_UPR",
                              'AXL_DC',
                              "Mono/Macro",
                              "B",
                              "NK",
                              "Proliferating",
                              "CD34Progenitor")
DotPlot(multiomic.object, assay = 'RNA', features = c('FLT3','CLEC9A','XCR1','BATF3',
                                                      'CD1C','ITGAX','CLEC10A',
                                                      'CCR7','IL7R','LAMP3',
                                                      'HSPA6','DNAJB1',
                                                      'AXL','SIGLEC6',
                                                      'CD14','LYZ','CD68','FCGR3A',
                                                      'CD19','MS4A1','CD79A',
                                                      'NKG7','NCAM1',
                                                      'MKI67','CD34'),
        group.by = 'celltype_hh')+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('Outputs/figures/Dotlpot_multi-omic_markers.pdf', 
       width = 10, height = 4)
rasterize(DimPlot(multiomic.object, 
                  cols = c("#C5E524", "#38C2E5","#384C94","#E679C5",
                           "#FF7F00","#FF6961","#FEDD00","#FFA07A",
                           "#009F6B","#B3B3B3","#DAA520"),pt.size = .1,
                  reduction = 'umap'),
          dpi = 300)
ggsave('Outputs/figures/Sup_fig_UMAP_peaks_all.pdf', width = 6, height = 4)
multiomic.object.cdcs <- multiomic.object[,Idents(multiomic.object) %in% c('cDC1', 
                                                                      'cDC2',
                                                                      'DC2pre-hm',
                                                                      'DC2hm')]
levels(multiomic.object.cdcs) <- c('cDC1', 'cDC2','DC2pre-hm', 'DC2hm')
multiomic.object.cdcs$transferred_label <- Idents(multiomic.object.cdcs)
multiomic.object.cdcs$transferred_label <- factor(multiomic.object.cdcs$transferred_label,
                                                  levels = c('cDC1', 'cDC2',
                                                             'DC2pre-hm', 'DC2hm'))
DefaultAssay(multiomic.object.cdcs) <- 'RNA'
multiomic.object.cdcs <- NormalizeData(multiomic.object.cdcs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
multiomic.object.cdcs <- RunUMAP(multiomic.object.cdcs, dims = 1:5)
multiomic.object.cdcs <- FindNeighbors(multiomic.object.cdcs) %>% FindClusters(resolution = .5)
multiomic.object.cdcs <- RenameIdents(multiomic.object.cdcs,
                                      '0'='cDC2','1'='cDC2','2'='DC2hm',
                                      '3'='DC2pre-hm','4'='cDC2','5'='cDC1',
                                      '6'='cDC2','7'='cDC2')
multiomic.object.cdcs$celltype_snRNA <- Idents(multiomic.object.cdcs)
levels(multiomic.object.cdcs) <- c('cDC1', 'cDC2','DC2pre-hm', 'DC2hm')
DotPlot(multiomic.object.cdcs,
        features = read.csv('../scRNA-seq analysis/Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Outputs/figures/Dotplot_snRNA-seq_markers.pdf', width = 10, height = 4)

DefaultAssay(multiomic.object.cdcs) <- 'peaks'
multiomic.object.cdcs <- RunTFIDF(multiomic.object.cdcs) %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction='lsi', dims=2:10)
rasterize(DimPlot(multiomic.object.cdcs, group.by = 'transferred_label',
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5'),
                  pt.size = 0.1),dpi = 300)
ggsave('Outputs/figures/UMAP_peaks_reduction.pdf', width = 4.5, height = 3)
DefaultAssay(multiomic.object.cdcs) <- 'RNA'
multiomic.object.cdcs <- NormalizeData(multiomic.object.cdcs) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
multiomic.object.cdcs <- RunUMAP(multiomic.object.cdcs, dims = 1:5)
rasterize(DimPlot(multiomic.object.cdcs, group.by = 'transferred_label',
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5'),
                  pt.size = 0.1),dpi = 300)
ggsave('Outputs/figures/UMAP_RNA_reduction.pdf', width = 4.5, height = 3)
multiomic.object.cdcs$celltype_snRNA <- factor(multiomic.object.cdcs$celltype_snRNA,
                                               levels = c('cDC1', 'cDC2',
                                                          'DC2pre-hm', 'DC2hm'))
rasterize(DimPlot(multiomic.object.cdcs, group.by = 'celltype_snRNA',
                    cols = c('#C5E524', '#38C2E5','#384C94','#E679C5'),
                    pt.size = 0.1),dpi = 300)
ggsave('Outputs/figures/UMAP_RNA_reduction_snRNA_clustering.pdf', width = 4.5, height = 3)
proportion.table <- as.data.frame(prop.table(table(multiomic.object.cdcs$celltype_snRNA,
                                     multiomic.object.cdcs$transferred_label),
                               margin = 2))
colnames(proportion.table) <- c('celltype_snRNA', 'transferred_label',
                               'Freq')
ggplot(proportion.table,aes(
           x=transferred_label,y=Freq,fill=celltype_snRNA))+
  geom_bar(stat = 'identity', position = 'stack')+
  theme_classic()+
  scale_fill_manual(values = c('#C5E524', '#38C2E5','#384C94','#E679C5'))
  # Figure 1N, feature peaks for each cDC cluster
DefaultAssay(multiomic.object.cdcs) <- 'peaks'
markers <- FindAllMarkers(multiomic.object.cdcs,
                          assay = 'peaks',test.use = 'LR',
                          latent.vars = 'nCount_peaks') 

# Correlation with scRNA-seq ----------------------------------------------
d6.scrna <- readRDS('Dnr6_cDCs.rds') # Could be found @ GEO
DotPlot(d6.scrna,
        features = read.csv('Highlighted_gene_selected.csv')$Gene)+
  scale_color_viridis_c()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Outputs/figures/Dotplot_scRNA-seq_markers.pdf', width = 10, height = 4)
d6.scrna$technique <- rep('Donor 6 scRNA-seq', ncol(d6.scrna))
Idents(multiomic.object.cdcs) <- multiomic.object.cdcs$transferred_label
multiomic.object.cdcs$technique <- rep('Donor 6 snRNA-seq', 
                                       ncol(multiomic.object.cdcs))
merged.object.snsc <- merge(d6.scrna,
                            multiomic.object.cdcs)
for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = FALSE)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], 
                                           selection.method = "vst", 
                                           nfeatures = 2000, verbose = FALSE)
}

integ_features <- SelectIntegrationFeatures(object.list = object.list)

integ_anchors <- FindIntegrationAnchors(
  object.list = object.list,
  anchor.features = integ_features,
  normalization.method = "LogNormalize" # Specify the method used
)

integrated_object <- IntegrateData(anchorset = integ_anchors)
DefaultAssay(integrated_object) <- "integrated"

integrated_object <- ScaleData(integrated_object, verbose = FALSE)
integrated_object <- RunPCA(integrated_object, verbose = FALSE)
integrated_object <- RunUMAP(integrated_object, dims = 1:10)
integrated_object$technique <- replace_na(integrated_object$technique,
                                          'Donor 6 snRNA-seq')
levels(integrated_object) <- c('cDC1', 'cDC2',
                               'DC2pre-hm', 'DC2hm')
# Visualize to confirm successful integration
p1 <- rasterize(DimPlot(integrated_object, reduction = "umap",
              group.by = "technique", pt.size = .1)+
  scale_color_aaas(),dpi = 300)
p2 <- rasterize(DimPlot(integrated_object, reduction = "umap",
              pt.size = .1,
              cols = c('#C5E524', '#38C2E5','#384C94','#E679C5')),dpi = 300)
p1 + p2
ggsave('Outputs/figures/Integrated_UMAP.pdf',
       width = 14, height = 6)

rasterize(DimPlot(integrated_object, reduction = "umap",
                  split.by = 'technique', pt.size = .1,
                  group.by = 'celltype_int',
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5')),
          dpi = 300)
ggsave('Outputs/figures/Integrated_UMAP_split_plot.pdf', width = 14, height = 6)
metadata <- integrated_object@meta.data

# Let's check the first few rows to confirm
head(metadata)
proportions_df <- metadata %>%
  # Group by the technology and then by the cell type
  group_by(technique, celltype_int) %>%
  # Count the number of cells in each group
  summarise(n = n(), .groups = 'drop') %>%
  # Ungroup and then group again only by technology to find the total cells per tech
  group_by(technique) %>%
  # Calculate the frequency (percentage)
  mutate(percentage = n / sum(n) * 100)

# View the resulting data frame
print(proportions_df)
ggplot(proportions_df, aes(x = celltype_int, y = percentage, fill = technique)) +
  # Use geom_col() or geom_bar(stat="identity") to plot the values directly
  geom_col(position = "dodge", alpha = 0.8, color = "black", linewidth = 0.5) +
  # Use a nice color scheme
  scale_fill_manual(values = c('#344184',
                               '#EB141B')) +
  labs(
    title = "Cell Type Proportions: scRNA-seq vs. snRNA-seq",
    x = "Cell Type",
    y = "Percentage of Total Cells (%)",
    fill = "Technology"
  ) +
  # Improve readability
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )
ggsave('Outputs/figures/Comparative_proportion_analysis.pdf',
       width = 6, height = 4)
integrated_object$celltype_technique <- paste(integrated_object$technique,
                                               Idents(integrated_object))
integrated_object$celltype_technique <- factor(integrated_object$celltype_technique,
                                               levels = c("Donor 6 scRNA-seq cDC1",
                                                          "Donor 6 snRNA-seq cDC1",
                                                          "Donor 6 scRNA-seq cDC2",
                                                          "Donor 6 snRNA-seq cDC2",
                                                          "Donor 6 scRNA-seq DC2pre-hm",
                                                          "Donor 6 snRNA-seq DC2pre-hm",
                                                          "Donor 6 scRNA-seq DC2hm",
                                                          "Donor 6 snRNA-seq DC2hm"))
DefaultAssay(integrated_object) <- 'RNA'
integrated_object <- NormalizeData(integrated_object) %>% ScaleData(features = rownames(integrated_object))
Idents(integrated_object) <- 'celltype_technique'
integrated_object$celltype_int <- Idents(integrated_object)
pseudobulk_expr <- AverageExpression(
  integrated_object,
  assays = "RNA",
  slot = 'data',
  group.by = c("technique", "celltype_int")
)

pseudobulk_df <- as.data.frame(pseudobulk_expr$RNA)
  #DC2hm comparison
cluster_to_compare <- "DC2hm"

# Create a data frame for plotting by pulling the relevant columns
comparison_df <- data.frame(
  scRNA = pseudobulk_df[['Donor 6 scRNA-seq_DC2hm']],
  snRNA = pseudobulk_df[['Donor 6 snRNA-seq_DC2hm']],
  gene = rownames(pseudobulk_df)
)

# Calculate the Pearson correlation coefficient
correlation <- cor(comparison_df$scRNA, comparison_df$snRNA, method = "pearson")
correlation_text <- paste0("Pearson's R = ", round(correlation, 3))
comparison_df$scRNA_log <- log1p(comparison_df$scRNA)
comparison_df$snRNA_log <- log1p(comparison_df$snRNA)


genes.to.annotated <- c(read.csv('Highlighted_gene_selected.csv')$Gene[7:31],
                        'MALAT1','MED13L','MT-CO1','RPL10',"RPL13" ,"RPS18")
# Filter the main dataframe to get the coordinates for  genes
annotation_df <- subset(comparison_df, gene %in% genes.to.annotated)
correlation_log <- cor(comparison_df$scRNA_log, comparison_df$snRNA_log, method = "pearson")
correlation_text_log <- paste0("Pearson's R = ", round(correlation_log, 3))
# Define an expression threshold to remove noisy and lowly expressed genes.
expression_threshold <- 0.2

comparison_df_filtered <- subset(comparison_df,
                                 scRNA_log > expression_threshold | snRNA_log > expression_threshold)

annotation_df_filtered <- subset(comparison_df_filtered, gene %in% genes.to.annotated)

correlation_filtered <- cor(comparison_df_filtered$scRNA_log,
                            comparison_df_filtered$snRNA_log,
                            method = "pearson")
correlation_text_filtered <- paste0("Pearson's R = ", round(correlation_filtered, 3))
rasterize(ggplot(comparison_df_filtered, aes(x = scRNA_log, y = snRNA_log)) +
  geom_point(alpha = 0.5, color = "#E679C5") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  geom_text_repel(data = annotation_df_filtered,
                  aes(label = gene),
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  labs(
    title = " Gene Expression Correlation: Cluster DC2hm",
    subtitle = correlation_text_filtered,
    x = "Log(Average Expression + 1) (scRNA-seq)",
    y = "Log(Average Expression + 1) (snRNA-seq)"
  ) +
  theme_bw(),dpi = 300)
ggsave('Outputs/figures/Correlation_analysis_scatter_DC2hm.pdf',
       width = 8, height = 8)


#cDC1 comparison
cluster_to_compare <- "cDC1"

comparison_df <- data.frame(
  scRNA = pseudobulk_df[['Donor 6 scRNA-seq_cDC1']],
  snRNA = pseudobulk_df[['Donor 6 snRNA-seq_cDC1']],
  gene = rownames(pseudobulk_df)
)

correlation <- cor(comparison_df$scRNA, comparison_df$snRNA, method = "pearson")
correlation_text <- paste0("Pearson's R = ", round(correlation, 3))
comparison_df$scRNA_log <- log1p(comparison_df$scRNA)
comparison_df$snRNA_log <- log1p(comparison_df$snRNA)


genes.to.annotated <- c(read.csv('Highlighted_gene_selected.csv')$Gene[1:3],
                        'MALAT1','MED13L','MT-CO1','RPL10',"RPL13" ,"RPS18")
annotation_df <- subset(comparison_df, gene %in% genes.to.annotated)

correlation_log <- cor(comparison_df$scRNA_log, comparison_df$snRNA_log, method = "pearson")
correlation_text_log <- paste0("Pearson's R = ", round(correlation_log, 3))

expression_threshold <- 0.2

comparison_df_filtered <- subset(comparison_df,
                                 scRNA_log > expression_threshold | snRNA_log > expression_threshold)

annotation_df_filtered <- subset(comparison_df_filtered, gene %in% genes.to.annotated)

correlation_filtered <- cor(comparison_df_filtered$scRNA_log,
                            comparison_df_filtered$snRNA_log,
                            method = "pearson")
correlation_text_filtered <- paste0("Pearson's R = ", round(correlation_filtered, 3))

rasterize(ggplot(comparison_df_filtered, aes(x = scRNA_log, y = snRNA_log)) +
  geom_point(alpha = 0.5, color = "#C5E524") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  geom_text_repel(data = annotation_df_filtered,
                  aes(label = gene),
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  labs(
    title = " Gene Expression Correlation: Cluster cDC1",
    subtitle = correlation_text_filtered,
    x = "Log(Average Expression + 1) (scRNA-seq)",
    y = "Log(Average Expression + 1) (snRNA-seq)"
  ) +
  theme_bw(),dpi = 300)
ggsave('Outputs/figures/Correlation_analysis_scatter_cDC1.pdf',
       width = 8, height = 8)

#cDC2 comparison
cluster_to_compare <- "cDC2"

comparison_df <- data.frame(
  scRNA = pseudobulk_df[['Donor 6 scRNA-seq_cDC2']],
  snRNA = pseudobulk_df[['Donor 6 snRNA-seq_cDC2']],
  gene = rownames(pseudobulk_df)
)

correlation <- cor(comparison_df$scRNA, comparison_df$snRNA, method = "pearson")
correlation_text <- paste0("Pearson's R = ", round(correlation, 3))
comparison_df$scRNA_log <- log1p(comparison_df$scRNA)
comparison_df$snRNA_log <- log1p(comparison_df$snRNA)


genes.to.annotated <- c(read.csv('Highlighted_gene_selected.csv')$Gene[4:31],
                        'MALAT1','MED13L','MT-CO1','RPL10',"RPL13" ,"RPS18")
annotation_df <- subset(comparison_df, gene %in% genes.to.annotated)

correlation_log <- cor(comparison_df$scRNA_log, comparison_df$snRNA_log, method = "pearson")
correlation_text_log <- paste0("Pearson's R = ", round(correlation_log, 3))

expression_threshold <- 0.2


comparison_df_filtered <- subset(comparison_df,
                                 scRNA_log > expression_threshold | snRNA_log > expression_threshold)


annotation_df_filtered <- subset(comparison_df_filtered, gene %in% genes.to.annotated)

correlation_filtered <- cor(comparison_df_filtered$scRNA_log,
                            comparison_df_filtered$snRNA_log,
                            method = "pearson")
correlation_text_filtered <- paste0("Pearson's R = ", round(correlation_filtered, 3))

rasterize(ggplot(comparison_df_filtered, aes(x = scRNA_log, y = snRNA_log)) +
  geom_point(alpha = 0.5, color = "#38C2E5") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  geom_text_repel(data = annotation_df_filtered,
                  aes(label = gene),
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  labs(
    title = " Gene Expression Correlation: Cluster cDC2",
    subtitle = correlation_text_filtered,
    x = "Log(Average Expression + 1) (scRNA-seq)",
    y = "Log(Average Expression + 1) (snRNA-seq)"
  ) +
  theme_bw(),dpi = 300)
# End of the line
ggsave('Outputs/figures/Correlation_analysis_scatter_cDC2.pdf',
       width = 8, height = 8)

#DC2pre-hm comparison
cluster_to_compare <- "DC2pre-hm"

comparison_df <- data.frame(
  scRNA = pseudobulk_df[['Donor 6 scRNA-seq_DC2pre-hm']],
  snRNA = pseudobulk_df[['Donor 6 snRNA-seq_DC2pre-hm']],
  gene = rownames(pseudobulk_df)
)

correlation <- cor(comparison_df$scRNA, comparison_df$snRNA, method = "pearson")
correlation_text <- paste0("Pearson's R = ", round(correlation, 3))
comparison_df$scRNA_log <- log1p(comparison_df$scRNA)
comparison_df$snRNA_log <- log1p(comparison_df$snRNA)


genes.to.annotated <- c(read.csv('../scRNA-seq analysis/Highlighted_gene_selected.csv')$Gene[4:31],
                        'MALAT1','MED13L','MT-CO1','RPL10',"RPL13" ,"RPS18")
annotation_df <- subset(comparison_df, gene %in% genes.to.annotated)

correlation_log <- cor(comparison_df$scRNA_log, comparison_df$snRNA_log, method = "pearson")
correlation_text_log <- paste0("Pearson's R = ", round(correlation_log, 3))
expression_threshold <- 0.2
comparison_df_filtered <- subset(comparison_df,
                                 scRNA_log > expression_threshold | snRNA_log > expression_threshold)
annotation_df_filtered <- subset(comparison_df_filtered, gene %in% genes.to.annotated)
correlation_filtered <- cor(comparison_df_filtered$scRNA_log,
                            comparison_df_filtered$snRNA_log,
                            method = "pearson")
correlation_text_filtered <- paste0("Pearson's R = ", round(correlation_filtered, 3))

rasterize(ggplot(comparison_df_filtered, aes(x = scRNA_log, y = snRNA_log)) +
  geom_point(alpha = 0.5, color = "#384C94") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick") +
  geom_text_repel(data = annotation_df_filtered,
                  aes(label = gene),
                  box.padding = 0.5,
                  max.overlaps = Inf) +
  labs(
    title = " Gene Expression Correlation: Cluster DC2pre-hm",
    subtitle = correlation_text_filtered,
    x = "Log(Average Expression + 1) (scRNA-seq)",
    y = "Log(Average Expression + 1) (snRNA-seq)"
  ) +
  theme_bw(),dpi = 300)
# End of the line
ggsave('Outputs/figures/Correlation_analysis_scatter_DC2pre-hm.pdf',
       width = 8, height = 8)


# Direct comparison for platform-sked getting platform-sked cell type markers
DimPlot(d6.scrna)
dc2hm.program.scrna <- FindAllMarkers(d6.scrna,only.pos = T)
DimPlot(multiomic.object.cdcs)
dc2hm.program.snrna  <- FindAllMarkers(multiomic.object.cdcs,only.pos = T)
# Filter for significant positive markers in the scRNA-seq dataset
significant.markers.scrna <- dc2hm.program.scrna %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) # Adjust thresholds if needed

# Create a named list of marker genes for each cell type (scRNA-seq)
marker.list.scrna <- significant.markers.scrna %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = 'drop') %>%
  tibble::deframe()

# Filter for significant positive markers in the multi-omic dataset
significant.markers.multi <- dc2hm.program.snrna  %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) # Adjust thresholds if needed

# Create a named list of marker genes for each cell type (multi-omic)
marker.list.multi <- significant.markers.multi %>%
  group_by(cluster) %>%
  summarise(genes = list(gene), .groups = 'drop') %>%
  tibble::deframe()

# Identify the common cell types (clusters) beten both datasets
common_clusters <- intersect(names(marker.list.scrna), names(marker.list.multi))
# Initialize an empty list to store the data frames for Excel
excel_output_list <- list()

for (cell_type in common_clusters) {
  
  genes_scrna <- marker.list.scrna[[cell_type]]
  genes_multi <- marker.list.multi[[cell_type]]
  
  comparison_list_venn <- list(
    scRNA = genes_scrna,
    snRNA = genes_multi
  )
  
  venn_plot <- ggvenn(
    comparison_list_venn, 
    fill_color = c("#0073C2FF", "#EFC000FF"),
    stroke_size = 0.5, 
    set_name_size = 5,
    text_size = 5
  ) +
    ggtitle(paste("Marker Gene Comparison for", cell_type)) +
    theme(plot.title = element_text(size = 16, face = "bold"))
  safe_filename <- gsub("[^A-Za-z0-9_.-]", "_", cell_type)
  
  ggsave(
    filename = paste0("Outputs/figures/Venn_Diagram_", safe_filename, ".pdf"), 
    plot = venn_plot, 
    width = 7, height = 6, # Standard PDF size in inches
    device = "pdf"
  )

  scrna_only <- setdiff(genes_scrna, genes_multi)
  multi_only <- setdiff(genes_multi, genes_scrna)
  overlap_genes <- intersect(genes_scrna, genes_multi)
  max_length <- max(length(scrna_only), length(multi_only), length(overlap_genes))
  
  comparison_df <- data.frame(
    scRNA_Only = c(scrna_only, rep(NA, max_length - length(scrna_only))),
    snRNA_Only = c(multi_only, rep(NA, max_length - length(multi_only))),
    Shared_Markers = c(overlap_genes, rep(NA, max_length - length(overlap_genes)))
  )
  
  safe_sheet_name <- substr(safe_filename, 1, 31)
  excel_output_list[[safe_sheet_name]] <- comparison_df
}

write.xlsx(excel_output_list,
           file = "Outputs/tables/cell_type_marker_comparison.xlsx")

# More focused DC2hm-related genes
dc2hm.program.scrna <- FindMarkers(d6.scrna, ident.1 = 'DC2hm', 'cDC2')
dc2hm.program.snrna <- FindMarkers(multiomic.object.cdcs,
                                   ident.1 = 'DC2hm', 'cDC2')
# --- 1. Process UP-regulated Genes (DC2hm > cDC2) ---

# Filter for significant UP-regulated genes
deg_scrna_up <- dc2hm.program.scrna %>%
  rownames_to_column(var = "gene") %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5)

deg_snrna_up <- dc2hm.program.snrna %>%
  rownames_to_column(var = "gene") %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5)

# Identify shared and unique UP-regulated genes and create data frames
shared_up_genes <- intersect(deg_scrna_up$gene, deg_snrna_up$gene)
scrna_only_up_df <- deg_scrna_up %>% filter(!gene %in% shared_up_genes)
snrna_only_up_df <- deg_snrna_up %>% filter(!gene %in% shared_up_genes)

shared_up_df <- deg_scrna_up %>%
  filter(gene %in% shared_up_genes) %>%
  select(gene, avg_log2FC.scRNA = avg_log2FC, p_val_adj.scRNA = p_val_adj) %>%
  inner_join(
    deg_snrna_up %>%
      filter(gene %in% shared_up_genes) %>%
      select(gene, avg_log2FC.snRNA = avg_log2FC, p_val_adj.snRNA = p_val_adj),
    by = "gene"
  )

# --- 2. Process DOWN-regulated Genes (DC2hm < cDC2) ---

# Filter for significant DOWN-regulated genes
deg_scrna_down <- dc2hm.program.scrna %>%
  rownames_to_column(var = "gene") %>%
  filter(p_val_adj < 0.05, avg_log2FC < -0.5)

deg_snrna_down <- dc2hm.program.snrna %>%
  rownames_to_column(var = "gene") %>%
  filter(p_val_adj < 0.05, avg_log2FC < -0.5)

# Identify shared and unique DOWN-regulated genes and create data frames
shared_down_genes <- intersect(deg_scrna_down$gene, deg_snrna_down$gene)
scrna_only_down_df <- deg_scrna_down %>% filter(!gene %in% shared_down_genes)
snrna_only_down_df <- deg_snrna_down %>% filter(!gene %in% shared_down_genes)

shared_down_df <- deg_scrna_down %>%
  filter(gene %in% shared_down_genes) %>%
  select(gene, avg_log2FC.scRNA = avg_log2FC, p_val_adj.scRNA = p_val_adj) %>%
  inner_join(
    deg_snrna_down %>%
      filter(gene %in% shared_down_genes) %>%
      select(gene, avg_log2FC.snRNA = avg_log2FC, p_val_adj.snRNA = p_val_adj),
    by = "gene"
  )

# --- Combine all data frames into a single list for export ---
combined_output <- list(
  "Shared_UP_DEGs" = shared_up_df,
  "scRNA_Only_UP" = scrna_only_up_df,
  "snRNA_Only_UP" = snrna_only_up_df,
  "Shared_DOWN_DEGs" = shared_down_df,
  "scRNA_Only_DOWN" = scrna_only_down_df,
  "snRNA_Only_DOWN" = snrna_only_down_df
)

write.xlsx(combined_output, 
           file = "Outputs/tables/DC2hm_vs_cDC2_Combined_Stats.xlsx")

# --- Create List for Upregulated Genes Venn Diagram ---
comparison_list_up <- list(
  scRNA = deg_scrna_up$gene,
  snRNA = deg_snrna_up$gene
)

# 2. Extract gene sets
set1 <- comparison_list_up[[1]]
set2 <- comparison_list_up[[2]]
set_names <- names(comparison_list_up)

# 3. Calculate the size of each segment
# n1 is the total in the first circle, n2 is the total in the second
# n12 is the intersection
area1_val <- length(set1)
area2_val <- length(set2)
cross_area_val <- length(intersect(set1, set2))

# 4. Draw the Venn diagram
# This function creates and saves the plot directly to a file (e.g., PNG)
# To prevent the function from creating a log file,  add this line:
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.plot <- draw.pairwise.venn(
  area1 = area1_val,
  area2 = area2_val,
  cross.area = cross_area_val,
  category = set_names,
  
  # --- Customizations ---
  scaled = TRUE,              # This is the crucial argument for proportionality
  fill = c("#CD534CFF", "#0073C2FF"), # Colors for the circles
  alpha = 0.7,                  # Color transparency
  lty = "blank",                # Removes the circle border lines
  
  # Text styling
  cex = 1.5,                    # Size for the numbers inside circles
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.8,                # Size for the category names
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.dist = c(0.05, 0.05),     # Distance of category names from circles
  
  # Title
  main = "Upregulated Genes: DC2hm vs cDC2",
  main.cex = 1.8,
  main.fontface = "bold"
)

grid.draw(venn.plot)


# 2. Extract DOWNREGULATED gene sets
#    Make sure to replace 'comparison_list_down' with  actual variable name
set1_down <- comparison_list_down[[1]]
set2_down <- comparison_list_down[[2]]
set_names_down <- names(comparison_list_down)

# 3. Calculate the size of each segment
area1_down <- length(set1_down)
area2_down <- length(set2_down)
cross_area_down <- length(intersect(set1_down, set2_down))

# 4. Suppress the log file creation
flog.threshold(ERROR, name = "VennDiagramLogger")

# 5. Draw the Venn diagram
venn.plot.down <- draw.pairwise.venn(
  area1 = area1_down,
  area2 = area2_down,
  cross.area = cross_area_down,
  category = set_names_down,
  
  # --- Customizations ---
  scaled = TRUE,              # Ensure circles are proportional to gene counts
  fill = c("#0073C2FF", "#56B4E9"), # Using shades of blue for downregulation
  alpha = 0.7,
  lty = "blank",
  
  # Text styling (kept consistent with the previous plot)
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.8,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.dist = c(0.05, 0.05),
  
  # Title for this plot
  main = "Downregulated Genes: DC2hm vs cDC2",
  main.cex = 1.8,
  main.fontface = "bold"
)

# --- Generate and Save Venn Diagram for Downregulated Genes ---
venn_plot_down <- ggvenn(
  comparison_list_down,
  fill_color = c("#00A087FF", "#3C5488FF"), # Distinct colors for DOWN
  stroke_size = 0.5,
  set_name_size = 5,
  text_size = 5
) +
  ggtitle("Downregulated Genes: DC2hm vs cDC2") +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave(
  filename = "Outputs/figures/Venn_DC2hm_vs_cDC2_Downregulated_Genes.pdf",
  plot = venn_plot_down,
  width = 7, height = 6
)

# Initialize an empty list to store the GO analysis results
go_results_list <- list()

# Loop through each data frame in  combined_output list
for (sheet_name in names(combined_output)) {
  
  # Get the current data frame
  current_df <- combined_output[[sheet_name]]
  
  # Extract the gene symbols and remove any NAs
  genes_to_test <- na.omit(current_df$gene)
  
  cat(paste("Processing:", sheet_name, "with", length(genes_to_test), "genes...\n"))
  
  # Only run if there are genes in the list
  if (length(genes_to_test) > 0) {
    
    # Run GO enrichment without the 'universe' argument
    ego <- enrichGO(gene          = genes_to_test,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2)
    
    # Check if any significant results re returned
    if (!is.null(ego) && nrow(ego) > 0) {
      # Add the simplified results to our list
      go_results_list[[sheet_name]] <- as.data.frame(ego) %>%
        dplyr::select(ONTOLOGY, ID, Description, p.adjust, GeneRatio, Count, geneID)
    } else {
      cat(paste("  -> No significant GO terms found for", sheet_name, "\n"))
    }
  }
}

# --- Save all results to a new Excel file ---
if (length(go_results_list) > 0) {
  write.xlsx(go_results_list, file = "Outputs/tables/GO_Enrichment_Results.xlsx")
  print("Analysis complete! Results saved to 'GO_Enrichment_Results.xlsx'")
} else {
  print("Analysis complete, but no significant enrichment results re found.")
}

# 3. Define the list of specific GO IDs for each category
designated_go_ids <- list(
  Shared_UP_DEGs = c("GO:0022407", "GO:0002696", "GO:0031295", "GO:0043122", "GO:0006568"),
  scRNA_Only_UP = c("GO:0019882", "GO:0002274", "GO:0030667", "GO:0003925", "GO:0061134"),
  snRNA_Only_UP = c("GO:0002181", "GO:0009612", "GO:0030055", "GO:0003714", "GO:0015629"),
  Shared_DOWN_DEGs = c("GO:0006909", "GO:0002431", "GO:0071356", "GO:0032612", "GO:0009615"),
  scRNA_Only_DOWN = c("GO:0050729", "GO:0002237", "GO:0009620", "GO:0045333", "GO:0035325"),
  snRNA_Only_DOWN = c("GO:0010506", "GO:0048538", "GO:0006638", "GO:0062197", "GO:0061726")
)

# 4. Loop through each category to create and save plots
for (result_name in names(designated_go_ids)) {
  
  full_go_df <- go_results_list[[result_name]]
  target_ids <- designated_go_ids[[result_name]]
  
  terms_to_plot <- full_go_df %>%
    filter(ID %in% target_ids)
  
  if (nrow(terms_to_plot) == 0) {
    warning(paste("For category '", result_name, "', none of the specified GO IDs re found. Skipping plot."))
    next
  }
  
  terms_to_plot <- terms_to_plot %>%
    arrange(Count) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  go_plot <- ggplot(terms_to_plot, aes(x = Description, y = Count, fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    coord_flip() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) +
    labs(
      title = paste("Designated GO Terms:", result_name),
      x = "",
      y = "Enriched Gene Count",
      fill = expression(-log[10]("P-adjust"))
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  if (grepl("UP", result_name)) {
    go_plot <- go_plot + scale_fill_gradient(low = "#FDE725FF", high = "#D95F02")
  } else {
    go_plot <- go_plot + scale_fill_gradient(low = "#85C1E9", high = "#1A5276")
  }
  
  print(go_plot)
  
  # 5. *** Save the plot to the specified folder ***
  ggsave(
    filename = file.path("Outputs/figures", paste0("GO_barplot_", result_name, ".pdf")),
    plot = go_plot,
    width = 8,
    height = 8,
    dpi = 300
  )
}

# Vlnplot comparison ------------------------------------------------------


VlnPlot(d6.scrna, features = c('CCR7','IL7R','CRLF2','CD80','BCL2','CCL22'),
        pt.size = 0, cols =  c('#C5E524', '#38C2E5','#384C94','#E679C5'))
ggsave('Outputs/figures/platform_common_up_scrna.pdf',
       width = 8, height = 6)
VlnPlot(d6.scrna, features = c('RPS23','RPL31','RPS12'),
        pt.size = 0, cols =  c('#C5E524', '#38C2E5','#384C94','#E679C5'))
ggsave('Outputs/figures/platform_snrna_up_scrna.pdf',
       width = 8, height = 3)
VlnPlot(d6.scrna, features = c('CCL17','CXCL9','IL32'),
        pt.size = 0, cols =  c('#C5E524', '#38C2E5','#384C94','#E679C5'))
ggsave('Outputs/figures/platform_scrna_up_scrna.pdf',
       width = 8, height = 3)
VlnPlot(multiomic.object.cdcs, 
        features = c('CCR7','IL7R','CRLF2','CD80','BCL2','CCL22'),
        pt.size = 0, cols =  c('#C5E524', '#38C2E5','#384C94','#E679C5'))
ggsave('Outputs/figures/platform_common_up_snrna.pdf',
       width = 8, height = 6)
VlnPlot(multiomic.object.cdcs, features = c('RPS23','RPL31','RPS12'),
        pt.size = 0, cols =  c('#C5E524', '#38C2E5','#384C94','#E679C5'))
ggsave('Outputs/figures/platform_snrna_up_snrna.pdf',
       width = 8, height = 3)
VlnPlot(multiomic.object.cdcs, features = c('CCL17','CXCL9','IL32'),
        pt.size = 0, cols =  c('#C5E524', '#38C2E5','#384C94','#E679C5'))
ggsave('Outputs/figures/platform_scrna_up_snrna.pdf',
       width = 8, height = 3)

                     
# Platform gene calling by direct comparison ---------------------------------------------------
Idents(integrated_object) <- 'celltype_int'
all_cell_types <- levels(Idents(integrated_object))

# Create an empty list to store the DE results for each cell type
de_results_list <- list()

# Loop through each cell type
for (cell_type in all_cell_types) {
  
  # Print progress
  message(paste("Processing:", cell_type))
  
  # Subset the Seurat object to only include the current cell type
  subset_object <- subset(integrated_object, idents = cell_type)
  
  # Run FindMarkers to compare the two technologies within this cell type
  #  use the 'RNA' assay for the original, non-integrated counts
  #  group by the 'tech' metadata column
  de_genes <- FindMarkers(
    subset_object,
    ident.1 = "Donor 6 scRNA-seq",
    ident.2 = "Donor 6 snRNA-seq",
    group.by = "technique",
    assay = "RNA",
    logfc.threshold = 0.25 # Set a minimum logFC threshold
  )
  
  # Add gene symbols as a column for plotting
  de_genes$gene <- rownames(de_genes)
  
  # Store the results in our list
  de_results_list[[cell_type]] <- de_genes
}

# Retrieve the results from the list
results_df <- de_results_list[[cell_type_to_plot]]

library(ggplot2)
library(ggrepel)
library(dplyr)

# (This code assumes 'results_df' and the threshold variables are already defined)

# --- Step 1: Define and count significant genes ---
log2fc_threshold <- 1
padj_threshold <- 0.05

up_regulated_genes <- results_df %>%
  filter(avg_log2FC > log2fc_threshold & p_val_adj < padj_threshold)

down_regulated_genes <- results_df %>%
  filter(avg_log2FC < -log2fc_threshold & p_val_adj < padj_threshold)

# --- Step 2: Prepare data for plotting ---
# (This section is the same as before, preparing 'plot_data')
plot_data <- results_df %>%
  arrange(avg_log2FC) %>%
  mutate(rank = 1:n()) %>%
  mutate(significance = case_when(
    avg_log2FC > log2fc_threshold & p_val_adj < padj_threshold ~ "Up-regulated",
    avg_log2FC < -log2fc_threshold & p_val_adj < padj_threshold ~ "Down-regulated",
    TRUE ~ "Not Significant"
  )) %>%
  mutate(gene_label = "")

top_up_genes <- tail(plot_data$gene, 10)
top_down_genes <- head(plot_data$gene, 10)
plot_data$gene_label[plot_data$gene %in% c(top_up_genes, top_down_genes)] <-
  plot_data$gene[plot_data$gene %in% c(top_up_genes, top_down_genes)]

# --- Step 3: Generate the updated plot with detailed labels ---
ggplot(plot_data, aes(x = rank, y = avg_log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(color = significance), alpha = 0.5, size = 1) +
  scale_color_manual(
    name = "Significance",
    values = c(
      "Up-regulated" = "#E41A1C",
      "Down-regulated" = "#377EB8",
      "Not Significant" = "grey"
    )
  ) +
  geom_text_repel(
    aes(label = gene_label),
    max.overlaps = Inf,
    box.padding = 0.5,
    color = "black"
  ) +
  # Updated labs() section with DEG counts and comparison info
  labs(
    title = paste("Ranked DE Genes in", cell_type_to_plot),
    subtitle = paste0(
      "Comparison: scRNA-seq vs. snRNA-seq\n",
      "Up: ", nrow(up_regulated_genes), " genes | ",
      "Down: ", nrow(down_regulated_genes), " genes (Log2FC > 1 & p-adj < 0.05)"
    ),
    x = "Gene Rank",
    y = "Log2 Fold Change"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "top"
  )
ggsave('Outputs/figures/Platform-gene calling DC2hm.pdf', width = 8,
       height = 8)

# --- Step 1: Initialize Lists to Store Results ---
# These lists will hold the GO enrichment data frames and the plots
go_up_results <- list()
go_down_results <- list()
go_plots <- list()

# Define DE thresholds
log2fc_threshold <- 1
padj_threshold <- 0.05

# Get the list of all cell types from  previous DE analysis
all_cell_types <- names(de_results_list)

# --- Step 2: Loop Through Each Cell Type ---
for (cell_type in all_cell_types) {
  
  message(paste("--- Processing GO Enrichment for:", cell_type, "---"))
  
  # Retrieve the DE results for the current cell type
  results_df <- de_results_list[[cell_type]]
  
  # --- Define gene sets ---
  # Up-regulated (enriched in scRNA-seq)
  up_genes <- results_df %>%
    filter(avg_log2FC > log2fc_threshold & p_val_adj < padj_threshold) %>%
    rownames()
  
  # Down-regulated (enriched in snRNA-seq)
  down_genes <- results_df %>%
    filter(avg_log2FC < -log2fc_threshold & p_val_adj < padj_threshold) %>%
    rownames()
  
  # Define the "universe" of all genes that re tested
  universe_genes <- rownames(results_df)
  
  # --- Run GO Enrichment for Up-regulated Genes ---
  #  check if there are enough genes to run the analysis
  if (length(up_genes) > 10) {
    go_up <- enrichGO(
      gene          = up_genes,
      universe      = universe_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = 'SYMBOL',
      ont           = "CC", # BP = Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    
    # Store the results
    go_up_results[[cell_type]] <- as.data.frame(go_up)
    
    # Create and store the plot
    p_up <- dotplot(go_up, showCategory = 10) + 
      ggtitle(paste(cell_type, "- Enriched in scRNA-seq")) +
      # Manually increase plot margins for better visualization
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    
    go_plots[[paste0(cell_type, "_up")]] <- p_up
  }
  
  # --- Run GO Enrichment for Down-regulated Genes ---
  if (length(down_genes) > 10) {
    go_down <- enrichGO(
      gene          = down_genes,
      universe      = universe_genes,
      OrgDb         = org.Hs.eg.db,
      keyType       = 'SYMBOL',
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    
    go_down_results[[cell_type]] <- as.data.frame(go_down)
    
    p_down <- dotplot(go_down, showCategory = 10) + 
      ggtitle(paste(cell_type, "- Enriched in snRNA-seq")) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    
    go_plots[[paste0(cell_type, "_down")]] <- p_down
  }
}

message("--- GO Enrichment Analysis Complete ---")

# --- Step 1: Create a Directory to Store the PDFs ---
# This keeps  main folder clean.
output_dir <- "Outputs/figures/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define DE thresholds
log2fc_threshold <- 1
padj_threshold <- 0.05

# Define how many top/bottom genes to annotate
num_genes_to_annotate <- 5 # Changed from 10 to 7 to reduce clutter

# Get the list of all cell types
all_cell_types <- names(de_results_list)

# --- Step 2: Loop Through Each Cell Type to Generate and Save Reports ---
for (cell_type in all_cell_types) {
  
  message(paste("--- Generating Report for:", cell_type, "---"))
  
  # === PART A: Generate the Ranked Log2FC Dot Plot ===
  
  results_df <- de_results_list[[cell_type]]
  
  # Define and count significant genes
  up_genes_df <- results_df %>% filter(avg_log2FC > log2fc_threshold & p_val_adj < padj_threshold)
  down_genes_df <- results_df %>% filter(avg_log2FC < -log2fc_threshold & p_val_adj < padj_threshold)
  unchanged_genes_count <- nrow(results_df) - nrow(up_genes_df) - nrow(down_genes_df)
  
  # Prepare data for ranked plot
  plot_data <- results_df %>%
    arrange(avg_log2FC) %>%
    mutate(rank = 1:n(),
           significance = case_when(
             avg_log2FC > log2fc_threshold & p_val_adj < padj_threshold ~ "Up-regulated",
             avg_log2FC < -log2fc_threshold & p_val_adj < padj_threshold ~ "Down-regulated",
             TRUE ~ "Not Significant"
           ),
           gene_label = "")
  
  # Identify the top N up-regulated (highest log2FC)
  top_up_genes_to_label <- plot_data %>% 
    filter(significance == "Up-regulated") %>%
    arrange(desc(avg_log2FC)) %>% 
    head(num_genes_to_annotate) %>% pull(gene)
  
  # Identify the top N down-regulated (lost log2FC)
  top_down_genes_to_label <- plot_data %>%
    filter(significance == "Down-regulated") %>%
    arrange(avg_log2FC) %>%
    head(num_genes_to_annotate) %>% pull(gene)
  
  # Assign labels only to the selected genes
  plot_data$gene_label[plot_data$gene %in% top_up_genes_to_label] <- 
    plot_data$gene[plot_data$gene %in% top_up_genes_to_label]
  plot_data$gene_label[plot_data$gene %in% top_down_genes_to_label] <- 
    plot_data$gene[plot_data$gene %in% top_down_genes_to_label]
  
  # MODIFIED PLOT CODE
  p_ranked <- ggplot(plot_data, aes(x = rank, y = avg_log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(aes(color = significance), alpha = 0.8, size = 2) +
    scale_color_manual(name = "Significance", 
                       values = c("Up-regulated" = "#E41A1C", 
                                  "Down-regulated" = "#377EB8", 
                                  "Not Significant" = "grey")) +
    geom_text_repel(
      aes(label = gene_label, color = significance), # Color labels to match points
      max.overlaps = Inf, # Allow some overlap to avoid lines going too far, increase if still too many lines
      box.padding = 0.6, # Slightly increase space around text
      point.padding = 0.3, # Slightly increase space from point
      segment.size = 0.3, # Thinner lines
      segment.alpha = 0.6, # More transparent lines
      force = 1, # Increase force to push labels further apart
      size = 3.5, # Slightly smaller text
      min.segment.length = 0.5 # Minimum length of segments (lines to points)
    ) +
    labs(
      title = "Ranked Differential Expression",
      subtitle = paste0("Up: ", nrow(up_genes_df), " | Down: ", nrow(down_genes_df), " | Unchanged: ", unchanged_genes_count),
      x = "Gene Rank", y = "Log2 Fold Change"
    ) +
    theme_bw() + theme(legend.position = "bottom")
  
  # === PART B: Generate GO Enrichment Bar Plots ===
  
  # Placeholder plots in case there are not enough genes for analysis
  p_up <- plot_spacer() + ggtitle("No significant up-regulated GO Terms")
  p_down <- plot_spacer() + ggtitle("No significant down-regulated GO Terms")
  
  # Plot for up-regulated genes (enriched in scRNA-seq)
  if (nrow(up_genes_df) > 10) { # Adjusted to use up_genes_df
    go_up <- enrichGO(gene = rownames(up_genes_df), universe = rownames(results_df), 
                      OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "CC")
    if (!is.null(go_up) && nrow(go_up) > 0) {
      go_up_df <- as.data.frame(go_up) %>% head(10)
      p_up <- ggplot(go_up_df, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_fill_gradient(low = "#E41A1C", high = "#FBB4AE", name = "p.adjust") +
        labs(title = "Enriched in scRNA-seq", x = "Gene Count", y = NULL) +
        theme_minimal() + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    }
  }
  
  # Plot for down-regulated genes (enriched in snRNA-seq)
  if (nrow(down_genes_df) > 10) { # Adjusted to use down_genes_df
    go_down <- enrichGO(gene = rownames(down_genes_df), universe = rownames(results_df),
                        OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "CC")
    if (!is.null(go_down) && nrow(go_down) > 0) {
      go_down_df <- as.data.frame(go_down) %>% head(10)
      p_down <- ggplot(go_down_df, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_fill_gradient(low = "#377EB8", high = "#B3CDE3", name = "p.adjust") +
        labs(title = "Enriched in snRNA-seq", x = "Gene Count", y = NULL) +
        theme_minimal() + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
    }
  }
  
  # === PART C: Combine and Save the Report ===
  
  combined_plot <- (p_ranked | (p_up / p_down)) +
    plot_layout(widths = c(2, 1))+
    plot_annotation(
      title = paste("Cross-Platform Analysis for:", cell_type),
      subtitle = "Comparison: scRNA-seq vs. snRNA-seq",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 12))
    )
  
  pdf_path <- file.path(output_dir, paste0("Report_", cell_type, ".pdf"))
  ggsave(pdf_path, plot = combined_plot, width = 12, height = 6, units = "in")
}

message(paste("--- All reports have been saved to the '", output_dir, "' folder. ---"))
go_term_list <- list(
  "cDC1" = c("GO:0022626", "GO:0098800", "GO:0005925", "GO:0030055", "GO:0098803", 
             "GO:0042611", "GO:0098553", "GO:0031975", "GO:0045202", "GO:1990351"),
  "cDC2" = c("GO:0022626", "GO:0005925", "GO:0042611", "GO:0005753", "GO:0070161", 
             "GO:0098553", "GO:0030669", "GO:1904813", "GO:0030139", "GO:0005764"),
  "DC2pre-hm" = c("GO:0022626", "GO:0030055", "GO:0042611", "GO:0034774", "GO:0060205", 
                  "GO:0101002", "GO:0045202", "GO:1990204", "GO:0016282", "GO:0005765"),
  "DC2hm" = c("GO:0022626", "GO:0042611", "GO:0070161", "GO:0098803", "GO:0045202", 
              "GO:0030139", "GO:0045334", "GO:0032588", "GO:0030136", "GO:0070820")
)
# Define thresholds and get cell type names
log2fc_threshold <- 1
padj_threshold <- 0.05
all_cell_types <- names(de_results_list)
output_dir_up_go <- output_dir
for (cell_type in all_cell_types) {
  
  message(paste("--- Generating UP GO plot for:", cell_type, "---"))
  
  # Get DE results and filter for up-regulated genes
  results_df <- de_results_list[[cell_type]]
  up_genes_df <- results_df %>% 
    filter(avg_log2FC > log2fc_threshold & p_val_adj < padj_threshold)
  
  # Proceed only if there are enough genes for analysis
  if (nrow(up_genes_df) > 10) {
    go_up <- enrichGO(gene = rownames(up_genes_df), universe = rownames(results_df),
                      OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "CC")
    
    if (!is.null(go_up) && nrow(go_up) > 0) {
      go_up_df_full <- as.data.frame(go_up)
      
      # Check if the current cell type has a custom GO list
      if (cell_type %in% names(go_term_list)) {
        go_up_df <- go_up_df_full %>% 
          filter(ID %in% go_term_list[[cell_type]])
      } else {
        go_up_df <- go_up_df_full %>% head(10)
      }
      
      # Generate and save the plot if any terms remain after filtering
      if (nrow(go_up_df) > 0) {
        p_up <- ggplot(go_up_df, aes(x = Count, y = reorder(Description, Count), fill = p.adjust)) +
          geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.3) +
          scale_fill_gradient(low = "#E41A1C", high = "#FBB4AE", name = "Adjusted p-value") +
          labs(
            title = paste(cell_type, "GO Enrichment"),
            subtitle = "Biological Processes Enriched in scRNA-seq",
            x = "Gene Count",
            y = NULL
          ) +
          theme_bw(base_size = 12) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right"
          )
        
        # Define the file path and save the plot
        pdf_path <- file.path(output_dir_up_go, paste0("GO_Up_", cell_type, ".pdf"))
        ggsave(pdf_path, plot = p_up, width = 8, height = 6, units = "in")
      }
    }
  }
}

message(paste("--- All plots have been saved to the '", output_dir_up_go, "' folder. ---"))
write.xlsx(go_down_results,'Outputs/figures/GO_down.xlsx')

# Thank god the tedious part is over.
                     
# ATAC markers ------------------------------------------------------------
marker.cluster.list <- list()
# add gene annotation to region
for (i in names(table(markers$cluster))) {
  marker.cluster <- markers[markers$cluster==i & markers$p_val_adj<0.05,]
  marker.cluster.list[[i]] <- marker.cluster[order(marker.cluster$avg_log2FC,
                                                   decreasing = T),]
}
write.xlsx(marker.cluster.list, 
           'Outputs/tables/Markers_by_clusters_unannotated.xlsx',
           rowNames=T)
  # Annotate gene closed to the feature accessible chromatin region
marker.cluster.list.anno <- marker.cluster.list
for (i in names(table(markers$cluster))) {
  marker.cluster.list.anno[[i]] <- marker.cluster.list.anno[[i]][marker.cluster.list.anno[[i]]$avg_log2FC>0,]
  gene.anno <- ClosestFeature(multiomic.object,
                              regions = rownames(marker.cluster.list.anno[[i]]))$gene_name
  
  marker.cluster.list.anno[[i]]$GeneAnno <- gene.anno
}
write.xlsx(marker.cluster.list.anno, 
           'Outputs/tables/Markers_by_clusters_annotated.xlsx',
           rowNames=T)
de.peaks.df <- rbind(marker.cluster.list.anno[[1]],
                     marker.cluster.list.anno[[2]],
                     marker.cluster.list.anno[[3]],
                     marker.cluster.list.anno[[4]])
write.csv(de.peaks.df, 'Outputs/tables/DE_peaks_dataframe.csv')


de.peaks.df <- read.csv('Outputs/tables/DE_peaks_dataframe.csv', row.names = 1)

de.region <- rownames(de.peaks.df)
rownames(de.peaks.df) <- de.region
DefaultAssay(multiomic.object) <- 'peaks'
de.peaks.mat <- AverageExpression(multiomic.object.cdcs,
                                  assays = 'peaks',
                                  features = unique(c(marker.cluster.list.anno[[1]]$gene,
                                                      marker.cluster.list.anno[[2]]$gene,
                                                      marker.cluster.list.anno[[3]]$gene,
                                                      marker.cluster.list.anno[[4]]$gene)))[[1]]
de.peaks.mat.scaled <-  t(apply(de.peaks.mat, 1, function(x) x / sum(x)))
write.csv(de.peaks.mat.scaled, 'Outputs/tables/DE_peaks_mat_scaled.csv')
de.peaks.mat.scaled <- read.csv('Outputs/tables/DE_peaks_mat_scaled.csv',
                                row.names = 1)
# Annotated gene
chr.region <- c('chr11-115379443-115380088',
                'chr2-112836308-112837090','chr1-247408079-247408501',
                'chr3-38151664-38152273','chr16-31350184-31350581',
                'chr1-158283542-158284402','chr21-41361462-41362886',
                'chr3-112364026-112364862','chr17-40577480-40578060',
                'chr19-10278125-10278589','chr19-44648615-44649244',
                'chr5-35854206-35854639','chr16-57413301-57413907',
                'chrX-1184058-1184825','chr3-122103690-122104432',
                'chr12-93573274-93574172','chr19-45016035-45017288',
                'chr11-35126120-35126898','chr9-5465763-5466362',
                'chr7-5582356-5583034','chr20-10676453-10676825',
                'chr3-183151602-183152325','chr6-113854565-113855586',
                'chr16-57340874-57341441','chr20-46117603-46119591',
                'chr19-4228843-4229963','chr17-42286729-42290256',
                'chr3-119562279-119563111')
de.peaks.df[chr.region, 'GeneAnno']
index <- match(chr.region, rownames(de.peaks.mat.scaled))
annotation <- de.peaks.df[chr.region, 'GeneAnno']
anno.markers <- rowAnnotation(ano = anno_mark(at = index, side = 'left',
                                              labels = annotation,
                                              labels_gp = gpar(fontsize = 5),
                                              lines_gp = gpar(),
                                              extend = unit(1, 'npc')))
# Let's create a nice palette                                
# Determine the maximum absolute value for symmetric scaling
max_abs_val <- max(abs(de.peaks.mat.scaled))

# Define the color scale breaks and corresponding viridis colors
color_breaks <- c(-max_abs_val, 0, max_abs_val)

# Create a color function using the Magma viridis option, spanning the breaks
viridis_color_fn <- colorRamp2(
  breaks = color_breaks,
  colors = viridis(3) # Creates three colors from the magma scale: low, middle (white/neutral), high
)

rasterize(as.ggplot(Heatmap(de.peaks.mat.scaled,
                            cluster_columns = FALSE, 
                            cluster_rows = FALSE,
                            show_row_names = FALSE,
                            width = unit(0.2, 'npc'),
                            left_annotation = anno.markers,
                            height = unit(0.6, 'npc'), 
                            border = TRUE, 
                            row_title_rot = 0,
                            # --- REVISED COLOR ARGUMENT ---
                            col = viridis_color_fn)),
          dpi = 300)
ggsave('Outputs/figures/Feature_peaks_heatmap.pdf', width = 3.5, height = 5)

# Wikipathway and KEGG enrichment
IDTransandenrichWP <- function(gene.list){
  result <- enrichWP(bitr(gene.list,
                          fromType = 'SYMBOL', toType = 'ENTREZID',
                          OrgDb = 'org.Hs.eg.db')$ENTREZID,
                     organism = 'Homo sapiens')
  return(result)
}
  # Figure 1O, bar plots for enrichment of DC2hm DE peaks
dc2hm.peaks.enrichwp <- IDTransandenrichWP(unique(marker.cluster.list.anno[[4]]$GeneAnno))
write.csv(dc2hm.peaks.enrichwp@result, 'Outputs/tables/DC2hm_peaks_wiki.csv')
enrichment.item.to.show <- c('Chemokine signaling',
                             'Apoptosis modulation and signaling',
                             'Modulators of TCR signaling and T cell activation',
                             'RANKL/RANK signaling',
                             'Thymic stromal lymphopoietin (TSLP) signaling')
dc2hm.peaks.enrichwp.show <- dc2hm.peaks.enrichwp@result[dc2hm.peaks.enrichwp@result$Description %in% enrichment.item.to.show,]
dc2hm.peaks.enrichwp.show$Description <- factor(dc2hm.peaks.enrichwp.show$Description,
                                              levels = rev(enrichment.item.to.show))
ggplot(dc2hm.peaks.enrichwp.show, aes(x=Description,y=Count,fill=p.adjust))+
  geom_bar(stat = 'identity')+
  coord_flip()+
  scale_fill_gradient(high = '#fef6b7', low = '#e15383')+
  theme_classic()
ggsave('Outputs/figures/DE_peaks_dc2hm_enrichment.pdf',
       width = 6, height = 1.5)

  # The fragments file is available @ GEO: GSE267255
UpdatePath(multiomic.object.cdcs@assays$ATAC@fragments[[1]],
           'atac_fragments.tsv.gz')
for (i in c('CLEC9A','CD1C','FSCN1', 'PVR', 'CRLF2','STAT5A','CD200')) {
  pdf(paste0('Outputs/figures/', i, '_coverage_plot.pdf'),
      width=3.5, height=3)
  plot(CoveragePlot(object = multiomic.object.cdcs, region = i, features = i, 
                    expression.assay = "SCT", assay = 'ATAC',
                    extend.upstream = 2000, extend.downstream = 2000)&
         scale_fill_manual(values =c('#C5E524', '#38C2E5','#384C94','#E679C5')))
  dev.off()
}

# Motif enrichment analysis -----------------------------------------------
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE",
              tax_group = 'vertebrates', all_versions = FALSE)
)
multiomic.object.cdcs <- AddMotifs(multiomic.object.cdcs, 
                              genome = BSgenome.Hsapiens.UCSC.hg38,
                              pfm = pfm, assay = 'peaks')
DefaultAssay(multiomic.object) <- 'peaks'
enrich.motifs.dc2hm <- FindMotifs(object = multiomic.object.cdcs, 
                                features = marker.cluster.list.anno[[4]]$gene)
  # Plot motif enrichment results, Figure 1P and S2M
# cDC2
enrich.motifs.cdc2 <- FindMotifs(object = multiomic.object.cdcs, 
                                 features = marker.cluster.list.anno[[2]]$gene)
# cDC1
enrich.motifs.cdc1 <- FindMotifs(object = multiomic.object.cdcs, 
                                 features = marker.cluster.list.anno[[1]]$gene)
# Save motif enrichment results
write.xlsx(list('cDC1'=enrich.motifs.cdc1,
                'cDC2'=enrich.motifs.cdc2,
                'DC2hm'=enrich.motifs.dc2hm),
           'Outputs/tables/cDC motif enrichment.xlsx')

# Plot Motif enrichment results 
# cDC1
enrich.motifs.cdc1.show <- enrich.motifs.cdc1[c('MA0652.1',
                                                'MA0517.1',
                                                'MA0835.2'),]
enrich.motifs.cdc1.show$motif.name <- factor(enrich.motifs.cdc1.show$motif.name,
                                             levels = c(
                                               'BATF3', 'STAT1::STAT2','IRF8'))
ggplot(enrich.motifs.cdc1.show,aes(y=motif.name,x=fold.enrichment))+
  geom_bar(stat = 'identity', fill='#C5E524')+
  theme_classic()
ggsave('Outputs/figures/Motif_enrich_barplot_cDC1.pdf',width = 2,height = 1 )
# cDC2
enrich.motifs.cdc2.show <- enrich.motifs.cdc2[c('MA0102.4',
                                                'MA0466.2',
                                                'MA1484.1'),]
enrich.motifs.cdc2.show$motif.name <- factor(enrich.motifs.cdc2.show$motif.name,
                                             levels = c(
                                               'ETS2', 'CEBPA','CEBPB'))
ggplot(enrich.motifs.cdc2.show,aes(y=motif.name,x=fold.enrichment))+
  geom_bar(stat = 'identity', fill='#38C2E5')+
  theme_classic()
ggsave('Outputs/figures/Motif_enrich_barplot_cDC2.pdf',width = 2,height = 1 )
# DC2hm
enrich.motifs.dc2hm.show <- enrich.motifs.dc2hm[c('MA0101.1',
                                              'MA1117.1',
                                              'MA0519.1'),]
enrich.motifs.dc2hm.show$motif.name <- factor(enrich.motifs.dc2hm.show$motif.name,
                                            levels = c(
                                              'Stat5a::Stat5b', 'RELB','REL'))
ggplot(enrich.motifs.dc2hm.show,aes(y=motif.name,x=fold.enrichment))+
  geom_bar(stat = 'identity', fill='#E679C5')+
  theme_classic()
ggsave('Outputs/figures/Motif_enrich_barplot_DC2hm.pdf',width = 2,height = 1 )

