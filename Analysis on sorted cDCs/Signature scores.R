library(DESeq2)
library(Seurat)
library(dplyr)
library(openxlsx)
library(ggsci)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
setwd("/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/Analysis on sorted cDCs")
# Read the counts and the Seurat d6.object.cdcs -----------------------------------
d6.object.cdcs <- readRDS('../scRNA-seq analysis/Outputs/rds/Dnr6_cDCs.rds')
d2.counts <- read.csv('Dnr2_DC_counts.csv',row.names = 1)
# Analysis and generate panels for Figure 2F----------------------------------------------------------------
fdata <- d2.counts
#IL7R_p_CCR7_p
group <- factor(c(rep("others", 3), rep("IL7R_p_CCR7_p", 1),rep("others", 1)), levels = c("others", "IL7R_p_CCR7_p"))
colData <- data.frame(row.names = colnames(fdata), group)
dds <- DESeqDataSetFromMatrix(fdata, colData, design = ~group)
dds <- DESeq(dds)
res = results(dds, contrast=c("group", "IL7R_p_CCR7_p", "others"))
IL7R_p_CCR7_pvsothers.genes <- as.data.frame(res)
IL7R_p_CCR7_pvsothers.genes$genes <- rownames(IL7R_p_CCR7_pvsothers.genes)
IL7R_p_CCR7_pvsothers.genes <- IL7R_p_CCR7_pvsothers.genes %>% arrange(desc(log2FoldChange)) %>% dplyr::select(genes,log2FoldChange,pvalue,padj) %>% dplyr::filter(pvalue < 0.05)
IL7R_p_CCR7_p_genes <- IL7R_p_CCR7_pvsothers.genes %>% top_n(n = 500, wt = log2FoldChange)
IL7R_p_CCR7_p_score <- list(rownames(IL7R_p_CCR7_p_genes))
d6.object.cdcs <- AddModuleScore(d6.object.cdcs,
                                 features = IL7R_p_CCR7_p_score,
                                 name = "IL7R_p_CCR7_p.score")
VlnPlot(d6.object.cdcs, features = "IL7R_p_CCR7_p.score1",
        cols = c('#C5E524', '#38C2E5','#384C94','#E679C5') ,pt.size = 0)
ggsave('Outputs/figures/IL7Rp_CCR7_p_score_vln.pdf', 
       width = 4.5, height = 4, dpi = 300)
FeaturePlot(d6.object.cdcs, "IL7R_p_CCR7_p.score1")+
  scale_color_viridis_c(option = 4)
ggsave('Outputs/figures/IL7Rp_CCR7_p_score_feature.pdf', 
       width =5 , height = 4, dpi = 300)
#IL7R_n_CCR7_p
group <- factor(c(rep("others", 4), rep("IL7R_n_CCR7_p", 1)), levels = c("others", "IL7R_n_CCR7_p"))
colData <- data.frame(row.names = colnames(fdata), group)
dds <- DESeqDataSetFromMatrix(fdata, colData, design = ~group)
dds <- DESeq(dds)
res = results(dds, contrast=c("group", "IL7R_n_CCR7_p", "others"))
IL7R_n_CCR7_pvsothers.genes <- as.data.frame(res)
IL7R_n_CCR7_pvsothers.genes$genes <- rownames(IL7R_n_CCR7_pvsothers.genes)
IL7R_n_CCR7_pvsothers.genes <- IL7R_n_CCR7_pvsothers.genes %>% arrange(desc(log2FoldChange)) %>% dplyr::select(genes,log2FoldChange,pvalue,padj) %>% dplyr::filter(pvalue < 0.05)
IL7R_n_CCR7_p_genes <- IL7R_n_CCR7_pvsothers.genes %>% top_n(n = 200, wt = log2FoldChange)
IL7R_n_CCR7_p_score <- list(rownames(IL7R_n_CCR7_p_genes))
d6.object.cdcs <- AddModuleScore(d6.object.cdcs,
                                 features = IL7R_n_CCR7_p_score,
                                 name = "IL7R_n_CCR7_p.score")
VlnPlot(d6.object.cdcs, features = "IL7R_n_CCR7_p.score1",
        cols = c('#C5E524', '#38C2E5','#384C94','#E679C5') ,pt.size = 0)
ggsave('Outputs/figures/IL7Rn_CCR7_p_score_vln.pdf', 
       width = 4.5, height = 4, dpi = 300)
FeaturePlot(d6.object.cdcs, "IL7R_n_CCR7_p.score1")+
  scale_color_viridis_c(option = 4)
ggsave('Outputs/figures/IL7Rn_CCR7_p_score_feature.pdf', 
       width = 5, height = 4, dpi = 300)
#IL7R_p_CCR7_n
group <- factor(c(rep("others", 2), rep("IL7R_p_CCR7_n", 1), rep("others", 2)),
                levels = c("others", "IL7R_p_CCR7_n"))
colData <- data.frame(row.names = colnames(fdata), group)
dds <- DESeqDataSetFromMatrix(fdata, colData, design = ~group)
dds <- DESeq(dds)
res = results(dds, contrast=c("group", "IL7R_p_CCR7_n", "others"))
IL7R_p_CCR7_nvsothers.genes <- as.data.frame(res)
IL7R_p_CCR7_nvsothers.genes$genes <- rownames(IL7R_p_CCR7_nvsothers.genes)
IL7R_p_CCR7_nvsothers.genes <- IL7R_p_CCR7_nvsothers.genes %>% arrange(desc(log2FoldChange)) %>% dplyr::select(genes,log2FoldChange,pvalue,padj) %>% dplyr::filter(pvalue < 0.05)
IL7R_p_CCR7_n_genes <- IL7R_p_CCR7_nvsothers.genes %>% top_n(n = 200, wt = log2FoldChange)
IL7R_p_CCR7_n_score <- list(rownames(IL7R_p_CCR7_n_genes))
d6.object.cdcs <- AddModuleScore(d6.object.cdcs, features = IL7R_p_CCR7_n_score, name = "IL7R_p_CCR7_n.score")
VlnPlot(d6.object.cdcs, features = "IL7R_p_CCR7_n.score1",
        cols = c('#C5E524', '#38C2E5','#384C94','#E679C5') ,pt.size = 0)
ggsave('Outputs/figures/IL7Rp_CCR7_n_score_vln.pdf', 
       width = 4.5, height = 4, dpi = 300)
FeaturePlot(d6.object.cdcs, "IL7R_p_CCR7_n.score1")+
  scale_color_viridis_c(option = 4)
ggsave('Outputs/figures/IL7Rp_CCR7_n_score_feature.pdf', 
       width = 5, height = 4, dpi = 300)
# cDC1 (CD141)
group <- factor(c(rep("mDC1", 1), rep("others", 4)), levels = c("others", "mDC1"))
colData <- data.frame(row.names = colnames(fdata), group)
dds <- DESeqDataSetFromMatrix(fdata, colData, design = ~group)
dds <- DESeq(dds)
res = results(dds, contrast=c("group", "mDC1", "others"))
mDC1vsothers.genes <- as.data.frame(res)
mDC1vsothers.genes$genes <- rownames(mDC1vsothers.genes)
mDC1vsothers.genes <- mDC1vsothers.genes %>% arrange(desc(log2FoldChange)) %>% dplyr::select(genes,log2FoldChange,pvalue,padj) %>% dplyr::filter(pvalue < 0.05)
mDC1_genes <- mDC1vsothers.genes %>% top_n(n = 200, wt = log2FoldChange)
mDC1_score <- list(rownames(mDC1_genes))
d6.object.cdcs <- AddModuleScore(d6.object.cdcs, 
                                 features = mDC1_score, 
                                 name = "cDC1.score")
VlnPlot(d6.object.cdcs, features = "cDC1.score1",
        cols = c('#C5E524', '#38C2E5','#384C94','#E679C5') ,pt.size = 0)
ggsave('Outputs/figures/cDC1.score_score_vln.pdf', 
       width = 4.5, height = 4, dpi = 300)
FeaturePlot(d6.object.cdcs, "cDC1.score1")+
  scale_color_viridis_c(option = 4)
ggsave('Outputs/figures/cDC1.score1_score_feature.pdf', 
       width = 5, height = 4, dpi = 300)
# cDC2 (CD1C)
group <- factor(c(rep("others", 1), rep("mDC2", 1), rep("others", 3)), levels = c("others", "mDC2"))
colData <- data.frame(row.names = colnames(fdata), group)
dds <- DESeqDataSetFromMatrix(fdata, colData, design = ~group)
dds <- DESeq(dds)
res = results(dds, contrast=c("group", "mDC2", "others"))
mDC2vsothers.genes <- as.data.frame(res)
mDC2vsothers.genes$genes <- rownames(mDC2vsothers.genes)
mDC2vsothers.genes <- mDC2vsothers.genes %>% arrange(desc(log2FoldChange)) %>% dplyr::select(genes,log2FoldChange,pvalue,padj) %>% dplyr::filter(pvalue < 0.05)
mDC2_genes <- mDC2vsothers.genes %>% top_n(n = 200, wt = log2FoldChange)
mDC2_score <- list(rownames(mDC2_genes))
d6.object.cdcs <- AddModuleScore(d6.object.cdcs, 
                                 features = mDC2_score, 
                                 name = "cDC2.score")
VlnPlot(d6.object.cdcs, features = "cDC2.score1",
        cols = c('#C5E524', '#38C2E5','#384C94','#E679C5') ,pt.size = 0)
ggsave('Outputs/figures/cDC2.score_score_vln.pdf', 
       width = 4.5, height = 4, dpi = 300)
FeaturePlot(d6.object.cdcs, "cDC2.score1")+ 
  scale_color_viridis_c(option = 4)
ggsave('Outputs/figures/cDC2.score1_score_feature.pdf', 
       width = 5, height = 4, dpi = 300)

# Save signatures
write.xlsx(list('IL7RpCCR7p'=IL7R_p_CCR7_p_genes,
                'IL7RnCCR7p'=IL7R_n_CCR7_p_genes,
                'IL7RpCCR7n-'=IL7R_p_CCR7_n_genes,
                'CD1cp'=mDC2_genes,
                'CD141p'=mDC1_genes),
           'Outputs/tables/Signature_genes.xlsx')

# Tissue enrichment -------------------------------------------------------
set2_palette <- colorRampPalette(brewer.pal(8, "Set2"))(35)
  # cDC1
gs <- GeneSet(
  geneIds = mDC1_genes$genes,
  organism = "Homo Sapiens",                 # Specify the organism
  geneIdType = SymbolIdentifier(),           # Define gene ID type as Symbol
  setName = "cDC1Genes"               # Optional: Give a name to the gene set
)
cdc1.tissue.enrich <- teEnrichment(gs)
seEnrichmentOutput<-cdc1.tissue.enrich[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),
                                      row.names = rowData(seEnrichmentOutput)[,1]),
                           colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
# Plotting the P-Values
ggplot(enrichmentOutput,aes(x=reorder(Tissue,Log10PValue),y=Log10PValue,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-Log10PValue')+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
          element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+
  ylim(c(0,2))+
  scale_fill_manual(values = set2_palette)
ggsave('Outputs/figures/cDC1_tissue_enrichment.pdf', width = 10, height = 4)
ggplot(enrichmentOutput,aes(x=reorder(Tissue,Tissue.Specific.Genes),
                            y=Tissue.Specific.Genes,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = 'Tissue.Specific.Genes')+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
          element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = set2_palette)+
  ylim(0,28)
ggsave('Outputs/figures/cDC1_tissue_enrichment_counts.pdf',
       width = 10, height = 4)
  # cDC2
gs <- GeneSet(
  geneIds = mDC2_genes$genes,
  organism = "Homo Sapiens",                 # Specify the organism
  geneIdType = SymbolIdentifier(),           # Define gene ID type as Symbol
  setName = "cDC2Genes"               # Optional: Give a name to the gene set
)
cdc2.tissue.enrich <- teEnrichment(gs)
seEnrichmentOutput<-cdc2.tissue.enrich[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),
                                      row.names = rowData(seEnrichmentOutput)[,1]),
                           colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
# Plotting the P-Values
ggplot(enrichmentOutput,aes(x=reorder(Tissue,Log10PValue),y=Log10PValue,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-Log10PValue')+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
          element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = set2_palette)
ggsave('Outputs/figures/cDC2_tissue_enrichment.pdf', width = 10, height = 4)
ggplot(enrichmentOutput,aes(x=reorder(Tissue,Tissue.Specific.Genes),
                            y=Tissue.Specific.Genes,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = 'Tissue.Specific.Genes')+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
          element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = set2_palette)+
  ylim(0,28)
ggsave('Outputs/figures/cDC2_tissue_enrichment_counts.pdf',
       width = 10, height = 4)
# DC2hm
gs <- GeneSet(
  geneIds = IL7R_p_CCR7_p_genes$genes,
  organism = "Homo Sapiens",                 # Specify the organism
  geneIdType = SymbolIdentifier(),           # Define gene ID type as Symbol
  setName = "DC2hmGenes"               # Optional: Give a name to the gene set
)
dc2hm.tissue.enrich <- teEnrichment(gs)
  # Fetch all hits
dc2hm.enrichment.hits <- c()
for (i in names(dc2hm.tissue.enrich[[2]])) {
  dc2hm.enrichment.hits <- c(dc2hm.enrichment.hits,
                             rownames(dc2hm.tissue.enrich[[2]][[i]]))
}
  # Cross ref with another datasets
genome.res.tra.list <- read.csv('Genome_Res_TRAlist.csv', row.names = 1)

seEnrichmentOutput<-dc2hm.tissue.enrich[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),
                                      row.names = rowData(seEnrichmentOutput)[,1]),
                           colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)

# Plotting the P-Values
ggplot(enrichmentOutput,aes(x=reorder(Tissue,Log10PValue),y=Log10PValue,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-Log10PValue')+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
          element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+
  ylim(c(0,2))+
  scale_fill_manual(values = set2_palette)
ggsave('Outputs/figures/DC2hm_tissue_enrichment.pdf', width = 10, height = 4)
ggplot(enrichmentOutput,aes(x=reorder(Tissue,Tissue.Specific.Genes),
                            y=Tissue.Specific.Genes,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = 'Tissue.Specific.Genes')+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
          element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major= element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = set2_palette)+
  ylim(0,28)
ggsave('Outputs/figures/DC2hm_tissue_enrichment_counts.pdf',
       width = 10, height = 4)

# Save results
write.xlsx(list('cDC1'=setNames(data.frame(assay(cdc1.tissue.enrich[[1]]),
                                           row.names = rowData(seEnrichmentOutput)[,1]),
                                colData(seEnrichmentOutput)[,1])[,1:3],
                'cDC2'=setNames(data.frame(assay(cdc2.tissue.enrich[[1]]),
                                           row.names = rowData(seEnrichmentOutput)[,1]),
                                colData(seEnrichmentOutput)[,1])[,1:3],
                'DC2hm'=setNames(data.frame(assay(dc2hm.tissue.enrich[[1]]),
                                           row.names = rowData(seEnrichmentOutput)[,1]),
                                colData(seEnrichmentOutput)[,1])[,1:3]),
           'Outputs/tables/Tissue_enrichment.xlsx', rowNames=T)


# Box plots for DC2hm upregulated TRAs
merge.tpm.tras <- na.omit(merged.tpm[read.csv('../TSLP-cDC2s vs cDC2s/DC2hm.up.tras.csv')$Gene,])
tissue_specificity_df <- data.frame(
  Gene = read.csv('../TSLP-cDC2s vs cDC2s/DC2hm.up.tras.csv')$Gene, 
  Tissue = read.csv('../TSLP-cDC2s vs cDC2s/DC2hm.up.tras.csv')$`RNA.TS.TPM` 
)
tissue_specificity_df$Tissue <- sub(":.*", "", tissue_specificity_df$Tissue)
tissue_specificity_df <- tissue_specificity_df[tissue_specificity_df$Gene %in% rownames(merge.tpm.tras),]
groups <- data.frame(
  Sample = colnames(merge.tpm.tras),
  Group = c(rep("cDC1", 4),
            rep("cDC2", 4),
            rep("DC2hm", 4))
)
long_df <- as.data.frame(merge.tpm.tras) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(groups, by = "Sample") %>%
  left_join(tissue_specificity_df, by = "Gene")

summary_stats <- long_df %>%
  group_by(Gene, Group) %>%
  summarise(
    Mean = mean(Expression),
    SD = sd(Expression),
    n = n(),
    SE = SD / sqrt(n),
    Tissue = unique(Tissue), 
    .groups = "drop"
  )
pdf('TRA_boxplots.pdf', width = 8, height = 8) 
# Generate and save each gene's plot in the PDF
unique_genes <- unique(long_df$Gene)
for (gene in unique_genes) {
  gene_data <- long_df %>% filter(Gene == gene)
  gene_stats <- summary_stats %>% filter(Gene == gene)
  tissue_label <- unique(gene_stats$Tissue)
  
  # Create the plot
  p <- ggplot() +
    geom_bar(data = gene_stats, aes(x = Group, y = Mean, fill = Group), 
             stat = "identity", color = "black", alpha = 0.7) +
    geom_errorbar(data = gene_stats, aes(x = Group, ymin = Mean - SE, ymax = Mean + SE),
                  width = 0.2, color = "black") +
    geom_jitter(data = gene_data, aes(x = Group, y = Expression), 
                width = 0.2, color = "black", size = 1.5) +
    labs(
      title = paste(tissue_label, ':', gene),
      x = "Group",
      y = "Expression Level"
    ) +
    theme_classic() +
    theme(legend.position = "none") # Remove legend for simplicity
  print(p)
}
dev.off()

