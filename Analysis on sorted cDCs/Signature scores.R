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
setwd("Analysis on sorted cDCs")
# Read the counts and the Seurat d6.object.cdcs -----------------------------------
d6.object.cdcs <- readRDS('Dnr6_cDCs.rds') # Available @ GSE267253
d2.counts <- read.csv('Dnr2_DC_counts.csv',row.names = 1) # Available @ GSE287912
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
