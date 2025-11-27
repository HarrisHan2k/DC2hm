library(pheatmap)
library(ggplotify)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(openxlsx)
setwd("Analysis on sorted cDCs")
# Plot 4 individual data analysis ------------------------
  # Read TPMs of different donors, available @ GSE287912
dnr2.tpm <- read.csv('Donor 2 5cDC TPM.csv',
                     row.names = 1, check.names = F)
dnr12.tpm <- read.csv('Donor 12 5cDC TPM.csv',
                      row.names = 1, check.names = F)
dnr1.tpm <- read.csv('Donor 1 4cDC TPM.csv',
                     row.names = 1, check.names = F)
dnr6.tpm <- read.csv('Donor 6 4cDC TPM.csv',
                     row.names = 1, check.names = F)
tpm.list <- list(dnr2.tpm,
                 dnr12.tpm,
                 dnr1.tpm,
                 dnr6.tpm)
donor.vec <- c('2','12','1','6')
genes.to.plot <- read.csv('Highlighted_gene_selected.csv')$Gene # Refer gene lists to the final version of the publication
for (i in seq(1,4)) {
  as.ggplot(pheatmap(tpm.list[[i]][genes.to.plot,],
                     cluster_cols = F, cluster_rows = F,
                     scale = 'row', cellwidth = 20, cellheight = 20))
  ggsave(paste0('Outputs/figures/Donor_',donor.vec[i],'_heatmap.pdf'),
         width = 10, height = 10)
}

# Plot TLR, inflammasome, and Treg induction genes ----------------------
  # Read batch-corrected TPM files
merged.tpm <- read.csv('cDC_adjusted_TPM.csv', row.names = 1,
                       check.names = F) # available @ GSE287912
  # Add cell group annotations
group.annotation.4dnrs <- data.frame('Celltype'=factor(c(rep('cDC2', 4),
                                                         rep('DC2hm', 4)),
                                                       levels = c('cDC2', 'DC2hm')),
                                     row.names = colnames(merged.tpm)[5:12])
  # TLR and adpaters
tlr.genes <- c('TLR1', 'TLR2', 'TLR3', 'TLR4', 'TLR6',
               'TLR7', 'TLR8', 'TLR9', 'MYD88', 'TICAM1')
as.ggplot(pheatmap(merged.tpm[tlr.genes,5:12],
                   cluster_rows = F, cluster_cols = F,
                   scale = 'row', cellwidth = 10, cellheight = 10, 
                   annotation_col = group.annotation.4dnrs,
                   annotation_colors = list('Celltype'=c('cDC2'='#38C2E5',
                                                         'DC2hm'='#E679C5'))))
ggsave('Outputs/figures/TLR_heatmap.pdf')
  # Inflammasome, IL-1beta secretion, and pyroptosis
inflammasome.genes <- c('NLRP1','NLRP3','CARD8','DDX3X',
                        'EIF2AK2','PYCARD','IL1B','IL18','CASP1','DHX33','CARD16',
                        'GSDMD')
as.ggplot(pheatmap(merged.tpm[inflammasome.genes,5:12], 
                   cluster_rows = F, cluster_cols = F,
                   scale = 'row', cellwidth = 10, cellheight = 10, 
                   annotation_col = group.annotation.4dnrs,
                   annotation_colors = list('Celltype'=c('cDC2'='#38C2E5',
                                                         'DC2hm'='#E679C5'))))
ggsave('Outputs/figures/Inflammasome_heatmap.pdf')
  # Treg induction
treg.induction <- c('ITGAV','ITGB8','IDO1', 'IDO2','PVR',
                    'DLL4','EBI3',
                    'ALDH1A1','ALDH1A2', 'ICOSLG', 'CCL22')
as.ggplot(pheatmap(merged.tpm[treg.induction,5:12],
                   cluster_rows = F, cluster_cols = F,
                   scale = 'row', cellwidth = 10, cellheight = 10, 
                   annotation_col = group.annotation.4dnrs,
                   annotation_colors = list('Celltype'=c('cDC2'='#38C2E5',
                                                         'DC2hm'='#E679C5'))))
ggsave('Outputs/figures/Treg_induction_heatmap.pdf')

 
