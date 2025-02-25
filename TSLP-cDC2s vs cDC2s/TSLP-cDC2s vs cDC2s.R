library(pacman)
p_load(sva, limma, pheatmap, ggplot2, ggplotify, dplyr, openxlsx,
       clusterProfiler,tidyverse, ggrepel, Seurat, ggsci, viridis,
       ggrastr, enrichplot) 
set.seed(0822)
setwd('/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/TSLP-cDC2s vs cDC2s')
# TSLP cDC2s vs fresh cDC2s -----------------------------------------------
  # Read count matrix, note that this experiments also include two inhibitor-treated groups
counts.tslp.vs.fresh <- read.csv('Counts_combined.csv', row.names = 1)
counts.tslp.vs.fresh <- counts.tslp.vs.fresh[rowSums(counts.tslp.vs.fresh)>60,]
counts.tslp.vs.fresh <- counts.tslp.vs.fresh[,c(1:6)]
  # Annotate groups and perform limma-voom DE analysis
group <- factor(c(rep('cDC2s',3),
                  rep('TSLP_cDC2s',3)))
mm <- model.matrix(~0+group)
y <- voom(counts.tslp.vs.fresh,
          mm, plot = T)
fit <- lmFit(y, mm)
contr.tslp.vs.fresh <- makeContrasts(groupTSLP_cDC2s-groupcDC2s,
                                   levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr.tslp.vs.fresh)
tmp <- eBayes(tmp)
de.result.tslp.vs.fresh <- topTable(tmp, sort.by = 'logFC', n = Inf)
de.result.tslp.vs.fresh <- de.result.tslp.vs.fresh[order(de.result.tslp.vs.fresh$logFC,
                                                         decreasing = T),]
write.xlsx(list('DE results'=de.result.tslp.vs.fresh[order(de.result.tslp.vs.fresh$logFC, decreasing = T),],
                'Upregulated'=de.result.tslp.vs.fresh[de.result.tslp.vs.fresh$logFC > 1 & de.result.tslp.vs.fresh$adj.P.Val < 0.05,],
                'Downregulated'=de.result.tslp.vs.fresh[de.result.tslp.vs.fresh$logFC < (-1) & de.result.tslp.vs.fresh$adj.P.Val < 0.05,]),
           'Outputs/tables/TSLP_vs_Fresh_cDC2_DEGs.xlsx', rowNames =T)
  # Figure 5G, barplot showing DC2hm marker and ISG induction by TSLP
gene.to.show <- read.csv('Highlighted_gene_selected.csv')
de.result.tslp.vs.fresh.show <- de.result.tslp.vs.fresh[gene.to.show$Gene,]
de.result.tslp.vs.fresh.show$Gene <- factor(rownames(de.result.tslp.vs.fresh.show),
                                              levels = rownames(de.result.tslp.vs.fresh.show))
ggplot(de.result.tslp.vs.fresh.show,
       aes(y=logFC,x=Gene, fill=-log10(adj.P.Val)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  scale_fill_gradient(low = '#F686BD', high = '#d930a6')+
  labs(y='Log2FC: TSLP-cDC2s vs Fresh cDC2s', x=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
  geom_hline(yintercept = c(0), color='black')
ggsave('Outputs/figures/TSLP-cDC2s_vs_fresh_cDC2_markers_barplot.pdf',
       width = 7.2, height = 2)
# TSLP signature score ----------------------------------------------------
  # Read Donor 6 Seurat object
d6.object <- readRDS('../scRNA-seq analysis/Outputs/rds/Dnr6_cDC2s.rds')
d6.object <- AddModuleScore(d6.object,
                            features = list(rownames(de.result.tslp.vs.fresh)[1:500]),
                            name = 'tslp.score')
  # Figure 5I,  of TSLP-cDC2 signatures
rasterize(FeaturePlot(d6.object, features = 'tslp.score1',pt.size = 0.1)+
            scale_color_viridis_c(option = 4),
          dpi = 300)
ggsave('Outputs/figures/TSLP_score_feature_plot.pdf',
       width = 4, height = 3)
VlnPlot(d6.object, features = 'tslp.score1',pt.size = 0.0,
        cols =c('cDC2'='#38C2E5','DC2pre-hm'='#384C94','DC2hm'='#E679C5'))
ggsave('Outputs/figures/TSLP_signature_violinplot.pdf', width = 4, height = 3)
# Med-cDC2s vs Fresh-cDC2s -------------------------------------------------------------
  # Read count matrix, note that this experiment also include two inhibitor-treated groups
counts.tslp.vs.med <- read.csv('Counts_with_Med-cDC2s.csv', row.names = 1)
counts.tslp.vs.med  <- counts.tslp.vs.med [,c('CT_1',  'CT_2',  'CT_3')]
colnames(counts.tslp.vs.med ) <- c('Med_1','Med_2','Med_3')
counts.tslp.vs.med  <- counts.tslp.vs.med [intersect(rownames(counts.tslp.vs.fresh),
                                                     rownames(counts.tslp.vs.med )),]
counts.med.vs.fresh <- cbind(counts.tslp.vs.med,
                             counts.tslp.vs.fresh[intersect(rownames(counts.tslp.vs.fresh),
                                                            rownames(counts.tslp.vs.med )),c(1,2,3)])
colnames(counts.med.vs.fresh) <- c('Med_1', 'Med_2', 'Med_3',
                                   'cDC2_1',  'cDC2_2',  'cDC2_3')
group <- factor(c(rep('Med_cDC2s',3),
                  rep('Fresh_cDC2',3)))
mm <- model.matrix(~0+group)
y <- voom(counts.med.vs.fresh, mm, plot = T)
fit <- lmFit(y, mm)
contr.tslp.vs.ctr <- makeContrasts(groupMed_cDC2s-groupFresh_cDC2,
                                   levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr.tslp.vs.ctr)
tmp <- eBayes(tmp)
de.result.med.vs.fresh <- topTable(tmp, sort.by = 'logFC', n = Inf)
write.xlsx(list('DE results'=de.result.med.vs.fresh,
                'Upregulated'=de.result.med.vs.fresh[de.result.med.vs.fresh$logFC>1 & de.result.med.vs.fresh$P.Value < 0.05,],
                'Downregulated'=de.result.med.vs.fresh[de.result.med.vs.fresh$logFC<(-1) & de.result.med.vs.fresh$P.Value < 0.05,]),
           'Outputs/tables/Med_cDC2s_vs_Fresh_cDC2_DEGs.xlsx',rowNames=T)
  # Figure S6B, barplot comparing the log2FCs
tslp.inducd.gene <- c('CD80', 'CD40', 'CD86', 'FSCN1', 'CD44', 'ICAM1',
                      'REL','CRLF2', 'IL7R', 'STAT5A',  'CCL17',
                      'SOCS2', 'CD274', 'CD200', 'ALDH1A2', 'PVR',
                      'DLL4', 'CCL22'
                      )
lfc.comparison.df <- rbind(de.result.med.vs.fresh[tslp.inducd.gene,],
                           de.result.tslp.vs.fresh[tslp.inducd.gene,])
lfc.comparison.df$Comparison <- c(rep('Med-cDC2s',length(tslp.inducd.gene)),
                                  rep('TSLP-cDC2s', length(tslp.inducd.gene)))
lfc.comparison.df$Gene <- factor(rep(tslp.inducd.gene,2),
                                 levels = tslp.inducd.gene)
ggplot(lfc.comparison.df, aes(x=Gene, y=logFC, fill=Comparison))+
  geom_bar(stat = 'identity', position = 'dodge',width = 0.5)+
  scale_fill_manual(values = c('#ffb3c6','#fb6f92'))+
  theme_classic()+
  geom_hline(yintercept = c(1), linetype='dashed',color='grey')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('Outputs/figures/Barplot_LFC_comparison_tslp_induced.pdf', 
       width = 6, height = 2.5) 



# TSLP-cDC2s vs Med-cDC2s -------------------------------------------------
counts.tslp.vs.med <- read.csv('Counts_with_Med-cDC2s.csv', row.names = 1)
counts.tslp.vs.med <- counts.tslp.vs.med[,c(1:6)]
colnames(counts.tslp.vs.med) <- c("Med_1",   "Med_2",   "Med_3",
                                  "TSLP_1", "TSLP_2", "TSLP_3")
group <- factor(c(rep('Med_cDC2s',3),
                  rep('TSLP_cDC2s',3)))
mm <- model.matrix(~0+group)
y <- voom(counts.tslp.vs.med, mm, plot = T)
fit <- lmFit(y, mm)
contr.tslp.vs.med <- makeContrasts(groupTSLP_cDC2s-groupMed_cDC2s,
                                   levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr.tslp.vs.med)
tmp <- eBayes(tmp)
de.result.tslp.vs.med <- topTable(tmp, sort.by = 'logFC', n = Inf)
write.xlsx(list('DE results'=de.result.tslp.vs.med,
                'Upregulated'=de.result.tslp.vs.med[de.result.tslp.vs.med$logFC>1 & de.result.tslp.vs.med$P.Value < 0.05,],
                'Downregulated'=de.result.tslp.vs.med[de.result.tslp.vs.med$logFC<(-1) & de.result.tslp.vs.med$P.Value < 0.05,]),
           'Outputs/tables/TSLP_cDC2s_vs_Med_cDC2_DEGs.xlsx',rowNames=T)
# GSEA --------------------------------------------------------------------
  # Figure 5H, GSEA
de.result.tslp.vs.fresh <- de.result.tslp.vs.fresh[order(de.result.tslp.vs.fresh$logFC,
                                                         decreasing = T),]
de.result.tslp.vs.fresh.ranked <-de.result.tslp.vs.fresh$logFC
names(de.result.tslp.vs.fresh.ranked) <- rownames(de.result.tslp.vs.fresh)
tslp.vs.fresh <- GSEA(de.result.tslp.vs.fresh.ranked,
                      TERM2GENE = rbind(read.gmt('LINDSTEDT_DENDRITIC_CELL_MATURATION_B.v2023.2.Hs.gmt'),
                                        read.gmt('WIERENGA_STAT5A_TARGETS_UP.v2023.2.Hs.gmt'),
                                        read.gmt('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING.v2023.2.Hs.gmt')))
gseaplot2(tslp.vs.fresh,geneSetID = 1:3,
          color = c('#FFBE0B','#00BBF9','#F15BB5'))
ggsave('Outputs/figures/TSLP-cDC2s_vs_cDC2s_GSEA.pdf',
       width = 4, height = 4)


# STAT5 inhibitor ---------------------------------------------------------
tslp.inhibitor.counts <- read.csv('Counts_matrix_TSLP_and_inhibitors.csv',
                                  row.names = 1)
group <- factor(c(rep('Ctr',3),
                  rep('TSLP',3),
                  rep('TSLP_c_rel_inhibitor',3),
                  rep('TSLP_STAT5_inhibitor',3)))
mm <- model.matrix(~0+group)
y <- voom(tslp.inhibitor.counts, mm, plot = T)
fit <- lmFit(y, mm)
contr.stat5.inhibitor.tslp <- makeContrasts(groupTSLP_STAT5_inhibitor-groupTSLP,
                                            levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr.stat5.inhibitor.tslp)
tmp <- eBayes(tmp)
de.result.stat5.inhibitor.tslp <- topTable(tmp, sort.by = 'logFC', n = Inf)
de.result.stat5.inhibitor.tslp <- de.result.stat5.inhibitor.tslp[order(de.result.stat5.inhibitor.tslp$logFC,
                                                                       decreasing = T),]
# Store DE analysis results
write.xlsx(list('DE results'=de.result.stat5.inhibitor.tslp[order(de.result.stat5.inhibitor.tslp$logFC, decreasing = T),],
                'Upregulated'=de.result.stat5.inhibitor.tslp[de.result.stat5.inhibitor.tslp$logFC > 0.25 & de.result.stat5.inhibitor.tslp$adj.P.Val < 0.05,],
                'Downregulated'=de.result.stat5.inhibitor.tslp[de.result.stat5.inhibitor.tslp$logFC < -(0.25) & de.result.stat5.inhibitor.tslp$adj.P.Val < 0.05,]),
           'Outputs/tables/STAT5_inhibitor_vs_TSLP.xlsx', rowNames =T)
  # plot genes to compare

de.result.stat5.inhibitor.tslp.show <- de.result.stat5.inhibitor.tslp[gene.to.show$Gene[4:26],]
de.result.stat5.inhibitor.tslp.show$Gene <- factor(rownames(de.result.stat5.inhibitor.tslp.show),
                                                   levels = rownames(de.result.stat5.inhibitor.tslp.show))
ggplot(de.result.stat5.inhibitor.tslp.show,
       aes(x=Gene, y=logFC, fill=-log10(adj.P.Val)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  scale_fill_gradient(low = '#B096ED',high = '#471ca8')
ggsave('Outputs/figures/STAT5_inhibiting_barplot.pdf', width = 6, height = 2.3)

