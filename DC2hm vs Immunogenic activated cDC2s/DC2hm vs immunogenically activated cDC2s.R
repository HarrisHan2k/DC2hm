library(pacman)
p_load(sva, limma, pheatmap, ggplot2, ggplotify, dplyr, openxlsx,
       clusterProfiler,tidyverse, ggrepel, Seurat, ggsci, ComplexHeatmap,
       circlize, enrichplot, TissueEnrich, RColorBrewer, ggplotify) 
set.seed(0822)
setwd('/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/DC2hm vs Immunogenic activated cDC2s')
# DE analysis on immunogenically activated cDC2s --------------------------
  # PAMP-cDC2s vs Ctr cDC2s
    # This is a batch-corrected (ComBat) matrix from 2 experiments including one with LPS+IFNg group that was anaylzed below
pamp.cdc2s.counts <- read.csv('Combined_PAMP_counts_batch_corrected.csv',
                              row.names = 1, check.names = F)
  
group <- c(rep('cDC2_Ctr',5),
           rep('cDC2_PAMP',6))
mm <- model.matrix(~0+group)
y <- voom(pamp.cdc2s.counts, mm, plot = T)
fit <- lmFit(y, mm)
pamp.vs.ctr.contrast <- makeContrasts(groupcDC2_PAMP-groupcDC2_Ctr,
                                      levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, pamp.vs.ctr.contrast)
tmp <- eBayes(tmp)
pamp.vs.ctr.result <- topTable(tmp, sort.by = 'logFC',n=Inf)
write.xlsx(list('DE results'=pamp.vs.ctr.result,
                'Upregulated'=pamp.vs.ctr.result[pamp.vs.ctr.result$logFC>1 & pamp.vs.ctr.result$adj.P.Val < 0.05,],
                'Downregulated'=pamp.vs.ctr.result[pamp.vs.ctr.result$logFC<(-1) & pamp.vs.ctr.result$adj.P.Val < 0.05,]),
           'Outputs/tables/PAMP-cDC2s_vs_Med-cDC2s.xlsx',rowNames=T)
  # LPS+IFN-g vs Ctr cDC2s
    # A PAMPs-stimulated group is also contained in this experiments
ifng.cdc2.counts <- read.csv('LPS+IFNg_or_PAMP_vs_Ctr_cDC2s_counts.csv',
                             row.names = 1, check.names = F)
group <- factor(c(rep('Ctr',3),
                  rep('IFNG_LPS',3),
                  rep('PAMP',3)))
mm <- model.matrix(~0+group)
y <- voom(ifng.cdc2.counts, mm, plot = T)
fit <- lmFit(y, mm)
# comparison
contr.ifng.lps.vs.ctr <- makeContrasts(groupIFNG_LPS-groupCtr,
                                       levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr.ifng.lps.vs.ctr)
tmp <- eBayes(tmp)
ifng.lps.result <- topTable(tmp, sort.by = 'P', n = Inf)
write.xlsx(list('DE results'=ifng.lps.result,
                'Upregulated'=ifng.lps.result[ifng.lps.result$logFC>1 & ifng.lps.result$adj.P.Val < 0.05,],
                'Downregulated'=ifng.lps.result[ifng.lps.result$logFC<(-1) & ifng.lps.result$adj.P.Val < 0.05,]),
           'Outputs/tables/LPS+IFNG-cDC2s_vs_Med-cDC2s.xlsx',
           rowNames=T)


# DE analysis on DC2hm ----------------------------------------------------
  # Single-cell results
    # Read rds files
d6.object <- readRDS('../scRNA-seq analysis/Outputs/rds/Dnr6_cDC2s.rds')
    # Perform DE analysis and save results
scrna.degs <- FindMarkers(d6.object, ident.1 = 'DC2hm', ident.2 = 'cDC2',
                          min.pct = 0, logfc.threshold = 0)
write.xlsx(list('DE results'=scrna.degs,
                'Upregulated'=scrna.degs[scrna.degs$avg_log2FC>1 & scrna.degs$p_val_adj < 0.05,],
                'Downregulated'=scrna.degs[scrna.degs$avg_log2FC<(-1) & scrna.degs$p_val_adj < 0.05,]),
           'Outputs/tables/DC2hm_markers_scRNA.xlsx',rowNames=T)
    # Assign and combine genes
pamp.up <- read.xlsx('Outputs/tables/PAMP-cDC2s_vs_Med-cDC2s.xlsx',
                     sheet = 'Upregulated', rowNames = T)
pamp.down <- read.xlsx('Outputs/tables/PAMP-cDC2s_vs_Med-cDC2s.xlsx',
                       sheet = 'Downregulated', rowNames = T)
pamp.all <- read.xlsx('Outputs/tables/PAMP-cDC2s_vs_Med-cDC2s.xlsx',
                      sheet = 'DE results', rowNames = T)
ifng.lps.up <- read.xlsx('Outputs/tables/LPS+IFNG-cDC2s_vs_Med-cDC2s.xlsx',
                         sheet = 'Upregulated', rowNames = T)
ifng.lps.down <- read.xlsx('Outputs/tables/LPS+IFNG-cDC2s_vs_Med-cDC2s.xlsx',
                           sheet = 'Downregulated', rowNames = T)
ifng.lps.all <- read.xlsx('Outputs/tables/LPS+IFNG-cDC2s_vs_Med-cDC2s.xlsx',
                          sheet = 'DE results', rowNames = T)
dc2hm.scrna <- read.xlsx('Outputs/tables/DC2hm_markers_scRNA.xlsx',
                         rowNames = T, sheet = 'DE results')
    # Assign genes to show
selected.genes <- c('ITGAX','CD1C','CCR7','CD40','CD80','LAMP3','CD274','HIVEP1',
                    'STAT5A','IL7R','CRLF2','REL','CCL22','ICOSLG','SREBF2',
                    'LPL','HMGCS1','SQLE','RELA','ISG15','OAS1','MX1','IFIT2',
                    'STAT1','IFITM2','IFI44L')
dc2hm.scrna.renamed <- dc2hm.scrna[,c("avg_log2FC", "p_val_adj")]
colnames(dc2hm.scrna.renamed) <- c("logFC", "adj.P.Val")
combined.df.scrna <- rbind(dc2hm.scrna.renamed,
                           pamp.all[,c("logFC", "adj.P.Val")],
                           ifng.lps.all [,c("logFC", "adj.P.Val")])
combined.df.scrna$Comparison <- c(rep('DC2hm vs cDC2s',
                                      nrow(dc2hm.scrna.renamed)),
                                  rep('PAMP-cDC2s vs cDC2s', 
                                      nrow(pamp.all)),
                                  rep('LPS+IFNG-cDC2s vs cDC2s',
                                      nrow(ifng.lps.all)))
combined.df.scrna$Gene <- c(rownames(dc2hm.scrna), rownames(pamp.all),
                            rownames(ifng.lps.all))
combined.df.scrna.vis <- combined.df.scrna[combined.df.scrna$Gene %in% selected.genes,]
combined.df.scrna.vis$Gene <- factor(combined.df.scrna.vis$Gene,
                                     levels = selected.genes)
combined.df.scrna.vis$Comparison <- factor(combined.df.scrna.vis$Comparison,
                                           levels = c('DC2hm vs cDC2s',
                                                      'PAMP-cDC2s vs cDC2s',
                                                      'LPS+IFNG-cDC2s vs cDC2s'))
ggplot(combined.df.scrna.vis,
       aes(x=Gene,y=logFC, fill=Comparison))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_npg()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust = .5))+
  geom_hline(yintercept = c(-1,1), linetype='dashed', color='lightgrey')
ggsave('Outputs/figures/Barplot_lfc_DC2hm_scRNA_vs_immunogenic.pdf', 
       width = 8, height = 3)
  # Sorted cDCs
sorted.dc2hm.counts <- read.csv('DC2hm_cDC2_batch_corrected_counts.csv',
                                row.names = 1, check.names = F)
    # Filter lowly expressed genes
sorted.dc2hm.counts <- sorted.dc2hm.counts[rowSums(sorted.dc2hm.counts)>40, ]
group <- c(rep('cDC2',4),
           rep('DC2hm',4))
mm <- model.matrix(~0+group)
y <- voom(sorted.dc2hm.counts, mm, plot = T)
fit <- lmFit(y, mm)
dc2hm.vs.cdc2.contrast <- makeContrasts(groupDC2hm-groupcDC2,
                                      levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, dc2hm.vs.cdc2.contrast)
tmp <- eBayes(tmp)
dc2hm.vs.cdc2.de.results <- topTable(tmp, sort.by = 'logFC',n=Inf)
write.xlsx(list('DE results'=dc2hm.vs.cdc2.de.results,
                'Upregulated'=dc2hm.vs.cdc2.de.results[dc2hm.vs.cdc2.de.results$logFC>1 & dc2hm.vs.cdc2.de.results$adj.P.Val < 0.05,],
                'Downregulated'=dc2hm.vs.cdc2.de.results[dc2hm.vs.cdc2.de.results$logFC<(-1) & dc2hm.vs.cdc2.de.results$adj.P.Val < 0.05,]),
           'Outputs/tables/DC2hm_vs_cDC2s.xlsx',rowNames=T)

dc2hm.all <- read.xlsx('Outputs/tables/DC2hm_vs_cDC2s.xlsx',
                                      sheet = 'DE results', rowNames = T)
dc2hm.up <- read.xlsx('Outputs/tables/DC2hm_vs_cDC2s.xlsx',
                       sheet = 'Upregulated', rowNames = T)
# Write all DE results
write.xlsx(list('DE results DC2hm scRNA'=scrna.degs,
                'Upregulated in DC2hm scRNA'=scrna.degs[scrna.degs$avg_log2FC>1 & scrna.degs$p_val_adj < 0.05,],
                'Downregulated in DC2hm scRNA'=scrna.degs[scrna.degs$avg_log2FC<(-1) & scrna.degs$p_val_adj < 0.05,],
                'DE results DC2hm'=dc2hm.vs.cdc2.de.results,
                'Upregulated in DC2hm'=dc2hm.vs.cdc2.de.results[dc2hm.vs.cdc2.de.results$logFC>1 & dc2hm.vs.cdc2.de.results$adj.P.Val < 0.05,],
                'Downregulated in DC2hm'=dc2hm.vs.cdc2.de.results[dc2hm.vs.cdc2.de.results$logFC<(-1) & dc2hm.vs.cdc2.de.results$adj.P.Val < 0.05,],
                'DE results PAMPs'=pamp.vs.ctr.result,
                'Upregulated in PAMPs'=pamp.vs.ctr.result[pamp.vs.ctr.result$logFC>1 & pamp.vs.ctr.result$adj.P.Val < 0.05,],
                'Downregulated in PAMPs'=pamp.vs.ctr.result[pamp.vs.ctr.result$logFC<(-1) & pamp.vs.ctr.result$adj.P.Val < 0.05,],
                'DE results LPS+IFNg'=ifng.lps.result,
                'Upregulated in LPS+IFNg'=ifng.lps.result[ifng.lps.result$logFC>1 & ifng.lps.result$adj.P.Val < 0.05,],
                'Downregulated in LPS+IFNg'=ifng.lps.result[ifng.lps.result$logFC<(-1) & ifng.lps.result$adj.P.Val < 0.05,]
                ), 'Outputs/tables/Homeostatic and immunogenic DE results.xlsx',
           rowNames=T)
# Comparison --------------------------------------------------------------
# Assign genes to show
selected.genes <- c('CD1C','FCER1A','CLEC10A','CD80','LAMP3','CD40',
                    'CCR7','FSCN1','MARCKS',
                    'CD274','CD200','FAS',
                    'CCL22','ICOSLG','EBI3',
                    'IL7R','CRLF2','STAT5A',
                    'REL','RELB','NFKB1','ISG15','OAS1','MX1',
                    'IFI16','IFIT2') 
  # Figure 2A
dc2hm.scrna.renamed <- dc2hm.scrna[,c("avg_log2FC", "p_val_adj")]
colnames(dc2hm.scrna.renamed) <- c("logFC", "adj.P.Val")
combined.df.scrna <- rbind(dc2hm.scrna.renamed,
                           pamp.all[,c("logFC", "adj.P.Val")],
                           ifng.lps.all [,c("logFC", "adj.P.Val")])
combined.df.scrna$Comparison <- c(rep('DC2hm vs cDC2s',
                                      nrow(dc2hm.scrna.renamed)),
                                  rep('PAMP-cDC2s vs cDC2s', 
                                      nrow(pamp.all)),
                                  rep('LPS+IFNG-cDC2s vs cDC2s',
                                      nrow(ifng.lps.all)))
combined.df.scrna$Gene <- c(rownames(dc2hm.scrna), rownames(pamp.all),
                            rownames(ifng.lps.all))
combined.df.scrna.vis <- combined.df.scrna[combined.df.scrna$Gene %in% selected.genes,]
combined.df.scrna.vis$Gene <- factor(combined.df.scrna.vis$Gene,
                                     levels = selected.genes)
combined.df.scrna.vis$Comparison <- factor(combined.df.scrna.vis$Comparison,
                                           levels = c('DC2hm vs cDC2s',
                                                      'PAMP-cDC2s vs cDC2s',
                                                      'LPS+IFNG-cDC2s vs cDC2s'))
ggplot(combined.df.scrna.vis,
       aes(x=Gene,y=logFC, fill=Comparison))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_npg()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust = .5))
ggsave('Outputs/figures/Barplot_lfc_DC2hm_scRNA_vs_immunogenic.pdf', 
       width = 8, height = 3)
  # Figure S5E
  # Bulk RNA-seq yield depth advantage, several lowly detected ones in scRNA-seq could be visualized here
selected.genes <- c('ITGAX','CD1C', 'CCR7','CD40','CD80','LAMP3','CD274',
                    'HIVEP1','STAT5A','IL7R','CRLF2', 'REL','CCL22',
                    'ICOSLG','SREBF2','LPL','HMGCS1','SQLE','RELA',
                    'ISG15','OAS1','MX1','STAT1','IFITM2','IFI44L') 
combined.df.sorted.dc2hm <- rbind(dc2hm.all[,c("logFC", "adj.P.Val")],
                           pamp.all[,c("logFC", "adj.P.Val")],
                           ifng.lps.all [,c("logFC", "adj.P.Val")])
combined.df.sorted.dc2hm$Comparison <- c(rep('DC2hm vs cDC2s',
                                      nrow(dc2hm.all)),
                                  rep('PAMP-cDC2s vs cDC2s', 
                                      nrow(pamp.all)),
                                  rep('LPS+IFNG-cDC2s vs cDC2s',
                                      nrow(ifng.lps.all)))
combined.df.sorted.dc2hm$Gene <- c(rownames(dc2hm.all),
                                   rownames(pamp.all),
                                   rownames(ifng.lps.all))
combined.df.sorted.dc2hm.vis <- combined.df.sorted.dc2hm[combined.df.sorted.dc2hm$Gene %in% selected.genes,]
combined.df.sorted.dc2hm.vis$Gene <- factor(combined.df.sorted.dc2hm.vis$Gene,
                                     levels = selected.genes)
combined.df.sorted.dc2hm.vis$Comparison <- factor(combined.df.sorted.dc2hm.vis$Comparison,
                                           levels = c('DC2hm vs cDC2s',
                                                      'PAMP-cDC2s vs cDC2s',
                                                      'LPS+IFNG-cDC2s vs cDC2s'))
ggplot(combined.df.sorted.dc2hm.vis,
       aes(x=Gene,y=logFC, fill=Comparison))+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_npg()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,  hjust = 1, vjust = .5))+
  geom_hline(yintercept = c(-1,1), linetype='dashed', color='lightgrey')
ggsave('Outputs/figures/Barplot_lfc_DC2hm_sorted_vs_immunogenic.pdf', 
       width = 8, height = 3)

# Scatter plot showing the fold change --------------------------------
    # Align the genes
genes <- intersect(rownames(dc2hm.all), c(rownames(pamp.all)))
genes <- intersect(genes,rownames(ifng.lps.all))
dc2hm.all <- dc2hm.all[genes,]
dc2hm.up <- dc2hm.up[genes,]
dc2hm.down <- dc2hm.down[genes,]
pamp.all <- pamp.all[genes,]
pamp.up <- pamp.up[genes,]
pamp.down <- pamp.down[genes,]
ifng.lps.all <- ifng.lps.all[genes,]
ifng.lps.up <- ifng.lps.up[genes,]
ifng.lps.down <- ifng.lps.down[genes,]
  # Figure S7B, scatter plot showing the different upregulated genes
df.show <- as.data.frame(cbind(dc2hm.all[genes,1],
                               pamp.all[genes,1]),row.names = genes)
colnames(df.show) <- c('logFC DC2hm', 'logFC PAMP-cDC2')
df.show$regulation <- ifelse(rownames(df.show)%in% intersect(rownames(dc2hm.up),
                                                             rownames(pamp.up)),
                             'Common up',
                             ifelse(rownames(df.show) %in% rownames(dc2hm.up),
                                    'Homeostatic up', 
                                    ifelse(rownames(df.show) %in% rownames(pamp.up),
                                           'Immunogenic up', 'Not upregulated')))
df.show$Gene <- rownames(df.show)
genes.to.annotate <- c('CCR7','CD80','CD40','IDO1','CD274','LAMP3','CD200','CCL19',
                       'CCL22','IL7R','JAG1','ICOSLG','SREBF2','CRLF2','MX2',
                       'ISG15','IFIT2','IFITM3','OASL','STAT1','OAS1',
                       'MX1')
rasterize(ggplot(df.show, aes(x=`logFC DC2hm`, `logFC PAMP-cDC2`, color=regulation))+
            geom_jitter(size=0.1)+
            theme_bw()+
            scale_color_manual(values = c('#FFBE0B','#F15BB5','#00BBF9','grey'))+
            ylim(-8,8)+
            xlim(-8,8)+
            ggrepel::geom_text_repel(
              data = subset(df.show, Gene %in% genes.to.annotate), # Filter to annotate specific genes
              aes(label = Gene), 
              size = 3, 
              box.padding = 0.5, 
              point.padding = 0.1, 
              segment.size = 0.3,
              max.overlaps = Inf, 
              force = 3,
              nudge_x = .2,
              nudge_y = .2,
              direction = 'both'
              ),
          dpi = 300)
ggsave('Outputs/figures/DC2hm_vs_PAMP_2-D_upregulation.pdf',
       width = 11, height = 8)

  # Figure S7C, scatter plot showing the different downregulated genes
df.show <- as.data.frame(cbind(dc2hm.all[genes,1],
                               pamp.all[genes,1]),row.names = genes)
colnames(df.show) <- c('logFC DC2hm', 'logFC PAMP-cDC2')
df.show$regulation <- ifelse(rownames(df.show)%in% intersect(rownames(dc2hm.down),
                                                             rownames(pamp.down)),
                             'Common down',
                             ifelse(rownames(df.show) %in% rownames(dc2hm.down),
                                    'Homeostatic down', 
                                    ifelse(rownames(df.show) %in% rownames(pamp.down),
                                           'Immunogenic down', 'Not downregulated')))
df.show$Gene <- rownames(df.show)
genes.to.annotate <- c('CD1C','NLRP3','TLR2','ITGAX','LYZ','IL1B','IFI16',
                       'IFITM2',
                       'LPL','CD1B','CD5')
rasterize(ggplot(df.show, aes(x=`logFC DC2hm`, `logFC PAMP-cDC2`, color=regulation))+
            geom_jitter(size=0.1)+
            theme_bw()+
            scale_color_manual(values = c('#FFBE0B','#F15BB5','#00BBF9','grey'))+
            ylim(-8,8)+
            xlim(-8,8)+
            ggrepel::geom_text_repel(
              data = subset(df.show, Gene %in% genes.to.annotate), # Filter to annotate specific genes
              aes(label = Gene), 
              size = 3, 
              box.padding = 0.5, 
              point.padding = 0.1, 
              segment.size = 0.3,
              max.overlaps = Inf, 
              force = 3,
              nudge_x = .2,
              nudge_y = .2,
              direction = 'both'
            ),
          dpi = 300)
ggsave('Outputs/figures/DC2hm_vs_PAMP_2-D_downregulation.pdf',
       width = 11, height = 8)

# Heatmap comparing DC2hm and PAMP-cDC2 up genes -------------------------------------------------------
  # Align the genes
genes <- intersect(rownames(dc2hm.all), c(rownames(pamp.all),
                                          rownames(ifng.lps.all)))
dc2hm.all <- dc2hm.all[genes,]
dc2hm.all <- dc2hm.all[order(dc2hm.all$adj.P.Val),]
dc2hm.up <- dc2hm.up[intersect(rownames(dc2hm.up),genes),]
dc2hm.up <- dc2hm.up[order(dc2hm.up$adj.P.Val),]
dc2hm.down <- dc2hm.down[intersect(rownames(dc2hm.down),genes),]
dc2hm.down <- dc2hm.down[order(dc2hm.down$adj.P.Val),]
pamp.all <- pamp.all[genes,]
pamp.all <- pamp.all[order(pamp.all$adj.P.Val),]
pamp.up <- pamp.up[intersect(rownames(pamp.up),genes),]
pamp.up <- pamp.up[order(pamp.up$adj.P.Val),]
pamp.down <- pamp.down[intersect(rownames(pamp.down),genes),]
pamp.down <- pamp.down[order(pamp.down$adj.P.Val),]
ifng.lps.all <- ifng.lps.all[genes,]
ifng.lps.all <- ifng.lps.all[order(ifng.lps.all$adj.P.Val),]
ifng.lps.up <- ifng.lps.up[intersect(rownames(ifng.lps.up),genes),]
ifng.lps.up <- ifng.lps.up[order(ifng.lps.up$adj.P.Val),]
ifng.lps.down <- ifng.lps.down[intersect(rownames(ifng.lps.down),genes),]
ifng.lps.down <- ifng.lps.down[order(ifng.lps.down$adj.P.Val),]
  # Create the dataframe
# Heatmap up
both.up <- intersect(rownames(dc2hm.up),
                     rownames(pamp.up))
homeostatic.up <- rownames(dc2hm.up[!rownames(dc2hm.up)%in% intersect(rownames(dc2hm.up),
                                                                      rownames(pamp.up)),])
immunogenic.up <- rownames(pamp.up[!rownames(pamp.up)%in% intersect(rownames(dc2hm.up),
                                                                    rownames(pamp.up)),])
log2fc.combined <- cbind(dc2hm.all[c(both.up,homeostatic.up, immunogenic.up),
                                   1:2],
                         pamp.all[c(both.up,homeostatic.up, immunogenic.up),1:2])
log2fc.combined <- log2fc.combined[,c(1,3)]
colnames(log2fc.combined) <- c('LFC in DC2hm',
                               'LFC in PAMP-cDC2')
# Enrichment
homeostatic.enrichment <- enrichGO(homeostatic.up, keyType = 'SYMBOL',
                                   ont = 'ALL', OrgDb = 'org.Hs.eg.db')
homeostatic.enrichment.show <- homeostatic.enrichment@result[homeostatic.enrichment@result$Description %in% c('axon development',
                                                                                                              'cell-cell junction organization', 
                                                                                                              'lipid storage',
                                                                                                              'cardiac muscle contraction',
                                                                                                              'neuron projection extension',
                                                                                                              'plasma membrane raft',
                                                                                                              'actin cytoskeleton'),]
homeostatic.enrichment.show$Description <- factor(homeostatic.enrichment.show$Description,
                                                  levels = homeostatic.enrichment.show[order(homeostatic.enrichment.show$Count),'Description'])
ggplot(homeostatic.enrichment.show, aes(x=Count,y=Description,fill=-log10(p.adjust)))+
  geom_bar(stat = 'identity')+
  scale_fill_gradient(high = '#fff7b7', low ='#e15383')+
  theme_classic()
ggsave('Outputs/figures/Homeostatic_up_enrichment.pdf',
       width = 5, height = 4)
immunogenic.enrichment <- enrichGO(immunogenic.up, keyType = 'SYMBOL',
                                   ont = 'ALL', OrgDb = 'org.Hs.eg.db')
immunogenic.enrichment.show <- immunogenic.enrichment@result[immunogenic.enrichment@result$Description %in% c('response to virus',
                                                                                                              'regulation of response to biotic stimulus',
                                                                                                              'regulation of innate immune response',
                                                                                                              'regulation of viral life cycle',
                                                                                                              'type II interferon production',
                                                                                                              'interleukin-1 production',
                                                                                                              'cytokine production involved in immune response',
                                                                                                              'interleukin-1 receptor binding'),]
immunogenic.enrichment.show$Description <- factor(immunogenic.enrichment.show$Description,
                                                  levels = immunogenic.enrichment.show[order(immunogenic.enrichment.show$Count),'Description'])
ggplot(immunogenic.enrichment.show, aes(x=Count,y=Description,fill=-log10(p.adjust)))+
  geom_bar(stat = 'identity')+
  scale_fill_gradient(high = '#fff7b7', low ='#e15383')+
  theme_classic()
ggsave('Outputs/figures/immunogenic_up_enrichment.pdf',
       width = 5, height = 5)
# Common up
common.enrichment <- enrichGO(both.up, keyType = 'SYMBOL',
                              ont = 'ALL', OrgDb = 'org.Hs.eg.db')
common.enrichment.show <- common.enrichment@result[common.enrichment@result$Description %in% c("positive regulation of lymphocyte activation",                                                                                                          
                                                                                               "positive regulation of leukocyte cell-cell adhesion"),]
common.enrichment.show$Description <- factor(common.enrichment.show$Description,
                                             levels = c("positive regulation of leukocyte cell-cell adhesion",
                                                        "positive regulation of lymphocyte activation"))
ggplot(common.enrichment.show, aes(x=Count,y=Description,fill=-log10(p.adjust)))+
  geom_bar(stat = 'identity')+
  scale_fill_gradient(high = '#fff7b7', low ='#e15383')+
  theme_classic()
ggsave('Outputs/figures/common_up_enrichment.pdf',
       width = 6, height = 1.7)

write.xlsx(list('Common up'=common.enrichment@result,
                'Homeostatic up'=homeostatic.enrichment@result,
                'Immunogenic up'=immunogenic.enrichment@result),
           'Outputs/tables/Upregulated_gene_module_enrichment.xlsx',rowNames=T)

  # Plot heatmap
up.show <- c('HIVEP1', 'HIVEP2', 'CD80', 'CD40', 'RELB', 'CD274',
             'CCR7','CCL22', 'ICOSLG', 'IL7R', 'CRLF2', 'REL', 'SREBF2',
             'MGLL', 'TNFSF4', 'LAMP3', 'OAS1', 'STAT1', 'ISG20', 'ISG15', 'IFIT2', 'MX1',
             'IFIT3', 'MX2', 'RELA', 'OASL', 'IL12B', 'IF116','BIRC3','LEF1','SQLE','LPL','HMGCS1',
             'TNNT2', 'AUTS2', 'SOX8')
group.annotation <- data.frame('Annotation'=c(rep('Both up', length(intersect(rownames(dc2hm.up),
                                                                              rownames(pamp.up)))),
                                              rep('Up only in DC2hm',
                                                  length(rownames(dc2hm.up[!rownames(dc2hm.up)%in% intersect(rownames(dc2hm.up),
                                                                                                             rownames(pamp.up)),]))),
                                              rep('Up only in PAMP sti',
                                                  length(rownames(pamp.up[!rownames(pamp.up)%in% intersect(rownames(dc2hm.up),
                                                                                                           rownames(pamp.up)),])))),
                               row.names = rownames(log2fc.combined))
index <- match(up.show,rownames(log2fc.combined))
anno.markers <- rowAnnotation(ano = anno_mark(at = index, 
                                              labels = up.show,
                                              labels_gp = gpar(fontsize = 6),
                                              lines_gp = gpar(),
                                              extend = unit(10, 'npc')))
pdf('Outputs/figures/heatmap_pamp_up.pdf', width = 4, height = 10)
Heatmap(log2fc.combined, 
        cluster_columns = F, cluster_rows = F,
        show_row_names = F, right_annotation = anno.markers,
        row_split = group.annotation$Annotation,gap = unit(0.01, 'npc'),
        width = unit(0.2, 'npc'),
        height = unit(0.45, 'npc'), border = T, row_title_rot = 0,
        col = colorRamp2(breaks = c(-5,0,10), 
                         colors = c('#2A9FF3','white' ,'#FF2B72')),
        raster_quality = 300, raster_device = 'tiff')
dev.off()



# GSEA --------------------------------------------------------------------
dc2hm.all <- dc2hm.all[order(dc2hm.all$logFC,
                             decreasing = T),]
dc2hm.all.ranked <- dc2hm.all
dc2hm.all.ranked <- dc2hm.all$logFC
names(dc2hm.all.ranked) <- rownames(dc2hm.all)
  # Figure 5D, DC2hm vs cDC2 TSLP/STAT5
dc2hm.vs.cdc2.tslp.stat5.gsea <- GSEA(dc2hm.all.ranked,
                                      TERM2GENE = rbind(read.gmt('WP_THYMIC_STROMAL_LYMPHOPOIETIN_TSLP_SIGNALING_PATHWAY.v2023.2.Hs.gmt'),
                                                        read.gmt('WIERENGA_STAT5A_TARGETS_UP.v2023.2.Hs.gmt')),
                                      seed = 0712)
gseaplot2(dc2hm.vs.cdc2.tslp.stat5.gsea, geneSetID = c(2,1),
          color = c('#e065b2','#4e7796'))
write.csv(dc2hm.vs.cdc2.tslp.stat5.gsea@result,
          'Outputs/tables/DC2hm_vs_cDC2_STLP_STAT5_gsea.csv')
ggsave('Outputs/figures/GSEA_DC2hm_vs_cDC2_STAT5_TSLP.pdf',
       width = 4, height = 4)

  # Figure S5E, GSEA
    # DC2hm
dc2hm.gsea <- GSEA(dc2hm.all.ranked,
                   TERM2GENE = rbind(read.gmt('LINDSTEDT_DENDRITIC_CELL_MATURATION_B.v2023.2.Hs.gmt'),
                                     read.gmt('WIERENGA_STAT5A_TARGETS_UP.v2023.2.Hs.gmt'),
                                     read.gmt('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING.v2023.2.Hs.gmt')),
                   seed = 0712, pvalueCutoff = 1)
gseaplot2(dc2hm.gsea, geneSetID = c(1,2,3),
          color = c('#FFBE0B','#00BBF9','#F15BB5'))
ggsave('Outputs/figures/GSEA_DC2hm_vs_cDC2_3pathway.pdf',
       width = 4, height = 4)
    # PAMP
pamp.all <- pamp.all[order(pamp.all$logFC,
                             decreasing = T),]
pamp.all.ranked <- pamp.all
pamp.all.ranked <- pamp.all$logFC
names(pamp.all.ranked) <- rownames(pamp.all)
pamp.vs.cdc2.gsea <- GSEA(pamp.all.ranked,
                          TERM2GENE = rbind(read.gmt('LINDSTEDT_DENDRITIC_CELL_MATURATION_B.v2023.2.Hs.gmt'),
                                            read.gmt('WIERENGA_STAT5A_TARGETS_UP.v2023.2.Hs.gmt'),
                                            read.gmt('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING.v2023.2.Hs.gmt')),
                          seed = 0712)
gseaplot2(pamp.vs.cdc2.gsea, geneSetID = c(1,2,3),
          color = c('#FFBE0B','#00BBF9','#F15BB5'))

ggsave('Outputs/figures/GSEA_PAMP_vs_cDC2_3pathway.pdf',
       width = 4, height = 4)
  # LPS+IFNg
ifng.lps.all <- ifng.lps.all[order(ifng.lps.all$logFC,
                                   decreasing = T),]
ifng.lps.all.ranked <- ifng.lps.all
ifng.lps.all.ranked <- ifng.lps.all$logFC
names(ifng.lps.all.ranked) <- rownames(ifng.lps.all)
ifng.lps.all.vs.cdc2.gsea <- GSEA(ifng.lps.all.ranked,
                                  TERM2GENE = rbind(read.gmt('LINDSTEDT_DENDRITIC_CELL_MATURATION_B.v2023.2.Hs.gmt'),
                                                    read.gmt('WIERENGA_STAT5A_TARGETS_UP.v2023.2.Hs.gmt'),
                                                    read.gmt('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING.v2023.2.Hs.gmt')),
                                  seed = 0712)
gseaplot2(ifng.lps.all.vs.cdc2.gsea, geneSetID = c(1,2,3),
          color = c('#FFBE0B','#00BBF9','#F15BB5'))
ggsave('Outputs/figures/GSEA_LPS+IFNG_vs_cDC2_3pathway.pdf',
       width = 4, height = 4)
write.xlsx(list('DC2hm vs cDC2s TSLP_STAT5'=dc2hm.vs.cdc2.tslp.stat5.gsea@result,
                'DC2hm vs cDC2s 3 pathways'=dc2hm.gsea@result,
                'PAMP-cDC2s 3 pathways'=pamp.vs.cdc2.gsea@result,
                'LPS+IFNG-cDC2s pathways'=ifng.lps.all.vs.cdc2.gsea@result),
           'Outputs/tables/GSEA_analysis.xlsx')
