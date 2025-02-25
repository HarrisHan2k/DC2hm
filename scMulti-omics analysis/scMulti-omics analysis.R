library(Seurat)
library(viridis)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(Signac)
library(ComplexHeatmap)
library(circlize)
library(ggplotify)
library(clusterProfiler)
library(openxlsx)
library(ggpubr)
library(future)
library(JASPAR2020)
library(TFBSTools)
library(ggrastr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)
library(ggsci)
set.seed(0822)
setwd('/media/dell/0E54E2B554E29EA9/HanRunpeng/mregDC_project/Lab_Sequencing_Data/Code deposit/scMulti-omics analysis')
# The preprocessing and clustering analysis was performed by Qiuchen Zhao
# Start with this processed object which could be obtained @GEO: GSE267255

# Clustering and marker analysis on cDCs ----------------------------------
multiomic.object <- readRDS('hu_spleen_multiomics.rds')
multiomic.object <- RenameIdents(multiomic.object, c('mDC1'='cDC1',
                                                     'mDC2'='cDC2',
                                                     'mDC2_stressed'='cDC2_UPR',
                                                     'mDC2_activated'='DC2pre-hm',
                                                     'DC3'='DC2hm',
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
  # Figure 1F, UMAP of scATAC peaks
multiomic.object.cdcs <- multiomic.object[,Idents(multiomic.object) %in% c('cDC1', 
                                                                      'cDC2',
                                                                      'DC2pre-hm',
                                                                      'DC2hm')]
levels(multiomic.object.cdcs) <- c('cDC1', 'cDC2','DC2pre-hm', 'DC2hm')
DefaultAssay(multiomic.object.cdcs) <- 'peaks'
multiomic.object.cdcs <- RunTFIDF(multiomic.object.cdcs) %>% FindTopFeatures(min.cutoff = 'q0') %>% RunSVD() %>% RunUMAP(reduction='lsi', dims=2:10)
rasterize(DimPlot(multiomic.object.cdcs, 
                  cols = c('#C5E524', '#38C2E5','#384C94','#E679C5'),
                  pt.size = 0.1),dpi = 300)
ggsave('Outputs/figures/UMAP_peaks_reduction.pdf', width = 4.5, height = 3)
  # Figure 1N, feature peaks for each cDC cluster
DefaultAssay(multiomic.object.cdcs) <- 'peaks'
markers <- FindAllMarkers(multiomic.object.cdcs,
                          assay = 'peaks',test.use = 'LR',
                          latent.vars = 'nCount_peaks') 
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
rasterize(as.ggplot(Heatmap(de.peaks.mat.scaled, 
                            cluster_columns = F, cluster_rows = F,
                            show_row_names = F, 
                            #row_split = group.annotation$Annotation,
                            width = unit(0.2, 'npc'),left_annotation = anno.markers,
                            height = unit(0.6, 'npc'), border = T, row_title_rot = 0,
                            col = colorRamp2(breaks = c(-0,0.5,1), 
                                             colors = c('#2A9FF3','white' ,'#FF2B72')))),
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

# Coverage plots for Figure 1M----------------------------------------------------------
  # The fragments file is available @ GEO: GSE267255
UpdatePath(multiomic.object.cdcs@assays$ATAC@fragments[[1]], '../../Donor_6_ATAC-seq/cellranger_arc/D6_ATAC/atac_fragments.tsv.gz')
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

