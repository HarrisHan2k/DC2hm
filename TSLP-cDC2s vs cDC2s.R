library(pacman)
p_load(sva, limma, pheatmap, ggplot2, ggplotify, dplyr, openxlsx,
       clusterProfiler,tidyverse, ggrepel, Seurat, ggsci, viridis,
       ggrastr, enrichplot, ggvenn,affy,gcrma, pheatmap,GEOquery, ggpubr,
       Biobase) 
set.seed(0822)
configure.args <- c(preprocessCore = "--disable-threading") # The multi-core GCRMA analysis seems to be troublesome
options(configure.args = configure.args)
setwd('TSLP-cDC2s vs cDC2s')
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
de.result.tslp.vs.fresh <- de.result.tslp.vs.fresh[order(de.result.tslp.vs.fresh$logFC,
                                                         decreasing = T),]
de.result.tslp.vs.fresh.ranked <-de.result.tslp.vs.fresh$logFC
names(de.result.tslp.vs.fresh.ranked) <- rownames(de.result.tslp.vs.fresh)
tslp.vs.fresh <- GSEA(de.result.tslp.vs.fresh.ranked,
                      TERM2GENE = rbind(read.gmt('LINDSTEDT_DENDRITIC_CELL_MATURATION_B.v2023.2.Hs.gmt'),
                                        read.gmt('WIERENGA_STAT5A_TARGETS_UP.v2023.2.Hs.gmt'),
                                        read.gmt('REACTOME_INTERFERON_ALPHA_BETA_SIGNALING.v2023.2.Hs.gmt')),
                      seed = 0712)
write.csv(tslp.vs.fresh@result,'Outputs/tables/TSLP_vs_Fresh_GSEA.csv')
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
                'Upregulated'=de.result.stat5.inhibitor.tslp[de.result.stat5.inhibitor.tslp$logFC > 1 & de.result.stat5.inhibitor.tslp$adj.P.Val < 0.05,],
                'Downregulated'=de.result.stat5.inhibitor.tslp[de.result.stat5.inhibitor.tslp$logFC < -(1) & de.result.stat5.inhibitor.tslp$adj.P.Val < 0.05,]),
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





# Process blood DC seq data (GSE59237) ---------------------------------------------

library(GEOquery)
library(limma)

getGEOSuppFiles("GSE59237")

# Unpack the .tar archive. This creates a directory with the .CEL.gz files
untar("GSE59237/GSE59237_RAW.tar", exdir = "GSE59237/CEL_files")

# List the files to confirm they are there
cel_files <- list.files("GSE59237/CEL_files", pattern = "\\.CEL\\.gz$", full.names = TRUE)
print(cel_files)

# Read the raw .CEL files into an AffyBatch object
raw_data <- ReadAffy(filenames = cel_files)
# Here GC-RMA normalization was performed according to the original paper
cat("Performing GC-RMA normalization (this may take a few minutes)...\n")
eset <- gcrma(raw_data)
cat("Normalization with GC-RMA complete. ✅\n")

# Step 4: Annotation, Filtering, and Scaling
# -----------------------------------------------------------------
cat("Filtering and scaling data...\n")
# Get feature annotation data from GEO
gse <- getGEO("GSE59237", GSEMatrix = TRUE)
feature_data_geo <- Biobase::fData(gse[[1]])

# Extract the normalized expression matrix (already log2-transformed by gcrma)
expression_matrix <- Biobase::exprs(eset)

# 4a. "Probes with no annotation were removed"
probes_to_keep <- rownames(feature_data_geo[feature_data_geo$`Gene Symbol` != "", ])
expression_matrix <- expression_matrix[intersect(rownames(expression_matrix), probes_to_keep), ]

# 4b. "Genes with small profile ranges...were filtered out"
gene_ranges <- apply(expression_matrix, 1, function(x) max(x) - min(x))
range_cutoff <- median(gene_ranges)
expression_matrix <- expression_matrix[gene_ranges > range_cutoff, ]

# 4c. "expression levels were centered and reduced" (Z-score scaling)
expression_matrix <- t(scale(t(expression_matrix)))

# Step 5: Differential Expression Analysis (TSLP vs. Fresh)
# -----------------------------------------------------------------
cat("Running differential expression analysis...\n")
# Define groups using the sample names from the normalized data
sample_names <- sampleNames(eset)
groups <- character(length(sample_names))
groups[grepl("Control_0h", sample_names)] <- "Control_0h"
groups[grepl("Control_6h", sample_names)] <- "Control_6h"
groups[grepl("TSLP_6h", sample_names)]    <- "TSLP_6h"
groups[grepl("TNF_6h", sample_names)]     <- "TNF_6h"
groups <- factor(groups, levels = c("Control_0h", "Control_6h", "TSLP_6h", "TNF_6h"))

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(expression_matrix, design)
cont.matrix <- makeContrasts(TSLP_vs_Fresh = TSLP_6h - Control_0h, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Step 6 & 7: Create and Clean the Final DE Table
# -----------------------------------------------------------------
cat("Generating final results...\n")
full_results <- topTable(fit2, coef = "TSLP_vs_Fresh", number = Inf, sort.by = "P")
full_results$Gene_Symbol <- feature_data_geo[rownames(full_results), "Gene Symbol"]

final_de_table <- full_results %>%
  filter(!is.na(Gene_Symbol) & Gene_Symbol != "") %>%
  group_by(Gene_Symbol) %>%
  slice_min(order_by = adj.P.Val, n = 1) %>%
  ungroup() %>%
  as.data.frame()

rownames(final_de_table) <- make.unique(as.character(final_de_table$Gene_Symbol))
final_de_table <- final_de_table %>% dplyr::select(-Gene_Symbol)

cat("\n## Preview of TSLP vs. Fresh DE results table (GC-RMA normalized):\n")
print(head(final_de_table))


sig_threshold_adj_p <- 0.05
microarray_sig_genes <- final_de_table[final_de_table$adj.P.Val < sig_threshold_adj_p, ]

microarray_up <- gsub("\\.\\d+$", 
                      "", 
                      rownames(microarray_sig_genes[microarray_sig_genes$logFC > 1, ]))
microarray_down <- gsub("\\.\\d+$",
                        "", 
                        rownames(microarray_sig_genes[microarray_sig_genes$logFC < (-1), ]))

rnaseq_up <- rownames(de.result.tslp.vs.fresh[de.result.tslp.vs.fresh$logFC > 1 & de.result.tslp.vs.fresh$adj.P.Val < 0.05, ])
rnaseq_down <- rownames(de.result.tslp.vs.fresh[de.result.tslp.vs.fresh$logFC < -1 & de.result.tslp.vs.fresh$adj.P.Val < 0.05, ])

venn_list_up <- list('Blood Microarray' = microarray_up, 'Spleen RNA-seq' = rnaseq_up)
ggvenn(venn_list_up, fill_color = c("#F8766D", "#00BFC4"), stroke_size = 0.5, set_name_size = 4) +
  labs(title = "Upregulated Genes in TSLP vs. Fresh")
ggsave("Outputs/figures/Venn_Upregulated_TSLP_vs_Fresh_GCRMA.pdf", width = 5, height = 4)

venn_list_down <- list('Blood Microarray' = microarray_down, 'Spleen RNA-seq' = rnaseq_down)
ggvenn(venn_list_down, fill_color = c("#F8766D", "#00BFC4"), stroke_size = 0.5, set_name_size = 4) +
  labs(title = "Downregulated Genes in TSLP vs. Fresh")
ggsave("Outputs/figures/Venn_Downregulated_TSLP_vs_Fresh_GCRMA.pdf", width = 5, height = 4)


# Your RNA-seq TSLP vs. Fresh results
rnaseq_tslp_fresh_df <- de.result.tslp.vs.fresh
rnaseq_tslp_fresh_df$Gene_Symbol <- rownames(rnaseq_tslp_fresh_df)
rnaseq_tslp_fresh_df <- rnaseq_tslp_fresh_df[, c("Gene_Symbol",
                                                 "logFC", "AveExpr",
                                                 "t", "P.Value", 
                                                 "adj.P.Val", "B")]

# Your RNA-seq TSLP vs. Medium results
rnaseq_tslp_med_df <- de.result.tslp.vs.med
rnaseq_tslp_med_df$Gene_Symbol <- rownames(rnaseq_tslp_med_df)
rnaseq_tslp_med_df <- rnaseq_tslp_med_df[, c("Gene_Symbol", 
                                             "logFC", "AveExpr", 
                                             "t", "P.Value", "adj.P.Val",
                                             "B")]

# The final microarray TSLP vs. Fresh results
microarray_final_df <- final_de_table
microarray_final_df$Gene_Symbol <- gsub("\\.\\d+$", "", rownames(microarray_final_df)) # Clean the unique names
microarray_final_df <- microarray_final_df[, c("Gene_Symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]


# Step 2: Get the lists of overlapping genes from the Venn diagram analysis
# -------------------------------------------------------------------------
# This entire section for creating the Venn diagram lists is UNCHANGED, as it is correct.
sig_threshold_adj_p <- 0.05
microarray_sig_genes <- final_de_table[final_de_table$adj.P.Val < sig_threshold_adj_p, ]
microarray_up <- gsub("\\.\\d+$", "", rownames(microarray_sig_genes[microarray_sig_genes$logFC > 1, ]))
microarray_down <- gsub("\\.\\d+$", "", rownames(microarray_sig_genes[microarray_sig_genes$logFC < (-1), ]))

rnaseq_up <- rownames(de.result.tslp.vs.fresh[de.result.tslp.vs.fresh$logFC > 1 & de.result.tslp.vs.fresh$adj.P.Val < 0.05, ])
rnaseq_down <- rownames(de.result.tslp.vs.fresh[de.result.tslp.vs.fresh$logFC < -1 & de.result.tslp.vs.fresh$adj.P.Val < 0.05, ])

# Find the intersection
overlap_upregulated_genes <- intersect(microarray_up, rnaseq_up)
overlap_downregulated_genes <- intersect(microarray_down, rnaseq_down)


# ############################################################################### #
# ### START OF CORRECTED SECTION TO BUILD EXCEL SHEETS ###
# ############################################################################### #

# Prepare the RNA-seq results table (this was already correct)
rnaseq_stats_for_merge <- de.result.tslp.vs.fresh %>%
  select(logFC, adj.P.Val) %>%
  rownames_to_column("Gene_Symbol")

# Prepare the microarray results table, ENSURING it is de-duplicated to prevent join errors
microarray_stats_for_merge <- final_de_table %>%
  rownames_to_column("Gene_Symbol_unique") %>%
  mutate(Gene_Symbol = gsub("\\.\\d+$", "", Gene_Symbol_unique)) %>%
  # THIS IS THE FIX: group by the clean gene symbol and take only the first entry in case of ties
  group_by(Gene_Symbol) %>%
  slice(1) %>%
  ungroup() %>%
  select(Gene_Symbol, logFC, adj.P.Val)

# Create the detailed dataframe for UPREGULATED overlap. This will now have the correct number of rows.
overlap_up_df <- data.frame(Gene_Symbol = overlap_upregulated_genes) %>%
  left_join(rnaseq_stats_for_merge, by = "Gene_Symbol") %>%
  rename(logFC_RNAseq = logFC, adj.P.Val_RNAseq = adj.P.Val) %>%
  left_join(microarray_stats_for_merge, by = "Gene_Symbol") %>%
  rename(logFC_Microarray = logFC, adj.P.Val_Microarray = adj.P.Val)

# Create the detailed dataframe for DOWNREGULATED overlap. This will also have the correct number of rows.
overlap_down_df <- data.frame(Gene_Symbol = overlap_downregulated_genes) %>%
  left_join(rnaseq_stats_for_merge, by = "Gene_Symbol") %>%
  rename(logFC_RNAseq = logFC, adj.P.Val_RNAseq = adj.P.Val) %>%
  left_join(microarray_stats_for_merge, by = "Gene_Symbol") %>%
  rename(logFC_Microarray = logFC, adj.P.Val_Microarray = adj.P.Val)



# --- NEW SECTION: Find unique genes and create dataframes for them ---
cat("Identifying unique genes for each dataset...\n")

# 1. Get the lists of unique genes using setdiff()
unique_up_rnaseq    <- setdiff(rnaseq_up, microarray_up)
unique_down_rnaseq  <- setdiff(rnaseq_down, microarray_down)
unique_up_microarray   <- setdiff(microarray_up, rnaseq_up)
unique_down_microarray <- setdiff(microarray_down, rnaseq_down)

# 2. Create dataframes with stats for the unique genes
unique_up_rnaseq_df <- data.frame(Gene_Symbol = unique_up_rnaseq) %>%
  left_join(rnaseq_stats_for_merge, by = "Gene_Symbol")

unique_down_rnaseq_df <- data.frame(Gene_Symbol = unique_down_rnaseq) %>%
  left_join(rnaseq_stats_for_merge, by = "Gene_Symbol")

unique_up_microarray_df <- data.frame(Gene_Symbol = unique_up_microarray) %>%
  left_join(microarray_stats_for_merge, by = "Gene_Symbol")

unique_down_microarray_df <- data.frame(Gene_Symbol = unique_down_microarray) %>%
  left_join(microarray_stats_for_merge, by = "Gene_Symbol")


# --- UPDATED: Create a named list of ALL data frames for export ---
list_of_datasets <- list(
  "RNAseq_TSLP_vs_Fresh" = rnaseq_tslp_fresh_df,
  "Microarray_TSLP_vs_Fresh" = microarray_final_df,
  "Overlap_Upregulated" = overlap_up_df,
  "Overlap_Downregulated" = overlap_down_df,
  "Unique_Up_RNAseq" = unique_up_rnaseq_df,
  "Unique_Down_RNAseq" = unique_down_rnaseq_df,
  "Unique_Up_Microarray" = unique_up_microarray_df,
  "Unique_Down_Microarray" = unique_down_microarray_df
)

# Write the list to a single .xlsx file
library(openxlsx)
write.xlsx(list_of_datasets, 
           file = "Outputs/tables/Comprehensive_Analysis_Summary.xlsx", 
           rowNames = FALSE)

# ############################################################################### #
cat("\n--- Starting Gene Ontology Analysis for Venn Segments ---\n")

# Step 1: Define the gene universe
universe_rnaseq <- rownames(de.result.tslp.vs.fresh)
universe_microarray <- gsub("\\.\\d+$", "", rownames(final_de_table))
common_universe <- intersect(universe_rnaseq, universe_microarray)
cat(paste("Using a common universe of", length(common_universe), "genes for GO analysis.\n"))

# Step 2: Loop through datasets and perform GO analysis
go_results_list <- list()

for (set_name in names(list_of_datasets)[3:8]) {
  
  current_df <- list_of_datasets[[set_name]]
  gene_list <- current_df$Gene_Symbol
  
  cat(paste("\nRunning GO analysis for:", set_name, " (", length(gene_list), " genes)\n"))
  
  if (length(gene_list) < 10) {
    cat("Skipping: Too few genes for GO analysis.\n")
    next
  }
  
  # Perform GO enrichment analysis
  ego <- enrichGO(gene          = gene_list,
                  universe      = common_universe,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.1)
  
  if (!is.null(ego) && nrow(ego@result) > 0) {
    cat(paste("Found", nrow(ego@result), "significant GO terms.\n"))
    
    # Store the full results
    go_results_list[[set_name]] <- as.data.frame(ego@result)
    
    # --- FINAL FIX APPLIED HERE ---
    # We will now create the dot plot manually with ggplot2 to avoid errors.
    
    # First, try to simplify if there are many terms.
    plot_object <- ego # Default to the original object
    title_suffix <- ""
    if (nrow(ego@result) > 30) {
      cat("Many terms found. Attempting to simplify for clarity...\n")
      # Use try() to catch potential errors during simplification
      ego_simplified_attempt <- try(
        clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min),
        silent = TRUE
      )
      
      # Check if simplification was successful and left any terms
      if (!inherits(ego_simplified_attempt, "try-error") && nrow(as.data.frame(ego_simplified_attempt)) > 0) {
        plot_object <- ego_simplified_attempt
        title_suffix <- " (Simplified)"
        cat("Simplification successful.\n")
      } else {
        cat("Simplification removed all terms or failed. Using original results.\n")
      }
    }
    
    # Manually create the dot plot using ggplot2 from the results dataframe
    # This is much more stable than the built-in dotplot function.
    plot_df <- as.data.frame(plot_object) %>%
      top_n(20, wt = -p.adjust)
    
    # Reorder factor levels for plotting
    plot_df$Description <- factor(plot_df$Description, levels = rev(plot_df$Description))
    
    p <- ggplot(plot_df, aes(x = GeneRatio, y = Description)) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_gradient(low = "red", high = "blue", name = "Adj. P-value") +
      theme_bw(base_size = 12) +
      ylab(NULL) +
      xlab("Gene Ratio") +
      ggtitle(paste("GO Enrichment for", gsub("_", " ", set_name), title_suffix))
    
    plot_filename <- paste0("Outputs/figures/GO_DotPlot_", set_name, ".pdf")
    ggsave(filename = plot_filename, plot = p, width = 10, height = 8) # Adjusted size for better readability
    cat(paste("Dot plot saved to:", plot_filename, "\n"))
    
  } else {
    cat("No significant GO terms found with the current cutoffs.\n")
  }
}

# Step 3: Save results to a single Excel file
if (length(go_results_list) > 0) {
  excel_filename <- "Outputs/tables/GO_Analysis_Venn_Diagram_Parts.xlsx"
  write.xlsx(go_results_list, file = excel_filename, rowNames = FALSE)
  cat(paste("\nAll GO analysis results have been compiled into:", excel_filename, "📄\n"))
}

cat("--- Gene Ontology Analysis Complete ---\n\n")
# LC markers --------------------------------------------------------------
de.result.tslp.vs.fresh.show <- de.result.tslp.vs.fresh[lc.markers,]
de.result.tslp.vs.fresh.show$Gene <- factor(rownames(de.result.tslp.vs.fresh.show),
                                            levels = rownames(de.result.tslp.vs.fresh.show))
p_rna_seq <- ggplot(de.result.tslp.vs.fresh.show,
       aes(y=logFC,x=Gene, fill=-log10(adj.P.Val)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  scale_fill_gradient(low = '#F686BD', high = '#d930a6')+
  labs(y='Log2FC: TSLP-cDC2s vs Fresh cDC2s', x=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
  geom_hline(yintercept = c(0), color='black')+
  geom_hline(yintercept = c(1), color='lightgrey')
p_rna_seq
final_de_table.show <- final_de_table[lc.markers,]
final_de_table.show$Gene <- factor(rownames(de.result.tslp.vs.fresh.show),
                                            levels = rownames(de.result.tslp.vs.fresh.show))
p_microarray <- ggplot(final_de_table.show,
       aes(y=logFC,x=Gene, fill=-log10(adj.P.Val)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  scale_fill_gradient(low = '#F686BD', high = '#d930a6')+
  labs(y='Log2FC: TSLP-cDC2s vs Fresh cDC2s', x=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
  geom_hline(yintercept = c(0), color='black')+
  geom_hline(yintercept = c(1), color='lightgrey')

  
ggarrange(p_rna_seq,p_microarray)





# tslp vs. med ------------------------------------------------------------
de.result.tslp.vs.med.show <- de.result.tslp.vs.med[lc.markers,]
de.result.tslp.vs.med.show$Gene <- factor(rownames(de.result.tslp.vs.med.show),
                                            levels = rownames(de.result.tslp.vs.med.show))
p_rna_seq <- ggplot(de.result.tslp.vs.med.show,
                    aes(y=logFC,x=Gene, fill=-log10(adj.P.Val)))+
  geom_bar(stat = 'identity')+
  theme_classic()+
  scale_fill_gradient(low = '#F686BD', high = '#d930a6')+
  labs(y='Log2FC: TSLP-cDC2s vs Med cDC2s', x=NULL)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5))+
  geom_hline(yintercept = c(0), color='black')+
  geom_hline(yintercept = c(1), color='lightgrey')
p_rna_seq



# Heatmap -----------------------------------------------------------------
gene.len <- read.csv('../TRA_vs_PAMP/Gene_length.csv')
# Select the 'Fresh cDC2s' columns from your data
fresh_counts <- counts.tslp.vs.fresh[, 1:3]
colnames(fresh_counts) <- c("Fresh_1", "Fresh_2", "Fresh_3")

# Select the 'TSLP cDC2s' columns
tslp_counts <- counts.tslp.vs.fresh[, 4:6]
colnames(tslp_counts) <- c("TSLP_1", "TSLP_2", "TSLP_3")

# Select the 'Medium cDC2s' columns
med_counts <- counts.med.vs.fresh[, 1:3] 
colnames(med_counts) <- c("Medium_1", "Medium_2", "Medium_3")

# Find the common genes shared across all datasets
common_genes <- intersect(rownames(fresh_counts), intersect(rownames(tslp_counts), rownames(med_counts)))

# Combine them into a single count matrix
combined_counts <- cbind(fresh_counts[common_genes, ], 
                         med_counts[common_genes, ], 
                         tslp_counts[common_genes, ])


# --- Step 2: Process Your Gene Length Data ---

# Rename the columns of your 'gene.len' object for clarity
colnames(gene.len) <- c("gene_name", "length")

# Set the row names to be the gene names for easy lookup
rownames(gene.len) <- gene.len$gene_name


# --- Step 3: Calculate TPM from Raw Counts ---

# Ensure the genes in your counts and length data match
common_genes_with_length <- intersect(rownames(combined_counts), rownames(gene.len))
counts_final <- combined_counts[common_genes_with_length, ]
lengths_final <- gene.len[common_genes_with_length, "length"]

# 1. Divide read counts by gene length in kilobases (RPK)
rpk <- counts_final / (lengths_final / 1000)

# 2. Calculate the "per million" scaling factor for each sample
scaling_factor <- colSums(rpk) / 1e6

# 3. Divide RPK values by the scaling factor to get TPM
tpm_matrix <- t(t(rpk) / scaling_factor)

cat("TPM matrix calculated successfully.\n")


# --- Step 4: Generate the Heatmap ---


# Subset the TPM matrix to only these genes
# We use intersect to avoid errors if a gene is not found
lc.markers <- c("CD1A", "CD207",'CCR6','CCL17','MMP12','CD209','CD14')
genes_found <- intersect(lc.markers, rownames(tpm_matrix))
tpm_subset <- tpm_matrix[genes_found, ]

# Scale the data by row (Z-score). This highlights relative changes
# for each gene across the conditions, which is best for visualization.
scaled_tpm <- t(scale(t(tpm_subset)))

# Create an annotation data frame for the columns
annotation_col <- data.frame(
  Condition = factor(rep(c("Fresh", "Medium", "TSLP"), each = 3))
)
rownames(annotation_col) <- colnames(scaled_tpm)
# Removing outliers
scaled_tpm <- scaled_tpm[,!colnames(scaled_tpm) %in% c('Fresh_3',
                                                      'Medium_2',
                                                      'TSLP_3')]
# Generate and save the heatmap
pheatmap(scaled_tpm,
         annotation_col = annotation_col,
         cluster_rows = F,          # Cluster genes based on similarity
         cluster_cols = FALSE,         # Keep samples in the Fresh-Med-TSLP order
         show_colnames = TRUE,
         main = "LC markers",
         fontsize_row = 10,         # Use the viridis color palette
         filename = "Outputs/figures/TPM_Heatmap_LC_Genes.pdf",
         width = 6, height = 5)

cat("Heatmap saved to Outputs/figures/TPM_Heatmap_TSLP_Genes.pdf\n")
