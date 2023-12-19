library(corrplot) # correlation viz
library(DESeq2) # RNA-seq analysis
library(heatmaply) # clustered corrplots with dendrograms
library(ComplexHeatmap) # heatmap viz
library(EnsDb.Hsapiens.v75) # data for MT genes
library(Homo.sapiens) # for MT genes
library(org.Hs.eg.db) # database for gene annotation
library(clusterProfiler) # gene annotation
library(AnnotationDbi) # gene annotation
library(EnhancedVolcano) # vulcano plot visualization
library(ggplot2) # data viz
library(purrr) # data manipulation
library(dplyr) # data manipulation

#-----------------READ AND CHECK DATA-----------------------------------
# open reads
expression.data <- read.table("/Users/suleymanov-ef/Desktop/Skoltech/Study/NGS/Project/readcounts.txt",
                              header = TRUE, row.names = 1,
                              sep = "\t", na.strings = "NA",
                              dec = ".")

expression.data <- expression.data %>%
  select(where(is.numeric)) %>%
  select(!Length)
  
expression.data <- expression.data[,c(1:8)]

colnames(expression.data) <- sub('trimmed_', '', colnames(expression.data))
colnames(expression.data) <- sub('control', 'C', colnames(expression.data))
colnames(expression.data) <- sub('treatment_', '', colnames(expression.data))

# correlation and correlation_matrix
correlation <- cor(expression.data)

# correlation plot
correlation_plot <- corrplot(correlation, method = 'pie', tl.cex=1)
correlation_plot

# dendrogram
correlation_plot_clustering <- heatmaply_cor(correlation, colors=magma)
correlation_plot_clustering

#-----------------FULL_DATA---------------------------------------------------------------------------------

# Differential expression analysis
expr.matrix <- as.matrix(expression.data)

# Provide description of our data to DESeq2
expr.design <- data.frame(row.names = colnames(expression.data),
                          condition = c(
                            "control", "control", 'control', 'control',
                            "myxothiazole", "myxothiazole", 'myxothiazole', 'myxothiazole'
                          ))


# Let's create a variable that summarizes all the data about our experiment
dds <- DESeqDataSetFromMatrix(countData = expr.matrix,
                              colData = expr.design,
                              design = ~ condition)

dds <- DESeq(dds) #Perform diff expr analysis
res <- results(dds) #save results of analysis

dds <- estimateSizeFactors(dds)
sizeFactors(dds) # Median (size factor) coefficients show
normalized.counts <- counts(dds, normalized = TRUE)

# Estimation of variance
plotDispEsts(dds)

# Results
plotMA(dds, ylim = c(-2,2), main = "DESeq2")

sum(res$padj < 0.05, na.rm = TRUE)
resSig <- subset(res, padj < 0.05)

# find genes with the biggest foldchange
the_biggest_fold_change <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]

# Calculate the number of up and down regulated genes
sum(resSig$log2FoldChange > 1, na.rm = TRUE)
sum(resSig$log2FoldChange < -1, na.rm = TRUE)

# Find upregulated and downregulated genes
resSigUp <- subset(resSig, log2FoldChange > 1)
resSigDown <- subset(resSig, log2FoldChange < -1)

#write.table(as.data.frame(resSigUp), "C:/Users/aisin/Desktop/Upreg_Leaf.txt",
#            sep = "\t", col.names = TRUE, row.names = TRUE,
#            na = "NA", quote = FALSE)

#writeLines(as.character(rownames(resSigUp)), "C:/Users/aisin/Desktop/Upreg_Mitochondria.txt")


#writeLines(as.character(rownames(resSigDown)), "C:/Users/aisin/Desktop/Downreg_Mitochondria.txt")


# visualize the results

names <- c(row.names(resSigUp), row.names(resSigDown)) #make list of DE genes
de.genes <- subset(normalized.counts, rownames(normalized.counts)
                   %in% names) # create table with norm_counts

# we can exclude 0 values, or add real low value to all our data
de.genes <- de.genes + 0.01
log.de.genes <- log10(de.genes) # log again

# Scale data
scaled.de.genes <- scale(log.de.genes)
scaled.de.genes <- t(scale(t(log.de.genes)))
Heatmap(scaled.de.genes, show_row_names = FALSE)

# Volcano plot
ens <- rownames(expression.data)

#symbols <- mapIds(org.Hs.eg.db, keys = ens,
#                  column = 'SYMBOL', keytype = 'ENSEMBL')

#symbols <- symbols[!is.na(symbols)]
#symbols <- symbols[match(rownames(res), names(symbols))]
#rownames(res) <- symbols

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=gsub("\\..*","",row.names(res)),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

EnhancedVolcano(res,
                lab = res$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 3.0,
                labSize = 5.0,
                encircleAlpha = 0.5,
                colAlpha = 1,
                legendLabels=c('Not sig.','FC','P value',
                               'Sign. Up/Down'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 2.0)



#-----------------5 HOURS-------------------------------------------------------------------------------

expression.data_5hours <- expression.data[,c(3,4,7,8)]

# Differential expression analysis
expr.matrix_5hours <- as.matrix(expression.data_5hours)

# Provide description of our data to DESeq2
expr.design_5hours <- data.frame(row.names = colnames(expression.data_5hours),
                          condition = c(
                            "control", "control",
                            "myxothiazole", "myxothiazole"
                          ))


# Let's create a variable that summarizes all the data about our experiment
dds_5hours <- DESeqDataSetFromMatrix(countData = expr.matrix_5hours,
                              colData = expr.design_5hours,
                              design = ~ condition)

dds_5hours <- DESeq(dds_5hours) #Perform diff expr analysis
res_5hours <- results(dds_5hours) #save results of analysis

dds_5hours <- estimateSizeFactors(dds_5hours)
sizeFactors(dds_5hours) # Median (size factor) coefficients show
normalized.counts_5hours <- counts(dds_5hours, normalized = TRUE)

# Estimation of variance
plotDispEsts(dds_5hours)

# Results
plotMA(dds_5hours, ylim = c(-2,2), main = "DESeq2")

sum(res_5hours$padj < 0.05, na.rm = TRUE)
resSig_5hours <- subset(res_5hours, padj < 0.05)

# find genes with the biggest foldchange
the_biggest_fold_change_5hours <- res_5hours[order(abs(res_5hours$log2FoldChange), decreasing = TRUE), ]

# Calculate the number of up and down regulated genes
sum(resSig_5hours$log2FoldChange > 1, na.rm = TRUE)
sum(resSig_5hours$log2FoldChange < -1, na.rm = TRUE)

# Find upregulated and downregulated genes
resSigUp_5hours <- subset(resSig_5hours, log2FoldChange > 1)
resSigDown_5hours <- subset(resSig_5hours, log2FoldChange < -1)

write.table(as.data.frame(resSigUp_5hours), "/Users/suleymanov-ef/Desktop/Skoltech/Study/NGS/Project/resSigUp_5hours.txt",
            sep = "\t", col.names = TRUE, row.names = TRUE,
            na = "NA", quote = FALSE)

write.table(as.data.frame(resSigDown_5hours), "/Users/suleymanov-ef/Desktop/Skoltech/Study/NGS/Project/resSigDown_5hours.txt",
            sep = "\t", col.names = TRUE, row.names = TRUE,
            na = "NA", quote = FALSE)

#writeLines(as.character(rownames(resSigUp)), "C:/Users/aisin/Desktop/Upreg_Mitochondria.txt")


#writeLines(as.character(rownames(resSigDown)), "C:/Users/aisin/Desktop/Downreg_Mitochondria.txt")


# visualize the results

names_5hours <- c(row.names(resSigUp_5hours), row.names(resSigDown_5hours)) #make list of DE genes
de.genes_5hours <- subset(normalized.counts_5hours, rownames(normalized.counts_5hours)
                   %in% names) # create table with norm_counts

# we can exclude 0 values, or add real low value to all our data
de.genes_5hours <- de.genes_5hours + 0.01
log.de.genes_5hours <- log10(de.genes_5hours) # log again

# Scale data
scaled.de.genes_5hours <- scale(log.de.genes_5hours)
scaled.de.genes_5hours <- t(scale(t(log.de.genes_5hours)))
Heatmap(scaled.de.genes_5hours, show_row_names = FALSE)

# Volcano plot
ens_5hours <- rownames(expression.data_5hours)

#symbols <- mapIds(org.Hs.eg.db, keys = ens,
#                  column = 'SYMBOL', keytype = 'ENSEMBL')

#symbols <- symbols[!is.na(symbols)]
#symbols <- symbols[match(rownames(res), names(symbols))]
#rownames(res) <- symbols

res_5hours$symbol <- mapIds(org.Hs.eg.db,
                     keys=gsub("\\..*","",row.names(res_5hours)),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

EnhancedVolcano(res_5hours,
                lab = res_5hours$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 3.0,
                labSize = 5.0,
                encircleAlpha = 0.5,
                colAlpha = 1,
                legendLabels=c('Not sig.','FC','P value',
                               'Sign. Up/Down'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 2.0)



#-----------------17 HOURS-------------------------------------------------------------------------------

expression.data_17hours <- expression.data[,c(1,2,5,6)]

# Differential expression analysis
expr.matrix_17hours <- as.matrix(expression.data_17hours)

# Provide description of our data to DESeq2
expr.design_17hours <- data.frame(row.names = colnames(expression.data_17hours),
                                 condition = c(
                                   "control", "control",
                                   "myxothiazole", "myxothiazole"
                                 ))


# Let's create a variable that summarizes all the data about our experiment
dds_17hours <- DESeqDataSetFromMatrix(countData = expr.matrix_17hours,
                                     colData = expr.design_17hours,
                                     design = ~ condition)

dds_17hours <- DESeq(dds_17hours) #Perform diff expr analysis
res_17hours <- results(dds_17hours) #save results of analysis

dds_17hours <- estimateSizeFactors(dds_17hours)
sizeFactors(dds_17hours) # Median (size factor) coefficients show
normalized.counts_17hours <- counts(dds_17hours, normalized = TRUE)

# Estimation of variance
plotDispEsts(dds_17hours)

# Results
plotMA(dds_17hours, ylim = c(-2,2), main = "DESeq2")

sum(res_17hours$padj < 0.05, na.rm = TRUE)
resSig_17hours <- subset(res_17hours, padj < 0.05)

# find genes with the biggest foldchange
the_biggest_fold_change_17hours <- res_17hours[order(abs(res_17hours$log2FoldChange), decreasing = TRUE), ]

# Calculate the number of up and down regulated genes
sum(resSig_17hours$log2FoldChange > 1, na.rm = TRUE)
sum(resSig_17hours$log2FoldChange < -1, na.rm = TRUE)

# Find upregulated and downregulated genes
resSigUp_17hours <- subset(resSig_17hours, log2FoldChange > 1)
resSigDown_17hours <- subset(resSig_17hours, log2FoldChange < -1)

write.table(as.data.frame(resSigUp_17hours), "/Users/suleymanov-ef/Desktop/Skoltech/Study/NGS/Project/resSigUp_17hours.txt",
            sep = "\t", col.names = TRUE, row.names = TRUE,
            na = "NA", quote = FALSE)

write.table(as.data.frame(resSigDown_17hours), "/Users/suleymanov-ef/Desktop/Skoltech/Study/NGS/Project/resSigDown_17hours.txt",
            sep = "\t", col.names = TRUE, row.names = TRUE,
            na = "NA", quote = FALSE)

#writeLines(as.character(rownames(resSigUp)), "C:/Users/aisin/Desktop/Upreg_Mitochondria.txt")


#writeLines(as.character(rownames(resSigDown)), "C:/Users/aisin/Desktop/Downreg_Mitochondria.txt")

# visualize the results

names_17hours <- c(row.names(resSigUp_17hours), row.names(resSigDown_17hours)) #make list of DE genes
de.genes_17hours <- subset(normalized.counts_17hours, rownames(normalized.counts_17hours)
                          %in% names) # create table with norm_counts

# we can exclude 0 values, or add real low value to all our data
de.genes_17hours <- de.genes_17hours + 0.01
log.de.genes_17hours <- log10(de.genes_17hours) # log again

# Scale data
scaled.de.genes_17hours <- scale(log.de.genes_17hours)
scaled.de.genes_17hours <- t(scale(t(log.de.genes_17hours)))
Heatmap(scaled.de.genes_17hours, show_row_names = FALSE)

# Volcano plot
ens_17hours <- rownames(expression.data_17hours)

#symbols <- mapIds(org.Hs.eg.db, keys = ens,
#                  column = 'SYMBOL', keytype = 'ENSEMBL')

#symbols <- symbols[!is.na(symbols)]
#symbols <- symbols[match(rownames(res), names(symbols))]
#rownames(res) <- symbols

res_17hours$symbol <- mapIds(org.Hs.eg.db,
                            keys=gsub("\\..*","",row.names(res_17hours)),
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")

EnhancedVolcano(res_17hours,
                lab = res_17hours$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 3.0,
                labSize = 5.0,
                encircleAlpha = 0.5,
                colAlpha = 1,
                legendLabels=c('Not sig.','FC','P value',
                               'Sign. Up/Down'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 2.0)


#-----------------GENE INTERSECTION-------------------------------------------------------------------------------

length(rownames(resSigUp_17hours))
length(rownames(resSigDown_17hours))
length(rownames(resSigUp_5hours))
length(rownames(resSigDown_5hours))

# Gene intersection
common_genes_upregulates <- intersect(rownames(resSigUp_5hours), rownames(resSigUp_17hours))
common_genes_downregulates <- intersect(rownames(resSigDown_5hours), rownames(resSigDown_17hours))

# Length of intersection of un- and down- regulated genes
length(common_genes_upregulates)
length(common_genes_downregulates)

# Upregulated genes GO enrichment
GO_results_upregulated_BP <- enrichGO(gene = common_genes_upregulates, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_results_upregulated_MF <- enrichGO(gene = common_genes_upregulates, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
GO_results_upregulated_CC <- enrichGO(gene = common_genes_upregulates, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")

upregulated_genes_plot_BP <- plot(barplot(GO_results_upregulated_BP, showCategory = 15))
png("upregulated_genes_plot_BP.png", res = 250, width = 1400, height = 1800)
print(upregulated_genes_plot_BP)
dev.off()

upregulated_genes_plot_MF <- plot(barplot(GO_results_upregulated_MF, showCategory = 15))
png("upregulated_genes_plot_MF.png", res = 250, width = 1400, height = 1800)
print(upregulated_genes_plot_MF)
dev.off()

upregulated_genes_plot_CC <- plot(barplot(GO_results_upregulated_CC, showCategory = 15))
png("upregulated_genes_plot_CC.png", res = 250, width = 1400, height = 1800)
print(upregulated_genes_plot_CC)
dev.off()


# Downregulated genes GO enrichment
GO_results_downregulated_BP <- enrichGO(gene = common_genes_downregulates, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_results_downregulated_MF <- enrichGO(gene = common_genes_downregulates, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
GO_results_downregulated_CC <- enrichGO(gene = common_genes_downregulates, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")

downregulated_genes_plot_BP <- plot(barplot(GO_results_downregulated_BP, showCategory = 15))
png("downregulated_genes_plot_BP.png", res = 250, width = 1400, height = 1800)
print(downregulated_genes_plot_BP)
dev.off()

downregulated_genes_plot_MF <- plot(barplot(GO_results_downregulated_MF, showCategory = 15))
png("downregulated_genes_plot_MF.png", res = 250, width = 1400, height = 1800)
print(downregulated_genes_plot_MF)
dev.off()

downregulated_genes_plot_CC <- plot(barplot(GO_results_downregulated_CC, showCategory = 15))
png("downregulated_genes_plot_CC.png", res = 250, width = 1400, height = 1800)
print(downregulated_genes_plot_CC)
dev.off()

#-----------------GENES FOR 17 HOURS-------------------------------------------------------------------------------


# Upregulated genes GO enrichment
GO_results_upregulated_17hours_BP <- enrichGO(gene = rownames(resSigUp_17hours), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_results_upregulated_17hours_MF <- enrichGO(gene = rownames(resSigUp_17hours), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
GO_results_upregulated_17hours_CC <- enrichGO(gene = rownames(resSigUp_17hours), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")

upregulated_genes_plot_17hours_BP <- plot(barplot(GO_results_upregulated_17hours_BP, showCategory = 15))
png("upregulated_genes_plot_17hours_BP.png", res = 250, width = 1400, height = 1800)
print(upregulated_genes_plot_17hours_BP)
dev.off()

upregulated_genes_plot_17hours_MF <- plot(barplot(GO_results_upregulated_17hours_MF, showCategory = 20))
png("upregulated_genes_plot_17hours_MF.png", res = 250, width = 1400, height = 1800)
print(upregulated_genes_plot_17hours_MF)
dev.off()

upregulated_genes_plot_17hours_CC <- plot(barplot(GO_results_upregulated_17hours_CC, showCategory = 20))
png("upregulated_genes_plot_17hours_CC.png", res = 250, width = 1400, height = 1800)
print(upregulated_genes_plot_17hours_CC)
dev.off()


# Downregulated genes GO enrichment
GO_results_downregulated_17hours_BP <- enrichGO(gene = rownames(resSigDown_17hours), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_results_downregulated_17hours_MF <- enrichGO(gene = rownames(resSigDown_17hours), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
GO_results_downregulated_17hours_CC <- enrichGO(gene = rownames(resSigDown_17hours), OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")

downregulated_genes_plot_17hours_BP <- plot(barplot(GO_results_downregulated_17hours_BP, showCategory = 15))
png("downregulated_genes_plot_17hours_BP.png", res = 250, width = 1400, height = 1800)
print(downregulated_genes_plot_17hours_BP)
dev.off()

downregulated_genes_plot_17hours_MF <- plot(barplot(GO_results_downregulated_17hours_MF, showCategory = 20))
png("downregulated_genes_plot_17hours_MF.png", res = 250, width = 1400, height = 1800)
print(downregulated_genes_plot_17hours_MF)
dev.off()

downregulated_genes_plot_17hours_CC <- plot(barplot(GO_results_downregulated_17hours_CC, showCategory = 15))
png("downregulated_genes_plot_17hours_CC.png", res = 250, width = 1400, height = 1800)
print(downregulated_genes_plot_17hours_CC)
dev.off()


#-----------------MITOCHONDRIAL GENES-------------------------------------------------------------------------------

mitochondrial_genes <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")

mitochondrial_genes_names <- mitochondrial_genes$gene_id

# Gene intersection with mitochondrial genes
common_genes_upregulates_5hours_mitochondrial <- intersect(rownames(resSigUp_5hours), mitochondrial_genes_names)
common_genes_upregulates_17hours_mitochondrial <- intersect(rownames(resSigUp_17hours), mitochondrial_genes_names)

common_genes_downregulates_5hours_mitochondrial <- intersect(rownames(resSigDown_5hours), mitochondrial_genes_names)
common_genes_downregulates_17hours_mitochondrial <- intersect(rownames(resSigDown_17hours), mitochondrial_genes_names)

length(common_genes_upregulates_5hours_mitochondrial)
length(common_genes_upregulates_17hours_mitochondrial)
length(common_genes_downregulates_5hours_mitochondrial)
length(common_genes_downregulates_17hours_mitochondrial)

# Upregulated genes GO enrichment
#GO_results_downregulated_mitochondrial_17hours_BP <- enrichGO(gene = common_genes_downregulates_17hours_mitochondrial, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
#GO_results_downregulated_mitochondrial_17hours_MF <- enrichGO(gene = common_genes_downregulates_17hours_mitochondrial, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "MF")
#GO_results_downregulated_mitochondrial_17hours_CC <- enrichGO(gene = common_genes_downregulates_17hours_mitochondrial, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "CC")

#downregulated_mitochondrial_genes_plot_17hours_BP <- plot(barplot(GO_results_downregulated_mitochondrial_17hours_BP, showCategory = 15))
#png("downregulated_mitochondrial_genes_plot_17hours_BP.png", res = 250, width = 1400, height = 1800)
#print(downregulated_mitochondrial_genes_plot_17hours_BP)
#dev.off()

#downregulated_mitochondrial_genes_plot_17hours_MF <- plot(barplot(GO_results_downregulated_mitochondrial_17hours_MF, showCategory = 20))
#png("downregulated_mitochondrial_genes_plot_17hours_MF.png", res = 250, width = 1400, height = 1800)
#print(downregulated_mitochondrial_genes_plot_17hours_MF)
#dev.off()

#downregulated_mitochondrial_genes_plot_17hours_CC <- plot(barplot(GO_results_downregulated_mitochondrial_17hours_CC, showCategory = 20))
#png("downregulated_mitochondrial_genes_plot_17hours_CC.png", res = 250, width = 1400, height = 1800)
#print(downregulated_mitochondrial_genes_plot_17hours_CC)
#dev.off()


#-----------------INTERSECTION PLOTS----------------------------------------------------------------------------

# Create DFs consising of gene names, fold_change type and time
resSigUp_17hours_df <- as.data.frame(rownames(resSigUp_17hours)) %>%
  mutate(Type = 'Upregulation') %>%
  mutate(Time = '17H')

resSigDown_17hours_df <- as.data.frame(rownames(resSigDown_17hours)) %>%
  mutate(Type = 'Downregulation') %>%
  mutate(Time = '17H')

resSigUp_5hours_df <- as.data.frame(rownames(resSigUp_5hours)) %>%
  mutate(Type = 'Upregulation') %>%
  mutate(Time = '5H')

resSigDown_5hours_df <- as.data.frame(rownames(resSigDown_5hours)) %>%
  mutate(Type = 'Downregulation') %>%
  mutate(Time = '5H')

# Change column names
colnames(resSigUp_17hours_df) <- c('Genes', 'Type', 'Time')
colnames(resSigDown_17hours_df) <- c('Genes', 'Type', 'Time')
colnames(resSigUp_5hours_df) <- c('Genes', 'Type', 'Time')
colnames(resSigDown_5hours_df) <- c('Genes', 'Type', 'Time')


intersected_genes <- bind_rows(resSigUp_17hours_df, resSigDown_17hours_df, resSigUp_5hours_df, resSigDown_5hours_df)

intersected_genes_upset <- intersected_genes %>%
  dplyr::select(c(Genes, c(Type, Time))) %>%
  mutate(Count = 1) %>%
  pivot_wider(names_from = Genes, values_from = Count) %>%
  replace(is.na(.), 0) %>%
  unite(Type_time, 'Type', 'Time')

genes <- colnames(intersected_genes_upset)[2:length(intersected_genes_upset)]

intersected_genes_upset[genes] = intersected_genes_upset[genes] == 1

intersected_genes_upset_matrix <- intersected_genes_upset %>%
  tibble::column_to_rownames(var="Type_time") %>%
  t() %>%
  as.data.frame()

genes_upset <- colnames(intersected_genes_upset_matrix)

intersection_plot <- upset(intersected_genes_upset_matrix, 
                           genes_upset, 
                           name = "Intersection", 
                           width_ratio=0.1,
                           labeller=ggplot2::as_labeller(c(
                             'Downregulation_17H'='DownReg_17',
                             'Downregulation_5H'='DownReg_5',
                             'Upregulation_17H'='UpReg_17',
                             'Upregulation_5H'='UpReg_5')))

intersection_plot



#-----------------BUBBLE PLOT MITOCHONDRIAL GENES------------------
res_mitochondrial <- res

res_mitochondrial$name <- row.names(res_mitochondrial)
res_mitochondrial <- subset(res_mitochondrial,res_mitochondrial$name %in% mitochondrial_genes_names)

EnhancedVolcano(res_mitochondrial,
                lab = rownames(res_mitochondrial),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-02,
                caption = bquote(~Log[2]~ "fold change cutoff, 2; p-value cutoff, 10e-4"),
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 3.0,
                labSize = 0.0,
                colAlpha = 1,
                legendLabels=c('Not sig.','FC','P value',
                               'Sign. Up/Down'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 2.0)

mitochondrial_upreg <- subset(res_mitochondrial, log2FoldChange > 1 & padj < 0.05)
mitochondrial_downreg <- subset(res_mitochondrial, log2FoldChange < -1 & padj < 0.05)








