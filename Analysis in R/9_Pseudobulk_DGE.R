

setwd("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/COSMX_colon_dataset/COSMX_Colonrectal_cancer_project/TMA_original/Seurat")

library(dplyr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggplot2)
library(dittoSeq)
library(Seurat)
library(Azimuth)
library(SeuratDisk)
library(SeuratData)
library(viridis)
library(RColorBrewer)
library(DoubletFinder)
library(DropletUtils)
library(dsb)
library(harmony)
library(EnhancedVolcano)
library(basilisk)
library(zellkonverter)
library(reticulate)
library(Matrix)
library(org.Hs.eg.db)
library(gridExtra)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(ggridges)
library(ExperimentHub)
library(DESeq2)
library(tidyverse)

options(ggrepel.max.overlaps = Inf)

fov_integrated <- readRDS("./fov_integrated.rds")

saveRDS(fov_integrated, "./fov_integrated.rds")

# Pseudo-bulk DGE using DESeq2

# Reference: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions
# Reference: https://www.sc-best-practices.org/conditions/differential_gene_expression.html#motivation
# Reference: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-are-some-p-values-set-to-na
# Reference: http://bioconductor.org/books/3.14/OSCA.multisample/multi-sample-comparisons.html#testing-for-between-label-differences

# Generally, both, pseudobulk methods with sum aggregation such as edgeR, DESeq2, or Limma[Ritchie et al., 2015] and mixed models such as MAST with random effect setting were found to be superior compared to naive methods, 
# such as the popular Wilcoxon rank-sum test or Seurat’s [Hao et al., 2021] latent models, which do not account for them[Junttila et al., 2022].

DefaultAssay(fov_integrated) <- "RNA"
Idents(fov_integrated) <- 'Epithelial_cluster'
levels(fov_integrated)

# 1- Aggregate expression based on CIPR cell types and Epithelial clusters

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample

# 1.1- counts matrix - sample level
# counts aggregate to sample level

View(fov_integrated@meta.data)
fov_integrated$DGE_epithelial_cluster_pseudobulk <- paste0(fov_integrated$sample_ID, sep = "_", fov_integrated$Epithelial_cluster)
View(fov_integrated@meta.data)

DefaultAssay(fov_integrated)

# Sum expression
cts <- AggregateExpression(fov_integrated, 
                           group.by = c("CIPR_ordered_clusters", "DGE_epithelial_cluster_pseudobulk"), #group by cell type and sample/condition
                           assays = 'RNA', #use default assay
                           slot = "counts", #DESeq2 uses raw counts.
                           return.seurat = FALSE)

# or average expression
cts <- AverageExpression(fov_integrated, 
                           group.by = c("CIPR_ordered_clusters", "DGE_epithelial_cluster_pseudobulk"), #group by cell type and sample/condition
                           assays = 'RNA', #use default assay
                           slot = "counts", #DESeq2 uses raw counts.
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)
write.csv(cts.t, "Epithelial_clusters_pseudo_bulk_expression.csv")

# get values where to split by cell types
splitRows <- gsub('_.*', '', rownames(cts.t))

# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))
# quick check
cts.split$'B cell'[1:10,1:10]

# fix colnames (remove cell type names) and transpose the matrix
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

cts.split.modified$'B cell'[1:10,1:10]

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101') test the regex to see if works

# 1.2- DE analysis with DESeq 2

# 1.2.1- Get counts matrix for each cell type of interest
counts_mac <- cts.split.modified$`SPP1+ macrophage`
counts_cd4 <- cts.split.modified$`CD8+ T cell`
counts_cd8 <- cts.split.modified$`CD4+ T cell`
counts_bcell <- cts.split.modified$`B cell`
counts_plasma <- cts.split.modified$`Plasma cell`
counts_fibro <- cts.split.modified$`Fibroblast`
counts_ep <- cts.split.modified$`Tumor epithelial cell`

write.csv(counts_mac, "SPP1_macrophage_epithelial_clusters_pseudo_bulk_expression.csv")
write.csv(counts_cd4, "CD4Tcell_epithelial_clusters_pseudo_bulk_expression.csv")
write.csv(counts_cd8, "CD8Tcell_epithelial_clusters_pseudo_bulk_expression.csv")
write.csv(counts_bcell, "Bcell_epithelial_clusters_pseudo_bulk_expression.csv")
write.csv(counts_plasma, "Plasmacell_epithelial_clusters_pseudo_bulk_expression.csv")
write.csv(counts_fibro, "Fibroblast_epithelial_clusters_pseudo_bulk_expression.csv")
write.csv(counts_ep, "Tumor_epithelialcell_epithelial_clusters_pseudo_bulk_expression.csv")

# 1.2.2. generate sample level metadata for each cell type of interest
colData <- data.frame(samples = colnames(counts_ep))

colData <- colData %>%
  mutate(condition = ifelse(grepl('-1', samples), 'Epithelial_cluster_1', ifelse(grepl('-2', samples),'Epithelial_cluster_2', 'Epithelial_cluster_3')),
         ) %>%
  column_to_rownames(var = 'samples')

#make a specific to cd4
colData2 <- data.frame(samples = colnames(counts_cd4))
colData2 <- colData2 %>%
  mutate(condition = ifelse(grepl('-1', samples), 'Epithelial_cluster_1', ifelse(grepl('-2', samples),'Epithelial_cluster_2', 'Epithelial_cluster_3')),
  ) %>%
  column_to_rownames(var = 'samples')

#make a specific to cd8
colData3 <- data.frame(samples = colnames(counts_cd8))
colData3 <- colData3 %>%
  mutate(condition = ifelse(grepl('-1', samples), 'Epithelial_cluster_1', ifelse(grepl('-2', samples),'Epithelial_cluster_2', 'Epithelial_cluster_3')),
  ) %>%
  column_to_rownames(var = 'samples')


# 1.2.3- perform DESeq2
# Create DESeq2 object   
dds_mac <- DESeqDataSetFromMatrix(countData = counts_mac, colData = colData, design = ~ condition)
dds_cd4 <- DESeqDataSetFromMatrix(countData = counts_cd4, colData = colData2, design = ~ condition)
dds_cd8 <- DESeqDataSetFromMatrix(countData = counts_cd8, colData = colData3, design = ~ condition)
dds_bcell <- DESeqDataSetFromMatrix(countData = counts_bcell, colData = colData, design = ~ condition)
dds_plasma <- DESeqDataSetFromMatrix(countData = counts_plasma, colData = colData, design = ~ condition)
dds_fibro <- DESeqDataSetFromMatrix(countData = counts_fibro, colData = colData, design = ~ condition)
dds_ep <- DESeqDataSetFromMatrix(countData = counts_ep, colData = colData, design = ~ condition)

# filter out genes that have lower than 10 reads
keep <- rowSums(counts(dds_mac)) >=10 
dds_mac <- dds_mac[keep,]

keep <- rowSums(counts(dds_cd4)) >=10 
dds_cd4 <- dds_cd4[keep,]

keep <- rowSums(counts(dds_cd8)) >=10 
dds_cd8 <- dds_cd8[keep,]

keep <- rowSums(counts(dds_bcell)) >=10 
dds_bcell <- dds_bcell[keep,]

keep <- rowSums(counts(dds_plasma)) >=10 
dds_plasma <- dds_plasma[keep,]

keep <- rowSums(counts(dds_fibro)) >=10 
dds_fibro <- dds_fibro[keep,]

keep <- rowSums(counts(dds_ep)) >=10 
dds_ep <- dds_ep[keep,]

# run DESeq2
dds_mac <- DESeq(dds_mac)
dds_cd4 <- DESeq(dds_cd4)
dds_cd8 <- DESeq(dds_cd8)
dds_bcell <- DESeq(dds_bcell)
dds_plasma <- DESeq(dds_plasma)
dds_fibro <- DESeq(dds_fibro)
dds_ep <- DESeq(dds_ep)

saveRDS(dds_mac, "dds_mac.rds")
saveRDS(dds_cd4, "dds_cd4.rds")
saveRDS(dds_cd8, "dds_cd8.rds")
saveRDS(dds_bcell, "dds_bcell.rds")
saveRDS(dds_plasma, "dds_plasma.rds")
saveRDS(dds_fibro, "dds_fibro.rds")
saveRDS(dds_ep, "dds_ep.rds")


# Check the coefficients for the comparison
resultsNames(dds_mac)
resultsNames(dds_cd4)
resultsNames(dds_cd8)
resultsNames(dds_bcell)
resultsNames(dds_plasma)
resultsNames(dds_fibro)
resultsNames(dds_ep)

# Generate results object

# Independent filtering and multiple testing
# Filtering criteria
# The goal of independent filtering is to filter out those tests from the procedure that have no, or little chance of showing significant evidence, without even looking at their test statistic. 
# Typically, this results in increased detection power at the same experiment-wide type I error. 
# Here, we measure experiment-wide type I error in terms of the false discovery rate.
# A good choice for a filtering criterion is one that
# is statistically independent from the test statistic under the null hypothesis,
# is correlated with the test statistic under the alternative, and
# does not notably change the dependence structure – if there is any – between the tests that pass the filter, compared to the dependence structure between the tests before filtering.
# The benefit from filtering relies on property (2), and we will explore it further below. 
# Its statistical validity relies on property (1) – which is simple to formally prove for many combinations of filter criteria with test statistics – 
# and (3), which is less easy to theoretically imply from first principles, but rarely a problem in practice. 
# We refer to (Bourgon, Gentleman, and Huber 2010) for further discussion of this topic.
# A simple filtering criterion readily available in the results object is the mean of normalized counts irrespective of biological condition, and so this is the criterion which is used automatically by the results function to perform independent filtering. 
# Genes with very low counts are not likely to see significant differences typically due to high dispersion. For example, we can plot the −log10 p values from all genes over the normalized mean counts:

# Genereta DFs comparing ep clusters 1 vs 3 (# 1 and 2 are similar, so no need to run 2 vs 3)
contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")

# SPP1 macrophage
res_mac <- results(dds_mac, contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")) 
res_mac2 <- as.data.frame(res_mac)
res_mac2 <- na.omit(res_mac2)
saveRDS(res_mac, "res_mac.rds")
write.csv(res_mac2, "SPP1Macro_EPcluster1_vs_EPcluster3.csv")

# CD4 T cell
res_cd4 <- results(dds_cd4, contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")) 
res_cd4_2 <- as.data.frame(res_cd4)
res_cd4_2 <- na.omit(res_cd4_2)
saveRDS(res_cd4, "res_cd4.rds")
write.csv(res_cd4_2, "CD4Tcell_EPcluster1_vs_EPcluster3.csv")

# CD8 T cell
res_cd8 <- results(dds_cd8, contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")) 
res_cd8_2 <- as.data.frame(res_cd8)
res_cd8_2 <- na.omit(res_cd8_2)
saveRDS(res_cd8, "res_cd8.rds")
write.csv(res_cd8_2, "CD8Tcell_EPcluster1_vs_EPcluster3.csv")

# B cell
res_bcell <- results(dds_bcell, contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")) 
res_bcell_2 <- as.data.frame(res_bcell)
res_bcell_2 <- na.omit(res_bcell_2)
saveRDS(res_bcell, "res_bcell.rds")
write.csv(res_bcell_2, "Bcell_EPcluster1_vs_EPcluster3.csv")

# Plasma cell
res_plasma <- results(dds_plasma, contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")) 
res_plasma_2 <- as.data.frame(res_plasma)
res_plasma_2 <- na.omit(res_plasma_2)
saveRDS(res_plasma, "res_plasma.rds")
write.csv(res_plasma_2, "Plasma_EPcluster1_vs_EPcluster3.csv")

# Fibroblast
res_fibro <- results(dds_fibro, contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")) 
res_fibro_2 <- as.data.frame(res_fibro)
res_fibro_2 <- na.omit(res_fibro_2)
saveRDS(res_fibro, "res_fibro.rds")
write.csv(res_fibro_2, "Fibro_EPcluster1_vs_EPcluster3.csv")

# Tumor ep cell
res_ep <- results(dds_ep, contrast=c("condition","Epithelial_cluster_1","Epithelial_cluster_3")) 
res_ep_2 <- as.data.frame(res_ep)
res_ep_2 <- na.omit(res_ep_2)
saveRDS(res_ep, "res_ep.rds")
write.csv(res_ep_2, "Tumor_EPcluster1_vs_EPcluster3.csv")

#res <- results(dds, name = "condition_Epithelial_cluster_2_vs_Epithelial_cluster_1") #use independentFiltering=FALSE if dont want NA results in padj values
#res
#res2 <- results(dds, name = "condition_Epithelial_cluster_3_vs_Epithelial_cluster_1") #use independentFiltering=FALSE if dont want NA results in padj values
#res2

# 1.3- DE analysis with edgeR

# Reference: https://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/single_cell_edger.html#:~:text=A%20pseudo%2Dbulk%20single%20cell,the%20minimal%20amount%20of%20code.
library(Glimma)
library(edgeR)

counts_mac <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Expression/SPP1_macrophage_epithelial_clusters_pseudo_bulk_expression.csv", row.names = 1)
counts_cd4 <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Expression/CD4Tcell_epithelial_clusters_pseudo_bulk_expression.csv", row.names = 1)
counts_cd8 <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Expression/CD8Tcell_epithelial_clusters_pseudo_bulk_expression.csv", row.names = 1)
counts_bcell <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Expression/Bcell_epithelial_clusters_pseudo_bulk_expression.csv", row.names = 1)
counts_plasma <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Expression/Plasmacell_epithelial_clusters_pseudo_bulk_expression.csv", row.names = 1)
counts_fibro <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Expression/Fibroblast_epithelial_clusters_pseudo_bulk_expression.csv", row.names = 1)

colData <- data.frame(samples = colnames(counts_mac))
colData <- colData %>%
  mutate(condition = ifelse(grepl('_1', samples), 'Epithelial_cluster_1', ifelse(grepl('_2', samples),'Epithelial_cluster_2', 'Epithelial_cluster_3')),
  ) %>%
  column_to_rownames(var = 'samples')

#make a specific to cd4
colData2 <- data.frame(samples = colnames(counts_cd4))
colData2 <- colData2 %>%
  mutate(condition = ifelse(grepl('.1', samples), 'Epithelial_cluster_1', ifelse(grepl('.2', samples),'Epithelial_cluster_2', 'Epithelial_cluster_3')),
  ) %>%
  column_to_rownames(var = 'samples')

#make a specific to cd8
colData3 <- data.frame(samples = colnames(counts_cd8))
colData3 <- colData3 %>%
  mutate(condition = ifelse(grepl('.1', samples), 'Epithelial_cluster_1', ifelse(grepl('.2', samples),'Epithelial_cluster_2', 'Epithelial_cluster_3')),
  ) %>%
  column_to_rownames(var = 'samples')

mac_dge <- DGEList(counts = counts_mac, samples = colData, group = colData$condition)
cd4_dge <- DGEList(counts = counts_cd4, samples = colData2, group = colData2$condition)
cd8_dge <- DGEList(counts = counts_cd8, samples = colData3, group = colData3$condition)
bcell_dge <- DGEList(counts = counts_bcell, samples = colData, group = colData$condition)
plasma_dge <- DGEList(counts = counts_plasma, samples = colData, group = colData$condition)
fibro_dge <- DGEList(counts = counts_fibro, samples = colData, group = colData$condition)

# With this we perform differential expression analysis using edgeR’s generalised linear models.
design <- model.matrix(~0 + condition, data = fibro_dge$samples)
colnames(design) <- make.names(gsub("condition", "", colnames(design)))
fibro_dge <- estimateDisp(fibro_dge, design)
contr <- makeContrasts("Epithelial_cluster_1 - Epithelial_cluster_3", levels = design)

mac_fit <- glmFit(fibro_dge, design)
mac_lrt <- glmLRT(mac_fit, contrast = contr)

# The results of this analysis can be visualised using glimmaMA() as it would be for bulk RNA-seq.
glimmaMA(mac_lrt, dge = fibro_dge)

# the DEA result for all the genes
mac_dea <- mac_lrt$table

# differentially expressed genes
# a data frame containing the elements 'logFC', the log-abundance ratio, i.e. fold change, for each tag in the two groups being compared, 
# 'logCPM', the log-average concentration/abundance for each tag in the two groups being compared, 'PValue', exact p-value for differential expression using the NB model. 
# When 'adjust.method' is not '"none"', there is an extra column of 'FDR' showing the adjusted p-value if 'adjust.method' is one of the '"BH"', '"BY"' and
# '"fdr"', or an extra column of 'FWER' if 'adjust.method' is one of the '"holm"', '"hochberg"', '"hommel"', and '"bonferroni"'.
toptag <- topTags(mac_lrt, n=958)
mac_deg <- toptag$table
write.csv(mac_deg, "SPP1Macro_EPcluster1_vs_EPcluster3_edgeR.csv")

# 1.4- Plot DEG - volcano plots and heatmaps

# Bcell, fibroblast, plasma cell, SPP1 macro, Tumor ep cell

options(ggrepel.max.overlaps = Inf)
DGE_bcell <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Lists_DEG_1vs3/Bcell_EPcluster1_vs_EPcluster3.csv")
DGE_plasma <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Lists_DEG_1vs3/Plasma_EPcluster1_vs_EPcluster3.csv")
DGE_fibro <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Lists_DEG_1vs3/Fibro_EPcluster1_vs_EPcluster3.csv")
DGE_mac <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Lists_DEG_1vs3/SPP1Macro_EPcluster1_vs_EPcluster3.csv")
DGE_ep <- read.csv("./Plots/10_DGE/2_Pseudo_bulk/1_Epithelial_clusters/Lists_DEG_1vs3/Tumor_EPcluster1_vs_EPcluster3.csv")

EnhancedVolcano(toptable = DGE_bcell, x = "log2FoldChange", y = "padj",
                lab = DGE_bcell$Features, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE B cell Epithelial cluster 1 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = DGE_plasma, x = "log2FoldChange", y = "padj",
                lab = DGE_plasma$Features, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Plasma cell Epithelial cluster 1 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = DGE_fibro, x = "log2FoldChange", y = "padj",
                lab = DGE_fibro$Features, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Fibroblast Epithelial cluster 1 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = DGE_mac, x = "log2FoldChange", y = "padj",
                lab = DGE_mac$Features, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1 macrophage Epithelial cluster 1 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = DGE_ep, x = "log2FoldChange", y = "padj",
                lab = DGE_ep$Features, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumor epithelial cell Epithelial cluster 1 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

###################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################


# 2- Aggregate expression based on CIPR cell types and Stromal clusters

DefaultAssay(fov_integrated) <- "RNA"
Idents(fov_integrated) <- 'aSMA_cluster'
levels(fov_integrated)

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample

# 1.1- counts matrix - sample level
# counts aggregate to sample level

View(fov_integrated@meta.data)
fov_integrated$DGE_stromal_cluster_pseudobulk <- paste0(fov_integrated$sample_ID, sep = "_", fov_integrated$aSMA_cluster)
View(fov_integrated@meta.data)

DefaultAssay(fov_integrated)

cts <- AggregateExpression(fov_integrated, 
                           group.by = c("CIPR_ordered_clusters", "DGE_stromal_cluster_pseudobulk"), #group by cell type and sample/condition
                           assays = 'RNA', #use default assay
                           slot = "counts", #DESeq2 uses raw counts.
                           return.seurat = FALSE)
cts <- cts$RNA

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)
write.csv(cts.t, "Stromal_clusters_pseudo_bulk_expression.csv")

# get values where to split by cell types
splitRows <- gsub('_.*', '', rownames(cts.t))

# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))
# quick check
cts.split$'B cell'[1:10,1:10]

# fix colnames (remove cell type names) and transpose the matrix
cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

cts.split.modified$'B cell'[1:10,1:10]

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101') test the regex to see if works

# 2.2- DE analysis

# 2.2.1- Get counts matrix for each cell type of interest
counts_mac <- cts.split.modified$`SPP1+ macrophage`
counts_cd4 <- cts.split.modified$`CD8+ T cell`
counts_cd8 <- cts.split.modified$`CD4+ T cell`
counts_bcell <- cts.split.modified$`B cell`
counts_plasma <- cts.split.modified$`Plasma cell`
counts_fibro <- cts.split.modified$`Fibroblast`

write.csv(counts_mac, "SPP1_macrophage_stromal_clusters_pseudo_bulk_expression.csv")
write.csv(counts_cd4, "CD4Tcell_stromal_clusters_pseudo_bulk_expression.csv")
write.csv(counts_cd8, "CD8Tcell_stromal_clusters_pseudo_bulk_expression.csv")
write.csv(counts_bcell, "Bcell_stromal_clusters_pseudo_bulk_expression.csv")
write.csv(counts_plasma, "Plasmacell_stromal_clusters_pseudo_bulk_expression.csv")
write.csv(counts_fibro, "Fibroblast_stromal_clusters_pseudo_bulk_expression.csv")

# 2.2.2. generate sample level metadata for each cell type of interest
colData <- data.frame(samples = colnames(counts_bcell))

colData <- colData %>%
  mutate(condition = ifelse(grepl('-1', samples), 'Stromal_cluster_1', ifelse(grepl('-2', samples),'Stromal_cluster_2', ifelse(grepl('-3', samples),'Stromal_cluster_3', 'Stromal_cluster_NA'))),
  ) %>%
  column_to_rownames(var = 'samples')

#make a specific to cd4
colData2 <- data.frame(samples = colnames(counts_cd4))
colData2 <- colData2 %>%
  mutate(condition = ifelse(grepl('-1', samples), 'Stromal_cluster_1', ifelse(grepl('-2', samples),'Stromal_cluster_2', ifelse(grepl('-3', samples),'Stromal_cluster_3', 'Stromal_cluster_NA'))),
  ) %>%
  column_to_rownames(var = 'samples')

#make a specific to cd8
colData3 <- data.frame(samples = colnames(counts_cd8))
colData3 <- colData3 %>%
  mutate(condition = ifelse(grepl('-1', samples), 'Stromal_cluster_1', ifelse(grepl('-2', samples),'Stromal_cluster_2', ifelse(grepl('-3', samples),'Stromal_cluster_3', 'Stromal_cluster_NA')))
         ) %>%
  column_to_rownames(var = 'samples')


# 2.2.3- perform DESeq2
# Create DESeq2 object   
dds_mac <- DESeqDataSetFromMatrix(countData = counts_mac, colData = colData, design = ~ condition)
dds_cd4 <- DESeqDataSetFromMatrix(countData = counts_cd4, colData = colData2, design = ~ condition)
dds_cd8 <- DESeqDataSetFromMatrix(countData = counts_cd8, colData = colData3, design = ~ condition)
dds_bcell <- DESeqDataSetFromMatrix(countData = counts_bcell, colData = colData, design = ~ condition)
dds_plasma <- DESeqDataSetFromMatrix(countData = counts_plasma, colData = colData, design = ~ condition)
dds_fibro <- DESeqDataSetFromMatrix(countData = counts_fibro, colData = colData, design = ~ condition)

# filter out genes that have lower than 10 reads
keep <- rowSums(counts(dds_mac)) >=10 
dds_mac <- dds_mac[keep,]

keep <- rowSums(counts(dds_cd4)) >=10 
dds_cd4 <- dds_cd4[keep,]

keep <- rowSums(counts(dds_cd8)) >=10 
dds_cd8 <- dds_cd8[keep,]

keep <- rowSums(counts(dds_bcell)) >=10 
dds_bcell <- dds_bcell[keep,]

keep <- rowSums(counts(dds_plasma)) >=10 
dds_plasma <- dds_plasma[keep,]

keep <- rowSums(counts(dds_fibro)) >=10 
dds_fibro <- dds_fibro[keep,]

# run DESeq2
dds_mac <- DESeq(dds_mac)
dds_cd4 <- DESeq(dds_cd4)
dds_cd8 <- DESeq(dds_cd8)
dds_bcell <- DESeq(dds_bcell)
dds_plasma <- DESeq(dds_plasma)
dds_fibro <- DESeq(dds_fibro)

saveRDS(dds_mac, "dds_mac_st.rds")
saveRDS(dds_cd4, "dds_cd4_st.rds")
saveRDS(dds_cd8, "dds_cd8_st.rds")
saveRDS(dds_bcell, "dds_bcell_st.rds")
saveRDS(dds_plasma, "dds_plasma_st.rds")
saveRDS(dds_fibro, "dds_fibro_st.rds")


# Check the coefficients for the comparison
resultsNames(dds_mac)
resultsNames(dds_cd4)
resultsNames(dds_cd8)
resultsNames(dds_bcell)
resultsNames(dds_plasma)
resultsNames(dds_fibro)

# Generate results object

# Independent filtering and multiple testing
# Filtering criteria
# The goal of independent filtering is to filter out those tests from the procedure that have no, or little chance of showing significant evidence, without even looking at their test statistic. 
# Typically, this results in increased detection power at the same experiment-wide type I error. 
# Here, we measure experiment-wide type I error in terms of the false discovery rate.
# A good choice for a filtering criterion is one that
# is statistically independent from the test statistic under the null hypothesis,
# is correlated with the test statistic under the alternative, and
# does not notably change the dependence structure – if there is any – between the tests that pass the filter, compared to the dependence structure between the tests before filtering.
# The benefit from filtering relies on property (2), and we will explore it further below. 
# Its statistical validity relies on property (1) – which is simple to formally prove for many combinations of filter criteria with test statistics – 
# and (3), which is less easy to theoretically imply from first principles, but rarely a problem in practice. 
# We refer to (Bourgon, Gentleman, and Huber 2010) for further discussion of this topic.
# A simple filtering criterion readily available in the results object is the mean of normalized counts irrespective of biological condition, and so this is the criterion which is used automatically by the results function to perform independent filtering. 
# Genes with very low counts are not likely to see significant differences typically due to high dispersion. For example, we can plot the −log10 p values from all genes over the normalized mean counts:

# Genereta DFs comparing ep clusters 1 vs 3 (# 1 and 2 are similar, so no need to run 2 vs 3)
contrast=c("condition","Stromal_cluster_1","Stromal_cluster_3")

# SPP1 macrophage
res_mac <- results(dds_mac, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_3")) 
res_mac2 <- as.data.frame(res_mac)
res_mac2 <- na.omit(res_mac2)
saveRDS(res_mac, "res_mac.rds")
write.csv(res_mac2, "SPP1Macro_STcluster1_vs_STcluster3.csv")

# CD4 T cell
res_cd4 <- results(dds_cd4, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_3")) 
res_cd4_2 <- as.data.frame(res_cd4)
res_cd4_2 <- na.omit(res_cd4_2)
saveRDS(res_cd4, "res_cd4.rds")
write.csv(res_cd4_2, "CD4Tcell_STcluster1_vs_STcluster3.csv")

# CD8 T cell
res_cd8 <- results(dds_cd8, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_3")) 
res_cd8_2 <- as.data.frame(res_cd8)
res_cd8_2 <- na.omit(res_cd8_2)
saveRDS(res_cd8, "res_cd8.rds")
write.csv(res_cd8_2, "CD8Tcell_STcluster1_vs_STcluster3.csv")

# B cell
res_bcell <- results(dds_bcell, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_3")) 
res_bcell_2 <- as.data.frame(res_bcell)
res_bcell_2 <- na.omit(res_bcell_2)
saveRDS(res_bcell, "res_bcell.rds")
write.csv(res_bcell_2, "Bcell_STcluster1_vs_STcluster3.csv")

# Plasma cell
res_plasma <- results(dds_plasma, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_3")) 
res_plasma_2 <- as.data.frame(res_plasma)
res_plasma_2 <- na.omit(res_plasma_2)
saveRDS(res_plasma, "res_plasma.rds")
write.csv(res_plasma_2, "Plasma_STcluster1_vs_STcluster3.csv")

# Fibroblast
res_fibro <- results(dds_fibro, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_3")) 
res_fibro_2 <- as.data.frame(res_fibro)
res_fibro_2 <- na.omit(res_fibro_2)
saveRDS(res_fibro, "res_fibro.rds")
write.csv(res_fibro_2, "Fibro_STcluster1_vs_STcluster3.csv")

contrast=c("condition","Stromal_cluster_1","Stromal_cluster_2")

# SPP1 macrophage
res_mac <- results(dds_mac, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_2")) 
res_mac2 <- as.data.frame(res_mac)
res_mac2 <- na.omit(res_mac2)
saveRDS(res_mac, "res_mac.rds")
write.csv(res_mac2, "SPP1Macro_STcluster1_vs_STcluster2.csv")

# CD4 T cell
res_cd4 <- results(dds_cd4, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_2")) 
res_cd4_2 <- as.data.frame(res_cd4)
res_cd4_2 <- na.omit(res_cd4_2)
saveRDS(res_cd4, "res_cd4.rds")
write.csv(res_cd4_2, "CD4Tcell_STcluster1_vs_STcluster2.csv")

# CD8 T cell
res_cd8 <- results(dds_cd8, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_2")) 
res_cd8_2 <- as.data.frame(res_cd8)
res_cd8_2 <- na.omit(res_cd8_2)
saveRDS(res_cd8, "res_cd8.rds")
write.csv(res_cd8_2, "CD8Tcell_STcluster1_vs_STcluster2.csv")

# B cell
res_bcell <- results(dds_bcell, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_2")) 
res_bcell_2 <- as.data.frame(res_bcell)
res_bcell_2 <- na.omit(res_bcell_2)
saveRDS(res_bcell, "res_bcell.rds")
write.csv(res_bcell_2, "Bcell_STcluster1_vs_STcluster2.csv")

# Plasma cell
res_plasma <- results(dds_plasma, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_2")) 
res_plasma_2 <- as.data.frame(res_plasma)
res_plasma_2 <- na.omit(res_plasma_2)
saveRDS(res_plasma, "res_plasma.rds")
write.csv(res_plasma_2, "Plasma_STcluster1_vs_STcluster2.csv")

# Fibroblast
res_fibro <- results(dds_fibro, contrast=c("condition","Stromal_cluster_1","Stromal_cluster_2")) 
res_fibro_2 <- as.data.frame(res_fibro)
res_fibro_2 <- na.omit(res_fibro_2)
saveRDS(res_fibro, "res_fibro.rds")
write.csv(res_fibro_2, "Fibro_STcluster1_vs_STcluster2.csv")

#res <- results(dds, name = "condition_Epithelial_cluster_2_vs_Epithelial_cluster_1") #use independentFiltering=FALSE if dont want NA results in padj values
#res
#res2 <- results(dds, name = "condition_Epithelial_cluster_3_vs_Epithelial_cluster_1") #use independentFiltering=FALSE if dont want NA results in padj values
#res2

# 2.2.4- Plot results

###################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 3- DGE analysis with MAST

# Reference: https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html
# Reference: https://github.com/RGLab/MAST/blob/devel/vignettes/MAITAnalysis.Rmd
# Reference: https://www.sc-best-practices.org/conditions/differential_gene_expression.html#single-cell-specific

library(ggplot2)
library(GGally)
library(GSEABase)
library(limma)
library(reshape2)
library(data.table)
library(knitr)
library(stringr)
library(NMF)
library(rsvd)
library(RColorBrewer)
library(MAST)

# 3.1- Loading and transforming data
freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)

data(maits, package='MAST') # large list object
dim(maits$expressionmat) # cells in rows and genes in columns 
head(maits$expressionmat)
head(maits$cdat) # cell metadata
head(maits$fdat) # genes metadata

scaRaw <- FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat) # convert large list to a LargeSingleCellAssay >> convert Seurat to SingleCellAssay

# 3.2- Convert Seurat to SCE to SCA
# Reference: https://www.bioconductor.org/packages/devel/bioc/vignettes/MAST/inst/doc/MAST-interoperability.html

# As a SingleCellExperiment-derived package, MAST can easily be inserted into workflows with packages such as scran, scater, zinbwave, SCnorm and others. 
# Moreover, subclassing SingleCellExperiment/SummarizedExperiment provides a flexible abstraction for the assay that contains the actual expression data. 
# It can use sparse Matrix and HDF5 as backends to save memory.
# To use MAST with such packages, you just need to upcast the SingleCellExperiment to MAST’s subclass SingleCellAssay with the function SceToSingleCellAssay that handles the coercion and checks the object for validity. 
# Going the other direction, generally SingleCellAssays should work in packages that use SingleCellExperiment, but if in doubt you could down-cast with as(sca, 'SingleCellExperiment').

fov_integrated <- readRDS("./fov_integrated.rds")
Idents(fov_integrated) <- 'Epithelial_cluster'

library(SingleCellExperiment)
library(scater)

# The only here consideration is making sure MAST can find log-like data, and possibly thresholding the data.
fov_sce <- as.SingleCellExperiment(fov_integrated, assay = "RNA")

sca = SceToSingleCellAssay(fov_sce)
saveRDS(sca, "./fov_inegrated_sca.rds")
sca <- readRDS("./fov_inegrated_sca.rds")

# Dimensions before subsetting
dim(sca)

# SKIP THIS - keep genes that are expressed in more than 2% of all cells
#sca2 <- sca[freq(sca)>0.02,]
#dim(sca2)

# add a column to the data which contains scaled number of genes that are expressed in each cell
# Recalculating the cellular detection rate (ngeneson)
# We and others have found that the number of genes detected in a sample, which we deemed the cellular detection rate is often the first principal component. 
# After we removed genes low expression, perhaps we want to recalculate it.
cdr2 <- colSums(assay(sca)>0)
colData(sca)$ngeneson <- scale(cdr2)
qplot(x=cdr2, y=colData(sca)$ngeneson) + xlab('New CDR') + ylab('Old CDR')

# PCA on filtered cells
set.seed(123)
plotPCA <- function(sca_obj){
  projection <- rpca(t(assay(sca_obj)), retx=TRUE, k=4)$x
  colnames(projection)=c("PC1","PC2","PC3","PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'ngeneson'),
                mapping=aes(color=Epithelial_cluster), upper=list(continuous='blank')))
  invisible(pca)
}

plotPCA(sca)

# Adaptive thresholding - SKIP THIS
# Single cell gene expression data are known to be zero-inflated and bimodal, which is feature we observe here as well.
scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x=logcounts))+geom_density() +facet_wrap(~primerid, scale='free_y')

# 3.3- Differential Expression using a Hurdle model - Epithelial clusters

# We’ll fit a hurdle model, modeling the condition and (centered) ngeneson factor, thus adjusting for the cellular detection rate.
# In order to have more interpretable coefficients, we’ll set the reference level of the factor to be the “unstimulated” cells.

# store the columns that we are interested in as factors - here the Epithelial clusters
# Comparing epithelial clusters 1 vs 3 (# 1 and 2 are similar, so no need to run 2 vs 3)
cond <-factor(colData(sca)$Epithelial_cluster)

# set the reference level
cond <- relevel(cond,"3")
colData(sca)$cond <- cond

celltype <- factor(colData(sca)$CIPR_ordered_clusters)
celltype <- make.names(celltype, unique=FALSE)
celltype <- factor(celltype)
colData(sca)$celltype <- celltype

# same for donors (which we need to model random effects)
replicate <- factor(colData(sca)$sample_ID)
colData(sca)$replicate <- replicate

# create a group per condition-celltype combination
colData(sca)$epithelial_group <- paste0(sca$celltype, sep = ".", sca$cond)
colData(sca)$epithelial_group <- factor(colData(sca)$epithelial_group)

# define and fit the model
zlmCond <- zlm(formula = ~ngeneson + epithelial_group + (1 | replicate), 
               sca=sca, 
               method='glmer', 
               ebayes=F, 
               strictConvergence=F,
               fitArgsD=list(nAGQ = 0)) # to speed up calculations

saveRDS(zlmCond, "./MAST_DGE_Epithelial_clusters.rds")

# perform likelihood-ratio test for the condition that we are interested in   
# SPP1..macrophage.1
# CD8..T.cell.1
# CD4..T.cell.1
# B.cell.1
# Plasma.cell
# Fibroblast

# Epithelial cluster 1 vs 3
summaryCond <- summary(zlmCond, doLRT='epithelial_groupSPP1..macrophage.1')
saveRDS(summaryCond, "./Epithelial_cluster1_SPP1_macrophage.rds")

summaryCond2 <- summary(zlmCond, doLRT='epithelial_groupCD8..T.cell.1')
saveRDS(summaryCond2, "./Epithelial_cluster1_CD8.rds")

summaryCond3 <- summary(zlmCond, doLRT='epithelial_groupCD4..T.cell.1')
saveRDS(summaryCond3, "./Epithelial_cluster1_CD4.rds")

summaryCond4 <- summary(zlmCond, doLRT='epithelial_groupB.cell.2')
saveRDS(summaryCond4, "./Epithelial_cluster1_Bcell.rds")

summaryCond5 <- summary(zlmCond, doLRT='epithelial_groupPlasma.cell.1')
saveRDS(summaryCond5, "./Epithelial_cluster1_Plasma.rds")

summaryCond6 <- summary(zlmCond, doLRT='epithelial_groupFibroblast.1')
saveRDS(summaryCond6, "./Epithelial_cluster1_Fibroblast.rds")

summaryCond7 <- summary(zlmCond, doLRT='epithelial_groupTumor.epithelial.cell.1')
saveRDS(summaryCond7, "./Epithelial_cluster1_Tumor_epithelial_cell.rds")

# Epithelial cluster 2 vs 3
summaryCond <- summary(zlmCond, doLRT='epithelial_groupSPP1..macrophage.2')
saveRDS(summaryCond, "./Epithelial_cluster2_SPP1_macrophage.rds")

summaryCond2 <- summary(zlmCond, doLRT='epithelial_groupCD8..T.cell.2')
saveRDS(summaryCond2, "./Epithelial_cluster2_CD8.rds")

summaryCond3 <- summary(zlmCond, doLRT='epithelial_groupCD4..T.cell.2')
saveRDS(summaryCond3, "./Epithelial_cluster2_CD4.rds")

summaryCond4 <- summary(zlmCond, doLRT='epithelial_groupB.cell.2')
saveRDS(summaryCond4, "./Epithelial_cluster2_Bcell.rds")

summaryCond5 <- summary(zlmCond, doLRT='epithelial_groupPlasma.cell.2')
saveRDS(summaryCond5, "./Epithelial_cluster2_Plasma.rds")

summaryCond6 <- summary(zlmCond, doLRT='epithelial_groupFibroblast.2')
saveRDS(summaryCond6, "./Epithelial_cluster2_Fibroblast.rds")

summaryCond7 <- summary(zlmCond, doLRT='epithelial_groupTumor.epithelial.cell.2')
saveRDS(summaryCond7, "./Epithelial_cluster2_Tumor_epithelial_cell.rds")

# get the table with log-fold changes and p-values
summaryDt <- summaryCond$datatable
summaryDt2 <- summaryCond2$datatable
summaryDt3 <- summaryCond3$datatable
summaryDt4 <- summaryCond4$datatable
summaryDt5 <- summaryCond5$datatable
summaryDt6 <- summaryCond6$datatable
summaryDt7 <- summaryCond7$datatable

# option 1
fcHurdle <- merge(summaryDt[contrast=='epithelial_groupSPP1..macrophage.1' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='epithelial_groupSPP1..macrophage.1' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)>FCTHRESHOLD], as.data.table(mcols(sca)), by='primerid')
setorder(fcHurdleSig, fdr)
 
#option 2 - use this to plot the volcano plots below

# SPP1 macrophages
result_macro <- merge(summaryDt[contrast=='epithelial_groupSPP1..macrophage.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                summaryDt[contrast=='epithelial_groupSPP1..macrophage.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_macro[,coef:=result_macro[,coef]/log(2)]
result_macro[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_macro = result_macro[result_macro$FDR<0.05,, drop=F]
result_macro <- stats::na.omit(as.data.frame(result_macro))
setorder(result_macro, FDR)
write.csv(result_macro, "DGE_MAST_Epithelial_cluster_2vs3_SPP1macrophage.csv")
result_macro <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_SPP1macrophage.csv")

# CD8 T cells
result_CD8 <- merge(summaryDt2[contrast=='epithelial_groupCD8..T.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                      summaryDt2[contrast=='epithelial_groupCD8..T.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_CD8[,coef:=result_CD8[,coef]/log(2)]
result_CD8[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_CD8 = result_CD8[result_CD8$FDR<0.05,, drop=F]
result_CD8<- stats::na.omit(as.data.frame(result_CD8))
setorder(result_CD8, FDR)
write.csv(result_CD8, "DGE_MAST_Epithelial_cluster_2vs3_CD8.csv")
result_CD8 <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_CD8.csv")

# CD4 T cells
result_CD4 <- merge(summaryDt3[contrast=='epithelial_groupCD4..T.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                    summaryDt3[contrast=='epithelial_groupCD4..T.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_CD4[,coef:=result_CD4[,coef]/log(2)]
result_CD4[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_CD4 = result_CD4[result_CD4$FDR<0.05,, drop=F]
result_CD4<- stats::na.omit(as.data.frame(result_CD4))
setorder(result_CD4, FDR)
write.csv(result_CD4, "DGE_MAST_Epithelial_cluster_2vs3_CD4.csv")
result_CD4 <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_CD4.csv")

# B cells
result_B <- merge(summaryDt4[contrast=='epithelial_groupB.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                    summaryDt4[contrast=='epithelial_groupB.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_B[,coef:=result_B[,coef]/log(2)]
result_B[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_B = result_B[result_B$FDR<0.05,, drop=F]
result_B <- stats::na.omit(as.data.frame(result_B))
setorder(result_B, FDR)
write.csv(result_B, "DGE_MAST_Epithelial_cluster_1vs3_Bcell.csv")
result_B <- read.csv("DGE_MAST_Epithelial_cluster_1vs3_Bcell.csv")

# Plasma cells
result_Plasma <- merge(summaryDt5[contrast=='epithelial_groupPlasma.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                  summaryDt5[contrast=='epithelial_groupPlasma.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_Plasma[,coef:=result_Plasma[,coef]/log(2)]
result_Plasma[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_Plasma = result_Plasma[result_Plasma$FDR<0.05,, drop=F]
result_Plasma <- stats::na.omit(as.data.frame(result_Plasma))
setorder(result_Plasma, FDR)
write.csv(result_Plasma, "DGE_MAST_Epithelial_cluster_2vs3_Plasma.csv")
result_Plasma <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_Plasma.csv")

# Fibroblasts
result_fibro <- merge(summaryDt6[contrast=='epithelial_groupFibroblast.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                       summaryDt6[contrast=='epithelial_groupFibroblast.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_fibro[,coef:=result_fibro[,coef]/log(2)]
result_fibro[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_fibro = result_fibro[result_fibro$FDR<0.05,, drop=F]
result_fibro <- stats::na.omit(as.data.frame(result_fibro))
setorder(result_fibro, FDR)
write.csv(result_fibro, "DGE_MAST_Epithelial_cluster_2vs3_Fibroblast.csv")
result_fibro <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_Fibroblast.csv")

# Tumor epithelial cell
result_tumor <- merge(summaryDt7[contrast=='epithelial_groupTumor.epithelial.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                      summaryDt7[contrast=='epithelial_groupTumor.epithelial.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_tumor[,coef:=result_tumor[,coef]/log(2)]
result_tumor[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_tumor = result_tumor[result_tumor$FDR<0.05,, drop=F]
result_tumor <- stats::na.omit(as.data.frame(result_tumor))
setorder(result_tumor, FDR)
write.csv(result_tumor, "DGE_MAST_Epithelial_cluster_2vs3_TEC.csv")
result_tumor <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_TEC.csv")

# 3.4- Visualization

# 3.4.1- Volcano plot
EnhancedVolcano(toptable = result_macro, x = "coef", y = "FDR",
                lab = result_macro$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1 macrophage Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_CD8, x = "coef", y = "FDR",
                lab = result_CD8$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 T cells Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_CD4, x = "coef", y = "FDR",
                lab = result_CD4$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 T cells Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_B, x = "coef", y = "FDR",
                lab = result_B$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE B cells Epithelial cluster 1 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_Plasma, x = "coef", y = "FDR",
                lab = result_Plasma$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Plasma cells Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_fibro, x = "coef", y = "FDR",
                lab = result_fibro$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Fibroblasts Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_tumor, x = "coef", y = "FDR",
                lab = result_tumor$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumor epithelial cell Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.4.2- Visualization of 50 most differentially expressed genes
entrez_to_plot <- fcHurdleSig[1:35,primerid]
symbols_to_plot <- fcHurdleSig[1:35,primerid]
flat_dat <- as(sca[entrez_to_plot,], 'data.table')

#thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 20, min_per_bin = 30)
#assays(sca, withDimnames = FALSE) <- list(thresh=thres$counts_threshold, tpm=assay(sca))

ggbase <- ggplot(flat_dat, aes(x=cond, y=counts, color=cond)) + geom_jitter()+facet_wrap(~primerid, scale='free_y')+ggtitle("DE Genes")
ggbase+geom_violin() 

flat_dat[,lmPred:=lm(counts~ngeneson + cond)$fitted, key=primerid]
ggbase +aes(x=ngeneson) + geom_line(aes(y=lmPred), lty=1) + xlab('Standardized Cellular Detection Rate')

# 3.4.3- Heatmap of most differentially expressed genes
mat_to_plot <- assay(sca[entrez_to_plot,])
rownames(mat_to_plot) <- symbols_to_plot
mat_to_plot2 <- as.matrix(mat_to_plot) #convert large dgCMatrix to a regular large matrix
aheatmap(mat_to_plot2, annCol=colData(sca)[,"Epithelial_cluster"],main="DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)))


# 3.5- Differential Expression using a Hurdle model - Stromal clusters

cond <-factor(colData(sca)$aSMA_cluster)

# set the reference level
cond <- relevel(cond,"2")
colData(sca)$cond <- cond

celltype <- factor(colData(sca)$CIPR_ordered_clusters)
celltype <- make.names(celltype, unique=FALSE)
celltype <- factor(celltype)
colData(sca)$celltype <- celltype

# same for donors (which we need to model random effects)
replicate <- factor(colData(sca)$sample_ID)
colData(sca)$replicate <- replicate

# create a group per condition-celltype combination
colData(sca)$stromal_group <- paste0(sca$celltype, sep = ".", sca$cond)
colData(sca)$stromal_group <- factor(colData(sca)$stromal_group)

# define and fit the model
zlmCond <- zlm(formula = ~ngeneson + stromal_group + (1 | replicate), 
               sca=sca, 
               method='glmer', 
               ebayes=F, 
               strictConvergence=F,
               fitArgsD=list(nAGQ = 0)) # to speed up calculations

saveRDS(zlmCond, "./MAST_DGE_Stromal_clusters.rds")

# perform likelihood-ratio test for the condition that we are interested in   
# SPP1..macrophage.1
# CD8..T.cell.1
# CD4..T.cell.1
# B.cell.1
# Plasma.cell
# Fibroblast

# Stromal cluster 1 vs 2
summaryCond <- summary(zlmCond, doLRT='stromal_groupSPP1..macrophage.1')
saveRDS(summaryCond, "./Stromal_cluster1_SPP1_macrophage.rds")

summaryCond2 <- summary(zlmCond, doLRT='stromal_groupCD8..T.cell.1')
saveRDS(summaryCond2, "./Stromal_cluster1_CD8.rds")

summaryCond3 <- summary(zlmCond, doLRT='stromal_groupCD4..T.cell.1')
saveRDS(summaryCond3, "./Stromal_cluster1_CD4.rds")

summaryCond4 <- summary(zlmCond, doLRT='stromal_groupB.cell.1')
saveRDS(summaryCond4, "./Stromal_cluster1_Bcell.rds")

summaryCond5 <- summary(zlmCond, doLRT='stromal_groupPlasma.cell.1')
saveRDS(summaryCond5, "./Stromal_cluster1_Plasma.rds")

summaryCond6 <- summary(zlmCond, doLRT='stromal_groupFibroblast.1')
saveRDS(summaryCond6, "./Stromal_cluster1_Fibroblast.rds")

summaryCond7 <- summary(zlmCond, doLRT='stromal_groupTumor.epithelial.cell.1')
saveRDS(summaryCond7, "./Stromal_cluster1_Tumor_epithelial_cell.rds")

# get the table with log-fold changes and p-values
summaryDt <- summaryCond$datatable
summaryDt2 <- summaryCond2$datatable
summaryDt3 <- summaryCond3$datatable
summaryDt4 <- summaryCond4$datatable
summaryDt5 <- summaryCond5$datatable
summaryDt6 <- summaryCond6$datatable
summaryDt7 <- summaryCond7$datatable

#option 2 - use this to plot the volcano plots below

# SPP1 macrophages
result_macro <- merge(summaryDt[contrast=='epithelial_groupSPP1..macrophage.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                      summaryDt[contrast=='epithelial_groupSPP1..macrophage.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_macro[,coef:=result_macro[,coef]/log(2)]
result_macro[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_macro = result_macro[result_macro$FDR<0.05,, drop=F]
result_macro <- stats::na.omit(as.data.frame(result_macro))
setorder(result_macro, FDR)
write.csv(result_macro, "DGE_MAST_Epithelial_cluster_2vs3_SPP1macrophage.csv")
result_macro <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_SPP1macrophage.csv")

# CD8 T cells
result_CD8 <- merge(summaryDt2[contrast=='epithelial_groupCD8..T.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                    summaryDt2[contrast=='epithelial_groupCD8..T.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_CD8[,coef:=result_CD8[,coef]/log(2)]
result_CD8[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_CD8 = result_CD8[result_CD8$FDR<0.05,, drop=F]
result_CD8<- stats::na.omit(as.data.frame(result_CD8))
setorder(result_CD8, FDR)
write.csv(result_CD8, "DGE_MAST_Epithelial_cluster_2vs3_CD8.csv")
result_CD8 <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_CD8.csv")

# CD4 T cells
result_CD4 <- merge(summaryDt3[contrast=='epithelial_groupCD4..T.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                    summaryDt3[contrast=='epithelial_groupCD4..T.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_CD4[,coef:=result_CD4[,coef]/log(2)]
result_CD4[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_CD4 = result_CD4[result_CD4$FDR<0.05,, drop=F]
result_CD4<- stats::na.omit(as.data.frame(result_CD4))
setorder(result_CD4, FDR)
write.csv(result_CD4, "DGE_MAST_Epithelial_cluster_2vs3_CD4.csv")
result_CD4 <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_CD4.csv")

# B cells
result_B <- merge(summaryDt4[contrast=='epithelial_groupB.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                  summaryDt4[contrast=='epithelial_groupB.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_B[,coef:=result_B[,coef]/log(2)]
result_B[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_B = result_B[result_B$FDR<0.05,, drop=F]
result_B <- stats::na.omit(as.data.frame(result_B))
setorder(result_B, FDR)
write.csv(result_B, "DGE_MAST_Epithelial_cluster_1vs3_Bcell.csv")
result_B <- read.csv("DGE_MAST_Epithelial_cluster_1vs3_Bcell.csv")

# Plasma cells
result_Plasma <- merge(summaryDt5[contrast=='epithelial_groupPlasma.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                       summaryDt5[contrast=='epithelial_groupPlasma.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_Plasma[,coef:=result_Plasma[,coef]/log(2)]
result_Plasma[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_Plasma = result_Plasma[result_Plasma$FDR<0.05,, drop=F]
result_Plasma <- stats::na.omit(as.data.frame(result_Plasma))
setorder(result_Plasma, FDR)
write.csv(result_Plasma, "DGE_MAST_Epithelial_cluster_2vs3_Plasma.csv")
result_Plasma <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_Plasma.csv")

# Fibroblasts
result_fibro <- merge(summaryDt6[contrast=='epithelial_groupFibroblast.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                      summaryDt6[contrast=='epithelial_groupFibroblast.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_fibro[,coef:=result_fibro[,coef]/log(2)]
result_fibro[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_fibro = result_fibro[result_fibro$FDR<0.05,, drop=F]
result_fibro <- stats::na.omit(as.data.frame(result_fibro))
setorder(result_fibro, FDR)
write.csv(result_fibro, "DGE_MAST_Epithelial_cluster_2vs3_Fibroblast.csv")
result_fibro <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_Fibroblast.csv")

# Tumor epithelial cell
result_tumor <- merge(summaryDt7[contrast=='epithelial_groupTumor.epithelial.cell.2' & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                      summaryDt7[contrast=='epithelial_groupTumor.epithelial.cell.2' & component=='logFC', .(primerid, coef)], by='primerid') # logFC coefficients
# MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
result_tumor[,coef:=result_tumor[,coef]/log(2)]
result_tumor[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')] # do multiple testing correction
result_tumor = result_tumor[result_tumor$FDR<0.05,, drop=F]
result_tumor <- stats::na.omit(as.data.frame(result_tumor))
setorder(result_tumor, FDR)
write.csv(result_tumor, "DGE_MAST_Epithelial_cluster_2vs3_TEC.csv")
result_tumor <- read.csv("DGE_MAST_Epithelial_cluster_2vs3_TEC.csv")

# 3.4- Visualization

# 3.4.1- Volcano plot
EnhancedVolcano(toptable = result_macro, x = "coef", y = "FDR",
                lab = result_macro$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1 macrophage Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_CD8, x = "coef", y = "FDR",
                lab = result_CD8$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 T cells Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_CD4, x = "coef", y = "FDR",
                lab = result_CD4$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 T cells Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_B, x = "coef", y = "FDR",
                lab = result_B$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE B cells Epithelial cluster 1 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_Plasma, x = "coef", y = "FDR",
                lab = result_Plasma$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Plasma cells Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_fibro, x = "coef", y = "FDR",
                lab = result_fibro$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Fibroblasts Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

EnhancedVolcano(toptable = result_tumor, x = "coef", y = "FDR",
                lab = result_tumor$primerid, pCutoff = 0.05, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumor epithelial cell Epithelial cluster 2 vs cluster 3", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()


# 4- Gene Set Enrichment Analysis

# bootstrap, resampling cells
# R should be set to >50 if you were doing this for real.
boots <- bootVcov1(zlmCond, R = 100)

module <- "BTM"
min_gene_in_module <- ?
packageExt <- system.file("extdata", package='MAST')
module_file <- list.files(packageExt, pattern = module, full.names = TRUE)
gene_set <- getGmt(module_file)
gene_ids <- geneIds(gene_set)
gene_ids <- gene_ids[!names(gene_ids)%like%"TBA"&!names(gene_ids)%like%"Macrophage"]
sets_indices <- limma::ids2indices(gene_ids, mcols(sca)$symbolid)
# Only keep modules with at least min_gene_in_module
sets_indices <- sets_indices[sapply(sets_indices, length) >= min_gene_in_module]





















