library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(EnhancedVolcano)

cts <- read.table("BM_PB_human_d0.txt", header=TRUE, row.names=1)
coldata <- read.table("BM_PB_coldata.txt", header=TRUE, row.names=1)

#1- Normalisation
#Check columns of the count matrix and the rows of the column data (information about samples) are in the same order. 
all(rownames(coldata) == colnames(cts))

# Creation of the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds

#Pre-filtering low count genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

#Differential expression analysis of BM_d0 vs PB_d0
dds <- DESeq(dds)
dds

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

res <- results(dds, contrast=c("condition","BM_d0","APB_d0"))
res
mcols(res)

resNoFilt <- results(dds, contrast=c("condition","BM_d42","ABM_d0"), independentFiltering=FALSE)
resNoFilt

## Getting better gene names

res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
#The following chunk of code uses the ENSEMBL mart, querying with the ENSEMBL transcript id and
#requesting the Entrez gene id and HGNC gene symbol.
library( "biomaRt" )
# this setting is needed, as it Rstudio does not pick up the proxy server of the university
#Sys.setenv("http_proxy" = "http://wwwcache.gla.ac.uk:8080")
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_transcript_id","ensembl_gene_id", "entrezgene_id","hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res$ensembl,
                  mart = ensembl )
idx <- match( res$ensembl, genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene_id[ idx ]
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res,20)

# save to disk
save(dds, file="dds.RData")
save(res, file="res.RData")

#To retrieve the normalized counts matrix from dds
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#2- Quality control
#Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="condition")
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 5:6], col="#00000020", pch=20, cex=0.3 )
cor(assay(rld)[, 1:6],method = "pearson")

sampleDists <- dist( t( assay(rld) ) ) # generated the distance
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Sample )
colnames(sampleDistMatrix) <- paste( rld$Sample )
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

### Extract the rlog matrix from the object
rld_mat <- assay(rld)    
### Compute pairwise correlation values
rld_cor <- cor(rld_mat)    ## cor() is a base R function
## check the output of cor(), make note of the row names and column names
head(rld_cor)

head(meta)
#Load pheatmap package
library(pheatmap)

### Plot heatmap using the correlation matrix and the metadata object
pheatmap(rld_cor, annotation = coldata)
#changing color
heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
pheatmap(rld_cor, annotation = coldata, color = heat.colors, border_color=NA, fontsize = 10, 
         fontsize_row = 10, height=20)

#3- DGE analysis itself
#Log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, res = res, coef="condition_BM_d0_vs_APB_d0", type="apeglm")
resLFC

## Getting better gene names
resLFC$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
#The following chunk of code uses the ENSEMBL mart, querying with the ENSEMBL transcript id and
#requesting the Entrez gene id and HGNC gene symbol.
library( "biomaRt" )
# this setting is needed, as it Rstudio does not pick up the proxy server of the university
#Sys.setenv("http_proxy" = "http://wwwcache.gla.ac.uk:8080")
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_transcript_id","ensembl_gene_id", "entrezgene_id","hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = resLFC$ensembl,
                  mart = ensembl )
idx <- match( resLFC$ensembl, genemap$ensembl_gene_id )
resLFC$entrez <- genemap$entrezgene_id[ idx ]
resLFC$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(resLFC,20)

#4- Summarizing results
#We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]
#We can summarize some basic tallies using the summary function.
summary(res)
#How many adjusted p-values were less than 0.01 or 0.05?
sum(res$padj < 0.01, na.rm=TRUE)


#Multiple hypothesis test correction.
library(dplyr)
# threshold of p = 0.01
res %>% 
  as.data.frame() %>% 
  filter(padj < 0.05) %>% 
  dim()

# threshold of p = 0.001
res %>% 
  as.data.frame() %>% 
  filter(padj < 0.01) %>% 
  dim()

# distribution of adjusted p-values
hist(res$padj, col="lightblue", main = "Adjusted p-value distribution")
# distribution of non-adjusted p-values
hist(res$pvalue, col="grey", main = "Non-adjusted p-value distribution")


#To obtain the list of significant DE genes:
diff = res %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(padj <= 0.01) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))
head(diff)
write.table(diff, file="deg_signpadj001.txt", sep="\t", quote=F, col.names=NA)

#To obtain the list of genes based on log2FC and padj values:
diff = res %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(abs(log2FoldChange) >= 1 & padj <= 0.05) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))
head(diff)
write.table(diff, file="SignificantDEG_log2FC1_padj005.txt", sep="\t", quote=F, col.names=NA)

#Summary for res without independent filtering
diff = resNoFilt %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(abs(log2FoldChange) >= 2) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))
head(diff)

#To obtain the most differentially up-regulated expressed genes do:
resSig <- res[which(res$padj < 0.05 ), ]
head( resSig[order(resSig$log2FoldChange), ] )
resSig <- res[which(res$padj < 0.05 & res$log2FoldChange>1),]
resSig
nrow(resSig)
rownames(resSig)

#5- Volcano Plots
# load the library if not done yet
library(EnhancedVolcano)

#Plotting withouth shirinking res
# The main function is named after the package
EnhancedVolcano(toptable = res,              # We use the shrunken log2 fold change as noise associated with low count genes is removed 
                x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
                y = "pvalue",                     # Name of the column in resLFC that contains the p-value
                lab = rownames(resLFC)
)
EnhancedVolcano(toptable = res,
                x = "log2FoldChange",
                y = "padj",
                lab = res$hgnc_symbol,
                xlim = c(-10, +10),
                ylim = c(0,7),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                col = c('black', 'blue', 'green', 'red'),
                title = "BM d0 vs PB d0\n (fold change cutoff = 1, padj cutoff = 0.05)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & Log2 fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black'
)

#selectLab = c('TLR3', 'TLR7', 'TLR8','ANGPT1', 'ANGPT2','VCAM1', 'ICAM1', 'SELE', 'TNF', 'IL1A', 'CXCL8', 'IFNA1','IFNA2',
             # 'IFNB1', 'IL7', 'IL7R','IL21', 'CXCR4','SPP1', 'CD34', 'TFRC'),

#Plotting after shirinking res
EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = resLFC$hgnc_symbol,
                xlim = c(-10, +10),
                ylim = c(0,7),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                col = c('black', 'blue', 'green', 'red'),
                title = "BM d0 vs PB d0\n (fold change cutoff = 2, padj cutoff = 0.01)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & Log2 fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = FALSE,
                widthConnectors = 0.15,
                colConnectors = ('black')
)


#6- HEATMAP of normalized counts of DEG 
contrast_oe=c("condition", "BM_d0", "APB_d0")
res_tableOE <- results(dds, contrast=contrast_oe, alpha = 0.05)
res_tableOE %>% 
  data.frame() %>% 
  View()
# Get information on each column in results
mcols(res_tableOE, use.names=T)

#Gene-level filtering
# Filter genes by zero expression
res_tableOE[which(res_tableOE$baseMean == 0),] %>% 
  data.frame() %>% 
  View()
# Filter genes that have an extreme outlier
res_tableOE[which(is.na(res_tableOE$pvalue) & 
                    is.na(res_tableOE$padj) &
                    res_tableOE$baseMean > 0),] %>% 
  data.frame() %>% 
  View()
# Filter genes below the low mean threshold
res_tableOE[which(!is.na(res_tableOE$pvalue) & 
                    is.na(res_tableOE$padj) & 
                    res_tableOE$baseMean > 0),] %>% 
  data.frame() %>% 
  View()
# Apply fold change shrinkage
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, type = "normal")
# MA plot using shrunken fold changes
plotMA(res_tableOE, ylim=c(-2,2))
# Create a tibble of results
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
# Subset the tibble to keep only significant genes
sigOE <- res_tableOE_tb %>%
  filter(padj < 0.01)
# Take a quick look at this tibble
sigOE
# DESeq2 creates a matrix when you use the counts() function
## First convert normalized_counts to a data frame and transfer the row names to a new column called "gene"
normalized_counts <- counts(dds, normalized=T) %>% 
  data.frame() %>%
  rownames_to_column(var="gene")
### Extract normalized expression for significant genes from the OE and control samples (2:4 and 7:9)
norm_OEsig <- normalized_counts %>% 
  filter(gene %in% sigOE$gene)  
rownames(norm_OEsig) <- norm_OEsig[,1]
norm_OEsig <- norm_OEsig[,-1]
### Set a color palette
heat_colors <- brewer.pal(6, "PuOr")

### Run pheatmap using the metadata data frame for the annotation
pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_rows = T,
         cluster_cols = F,
         show_rownames = F,
         annotation = coldata, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

#7- Clustering DEG
# Obtain rlog values for those significant genes
rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)  
cluster_rlog <- rld_mat[diff$gene, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
#metadata: the metadata dataframe that corresponds to samples
#time: character column name in metadata that will be used as variable that changes
#col: character column name in metadata to separate samples

meta = data.frame(coldata, row.names = colnames(cts))
clusters <- degPatterns(cluster_rlog, metadata = coldata, col = "col")

class(clusters)
names(clusters)
head(clusters&df)

# Extract the Group 1 genes
cluster_groups <- clusters$df
group1 <- clusters$df %>%
  filter(cluster == 1)

cluster_groups <- clusters$df
group2 <- clusters$df %>%
  filter(cluster == 2)

#After extracting a group of genes, we can use annotation packages to obtain additional information. 
#We can also use these lists of genes as input to downstream functional analysis tools to obtain more biological insight 
#and see whether the groups of genes share a specific function.

group1$genes <- sapply( strsplit( rownames(group1), split="\\+" ), "[", 1 )
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_transcript_id","ensembl_gene_id", "entrezgene_id","hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = group1$genes,
                  mart = ensembl )
idx <- match( group1$genes , genemap$ensembl_gene_id )
group1$genes <- genemap$hgnc_symbol[ idx ]
write.table(group1, file="cluster_genes_1.txt", sep="\t", quote=F, col.names=NA)



group2$genes <- sapply( strsplit( rownames(group2), split="\\+" ), "[", 1 )
genemap <- getBM( attributes = c("ensembl_transcript_id","ensembl_gene_id", "entrezgene_id","hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = group2$genes,
                  mart = ensembl )
idx <- match( group2$genes , genemap$ensembl_gene_id )
group2$genes <- genemap$hgnc_symbol[ idx ]
write.table(group2, file="cluster_genes_2.txt", sep="\t", quote=F, col.names=NA)
