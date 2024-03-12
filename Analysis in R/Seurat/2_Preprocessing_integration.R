
setwd("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/COSMX_colon_dataset/COSMX_Colonrectal_cancer_project/TMA_original/Seurat")

library(Seurat)
library(SeuratWrappers)
library(SeuratData)
library(Azimuth)
library(data.table)
library(dplyr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggplot2)
library(dittoSeq)
library(SeuratDisk)
library(viridis)
library(RColorBrewer)
library(DoubletFinder)
library(DropletUtils)
library(dsb)
library(harmony)

fov_merged <- readRDS("./RDS_files/3_Merged/fov_merged.rds")
saveRDS(fov_merged, "fov_merged.rds")

saveRDS(fov_subset, "fov_subset.rds")
fov_subset <- readRDS("./RDS_files/4_Subset/fov_subset.rds")

fov_integrated <- readRDS("./RDS_files/5_Integration/All/fov_integrated.rds")
saveRDS(fov_integrated, "fov_integrated.rds")

# 1- Standard pre-processing workflow

# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
# These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.

# 1.1- Calculate QC

# RENAME UMAP TO RNA counts
colnames(fov_merged@UMI)[index.of.column] <- "new.name"

fov_merged$RNA_counts <- fov_merged$UMI
fov_merged@meta.data[["UMI"]] <- NULL

# Visualize QC metrics as a violin plot before filtering
Idents(fov_merged) <- 'sample_ID'
levels(fov_merged)

fov_merged[["percent.neg"]] <- PercentageFeatureSet(fov_merged, pattern = "^NegPrb") 
fov_merged[["percent.mt"]] <- PercentageFeatureSet(fov_merged, pattern = "^MT") 
fov_merged[["percent.ribo"]] <- PercentageFeatureSet(fov_merged, pattern = "^RPL") 
fov_merged[["percent.hb"]] <- PercentageFeatureSet(fov_merged, pattern = "^HBB") 

feats <- c("Genes", "nCount_RNA", "percent.neg", "percent.mt", "percent.ribo", "percent.hb")
VlnPlot(fov_merged, features = feats, pt.size = 0.1, ncol = 2) +
  NoLegend()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# plot1 <- FeatureScatter(fov_merged, feature1 = "UMI", feature2 = "percent.mt")
plot2 <- FeatureScatter(fov_merged, feature1 = "nCount_RNA", feature2 = "Genes")
plot2

# Convert fov_merged to anndata object to plot QC features before filtering
library(reticulate)
use_condaenv("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate")
use_python("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/bin/python")
library(sceasy)
sceasy::convertFormat(fov_merged, from="seurat", to="anndata",
                      outFile='adata_fov_merged.h5ad')

# 1.2- Filtering

# 1.2.1- Detection-based filtering and removal of cells with high % of mitochondrial genes (optional: high % of ribossomal genes)

# A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. 
# Here we will only consider cells with at least 5 detected genes. 
# Please note that those values are highly dependent on the library preparation method used.
fov_subset <- subset(fov_merged, subset = nFeature_RNA > 5 & percent.neg < 10)
head(fov_subset)
dim(fov_subset) # Gives you the numbers of genes (first number) and cells (second number) in your seurat object
head(fov_subset) 
#saveRDS(fov_subset, "fov_subset_neg_probes.rds")

unique(sapply(X = strsplit(colnames(fov_subset), split = "_"), FUN = "[", 1))
table(fov_subset$orig.ident)

saveRDS(fov_merged, "fov_merged.rds")
#fov_merged <- readRDS("./RDS files/fov_merged.rds")

saveRDS(fov_subset, "fov_subset.rds")
fov_subset <- readRDS("./fov_subset.rds")

# 1.2.2- Extract the count matrix with the negative probe expression and calculate cells background
counts <- LayerData(fov_subset, assay = "RNA")
saveRDS(counts, "total_counts.rds")
neg <- counts[(which(rownames(counts) %in% c('NegPrb3',	'NegPrb5', 'NegPrb6',	'NegPrb7', 'NegPrb8',	'NegPrb9', 'NegPrb10',	'NegPrb11',	'NegPrb12',	'NegPrb13',	'NegPrb14',	'NegPrb15',	'NegPrb16', 'NegPrb18',	'NegPrb19',	'NegPrb20',	'NegPrb21',	'NegPrb22',	'NegPrb23'))),]
saveRDS(neg, "negative_probes_counts.rds")

# cells should be in rows and genes in columns
counts <- as.matrix(counts)
str(counts)
counts <- t(counts)
str(counts)

neg <- as.matrix(neg)
neg <- t(neg)
str(neg)

# Vector giving each cell’s mean negative control value:
negmean <- Matrix::rowMeans(neg)
head(negmean)
negmean2 <- as.matrix(negmean)

# estimate per-cell bg as a fraction of total counts:
negmean.per.totcount <- mean(negmean) / mean(rowSums(counts))
per.cell.bg <- rowSums(counts) * negmean.per.totcount
per.cell.bg2 <- as.matrix(per.cell.bg)

fov_subset$meanneg_counts <- negmean2
fov_subset$per_cell_bg <- per.cell.bg2

VlnPlot(fov_subset, features = c("meanneg_counts", "per_cell_bg"), pt.size = 0.1, ncol = 1) +
  NoLegend()

saveRDS(fov_subset, "fov_subset.rds")

sceasy::convertFormat(fov_subset, from="seurat", to="anndata",
                      outFile='adata_fov_subset.h5ad')

# 1.2.3- Filter genes

# First see which genes contribute the most to such reads. We can for instance plot the percentage of counts per gene.
# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way will take a very long time.
par(mar = c(4, 8, 2, 1))
C <- fov_subset@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE) #As the level of expression of mitochondrial and MALAT1 genes are judged as mainly technical, it can be wise to remove them from the dataset bofore any further analysis.

# Filter MALAT1 and negative probe counts - careful this removes the protein assay slot
counts <- GetAssayData(fov_subset, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c('MALAT1', 'NegPrb3',	'NegPrb5', 'NegPrb6',	'NegPrb7', 'NegPrb8',	'NegPrb9', 'NegPrb10',	'NegPrb11',	'NegPrb12',	'NegPrb13',	'NegPrb14',	'NegPrb15',	'NegPrb16', 'NegPrb18',	'NegPrb19',	'NegPrb20',	'NegPrb21',	'NegPrb22',	'NegPrb23'))),]

# Re add the protein assay data
counts2 <- GetAssayData(fov_subset, assay = "Protein") 
counts2 <- counts2[-(which(rownames(counts2) %in% c('panck-positive'))),]
fov_subset <- subset(fov_subset, features = rownames(counts))
ptn_assay <- CreateAssayObject(counts = counts2)
# add this assay to the previously created Seurat object
fov_subset[["Protein"]] <- ptn_assay

# Filter Mitocondrial
# BM_subset2 <- BM_subset2[!grepl("^MT-", rownames(BM_subset)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
# <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
#data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]
dim(fov_subset)

# Plot filtered QC metrics

feats <- c("Genes", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb")
VlnPlot(fov_subset, features = feats, pt.size = 0.1, ncol = 2) + NoLegend()

FeatureScatter(fov_subset, feature1 = "nCount_RNA", feature2 = "Genes")

# Most expressed gebnes after filtering MALAT1
par(mar = c(4, 8, 2, 1))
C <- fov_subset@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE) #As the level of expression of mitochondrial and MALAT1 genes are judged as mainly technical, it can be wise to remove them from the dataset bofore any further analysis.

saveRDS(fov_subset, "fov_subset.rds")

# 1.2.3- Sample sex

genes.file = "./genes.table.csv"

if (!file.exists(genes.file)) {
  suppressMessages(require(biomaRt))
  
  # initialize connection to mart, may take some time if the sites are
  # unresponsive.
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  
  # fetch chromosome info plus some other annotations
  genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                   "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
                                    useCache = F))
  
  if (!dir.exists("data/results")) {
    dir.create("data/results")
  }
  if (is.data.frame(genes.table)) {
    write.csv(genes.table, file = genes.file)
  }
  
  if (!file.exists(genes.file)) {
    download.file("https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/misc/genes.table.csv",
                  destfile = "data/results/genes.table.csv")
    genes.table = read.csv(genes.file)
  }
  
} else {
  genes.table = read.csv(genes.file)
}

genes.table <- genes.table[genes.table$external_gene_name %in% rownames(fov_subset),
]

# Now that we have the chromosome information, we can calculate per cell the proportion of reads that comes from chromosome Y.

chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
fov_subset$pct_chrY = colSums(fov_subset@assays$RNA@counts[chrY.gene, ])/colSums(fov_subset@assays$RNA@counts)
# Then plot XIST expression vs chrY proportion. 
# As you can see, the samples are clearly on either side, even if some cells do not have detection of either.
FeatureScatter(BM_subset, feature1 = "XIST", feature2 = "pct_chrY")
VlnPlot(BM_subset, features = c("XIST", "pct_chrY"))


########################################################################################################################################################################################################################################################################################################################

# 2- Normalize the data

# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in BM1_subset[["RNA"]]@data.

fov_subset <- NormalizeData(fov_subset, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

# 2.1- Identification of highly variable features (feature selection)

# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
# We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
# By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

fov_subset <- FindVariableFeatures(fov_subset, selection.method = "vst", nfeatures = 500)

# Identify the 10 most highly variable genes
top500 <- head(VariableFeatures(fov_subset), 500)

# plot variable features with and without labelsoptions(ggrepel.max.overlaps = 50)
options(ggrepel.max.overlaps = 100)
display.brewer.all(n=10, exact.n=FALSE)

plot1 <- VariableFeaturePlot(fov_subset)
plot2 <- LabelPoints(plot = plot1, points = top500, repel = TRUE)
plot1
plot2

# 2.2- Scaling the data

# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(fov_subset)
fov_subset <- ScaleData(fov_subset, features = all.genes)

# 3- SCTransform

# Different algorithm to normalise and scale the data, the results are stored in a different assay, so our previous normalisation is not overwritten
# Biological heterogeneity in single-cell RNA-seq data is often confounded by technical factors including sequencing depth.
# This procedure omits the need for heuristic steps including pseudocount addition or log-transformation and improves common downstream analytical tasks 
# such as variable gene selection, dimensional reduction, and differential expression.
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage
# The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure. It can be invoked by specifying method="glmGamPoi".

library(sctransform)
fov_subset <- SCTransform(fov_subset, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

saveRDS(fov_subset, "fov_subset.rds")
fov_subset <- readRDS("./fov_subset.rds")

sceasy::convertFormat(fov_merged, from="seurat", to="anndata",
                      outFile='adata_fov_merged.h5ad')

# 3.1- Perform linear dimensional reduction

# Calculate PCs using variable features determined by SCTransform (3000 by default)
set.seed(22022)
fov_subset <- RunPCA(fov_subset, assay = "SCT", npcs = 50)

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
VizDimLoadings(fov_subset, dims = 1:2, reduction = "pca")
DimPlot(fov_subset, reduction = "pca")

# In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and features are ordered according to their PCA scores. 
# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. 
# Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.
DimHeatmap(fov_subset, dims = 1, cells = 1000, balanced = TRUE)
DimHeatmap(fov_subset, dims = 1:15, cells = 500, balanced = TRUE)

# 3.2- Determine the ‘dimensionality’ of the dataset
ElbowPlot(fov_subset, ndims = 50)

# 3.3- Run non-linear dimensional reduction (UMAP/tSNE) before integration

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
set.seed(22022)
fov_subset <- RunUMAP(fov_subset, dims = 1:50, reduction = "pca", reduction.name = "umap")

# 3.4- Visualise un-corrected (before integration) metadata in UMAP

DimPlot(fov_subset, reduction = "umap", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(fov_subset, var = "sample_ID",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "CMS",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "KM_Score",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "MMR_Status_Jen",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "IF_type",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "IF_type", split.by = "sample_ID",
             reduction.use = "umap", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_subset, var = "sample_ID", split.by = "sample_ID",
             reduction.use = "umap", size = 0.75,
             do.label = FALSE)

saveRDS(fov_subset, "fov_subset.rds")

sceasy::convertFormat(fov_subset, from="seurat", to="anndata",
                      outFile='adata_fov_subset.h5ad')

########################################################################################################################################################################################################################################################################################

# 4- Harmony integration
DefaultAssay(fov_subset) <- "SCT" 

# 4.1- Run Harmony
set.seed(22022)
fov_integrated <- RunHarmony(fov_subset, 
                            group.by.vars = c("sample_ID"), 
                            reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

# 4.2- Run UMAP with the corrected Harmony reduced dimensions
set.seed(22022)
fov_integrated <- RunUMAP(fov_integrated, reduction = "harmony", assay = "SCT", dims = 1:50, reduction.name = "umap.harmony")

# 4.3- Visualise UMAPs of the Harmony corrected (integrated) data
#options(ggrepel.max.overlaps = Inf)
DimPlot(fov_integrated, reduction = "umap.harmony", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(fov_integrated, var = "sample_ID",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_integrated, var = "IF_type",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell types based on IF on UMAP")

dittoDimPlot(fov_integrated, var = "IF_type", split.by = "sample_ID",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_integrated, var = "sample_ID", split.by = "sample_ID",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = FALSE)

# 4.4- Visualise epithelial and stromal canonical genes in the harmony corrected UMAP

DefaultAssay(fov_integrated) <- 'RNA'

features <- c("KRT8", "KRT18", "PTPRC", "CD3E")

FeaturePlot(fov_integrated, features = features, ncol = 2, reduction = "umap.harmony") # run UMAP again when changing the correction embedding

# 4.5- Visualise epithelial and stromal canonical protein markers and IF type in the harmony-corrected UMAP

DefaultAssay(fov_integrated) <- 'Protein'

features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

FeaturePlot(fov_integrated, features = features, ncol = 2, reduction = "umap.harmony") # run UMAP again when changing the correction embedding

saveRDS(fov_integrated, "fov_integrated.rds")
fov_integrated <- readRDS("./RDS_files/5_Integration/All/fov_integrated.rds")

########################################################################################################################################################################################################################################################################################

# 5- fastMNN integration 

# 5.1- Test fastMNN integration with the new function and compare outputs with the previous runs - (install version 5 to run functions below, then reinstall v4.3.0)
# Perform streamlined (one-line) integrative analysis
# Reference: https://satijalab.org/seurat/articles/seurat5_integration.html
# Each of these methods performs integration in low-dimensional space, and returns a dimensional reduction that aims to co-embed shared cell types across batches.
DefaultAssay(fov_integrated) <- 'SCT'

# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
set.seed(22022)
fov_integrated <- IntegrateLayers(
  object = fov_integrated, method = FastMNNIntegration, batch = fov_integrated$sample_ID,
  orig.reduction = "pca", new.reduction = "integrated.mnn", assay = "SCT", auto.merge = TRUE, #changed to SCT as previous version was run based on logcounts of RNA slot.
  verbose = FALSE
)

set.seed(22022)
fov_integrated <- RunUMAP(fov_integrated, reduction = "integrated.mnn", assay = "SCT", dims = 1:50, reduction.name = "umap.mnn")

DimPlot(fov_integrated, reduction = "umap.mnn", label = FALSE)

dittoDimPlot(fov_integrated, var = "sample_ID",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_integrated, var = "IF_type",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_integrated, var = "IF_type", split.by = "sample_ID",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_integrated, var = "sample_ID", split.by = "sample_ID",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE)

# Visualise epithelial and stromal canonical genes in the harmony corrected UMAP

DefaultAssay(fov_integrated) <- 'RNA'

features <- c("KRT8", "KRT18", "PTPRC", "CD3E")

FeaturePlot(fov_integrated, features = features, ncol = 2, reduction = "umap.mnn") # run UMAP again when changing the correction embedding

# Visualise epithelial and stromal canonical protein markers and IF type in the harmony-corrected UMAP

DefaultAssay(fov_integrated) <- 'Protein'

features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

FeaturePlot(fov_integrated, features = features, ncol = 2, reduction = "umap.mnn") # run UMAP again when changing the correction embedding


# 5.2- Test scVI integration with the new function and compare outputs with the previous runs

# Reference: https://docs.scvi-tools.org/en/stable/installation.html
# Reference: https://github.com/satijalab/seurat/issues/7164

library(reticulate)
#py_install("scvi-tools", envname = "/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/", method="conda", python_version=3.9, pip=TRUE)
#py_install("scanpy", envname = "/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/", method="auto", python_version=3.9, pip=TRUE)
use_condaenv(condaenv="/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/", required=FALSE)
scvi <- import("scvi")
scanpy <- import("scanpy")

DefaultAssay(fov_integrated) <- "RNA"
options(Seurat.object.assay.version = "v5")
# create a new seurat object, split by sample_ID and convert to seurat v5 assay
fov_integrated2 <- fov_integrated
fov_integrated2[["RNA"]] <- split(fov_integrated2[["RNA"]], f = fov_integrated2$sample_ID)
fov_integrated2

# Run scVI integration
fov_integrated2 <- IntegrateLayers(
  object = fov_integrated2, orig.reduction = "pca", method = scVIIntegration, assay = "RNA",
  new.reduction = "integrated.scvi",
  conda_env = "/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/", verbose = FALSE
)

########################################################################################################################################################################################################################################################################################

# 6- Apply BBKNN batch effect correction method directly in the seurat object

# Reference: https://cran.r-project.org/web/packages/bbknnR/vignettes/bbknnR-tutorial.html
library(bbknnR)

DefaultAssay(fov_integrated) <- "SCT" 
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

set.seed(220228)
fov_integrated <- RunBBKNN(fov_integrated, batch_key = "sample_ID", reduction = "pca", n_pcs = 50, graph_name = "BBKNN",
                           run_TSNE = FALSE, run_UMAP = TRUE,
                           UMAP_name = "umap.bbknn",
                           UMAP_key = "UMAP_")

saveRDS(fov_integrated, "fov_integrated.rds")

library(cowplot)
library(dittoSeq)
library(viridis)

# visualize sample id 
dittoDimPlot(fov_integrated, var = "sample_ID",
             reduction.use = "umap.bbknn", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_integrated, var = "IF_type",
             reduction.use = "umap.bbknn", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell types based on IF on UMAP")

dittoDimPlot(fov_integrated, var = "IF_type", split.by = "sample_ID",
             reduction.use = "umap.bbknn", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_integrated, var = "sample_ID", split.by = "sample_ID",
             reduction.use = "umap.bbknn", size = 0.75,
             do.label = FALSE)

# Visualise epithelial and stromal canonical genes in the BBKNN-corrected UMAP

DefaultAssay(fov_integrated) <- 'RNA'

features <- c("KRT8", "KRT18", "PTPRC", "CD3E")

FeaturePlot(fov_integrated, features = features, ncol = 2, reduction = "umap.bbknn") # run UMAP again when changing the correction embedding

# Visualise epithelial and stromal canonical protein markers and IF type in the harmony and fastMNN-corrected UMAP

DefaultAssay(fov_integrated) <- 'Protein'

features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

FeaturePlot(fov_integrated, features = features, ncol = 2, reduction = "umap.bbknn") # run UMAP again when changing the correction embedding


########################################################################################################################################################################################################################################################################################

# 7- Normalisation of protein data from IF using CLR from Seurat

# The resulting values can be interpreted as either a natural log ratio of the count for a given protein relative to the other proteins in the cell (CLR “across proteins”, as implemented in the original report of CITEseq 2 ) 
# or relative to other cells (CLR “across cells”, a modiﬁcation used in later work by the authors 19 , which renders CLR less dependent on the composition of the antibody panel). 
# The CLR transformation helps to better separate cell populations, but it does not directly estimate and correct for speciﬁc sources of technical noise including the apparent background noise.
# keep the raw IF intensity data as well.

fov_integrated <- readRDS("./RDS_files/5_Integration/All/fov_integrated.rds")

DefaultAssay(fov_integrated) <- "Protein"
fov_integrated <- NormalizeData(fov_integrated, normalization.method = "CLR", margin = 2) # If performing CLR normalization, normalize across features (1) or cells (2)
fov_integrated@assays[["Protein"]]@data <- as(fov_integrated@assays[["Protein"]]@data, "dgCMatrix")

saveRDS(fov_integrated, "fov_integrated.rds")

p1 <- FeaturePlot(fov_integrated, "MeanPanCK", slot = "data", cols = c("lightgrey", "darkgreen")) + ggtitle("MeanPanCK")
p1

# save.image(file = "COSMX_Seurat.RData")

########################################################################################################################################################################################################################################################################################

# 8- Run Scanpy in R to run BBKNN integration

fov_integrated <- readRDS("./RDS_files/5_Integration/Harmony/fov_integrated.rds")

# Reference: https://theislab.github.io/scanpy-in-R/

# 8.1-  Creating the R environment
renv::init()

# 8.2- Creating the Python environment
renv::install("reticulate")
renv::use_python()

# 8.3- Installing R packages
pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)

bioc_pkgs <- c(
  "SingleCellExperiment",
  "scater",
  "multtest"
)

# If you are using an {renv} environment
renv::install(pkgs)
renv::install(bioc_pkgs)

# 8.4- Installing Python packages
py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)
reticulate::py_install(py_pkgs)
reticulate::py_install('numpy==1.6.1')

# 8.5- Snapchot

# If you have chosen to use an {renv} environment you should run the renv::snapshot() command after installing packages. 
# This records the changes you have made so they can be restored later if necessary.
renv::snapshot()

# 8.6- The {reticulate} approach

# The {reticulate} approach involves calling scanpy functions from an R session. 
# This means that we keep everything inside a single R session and only need to write R code but the downside is we cannot directly copy examples that have been written in Python.
# We first need to load the R packages we need.
suppressPackageStartupMessages({
  library("reticulate")
  library("ggplot2")
  library("SingleCellExperiment")
  library("scater")
  library("Seurat")
})

# This time to load a Python package we don’t use a Python chunk but we use the reticulate::import() function.
sc <- import("scanpy")
np <- reticulate::import("numpy")
print(np$version$full_version)

config <- py_config()
config$numpy

# 8.6.1- Creating AnnData from Seurat

adata_seurat <- sc$AnnData(
  X   = t(GetAssayData(fov_integrated)),
  obs = fov_integrated[[]],
  var = GetAssay(fov_integrated)[[]]
)
#adata_seurat$obsm$update(umap = Embeddings(fov_integrated, "umap")) #skip this...
adata_seurat

sc$AnnData$write_h5ad(adata_seurat, 'adata_fov_integrated.h5ad')

# 8.6.2- BBKNN correction - run it in python with the saved anndata object created above, then load it back here to add the BBKNN-corrected reduced dim in the seurat object

#batch_correction_obs = 'sample_ID'
# Calculate PCA, this must be done before BBKNN
# sc$tl$pca(adata_seurat, n_comps=957)
# BBKNN - it is used in place of the scanpy 'neighbors' command that calculates nearest neighbours in the feature space
# py_pkgs <- c( "bbknn")
# reticulate::py_install(py_pkgs)# it does not work cuz there is an error when installing bbknn with pip install
# conda_install('r-reticulate',"bbknn") #it does not work
# sc$external$pp$bbknn(adata_seurat, batch_key='sample_ID', n_pcs=957)# it can also be run in the Scanpy Jupyter Notebook

# 8.6.3- Creating a Seurat object from AnnData

adata_fov_integrated = sc$read("./adata_fov_integrated.h5ad")
 
# Get the expression matrix
exprs <- t(adata_fov_integrated$X)
colnames(exprs) <- adata_fov_integrated$obs_names$to_list()
rownames(exprs) <- adata_fov_integrated$var_names$to_list()
# Create the Seurat object
seurat <- CreateSeuratObject(exprs)
# Set the expression assay
seurat <- SetAssayData(seurat, "data", exprs)
# Add observation metadata
seurat <- AddMetaData(seurat, adata_fov_integrated$obs)
# Add fetaure metadata
seurat[["RNA"]][["n_cells"]] <- adata_fov_integrated$var["n_cells"]
# Add embedding
embedding <- adata_fov_integrated$obsm["X_umap"]
rownames(embedding) <- adata_fov_integrated$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

########################################################################################################################################################################################################################################################################################

# 9- fastMNN integration in the SingleCellExperiment object
library(batchelor)
library(SingleCellExperiment)
library(scater)
library(SingleR)

DefaultAssay(fov_integrated) <- "RNA" 
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)
fov_sce <- as.SingleCellExperiment(fov_integrated, assay = "RNA")
saveRDS(fov_sce, "fov_integrated_sce.rds")

set.seed(220228)
out <- fastMNN(fov_sce, batch = fov_sce$sample_ID,
               auto.merge = TRUE,
               assay.type = "logcounts")
saveRDS(out, "sce_fastMNN_corrected.rds")

# 9.1- Transfer the correction results to the main sce object
reducedDim(fov_sce, "fastMNN") <- reducedDim(out, "corrected")

# 5.2- Transfer the fastMNN correction results to the harmony-integrated seurat object
fMNN <- reducedDim(out, "corrected")
fov_integrated[["fastMNN"]] <- CreateDimReducObject(embeddings = fMNN, key = "fastMNN_", assay = DefaultAssay(fov_integrated))

# 9.3- Quality control and visualisation of correction results

# The fastMNN function further returns outputs that can be used to assess the quality of the batch correction.
# The metadata(out)$merge.info entry collects diagnostics for each individual merging step.
# Here, the batch.size and lost.var entries are important. 
# The batch.size entry reports the relative magnitude of the batch effect and the lost.var entry represents the percentage of lost variance per merging step. 
# A large batch.size and low lost.var indicate sufficient batch correction.

merge_info <- metadata(out)$merge.info 

DataFrame(left = merge_info$left,
          right = merge_info$right,
          batch.size = merge_info$batch.size,
          max_lost_var = rowMax(merge_info$lost.var))

library(scater)
set.seed(220228)

# 9.4- Embedding corrected fastMNN in UMAP in the SCE and Seurat objects
fov_sce <- runUMAP(fov_sce, dimred= "fastMNN", name = "UMAP_mnnCorrected") 
fov_integrated <- RunUMAP(fov_integrated, reduction = "fastMNN", assay = "SCT", dims = 1:50, reduction.name = "umap.fastMNN")
DimPlot(fov_integrated, reduction = "umap.fastMNN", label = FALSE)

library(cowplot)
library(dittoSeq)
library(viridis)

# 9.5- Visualise UMAPs of the fastMNN corrected (integrated) data
dittoDimPlot(fov_sce, var = "sample_ID",
             reduction.use = "UMAP_mnnCorrected", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_sce, var = "CMS",
             reduction.use = "UMAP_mnnCorrected", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_sce, var = "sample_ID",
             reduction.use = "UMAP", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_integrated, var = "CMS",
             reduction.use = "umap.fastMNN", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")


########################################################################################################################################################################################################################################################################################


# 10- Adding embeddings from different batch-correction methods

fov_subset <- readRDS("./RDS_files/5_Integration/fov_subset_integration.rds")

# 10.1- Adding Harmony integration (from previous run in fov_integrated)

DefaultAssay(fov_subset) <- "RNA" 
fov_subset@reductions[["harmony"]] <- fov_integrated@reductions[["harmony"]]
fov_subset@reductions[["umap.harmony"]] <- fov_integrated@reductions[["umap.harmony"]]


# 10.2- Adding integrations from anndata
library(reticulate)
use_condaenv("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate")
use_python("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/bin/python")
library(sceasy)
library(rhdf5)
library(zellkonverter)
library(scater)
library(loomR)
library(Seurat)
library(patchwork)

adata_fov_subset_integration <- readH5AD('./h5ad_files/adata_fov_subset_integration.h5ad') 
adata_fov_subset_integration <- logNormCounts(adata_fov_subset_integration)
adata_fov_subset_integration2 <- as.Seurat(adata_fov_subset_integration)

fov_subset@reductions[["STAligner"]] <- adata_fov_subset_integration2@reductions[["STAligner"]]
fov_subset@reductions[["Scanorama"]] <- adata_fov_subset_integration2@reductions[["Scanorama"]]
fov_subset@reductions[["scVI"]] <- adata_fov_subset_integration2@reductions[["scVI"]]
fov_subset@reductions[["scANVI"]] <- adata_fov_subset_integration2@reductions[["X_scANVI"]]
fov_subset@reductions[["CellHint"]] <- adata_fov_subset_integration2@reductions[["CellHint"]]

fov_subset@meta.data[["majority_voting"]] <- adata_fov_subset_integration2@meta.data[["majority_voting"]]
fov_subset@meta.data[["conf_score"]] <- adata_fov_subset_integration2@meta.data[["conf_score"]]
fov_subset@meta.data[["CIPR_annotated_clusters"]] <- adata_fov_subset_integration2@meta.data[["CIPR_annotated_clusters"]]

fov_subset@meta.data[["majority_voting2"]] <- adata_fov_subset_integration2@meta.data[["majority_voting2"]]
fov_subset@meta.data[["conf_score2"]] <- adata_fov_subset_integration2@meta.data[["conf_score2"]]

# 10.3- Run UMAP with new integrations

set.seed(22022)
fov_subset <- RunUMAP(fov_subset, reduction = "STAligner", assay = "SCT", dims = 1:30, reduction.name = "umap.staligner") 
fov_subset <- RunUMAP(fov_subset, reduction = "Scanorama", assay = "SCT", dims = 1:50, reduction.name = "umap.scanorama")
fov_subset <- RunUMAP(fov_subset, reduction = "scVI", assay = "SCT", dims = 1:10, reduction.name = "umap.scvi")
fov_subset <- RunUMAP(fov_subset, reduction = "scANVI", assay = "SCT", dims = 1:10, reduction.name = "umap.scanvi")

DimPlot(fov_subset, reduction = "umap.scanvi", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(fov_subset, var = "sample_ID",
             reduction.use = "umap.scanvi", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "sample_ID", split.by = "sample_ID",
             reduction.use = "umap.scanvi", size = 0.75,
             do.label = FALSE)

dittoDimPlot(fov_subset, var = "IF_type",
             reduction.use = "umap.scanvi", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell types based on IF on UMAP")

dittoDimPlot(fov_subset, var = "IF_type", split.by = "sample_ID",
             reduction.use = "umap.scanvi", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_subset, var = "majority_voting",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "majority_voting2",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

dittoDimPlot(fov_subset, var = "conf_score",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = FALSE, labels.size = 3, labels.highlight = FALSE) +
  theme(legend.title = element_blank()) +
  ggtitle("Confidence score on UMAP")

FeaturePlot(fov_subset, features = "conf_score", reduction = "umap.harmony", cols = c("#00924C", "#FEFA35"))
FeaturePlot(fov_subset, features = "conf_score2", reduction = "umap.harmony", cols = c("#00924C", "#FEFA35"))

saveRDS(fov_subset, "fov_subset_integration.rds")
saveRDS(adata_fov_subset_integration, "adata_fov_subset_integration_sce.rds")
saveRDS(adata_fov_subset_integration2, "adata_fov_subset_integration_seurat.rds")


# 10.4- Plot cells types per fov - comapre previous annotated clusters with CIPR vs CellTypist annotation

# following Colin's suggestions
cell_type_color <- setNames(c("#00924C","#00EA7B","#D5FFEB","#783D0D",
                              "#BC9201",'#93B220',"#6FB19B", "#D4E8E1",
                              "#AD0791","#FCCDE5","#BBABFF","#7853FF",
                              "#F62ED5","#FEFA35","#E42520","#FB9797",
                              "#FDC77F","#013F89","#01A6C7","#FF7F00",
                              "#5DA5FD","#7AE8FE"
),
c("Tumor epithelial cell", "ChemokinesHigh Tumor epithelial cell", "IGF2+AREG+ Tumor epithelial cell", "Stem cell", 
  "IGF2+AREG+ stem cell", "Goblet cell", "Tuft cell", "ChemokinesHigh tuft cell",
  "Plasma cell", 'B cell', "CD8+ T cell", "CD4+ T cell", "Treg",
  'Fibroblast', 'Pericyte', "Endothelial cell", 'Myofibroblast',
  'SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage',
  'ChemokinesHigh monocyte'))

Idents(fov_subset) <- 'CIPR_annotated_clusters'
levels(fov_subset)

# 19.4- Visualise clusters in a spatial context
ImageDimPlot(fov_subset, fov = "fov01", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov02", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov03", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov04", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov05", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov06", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov07", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov08", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov09", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov10", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov11", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov12", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov13", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov14", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov15", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov16", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov17", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov18", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov19", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov20", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov21", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov22", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov23", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov24", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
ImageDimPlot(fov_subset, fov = "fov25", axes = FALSE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)


cell_type_color <- setNames(c("#00924C","#00EA7B","#D5FFEB","#783D0D",
                              "#BC9201",'#93B220',"#6FB19B", "#D4E8E1",
                              "#AD0791","#FCCDE5","#BBABFF","#7853FF",
                              "#F62ED5","#FEFA35","#E42520","#FB9797",
                              "#FDC77F","#013F89","#01A6C7","#FF7F00",
                              "#5DA5FD","#7AE8FE"
),
c("Tumor epithelial cell", "ChemokinesHigh Tumor epithelial cell", "IGF2+AREG+ Tumor epithelial cell", "Stem cell", 
  "IGF2+AREG+ stem cell", "Goblet cell", "Tuft cell", "ChemokinesHigh tuft cell",
  "Plasma cell", 'B cell', "CD8+ T cell", "CD4+ T cell", "Treg",
  'Fibroblast', 'Pericyte', "Endothelial cell", 'Myofibroblast',
  'SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage',
  'ChemokinesHigh monocyte'))

ImageDimPlot(fov_subset, fov = "fov01", axes = FALSE, group.by = 'majority_voting', cols = "alphabet", size = 1, boundaries = "segmentation", dark.background = TRUE, coord.fixed = TRUE)
 


# Extract global position of cells per fov:
df1 <- fov_subset@images[["fov01"]]@boundaries[["centroids"]]@coords
df2 <- fov_subset@images[["fov02"]]@boundaries[["centroids"]]@coords
df3 <- fov_subset@images[["fov03"]]@boundaries[["centroids"]]@coords
df4 <- fov_subset@images[["fov04"]]@boundaries[["centroids"]]@coords
df5 <- fov_subset@images[["fov05"]]@boundaries[["centroids"]]@coords
df6 <- fov_subset@images[["fov06"]]@boundaries[["centroids"]]@coords
df7 <- fov_subset@images[["fov07"]]@boundaries[["centroids"]]@coords
df8 <- fov_subset@images[["fov08"]]@boundaries[["centroids"]]@coords
df9 <- fov_subset@images[["fov09"]]@boundaries[["centroids"]]@coords
df10 <- fov_subset@images[["fov10"]]@boundaries[["centroids"]]@coords
df11 <- fov_subset@images[["fov11"]]@boundaries[["centroids"]]@coords
df12 <- fov_subset@images[["fov12"]]@boundaries[["centroids"]]@coords
df13 <- fov_subset@images[["fov13"]]@boundaries[["centroids"]]@coords
df14 <- fov_subset@images[["fov14"]]@boundaries[["centroids"]]@coords
df15 <- fov_subset@images[["fov15"]]@boundaries[["centroids"]]@coords
df16 <- fov_subset@images[["fov16"]]@boundaries[["centroids"]]@coords
df17 <- fov_subset@images[["fov17"]]@boundaries[["centroids"]]@coords
df18 <- fov_subset@images[["fov18"]]@boundaries[["centroids"]]@coords
df19 <- fov_subset@images[["fov19"]]@boundaries[["centroids"]]@coords
df20 <- fov_subset@images[["fov20"]]@boundaries[["centroids"]]@coords
df21 <- fov_subset@images[["fov21"]]@boundaries[["centroids"]]@coords
df22 <- fov_subset@images[["fov22"]]@boundaries[["centroids"]]@coords
df23 <- fov_subset@images[["fov23"]]@boundaries[["centroids"]]@coords
df24 <- fov_subset@images[["fov24"]]@boundaries[["centroids"]]@coords
df25 <- fov_subset@images[["fov25"]]@boundaries[["centroids"]]@coords

df <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19,
            df20, df21, df22, df23, df24, df25)


df26 <- as.data.frame(fov_subset@meta.data[["sample_ID"]])
df <- cbind(df, df26)
write.csv(df,'spatial_fov_fov_subset_integrated.csv')

# Edit spatial and spatial_fov files based on the retained cells after QC in Seurat

df1 <- read.csv("S1_TMA_metadata_file_seurat_scanpy.csv")
df2 <- read.csv("S1_TMA_metadata_file.csv")

idx <- match(df1$CELL_ID2, df2$CELL_ID2)
df1$CenterX_local_px <- df2$CenterX_local_px[ idx ]
df1$CenterY_local_px <- df2$CenterY_local_px[ idx ]
write.csv(df1, "S1_TMA_metadata_file_seurat_scanpy_corrected.csv")




