
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

options(ggrepel.max.overlaps = Inf)

fov_integrated <- readRDS("./RDS_files/5_Integration/All/fov_integrated.rds")
saveRDS(fov_integrated, "./fov_integrated.rds")

sce <- readRDS("./RDS_files/6_Gated/fov_epithelial_sce.rds")
saveRDS(sce, "fov_integrated_sce.rds")

fov_epithelial <- readRDS("./RDS_files/6_Gated/fov_epithelial.rds")

saveRDS(fov_epithelial, "./fov_epithelial.rds")


# 14- Cluster annotation with marker gene databases

# 14.1- Cell Annotation with scCATCH

# Reference: https://github.com/ZJUFanLab/scCATCH
# Reference: https://raw.githack.com/ZJUFanLab/scCATCH_performance_comparison/master/scCATCH/tutorial.html
# This package is a single cell Cluster-based auto-Annotation Toolkit for Cellular Heterogeneity (scCATCH) from cluster potential marker genes identification to cluster annotation based on evidence-based score 
# by matching the potential marker genes with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).
library(scCATCH)

set.seed(1234)
DefaultAssay(fov_epithelial) <- "RNA"
Idents(fov_epithelial) <- 'ep_harmony_clusters'
levels(fov_epithelial)

# Create scCATCH object from Seurat object:
obj <- createscCATCH(data = fov_epithelial[['RNA']]@data, cluster = as.character(Idents(fov_epithelial)))
saveRDS(obj, "fov_epithelial_scCATCH_object.rds")

# 14.1.1- Find marker gene for each cluster
# Users need to provided the speices, tissue, or cancer information. 
# Use similar parmaters to those used with the FindAllMarkers function from Seurat
# Available tissues and cancers at https://github.com/ZJUFanLab/scCATCH/wiki
# Select different combination of tissues or cancers for annotation
#cellmatch$tissue
# filter cellmatch
#cellmatch <- cellmatch[cellmatch$species == "Human", ]
#cellmatch <- cellmatch[cellmatch$tissue %in% c("Blood", "Peripheral blood", "Serum", "Colon", "Colorectum", "Intestine"), ]
#cellmatch <- cellmatch[cellmatch$cancer %in% c("Colon Cancer", "Colorectal Cancer"), ]

obj <- findmarkergene(obj, species = 'Human', cluster = "All", if_use_custom_marker = FALSE, marker = cellmatch, 
                      tissue = c("Colon", "Colorectum", "Intestine"),
                      cancer = c("Colorectal Cancer"),
                      use_method = "2", cell_min_pct = 0.1,logfc = 0.25, pvalue = 0.05, verbose = TRUE)

# 14.1.2- Find cell type for each cluster
obj <- findcelltype(obj)

list <- obj@celltype
write.csv(list, "scCatch_annotated_epithelial_clusters")

# scCatch cluster annotation and visualisation

SCATCH_clusters <- c("T cell", "T cell", "Tumor epithelial cell", "Tumor epithelial cell", "S100P+ epithelial cell",
                     "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell",
                     "T cell", "T cell", "Cancer-associated fibroblast", "SPP1+CXCL8+CXCL1+ monocyte", "SPP1+MMP12+CXCL10+ macrophage",
                     "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell", "Tumor epithelial cell",
                     'B cell', "B cell", "B cell")
names(SCATCH_clusters) <- levels(fov_epithelial)
fov_epithelial <- RenameIdents(fov_epithelial, SCATCH_clusters)

# Adding object metadata with
fov_epithelial$SCATCH_clusters <- factor(fov_epithelial@active.ident)

dittoDimPlot(fov_epithelial, var = "SCATCH_clusters",
                   reduction.use = "umap.harmony", size = 1,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("scCATCH epithelial clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
Idents(fov_epithelial) <- "SCATCH_clusters"
ImageDimPlot(fov_epithelial, fov = "fov02", group.by = "SCATCH_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", group.by = "SCATCH_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_epithelial, var = "SCATCH_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "SCATCH_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "SCATCH_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "SCATCH_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_epithelial, features = features, sort = TRUE, ncol = 2)

# 14.2- Cell Annotation with SCSA

# Reference: https://github.com/bioinfo-ibms-pumc/SCSA
# To annotate a human scRNA-seq sets generated by 'FindAllMarkers' function of Seurat with ensemblIDs, use the following code:
# -d DB, --db DB        Database for annotation. (whole.db)
# -s SOURCE, --source Source of marker genes. (cellranger,[seurat],[scanpy], [scran])
# -i INPUT, --input INPUT Input file for marker annotation(Only CSV format supported).
# -k TISSUE, --tissue Tissue for annotation. Only used for cellmarker database. Multiple tissues should be seperated by commas.Run '-l' option to see all tissues.
# In linux platform:('All',['Bone marrow'],['Bone marrow,Brain,Blood'][...])
# python3 SCSA.py -d whole_v2.db -s seurat -i Markers_epithelial_hamorny_clusters.csv -k All -E -g Human -p 0.05 -f 1.5

# Add SCSA cluster lables
Idents(fov_epithelial) <- "ep_harmony_clusters"
levels(fov_epithelial)
SCSA_clusters <- c("T cell", "Mast cell", "Epithelial cell", "Epithelial cell", "Epithelial cell", "Monocyte",
                   "Epithelial cell", "Fibroblast", "Epithelial cell", "T cell", "Plasma cell", "Fibroblast",
                   "Monocyte_macrophage", "Macrophage", "Epithelial cell", "Epithelial cell", "Epithelial cell", "Epithelial cell", "Stem cell",
                   "Epithelial cell", "Epithelial cell", "T cell", "Plasma cell", "B cell", "B cell") 
names(SCSA_clusters) <- levels(fov_epithelial)
fov_epithelial <- RenameIdents(fov_epithelial, SCSA_clusters)

# Adding object metadata with
fov_epithelial$SCSA_clusters <- factor(fov_epithelial@active.ident)

multi_dittoDimPlot(fov_epithelial, vars = c("SCATCH_clusters","SCSA_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
Idents(fov_epithelial) <- "SCSA_clusters"
ImageDimPlot(fov_epithelial, fov = "fov02", group.by = c("SCSA_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", group.by = c("SCSA_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_epithelial, var = "SCSA_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "SCSA_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "SCSA_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "IF_marker", group.by = "SCSA_clusters",  scale = c("percent"))
dittoBarPlot(fov_epithelial, var = "IF_type2", group.by = "SCSA_clusters",  scale = c("percent"))

# Plot normalised protein expression per cluster
DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "SCSA_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_epithelial, features = features, sort = TRUE, ncol = 2)

########################################################################################################################################################################################################################################################################################

# 15- Cluster annotation with correlation-based annotation

# 15.1- Cell Annotation with SingleR

# Reference: https://github.com/dviraran/SingleR
# Reference: https://github.com/LTLA/SingleR/blob/master/README.md
# Reference: Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.
# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html
# Reference: http://bioconductor.org/books/release/SingleRBook/
# Reference: https://nbisweden.github.io/excelerate-scRNAseq/session-celltypeid/celltypeid.html

# Example of the use of SingleR() with the Human Primary Cell Atlas dataset (Mabbott et al. 2013) as the reference.
# Loading reference data with Ensembl annotations.
# Each celldex dataset actually has three sets of labels that primarily differ in their resolution. 
# For the purposes of this demonstration, we will use the “fine” labels in the label.fine metadata field, which represents the highest resolution of annotation available for this dataset.
ref.data <- readRDS("./RDS_files/7_Reference_datasets/HumanPrimaryCellAtlasData/HumanPrimaryCellAtlasData.rds")
ref.data2 <- readRDS("./RDS_files/7_Reference_datasets/BlueprintEncodeData/BlueprintEncodeData.rds")

# Human databases in the library celldex:
# HumanPrimaryCellAtlasData - ownload and cache the normalized expression values of the data stored in the Human Primary Cell Atlas.
# BlueprintEncodeData - Download and cache the normalized expression values of 259 RNA-seq samples of pure stroma and immune cells as generated and supplied by Blueprint and ENCODE.
# DatabaseImmuneCellExpressionData - Download and cache the normalized expression values of 1561 bulk RNA-seq samples of sorted cell populations from the Database of Immune Cell Expression (DICE).
# NovershternHematopoieticData - ownload and cache the normalized expression values of 211 bulk human microarray samples of sorted hematopoietic cell populations that can be found in GSE24759.
# Samples were additionally annotated to 38 fine cell types ("label.fine").
# MonacoImmuneData - Download and cache the normalized expression values of 114 bulk RNA-seq samples of sorted immune cell populations

# 15.1.1- Performing predictions. Convert Seurat object to SingleCellExperiment object
library(SingleCellExperiment)
library(scater)
library(SingleR)

DefaultAssay(fov_epithelial) <- "RNA" 
Idents(fov_epithelial) <- 'ep_harmony_clusters'
levels(fov_epithelial)
sce <- as.SingleCellExperiment(fov_epithelial, assay = "RNA")
saveRDS(sce, "fov_epithelial_sce.rds")

# We perform annotation by calling SingleR() on our test dataset and the reference (ImmGen) dataset, 
# leaving the default of de.method="classic" to use the original marker detection scheme. 
# This applies the algorithm described in Section 1.2, returning a DataFrame where each row contains prediction results for a single cell in the sce object. 
# Labels are provided in the labels column.
# See 'Choices of assay data' for 'assay.type.test=' explanation.

# 15.1.2- Using the HumanPrimaryCellAtlasData ref 

# main labels
pred <- SingleR(test = sce, ref = ref.data, 
                labels = ref.data$label.main, assay.type.test=1)
colnames(pred)
saveRDS(pred, "./fov_epithelial_SingleR.rds")
pred <- readRDS("./fov_epithelial_SingleR.rds")

# fine labels
pred2 <- SingleR(test = sce, ref = ref.data, 
                 labels = ref.data$label.fine, assay.type.test=1)
colnames(pred2)
saveRDS(pred2, "./fov_integrated_SingleR_pred2.rds")
pred2 <- readRDS("./fov_integrated_SingleR_pred2.rds")

# 15.1.3- Using the Blueprintenconde ref

# main labels
pred3 <- SingleR(test = sce, ref = ref.data2, 
                 labels = ref.data2$label.main, assay.type.test=1)
colnames(pred3)
saveRDS(pred3, "./fov_integrated_SingleR_pred3.rds")
pred3 <- readRDS("./fov_integrated_SingleR_pred3.rds")

# fine labels
pred4 <- SingleR(test = sce, ref = ref.data2, 
                 labels = ref.data2$label.fine, assay.type.test=1)
colnames(pred4)
saveRDS(pred4, "./fov_integrated_SingleR_pred4.rds")
pred4 <- readRDS("./fov_integrated_SingleR_pred4.rds")

# 15.1.4- Adding SingleR annotated epithelial clusters and visualisation
SINGLER_clusters <- pred@listData[["labels"]]
fov_epithelial@meta.data <-cbind(fov_epithelial@meta.data, SINGLER_clusters)

dittoDimPlot(fov_epithelial, var = "SINGLER_clusters",
                   reduction.use = "umap.harmony", size = 1,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

multi_dittoDimPlot(fov_epithelial, vars = c("SCATCH_clusters","SCSA_clusters", "SINGLER_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
Idents(fov_epithelial) <- "SINGLER_clusters"
ImageDimPlot(fov_epithelial, fov = "fov02", group.by = c("SINGLER_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", group.by = c("SINGLER_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_epithelial, var = "SINGLER_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "SINGLER_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "SINGLER_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "SINGLER_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_epithelial, features = features, sort = TRUE, ncol = 2)


# 15.2- Cell Annotation with CIPR

# Reference: https://github.com/atakanekiz/CIPR-Shiny
library(CIPR)
DefaultAssay(fov_epithelial) <- "RNA" 
Idents(fov_epithelial) <- 'ep_harmony_clusters'
levels(fov_epithelial)

# 15.2.1- Calculate average gene expression per cluster
# This is the step where we generate the input for CIPR's all-genes correlation methods
avgexp <- AverageExpression(fov_epithelial)
avgexp <- as.data.frame(avgexp$RNA)
avgexp$gene <- rownames(avgexp)
write.csv(avgexp, "avgexp_epithelial.csv")

# CIPR analysis

# 15.2.2- Standard logFC comparison method
# Immunological Genome Project (ImmGen) microarray data from sorted mouse immune cells. This dataset is prepared by using both V1 and V2 ImmGen releases and it contains 296 samples from 20 different cell types (253 subtypes).
# Mouse RNAseq data from sorted cells reported in Benayoun et al. (2019). This dataset contains 358 sorted immune and nonimmune samples from 18 different lineages (28 subtypes).
# Blueprint/Encode RNAseq data that contains 259 sorted human immune and nonimmune samples from 24 different lineages (43 subtypes).
# Human Primary Cell Atlas (hpca) that contains microarray data from 713 sorted immune and nonimmune cells (37 main cell types and 157 subtypes).
# DICE (Database for Immune Cell Expression(/eQTLs/Epigenomics) that contains 1561 human immune samples from 5 main cell types (15 subtypes). To reduce object sizes, mean TPM values per cell type is used.
# Human microarray data from sorted hematopoietic cells reported in Novershtern et al. (2011). This dataset contains data from 211 samples and 17 main cell types (38 subtypes)
# Human RNAseq (hsrna) data from sorted cells reported in Monaco et al. (2019). This dataset contains 114 samples originating from 11 main cell types (29 subtypes)
# A custom reference dataset provided by the user. This dataset can be obtained from a number of high througput methods including microarray and bulk RNAseq. For details about how to prepare custom reference, please see the How-to tab on the Shiny website.
        
allmarkers <- read.csv("./Plots/7_Clustering/Option 2/Epithelial/Harmony/Markers_epithelial_hamorny_clusters_RNA.csv")

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hsrnaseq", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T,
     # axis.text.x=element_text(color="red") # arguments to pass to ggplot2::theme() to change plotting parameters
)

# Plot identity scores for a select cluster
library(ggplot2)

ind_clu_plots$cluster0 +
  theme(axis.text.y = element_text(color="red"),
        axis.text.x = element_text(color="blue")) +
  labs(fill="Reference")+
  ggtitle("Automated cluster annotation results are shown for cluster 0") +
  annotate("text", label="2 sd range", x=10, y= 500, size=8, color = "steelblue")+
  annotate("text", label= "1 sd range", x=10, y=175, size=8, color ="orange2")+
  geom_rect(aes(xmin=94, xmax=99, ymin=550, ymax=900), fill=NA, linewidth=3, color="red")

# Tabulate CIPR results
# CIPR results (both top 5 scoring reference types per cluster and the entire analysis) are saved as global objects (`CIPR_top_results` and `CIPR_all_results` respectively) to allow users to explore the outputs and generate specific plots and tables.
DT::datatable(CIPR_top_results)
DT::datatable(head(CIPR_all_results))

write.csv(CIPR_top_results, "CIPR_top_results_hpca.csv")
write.csv(CIPR_all_results, "CIPR_all_results_hpca.csv")
write.csv(CIPR_top_results, "CIPR_top_results_avgexp_correlation.csv")

# 15.2.3- Standard all-genes correlation method
# CIPR also implements a simple correlation approach in which overall correlation in gene expression is calculated for the pairs of unknown clusters and the reference samples (regardless of the differential expression status of the gene). 
# This approach is conceptually similiar to some other automated identity prediction pipelines such as [SingleR](https://www.ncbi.nlm.nih.gov/pubmed/30643263) and [scMCA](https://www.ncbi.nlm.nih.gov/pubmed/30758821).
# Users can select one of the following methods:
# Spearman's correlation:__ It calculates correlation based on ranked gene expression. It can be suitable for comparing experimental and reference data which were obtained using different technologies.
# Pearson's correlation:__ It calculates linear correlations. This can be useful when the user would like to provide a custom reference dataset to CIPR.

CIPR(input_dat = avgexp,
     comp_method = "all_genes_spearman", 
     reference = "hpca", 
     plot_ind = T,
     plot_top = F, 
     global_results_obj = T, 
     global_plot_obj = T)


# 15.2.4- Add CIPR cluster labels
Idents(fov_epithelial) <- "ep_harmony_clusters"
levels(fov_epithelial)
CIPR_hpca_clusters <- c("T cell", "T cell", "Epithelial cell",  "Epithelial cell",  "Epithelial cell", "Monocyte",  "Epithelial cell",  "Epithelial cell",
                   "Monocyte-derived macrophage", "Neuronal cell", "T cell", "Tissue stem cell", "SPP1+ Monocyte-derived DC", "Macrophage", "Monocyte-derived DC", 
                   "Epithelial cell",  "Epithelial cell",  "Epithelial cell", "NK cell", "NK cell", "Monocyte", "T cell", "Neuronal cell", "B cell", "B cell") 
names(CIPR_hpca_clusters) <- levels(fov_epithelial)
fov_epithelial <- RenameIdents(fov_epithelial, CIPR_hpca_clusters)

# Adding object metadata with
fov_epithelial$CIPR_hpca_clusters <- factor(fov_epithelial@active.ident)

dittoDimPlot(fov_epithelial, var = "CIPR_hpca_clusters",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

multi_dittoDimPlot(fov_epithelial, vars = c("SCATCH_clusters","SCSA_clusters", "SINGLER_clusters", "CIPR_hpca_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

# Visualise labelled clusters in spatial context
Idents(fov_epithelial) <- "CIPR_hpca_clusters"
ImageDimPlot(fov_epithelial, fov = "fov02", group.by = "CIPR_hpca_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", group.by = "CIPR_hpca_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_epithelial, var = "CIPR_hpca_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "CIPR_hpca_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "CIPR_hpca_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "CIPR_hpca_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_epithelial, features = features, sort = TRUE, ncol = 2)

# 15.3- Making custom references with the Nanostring datasets

# Reference: https://github.com/Nanostring-Biostats/CellProfileLibrary/tree/master/Human/Adult
# Reference: https://aekiz.shinyapps.io/CIPR/
# • See code in the file "cipr_ref_prepation.rmd"

# This repo contains a library of "cell profile matrices" - matrices giving the average expression profiles of all cell types found in a tissue. Each matrix in the library was derived from a single scRNA-seq experiment. These matrices can be used with cell type deconvolution packages like SpatialDecon to get cell type proportions or identities.
# Each RData file contains 3 file types:#
# Cell Profile Matrix: average expression profile of all cell types found in a tissue.
# Cell Groups
# Dataset Metadata: Includes database, tissue, age, and paper info
# Matrices were generated using published datasets that had annotated cell types in human and mouse. Datasets were normalized, by total gene count, if raw data was used otherwise the publication’s normalization used. These datasets were filtered for cells expressing (count > 0) at least 100 genes and only calculated cell type profiles for cell types with 15+ viable cells. Profiles were created by taking the average expression of each gene across all viable cells of each cell type.

#log normalised references
gut_ref_epithelial <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/gut_ref_epithelial.csv")
rectum_ref <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/rectum_ref.csv")

merged <- merge(rectum_ref, gut_ref_epithelial, by.x = "Gene", by.y = "Gene") # gave the best performance for the annotation
write.csv(merged,"./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/rectum_gut_ref_epithelial.csv")
rectum_gut_ref_epithelial <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/rectum_gut_ref_epithelial.csv")

# Annotation with CIPR
allmarkers <- read.csv("./Plots/7_Clustering/Option 2/Epithelial/Harmony/Markers_epithelial_hamorny_clusters_RNA.csv")
allmarkers <- as.data.frame(allmarkers)

avgexp <- read.csv("./Plots/7_Clustering/Option 2/Epithelial/Harmony/avgexp_epithelial.csv", sep=",", row.names=1, header = TRUE)

# correlation method
CIPR(input_dat = avgexp,
     comp_method = "all_genes_spearman", 
     reference = "custom",
     custom_reference = rectum_gut_ref_epithelial,  
     plot_ind = F,
     plot_top = T, 
     global_results_obj = T, 
     global_plot_obj = F)

write.csv(CIPR_top_results, "CIPR_top_results_gut_ref_epithelial_ref.csv")


# 15.2.4- Add CIPR cluster lables
Idents(fov_epithelial) <- "ep_harmony_clusters"
levels(fov_epithelial)
#fov_epithelial$CIPR_rectum_gut_clusters <- NULL
CIPR_rectum_gut_clusters <- c("Tumor epthelial cell", 'Goblet cell', "Tumor epthelial cell", "Tumor epthelial cell", "Goblet cell",
                          "ChemokinesHigh Tumor epthelial cell", 'Tumor epthelial cell', "Tuft cell", "IGF2+AREG+ stem cell", "Tumor epthelial cell", "Tumor epthelial cell",
                          "Fibroblast","ChemokinesHigh monocyte", "SPP1+ macrophage", "IGF2+AREG+ Tumor epthelial cell", "Goblet cell", "Tuft cell", "Goblet cell", "Stem cell",
                          "IGF2+AREG+ stem cell", "ChemokinesHigh tuft cell", "Tumor epthelial cell", "B cell", "B cell", "Plasma cell") 
names(CIPR_rectum_gut_clusters) <- levels(fov_epithelial)
fov_epithelial <- RenameIdents(fov_epithelial, CIPR_rectum_gut_clusters)

# Adding object metadata with
fov_epithelial$CIPR_rectum_gut_clusters <- factor(fov_epithelial@active.ident)

dittoDimPlot(fov_epithelial, var = "CIPR_rectum_gut_clusters",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

multi_dittoDimPlot(fov_epithelial, vars = c("SCATCH_clusters","SCSA_clusters", "CIPR_rectum_clusters", "CIPR_rectum_gut_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
Idents(fov_epithelial) <- "CIPR_rectum_gut_clusters"
ImageDimPlot(fov_epithelial, fov = "fov02", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_epithelial, var = "CIPR_rectum_gut_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "CIPR_rectum_gut_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "CIPR_rectum_gut_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "CIPR_rectum_gut_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_epithelial, features = features, sort = TRUE, ncol = 2)

DefaultAssay(fov_epithelial) <- "RNA" 
fov_epithelial.markers <- FindAllMarkers(fov_epithelial, only.pos = TRUE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_epithelial.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_epithelial.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap((subset(fov_epithelial, downsample = 1000)), features = top10$gene, size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 8))

dittoHeatmap(fov_epithelial, top10$gene, 
             cluster_rows = FALSE, 
             annot.by = "CIPR_rectum_gut_clusters", 
             fontsize_row = 6)

avgexp <- AverageExpression(subset(fov_epithelial, features = top10$gene))
avgexp <- as.data.frame(avgexp$RNA)
write.csv(avgexp, "Average_expression_Top_expressed_genes_CIPR_epithelial_clusters_RNA.csv")
#avgexp$gene <- rownames(avgexp)
#avgexp <- avgexp[,1:13]

#avgexp2 <- scale(avgexp)
avgexp3 <- t(avgexp)
avgexp4 <- scale(avgexp3)
avgexp5 <- t(avgexp4)

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
heat <- Heatmap(avgexp5, cluster_columns = TRUE, cluster_rows = FALSE,
                col = colorRamp2(c(1, 0.6, 0, -0.6, -1), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                heatmap_legend_param = list(color_bar = "continuous"),
                show_row_dend = TRUE, show_column_dend = TRUE,
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 10))

print(heat)

########################################################################################################################################################################################################################################################################################

# 16- Cluster annotation with supervised classiﬁcation

# 16.1- Cell Annotation with CHETAH (CHaracterization of cEll Types Aided by Hierarchical classification)

# Reference: https://github.com/jdekanter/CHETAH
# http://www.bioconductor.org/packages/devel/bioc/vignettes/CHETAH/inst/doc/CHETAH_introduction.html
# Cell types are assigned by correlating the input data to a reference in a hierarchical manner. 
# This creates the possibility of assignment to intermediate types if the data does not allow to fully classify to one of the types in the reference.
# CHETAH is built to work with scRNA-seq references, but will also work (with limited capabilities) with RNA-seq or micro-array reference datasets.
# to run CHETAH, you will only need: your input data and a reference dataset, annotated with cell types. Both as a SingleCellExperiment.
# CHETAH constructs a classification tree by hierarchically clustering the reference data. The classification is guided by this tree. In each node of the tree, input cells are either assigned to the right, or the left branch. A confidence score is calculated for each of these assignments. When the confidence score for an assignment is lower than the threshold (default = 0.1), the classification for that cell stops in that node.
library(CHETAH)

# 16.1.1- Preparing the input data
# For the input we define a "counts" assay and "TSNE" or "UMAP" for reduced dimensions
library(SingleCellExperiment)
library(scater)

DefaultAssay(fov_epithelial) <- "RNA" 
Idents(fov_epithelial) <- 'ep_harmony_clusters'
levels(fov_epithelial)
sce <- as.SingleCellExperiment(fov_epithelial, assay = "RNA")
saveRDS(sce, "fov_epithelial_sce.rds")

# Load CHETAH reference
chetah_ref <- get(load('./RDS_files/7_Reference_datasets/Chetah_reference/CHETAH_TME_reference.Rdata'))

# Customised reference data:
colon_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/colon_ref_chetah.rds")
gut_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/gut_ref_chetah.rds")
ileum_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/ileum_ref_chetah.rds")
rectum_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/rectum_ref_chetah.rds")
colon_gut_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/colon_gut_ref_chetah.rds")

# Running CHETAH
sce <- CHETAHclassifier(input = sce, ref_cells = gut_ref, thresh = 0.1) # it took 2h to run...

# output
# CHETAH returns the input object, but added: input$celltype_CHETAH - a named character vector that can directly be used in any other workflow/method.
# “hidden” int_colData and int_metadata, not meant for direct interaction, but which can all be viewed and interacted with using: PlotCHETAH and CHETAHshiny
PlotCHETAH(input = sce, redD = "UMAP.HARMONY")

########################################################################################################################################################################################################################################################################################################################################################################################################

# 16.2- Cell Annotation with InSituType

# Reference: file:///Users/joaoluizsfilho/Packages/scRNAseq/InSituType-main/vignettes/NSCLC-semi-supervised-cell-typing-vignette.html
library(InSituType)

# 16.2.1- Preparing reference matrixes

# Use the original references from the CellProfileLibrary
# Necessary inputs: matrix of counts data, cells x genes in the SingleCelllExperiment format
# A “reference matrix” giving the expected expression profile of each cell type, with genes in rows and cell types in columns. 
# The reference matrix must be in linear-scale, not log-scale. Insitutype can handle aligning its genes (rows) to your counts matrix. 
# For supervised cell typing, the reference matrix should contain every cell type present in the tissue. 
# If you suspect the tissue has cell types outside your reference matrix, use insitutype’s semi-supervised clustering capability.

# Customised rectum reference data:
rectum_ref2 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_ref.csv", sep=",", row.names=1)
rectum_ref3  <- as.matrix(rectum_ref2)
class(rectum_ref3)
str(rectum_ref3)
saveRDS(rectum_ref3, "rectum_ref_insitutype.rds")
rectum_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_ref_insitutype.rds")

# Customised gut refs data with non-epithelial cells removed:
gut_ref2 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_epithelial.csv", sep=",", row.names=1)
gut_ref3  <- as.matrix(gut_ref2)
str(gut_ref3)
saveRDS(gut_ref3, "./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_epithelial_insitutype.rds")
gut_ref_epithelial <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_epithelial_insitutype.rds")

# Customised rectum+gut refs data:
gut_ref2 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_epithelial.csv")
rectum_ref2 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_ref.csv")
merged <- merge(rectum_ref2, gut_ref2, by.x = "Gene", by.y = "Gene") # gave the best performance for the annotation
write.csv(x = merged, file = "./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_gut_ref_epithelial.csv", row.names = TRUE, quote = FALSE)

rectum_gut_ref_epithelial <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_gut_ref_epithelial.csv", sep=",", row.names=1)
rectum_gut_ref_epithelial2  <- as.matrix(rectum_gut_ref_epithelial)
str(rectum_gut_ref_epithelial2)
saveRDS(rectum_gut_ref_epithelial2, "./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_gut_ref_epithelial_insitutype.rds")
rectum_gut_ref_epithelial <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_gut_ref_epithelial_insitutype.rds")

# Customised colon+gut and epithelial+gut refs data:
colon_rectum_gut_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/colon_rectum_gut_ref_insitutype.rds")

# 16.2.2- Extract the count matrix with the negative probe expression and calculate cells background

fov_subset <- readRDS("./RDS_files/4_Subset/fov_subset_neg_probes.rds")
Idents(fov_subset) <- "IF_type"
levels(fov_subset)
fov_subset <- subset(fov_subset, idents = c("Epithelial"))

counts <- LayerData(fov_subset, assay = "RNA")
saveRDS(counts, "epithelial_all_counts.rds")
neg <- counts[(which(rownames(counts) %in% c('NegPrb3',	'NegPrb5', 'NegPrb6',	'NegPrb7', 'NegPrb8',	'NegPrb9', 'NegPrb10',	'NegPrb11',	'NegPrb12',	'NegPrb13',	'NegPrb14',	'NegPrb15',	'NegPrb16', 'NegPrb18',	'NegPrb19',	'NegPrb20',	'NegPrb21',	'NegPrb22',	'NegPrb23'))),]
saveRDS(neg, "epithelial_negative_counts.rds")

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

# estimate per-cell bg as a fraction of total counts:
negmean.per.totcount <- mean(negmean) / mean(rowSums(counts))
per.cell.bg <- rowSums(counts) * negmean.per.totcount

# 16.2.3- Updating the reference profiles to account for platform effects

# As in most cases, our reference matrix was derived using a different platform. As such, it will be a noisy representation of what cells will look like in CosMx data.
# To avoid disadvantaging our reference cell types compared to newly-derived clusters, we’ll have to update the reference profiles to better fit CosMx data.
# 1. Identify “anchor cells”: cells that can confidently be assigned to a reference cell type. Anchor cells can be hand-selected based on marker genes or morphology. 
# But for most cases, they will be selected using the data-driven approach implemented below. 
# 2. Estimate the mean profile of each cell type’s anchors. These mean profiles form the updated reference matrix.
# This step is slow for big datasets. It's recommended to run this once and save the results. Then you can iteratively adjust your choices when selecting anchor cells)

astats <- get_anchor_stats(counts = counts,
                           neg = negmean,
                           profiles = rectum_ref)
saveRDS(astats, "anchor_stats_insitutype.rds")

# now choose anchors:
anchors <- choose_anchors_from_stats(counts = counts, 
                                     neg = negmean, 
                                     bg = per.cell.bg,
                                     anchorstats = astats, 
                                     # a very low value chosen for the mini
                                     # dataset. Typically hundreds of cells
                                     # would be better.
                                     n_cells = 50, 
                                     min_cosine = 0.4, 
                                     min_scaled_llr = 0.03, 
                                     insufficient_anchors_thresh = 5)

# plot the anchors atop the UMAP:
par(mfrow = c(1, 1))
plot(fov_$UMAP.HARMONY, pch = 16, cex = 0.1, col = "peachpuff1", xaxt = "n",  yaxt = "n", xlab = "", ylab = "",
     main = "Selected anchor cells")
points(sce$UMAP.HARMONY[!is.na(anchors), ], col = iocolors[anchors[!is.na(anchors)]], pch = 16, cex = 0.6)
legend("topright", pch = 16, col = iocolors[setdiff(unique(anchors), NA)], legend = setdiff(unique(anchors), NA), cex = 0.65)


# 16.2.4- Semi-supervised cell annotation
updatedprofiles <- updateReferenceProfiles(reference_profiles = rectum_ref, 
                                           counts = counts, 
                                           neg = neg, 
                                           bg = per.cell.bg,
                                           anchors = anchors) 
str(updatedprofiles)

# OTHER OPTION TO RUN WITHOUT ANCHORS
updatedprofiles <- updateReferenceProfiles(reference_profiles = rectum_ref, 
                                           counts = counts, 
                                           neg = neg, 
                                           bg = per.cell.bg,
                                           anchors = NULL) 

# 16.2.5- Supervised cell annotation
# Run with rectum ref
sup <- insitutypeML(x = counts,
                    neg = negmean,
                    reference_profiles = rectum_ref)   

Insitutype_rectum <- sup$clust
write.csv(Insitutype_rectum, "insituclusters_supervised_rectum_ref.csv")

# Run with ileum ref
sup2 <- insitutypeML(x = counts,
                    neg = negmean,
                    reference_profiles = gut_ref_epithelial)   

Insitutype_gut <- sup2$clust
write.csv(Insitutype_gut, "insituclusters_supervised_gut_ref_epithelial.csv")

# Run with colon ref
sup3 <- insitutypeML(x = counts,
                     neg = negmean,
                     reference_profiles = rectum_gut_ref_epithelial)   

Insitutype_rectum_gut <- sup3$clust
write.csv(Insitutype_rectum_gut, "insituclusters_supervised_rectum_gut_ref_epithelial.csv")

sup4 <- insitutypeML(x = counts,
                     neg = negmean,
                     reference_profiles = colon_rectum_gut_ref)   

Insitutype_colon_rectum_gut <- sup4$clust
write.csv(Insitutype_colon_rectum_gut, "insituclusters_supervised_colon_rectum_gut_ref_epithelial.csv")

#fov_epithelial$Insitutype_gut <- NULL
#fov_epithelial$Insitutype_rectum <- NULL
#fov_epithelial$Insitutype_colon_rectum_gut <- NULL

fov_epithelial@meta.data <-cbind(fov_epithelial@meta.data, Insitutype_rectum_gut)
fov_epithelial@meta.data <-cbind(fov_epithelial@meta.data, Insitutype_colon_rectum_gut)

Idents(fov_epithelial) <- "Insitutype_colon_rectum_gut"
levels(fov_epithelial)
library(dittoSeq)
dittoDimPlot(fov_epithelial, var = "Insitutype_colon_rectum_gut",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
ImageDimPlot(fov_epithelial, fov = "fov02", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = sup3)
class(fp)
print(fp)

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_epithelial, var = "Insitutype_rectum_gut", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "Insitutype_rectum_gut", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "Insitutype_rectum_gut", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "Insitutype_rectum_gut"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_epithelial, features = features, sort = TRUE, ncol = 2)

# Updating clustering results
# After close review, you might want to update your clustering results.
# The function refineClusters can modify your results in the following ways:
# merge closely-related or frequently-confused clusters
# delete clusters (Sometimes a “catch-all” cluster arises that is frequently confused wiht multiple other clusters. It’s best to delete these clusters and let their cells get reassigned to meaningful clusters.)
# sub-cluster large clusters
# Here we use the refineClusters function to perform all these operations. We also use its “merges” argument to rename cluster c. 

# refine the clusters:
newclusts <- refineClusters(
  logliks = semisup$logliks,
  merges = c("fibroblast" = "stroma", "endothelial" = "stroma"),
  to_delete = c("a"),
  # subclustering via refineClusters is not recommended for semi-supervised
  # results
  subcluster = NULL,
  counts = counts,
  neg = negmean
) 
str(newclusts)


# 16.3- Cell Annotation with Garnett
library(garnett)
library(org.Hs.eg.db)

# 16.3.1- Build the CDS object - see Monocle documentation for help

# http://cole-trapnell-lab.github.io/monocle-release/docs/#the-celldataset-class
# https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/#cell_data_set


# The CellDataSet class:
# Because Garnett builds on Monocle 3, data for Garnett is held in objects of the cell_data_set (CDS) class. 
# This class is derived from the Bioconductor SingleCellExperiment class, which provides a common interface familiar to those who have analyzed single-cell data with Bioconductor. 
#Monocle 3 provides detailed documentation about how to generate an input CDS here.

# Monocle holds single cell expression data in objects of the CellDataSet class. 
# The class requires three input files:
# exprs, a numeric matrix of expression values, where rows are genes, and columns are cells
# phenoData, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
# featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.

# The expression value matrix must:
# have the same number of columns as the phenoData has rows.
# have the same number of rows as the featureData data frame has rows.

# Additionally:
# row names of the phenoData object should match the column names of the expression matrix.
# row names of the featureData object should match row names of the expression matrix.
# one of the columns of the featureData should be named "gene_short_name".
# this function was built for gene names, so convert protein IDs to gene Symbols from ENSEMBL

# Option 1 - convert Seurat object to SCE object
sce <- readRDS("./RDS_files/6_Gated/fov_epithelial_sce.rds")
saveRDS(sce, "./sce.rds")

expr_matrix <- counts(sce) #exprs matrix based on the raw RNA (UMI) counts. THIS IS IMPORTANT! DO NOT USE NORMALISED DATA!

# convert to sparse matrix class dgMatrix
# expr_matrix <- Matrix(expr_matrix, sparse = TRUE)
# writeMM(expr_matrix, "expr_matrix.mtx")
# write.table(expr_matrix,"expr_matrix.txt",sep="\t",row.names=FALSE)

pData <- colData(sce) # cells and their metadata
#write.table(pData,"pdata.txt",sep="\t",row.names=FALSE)
#pdata <- read.delim("pdata.txt")
#pdata <- as.data.frame(pdata)
#pdata <- new("AnnotatedDataFrame", data = pdata)

fData <- rownames(sce)
write.table(fData,"fdata.txt",sep="\t",row.names=FALSE)
fdata <- read.delim("fdata.txt") # metadata of proteins - change column name to gene_short_name in the original panel.csv file
fdata <- as.data.frame(fdata)
#fdata <- new("AnnotatedDataFrame", data = fdata)

colnames(expr_matrix) <- rownames(pData) #this is important
rownames(expr_matrix) <- fdata$gene_short_name
rownames(fdata) <- rownames(expr_matrix) #this is important

monocle_cds <- new_cell_data_set(expression_data = expr_matrix,
                                 cell_metadata = pData,
                                 gene_metadata = fdata)
saveRDS(monocle_cds, "fov_epithelial_cds.rds")
monocle_cds <- readRDS("./fov_epithelial_cds.rds")

# Option 2 - Extract data, phenotype data, and feature data from the SeuratObject - USE THIS APPROACH!!
# https://github.com/cole-trapnell-lab/monocle-release/issues/262

data <- as(as.matrix(fov_epithelial@assays$RNA@counts), 'sparseMatrix') # if you do have UMI data, you should not normalize it yourself prior to creating your CellDataSet. Use the raw RNA counts.
pd <- new('AnnotatedDataFrame', data = fov_epithelial@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds2 <- newCellDataSet(as(data, "dgCMatrix"),
                              phenoData = pd,
                              featureData = fd)
str(monocle_cds2) 

# generate size factors for normalization later
monocle_cds <- estimateSizeFactors(monocle_cds)

saveRDS(monocle_cds, "fov_epithelial_cds.rds")
monocle_cds <- readRDS("./fov_epithelial_cds.rds")


# 16.3.2- Automated annotation with Garnett 
# Reference: https://cole-trapnell-lab.github.io/monocle3/docs/clustering/#garnett
# https://cole-trapnell-lab.github.io/garnett/docs_m3/

# Constructing a marker file
# Create a Colon_markers.txt file based on the scRNAseq colon paper. 
# As done in Giotto, create a file resulted from the matching of the top expressed genes per cluster from Seurat with the top expressed genens per cluster listed by the paper.
# List of marker genes from the paper: E-MTAB-8410.marker_genes_inferred_cell_type_-_ontology_labels
# Create a list of marker genes per cluster from Seurat. Use assay = RNA and slot = data to calculate the expressed genes and select only the positively expressed genes. 
# Seurat marker list file: 
DefaultAssay(fov_epithelial) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_epithelial) <- 'ep_harmony_clusters'
levels(fov_epithelial)
fov_epithelial.marker <- FindAllMarkers(object = fov_epithelial, only.pos = TRUE, logfc.threshold = 0.25)
clipr::write_clip(fov_epithelial.marker) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# Match features from epithelial clusters to the 2020 reference marker gene list (2020_Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer)
markers_ref <- read.csv("./RDS_files/7_Reference_datasets/Garnett/2020reference/E-MTAB-8410.marker_genes_inferred_cell_type_-_ontology_labels.csv")
fov_epithelial.marker

ortho <- markers_ref
cluster <- fov_epithelial.marker
idx <- match(ortho$geneID, cluster$gene)
ortho$gene <- cluster$gene[ idx ]
ortho <- na.omit(ortho) 
write.csv(ortho, "marker_genes_2020reference_matched_epithelial_clusters.csv")

# Exclude genes that are good markers for more than one cell type:

garnett_markers <- ortho %>% 
  group_by(gene) %>%
  filter(n() == 1)

write.csv(garnett_markers, "good_marker_genes_2020reference_matched_epithelial_clusters.csv")

# 16.3.3- Defining and checking your markers for the classifier 

system.file(package = "garnett")
marker_file_path <- system.file("extdata", "Epithelial_2020ref_top5_markers.txt", #save marker files in this directory
                                package = "garnett")

# Because defining the marker file is often the hardest part of the process, 
# Garnett includes functions to check whether your markers are likely to work well. 
# The two functions relevant are check_markers and plot_markers. 
# Check_markers generates a table of information about your markers and plot_markers plots the most relevant information.
# In addition to the small included dataset, we have included two example marker files with the package. 

marker_check <- check_markers(cds = monocle_cds, marker_file_path, #check markers file in the Visual Studio Code if there are errors with typing
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
plot_markers(marker_check)

write.csv(marker_check, "Check_colon_cluster_v1.csv")

# Several genes showed high ambiguity. 
# Remove genes with high ambiguity or check in the matched reference which cell type its expression based on logFC is higher and use to define that particular cell type.  

# 16.3.4- Train the classifier
# Now it's time to train the classifier. The arguments should be pretty close to those for check_markers. 
# The one parameter I am changing from default below is the num_unknown argument. 
# This tells Garnett how many outgroup cells it should compare against. 
# The default is 500, but in this toy dataset with so few cells we want fewer.

set.seed(260)
classifier <- train_cell_classifier(cds = monocle_cds,
                                    marker_file = marker_file_path,
                                    db=org.Hs.eg.db,
                                    cds_gene_id_type = "SYMBOL", #convert protein symbols to geneID symbols in the original files generating the cds object
                                    num_unknown = 50,
                                    marker_file_gene_id_type = "SYMBOL", cores = 4)

saveRDS(classifier, "Garnett_classifier_epithelial_all_markers.rds")
#classifier <- readRDS("./Garnett_classifier_allcells_top5markers.rds")

# 16.3.5- Viewing the classification genes
#Garnett classification is trained using a multinomial elastic-net regression. 
#This means that certain genes are chosen as the relevant genes for distinguishing between cell types. 
#Which genes are chosen may be of interest, so Garnett includes a function to access the chosen genes. 
#Note: Garnett does not regularize the input markers, so they will be included in the classifier regardless.
#The function we use to see the relevant genes is get_feature_genes. 
#The arguments are the classifier, which node you'd like to view (if your tree is hierarchical) - use "root" for the top node and the parent cell type name for other nodes, and the db for your species.
#The function will automatically convert the gene IDs to SYMBOL if you set convert_ids = TRUE.

feature_genes <- get_feature_genes(classifier, 
                                   node = "root",
                                   db = org.Hs.eg.db, convert_ids = TRUE)
head(feature_genes)

# 16.3.6- Classifying cells

monocle_cds <- classify_cells(monocle_cds, classifier,
                              db = org.Hs.eg.db,
                              cluster_extend = TRUE,
                              cds_gene_id_type = "SYMBOL")

head(pData(monocle_cds))
table(pData(monocle_cds)$cell_type)
table(pData(monocle_cds)$cluster_ext_type)

write.csv(table(pData(monocle_cds)$cell_type), "Garnett_epithelial_clusters_top5markers.csv")
write.csv(table(pData(monocle_cds)$cell_type), "Garnett_epithelial_clusters_allmarkers.csv")


# 16.3.6-Store labels of classified cells in Spatial Experiment object

fov_epithelial$Garnett_epithelial_v1 <- pData(monocle_cds)$cell_type #clusters from top5 epithelial markers
table(fov_epithelial$Garnett_epithelial_v1)
fov_epithelial$Garnett_epithelial_v2 <- pData(monocle_cds)$cell_type #clusters from all epithelial markers
table(fov_epithelial$Garnett_epithelial_v2)

# Save the updated Seurat object.
saveRDS(fov_epithelial, "fov_epithelial.rds")

# 16.3.7- Visualisation of annotated clusters using Garnett
Idents(fov_epithelial) <- "Garnett_epithelial_v1"
levels(fov_epithelial)

dittoDimPlot(fov_epithelial, var = "Garnett_epithelial_v1",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
ImageDimPlot(fov_epithelial, fov = "fov02", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# 16.3.7- Constructing a marker file from the 2022_Spatially organized multicellular immune hubs in human colorectal cancer

monocle_cds <- readRDS("./fov_epithelial_cds.rds")

markers_ref <- read.csv("./RDS_files/7_Reference_datasets/Garnett/2022reference/2022_reference_epithelial.csv")
fov_epithelial.marker <- read.csv("./Plots/7_Clustering/Option 2/Epithelial/Harmony/Top_positive_marker_genes_epithelial_clusters.csv")

ortho <- markers_ref
cluster <- fov_epithelial.marker
idx <- match(ortho$geneID, cluster$gene)
ortho$gene <- cluster$gene[ idx ]
ortho <- na.omit(ortho) 
write.csv(ortho, "marker_genes_2022reference_matched_epithelial_clusters.csv")

# Exclude genes that are good markers for more than one cell type:

garnett_markers <- ortho %>% 
  group_by(gene) %>%
  filter(n() == 1)

write.csv(garnett_markers, "good_marker_genes_2022reference_matched_epithelial_clusters.csv")

# 16.3.7.1- Defining and checking your markers for the classifier 

system.file(package = "garnett")
marker_file_path <- system.file("extdata", "Epithelial_2022ref_good_markers.txt", #save marker files in this directory
                                package = "garnett")

marker_check <- check_markers(cds = monocle_cds, marker_file_path, #check markers file in the Visual Studio Code if there are errors with typing
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
plot_markers(marker_check)

set.seed(260)
classifier <- train_cell_classifier(cds = monocle_cds,
                                    marker_file = marker_file_path,
                                    db=org.Hs.eg.db,
                                    cds_gene_id_type = "SYMBOL", #convert protein symbols to geneID symbols in the original files generating the cds object
                                    num_unknown = 50,
                                    marker_file_gene_id_type = "SYMBOL", cores = 4)

saveRDS(classifier, "Garnett_classifier_epithelial_good_markers_2022ref.rds")
classifier <- readRDS("./Garnett_classifier_epithelial_all_markers_2022ref.rds")

# 16.3.7.2- Classifying cells

monocle_cds <- classify_cells(monocle_cds, classifier,
                              db = org.Hs.eg.db,
                              cluster_extend = TRUE,
                              cds_gene_id_type = "SYMBOL")

head(pData(monocle_cds))
table(pData(monocle_cds)$cell_type)
table(pData(monocle_cds)$cluster_ext_type)

write.csv(table(pData(monocle_cds)$cell_type), "Garnett_epithelial_clusters_allmarkers_2022ref.csv")
write.csv(table(pData(monocle_cds)$cell_type), "Garnett_epithelial_clusters_goodmarkers_2022ref.csv")

# 16.3.6-Store labels of classified cells in Spatial Experiment object

fov_epithelial$Garnett_epithelial_v3 <- pData(monocle_cds)$cell_type #clusters from all epithelial markers 2022ref

# Save the updated Seurat object.
saveRDS(fov_epithelial, "fov_epithelial.rds")

# 16.3.7- Visualisation of annotated clusters using Garnett
Idents(fov_epithelial) <- "Garnett_epithelial_v3"
levels(fov_epithelial)

dittoDimPlot(fov_epithelial, var = "Garnett_epithelial_v3",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
ImageDimPlot(fov_epithelial, fov = "fov02", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_epithelial, var = "Garnett_epithelial_v3", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "Garnett_epithelial_v3", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "Garnett_epithelial_v3", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "Garnett_epithelial_v3"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_epithelial, features = features, sort = TRUE, ncol = 2)


######################################################################################################################################################################################################################################################################################################################################################################################################################

# 17- Manual cell annotation of Seurat clusters based on the marker gene list of the reference dataset

Idents(fov_integrated) <- 'seurat_clusters'
levels(fov_integrated)

# 9.1- Rename idents of dsb clusters 
cluster_ID <- c("GNLY+ T cell", "ChemokineHighVEGFA+KRT15+ Epithelial cell", "Epithelial cell", "PIGR+ Epithelial cell", 
                "VEGFA+KRT15+ Epithelial cell", "Stromal cell", "Activated CD16+CD68+ Macrophage", "Tuft cell",
                "GNLY+ T cell", "SPP1+CCL2+ CD16+CD68+ Macrophage", "CXCL8+SOD2+ Neutrophil", "Plasma cell", "Plasma cell", "Plasma cell",
                "CXCL8+SOD2+ Neutrophil")

names(cluster_ID) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, cluster_ID)

# 9.2- Adding object metadata with cluster names in order - to be used in the bar plots
fov_integrated$cluster_ID <- factor(fov_integrated@active.ident)

# 9.3- Visualisation of manually annotated clusters

ImageDimPlot(fov_integrated, fov = "fov01", group.by = "cluster_ID", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov22", group.by = "cluster_ID", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
#Options for cols are "alphabet", "alphabet2", "glasbey", "polychrome", "stepped", and "parade".

# UMAP of the identified clusters
DimPlot(fov_integrated, reduction = "umap", label = FALSE, cols = "polychrome", pt.size = 2) #try other ploting packages and explore the function

dittoDimPlot(fov_integrated, var = "cluster_ID",
             reduction.use = "umap", size = 2,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Manually annotated clusters on UMAP")

# Plot cell numbers and frequencies by sample
p <- dittoBarPlot(fov_integrated, var = "cluster_ID", group.by = "sample_ID",  scale = c("percent")) 
p
dittoBarPlot(fov_integrated, var = "cluster_ID", group.by = "sample_ID",  scale = c("count")) 

# Manually define color code for each cluster label
cluster_ID_color <- setNames(c( "#023FA5", "grey",  "#FF00FF", "#F71F0F", "#9933FF",
                                "#FFD966", "#66D9FF", "#006633",  "#00FF00",
                                "#FF8000", "#606060"),
                             c("Epithelial cell", "PIGR+ Epithelial cell", "VEGFA+KRT15+ Epithelial cell", "ChemokineHighVEGFA+KRT15+ Epithelial cell", "Tuft cell",
                               "Stromal cell", "CXCL8+SOD2+ Neutrophil", "Activated CD16+CD68+ Macrophage", "SPP1+CCL2+ CD16+CD68+ Macrophage", 
                               "Plasma cell", "GNLY+ T cell"))

# Plot using the customised vector of colors above
ImageDimPlot(fov_integrated, fov = "fov01", group.by = "cluster_ID", axes = TRUE, cols = cluster_ID_color, size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov22", group.by = "cluster_ID", axes = TRUE, cols = cluster_ID_color, size = 1, boundaries = "segmentation")
# UMAP of the identified clusters
DimPlot(fov_integrated, reduction = "umap", label = FALSE, cols = cluster_ID_color , pt.size = 2)

# Plot specific cell types
# Epithelial compartment
p1 <- ImageDimPlot(fov_integrated, fov = "Primary_fov01", group.by = "cluster_ID", axes = TRUE, cols = cluster_ID_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("Epithelial cell", "PIGR+ Epithelial cell", "VEGFA+KRT15+ Epithelial cell", "ChemokineHighVEGFA+KRT15+ Epithelial cell", "Tuft cell")))
p2 <- ImageDimPlot(fov_integrated, fov = "Metastasis_fov22", group.by = "cluster_ID", axes = TRUE, cols = cluster_ID_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("Epithelial cell", "PIGR+ Epithelial cell", "VEGFA+KRT15+ Epithelial cell", "ChemokineHighVEGFA+KRT15+ Epithelial cell", "Tuft cell")))
p1
p2

# Stromal compartment
p3 <- ImageDimPlot(fov_integrated, fov = "Primary_fov01", group.by = "cluster_ID", axes = TRUE, 
                   cols = cluster_ID_color, size = 1, boundaries = "segmentation", 
                   cells = WhichCells(fov_integrated, idents = c("Stromal cell")))

p4 <- ImageDimPlot(fov_integrated, fov = "Metastasis_fov22", group.by = "cluster_ID", axes = TRUE, 
                   cols = cluster_ID_color, size = 1, boundaries = "segmentation", 
                   cells = WhichCells(fov_integrated, idents = c("Stromal cell")))
p3
p4

# Immune compartment
p5 <- ImageDimPlot(fov_integrated, fov = "Primary_fov01", group.by = "cluster_ID", axes = TRUE, 
                   cols = cluster_ID_color, size = 1, boundaries = "segmentation", 
                   cells = WhichCells(fov_integrated, idents = c("CXCL8+SOD2+ Neutrophil", "Activated CD16+CD68+ Macrophage",
                                                                 "SPP1+CCL2+ CD16+CD68+ Macrophage", "Plasma cell", "GNLY+ T cell")))

p6 <- ImageDimPlot(fov_integrated, fov = "Metastasis_fov22", group.by = "cluster_ID", axes = TRUE, 
                   cols = cluster_ID_color, size = 1, boundaries = "segmentation", 
                   cells = WhichCells(fov_integrated, idents = c("CXCL8+SOD2+ Neutrophil", "Activated CD16+CD68+ Macrophage",
                                                                 "SPP1+CCL2+ CD16+CD68+ Macrophage", "Plasma cell", "GNLY+ T cell")))
p5
p6


# See error when loading data and images to create Seurat object...
p7 <- ImageDimPlot(fov_integrated, fov = "fov22", group.by = "cluster_ID", axes = TRUE, 
                   cols = cluster_ID_color, size = 1, boundaries = "segmentation", 
                   cells = WhichCells(fov_integrated, idents = c("Epithelial cell", "PIGR+ Epithelial cell", "VEGFA+KRT15+ Epithelial cell", "ChemokineHighVEGFA+KRT15+ Epithelial cell", "Tuft cell")), 
                   alpha = 0.3, molecules = c("PIGR", "CEACAM6"), mols.size = 0.01, mols.alpha = 1, nmols = 50000, border.color = "black", coord.fixed = FALSE)

p7

################################################################################################################################################################################


