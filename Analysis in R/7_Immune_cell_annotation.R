
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

fov_integrated <- readRDS("./fov_integrated.rds")

saveRDS(fov_integrated, "./fov_integrated.rds")

sce <- readRDS("./RDS_files/6_Gated/fov_stromal_sce.rds")
saveRDS(sce, "fov_stromal_sce.rds")

fov_immune <- readRDS("./RDS_files/6_Gated/fov_immune.rds")

saveRDS(fov_immune, "./fov_immune.rds")


# 14- Cluster annotation with marker gene databases

# 14.1- Cell Annotation with scCATCH

# Reference: https://github.com/ZJUFanLab/scCATCH
# Reference: https://raw.githack.com/ZJUFanLab/scCATCH_performance_comparison/master/scCATCH/tutorial.html
# This package is a single cell Cluster-based auto-Annotation Toolkit for Cellular Heterogeneity (scCATCH) from cluster potential marker genes identification to cluster annotation based on evidence-based score 
# by matching the potential marker genes with known cell markers in tissue-specific cell taxonomy reference database (CellMatch).
library(scCATCH)

set.seed(1234)
DefaultAssay(fov_stromal) <- "RNA"
Idents(fov_stromal) <- 'st_harmony_clusters'
levels(fov_stromal)

# Create scCATCH object from Seurat object:
obj <- createscCATCH(data = fov_stromal[['RNA']]@data, cluster = as.character(Idents(fov_stromal)))
saveRDS(obj, "fov_stromal_scCATCH_object.rds")

# 14.1.1- Find marker gene for each cluster

# Users need to provided the species, tissue, or cancer information. 
# Use similar paramaters to those used with the FindAllMarkers function from Seurat
# Available tissues and cancers at https://github.com/ZJUFanLab/scCATCH/wiki
# Select different combination of tissues or cancers for annotation
# cellmatch$tissue
# filter cellmatch
#cellmatch <- cellmatch[cellmatch$species == "Human", ]
#cellmatch <- cellmatch[cellmatch$tissue %in% c("Blood", "Peripheral blood", "Serum", "Colon", "Colorectum", "Intestine"), ]
#cellmatch <- cellmatch[cellmatch$cancer %in% c("Colon Cancer", "Colorectal Cancer"), ]

obj <- findmarkergene(obj, species = 'Human', cluster = "All", if_use_custom_marker = FALSE, marker = cellmatch, 
                      tissue = c("Colon", "Colorectum", "Intestine", ),
                      cancer = c("Colorectal Cancer"),
                      use_method = "2", cell_min_pct = 0.1,logfc = 0.25, pvalue = 0.05, verbose = TRUE)

# 14.1.2- Find cell type for each cluster
obj <- findcelltype(obj)

list <- obj@celltype
write.csv(list, "scCatch_annotated_stromal_clusters.csv")

# scCatch cluster annotation and visualisation

SCATCH_clusters <- c("Cancer-associated fibroblast", "T cell", "Cancer stem cell", "SPP1+MMP12+ macrophage", "Endothelial cell",
                     "Cancer-associated fibroblast", "Cancer stem cell", "SPP1+CXCL8+CXCL1+ monocyte", "B cell",
                     "Cancer-associated fibroblast", "Cancer-associated fibroblast", "Cancer-associated fibroblast", "Cancer stem cell", "S100P+ epithelial cell",
                     "B cell", "Cancer-associated fibroblast", "B cell", "SPP1+ monocyte", "Cancer stem cell", "Macrophage", "Cancer stem cell", "CXCL8+CXCL1+ monocyte",
                     'Cancer stem cell')
names(SCATCH_clusters) <- levels(fov_stromal)
fov_stromal <- RenameIdents(fov_stromal, SCATCH_clusters)

# Adding object metadata with
fov_stromal$SCATCH_clusters <- factor(fov_stromal@active.ident)

dittoDimPlot(fov_stromal, var = "SCATCH_clusters",
                   reduction.use = "umap.harmony", size = 1,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("scCATCH epithelial clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
Idents(fov_stromal) <- "SCATCH_clusters"
ImageDimPlot(fov_stromal, fov = "fov02", group.by = "SCATCH_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov23", group.by = "SCATCH_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_stromal, var = "SCATCH_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "SCATCH_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "SCATCH_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_stromal) <- "Protein"
Idents(fov_stromal) <- "SCATCH_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_stromal, features = features, sort = TRUE, ncol = 2)

# 14.2- Cell Annotation with SCSA

# Reference: https://github.com/bioinfo-ibms-pumc/SCSA
# To annotate a human scRNA-seq sets generated by 'FindAllMarkers' function of Seurat with ensemblIDs, use the following code:
# -d DB, --db DB        Database for annotation. (whole.db)
# -s SOURCE, --source Source of marker genes. (cellranger,[seurat],[scanpy], [scran])
# -i INPUT, --input INPUT Input file for marker annotation(Only CSV format supported).
# -k TISSUE, --tissue Tissue for annotation. Only used for cellmarker database. Multiple tissues should be seperated by commas.Run '-l' option to see all tissues.
# In linux platform:('All',['Bone marrow'],['Bone marrow,Brain,Blood'][...])
# python3 SCSA.py -d whole_v2.db -s seurat -i Markers_stromal_hamorny_clusters.csv -k All -E -g Human -p 0.05 -f 1.5

# Add SCSA cluster lables
Idents(fov_stromal) <- "st_harmony_clusters"
levels(fov_stromal)

fov_stromal$SCSA_clusters <- NULL

SCSA_clusters <- c("Fibroblast", "Epithelial cell", "Epithelial cell", "SPP1+MMP12+ Macrophage", "Endothelial cell",
                   "Fibroblast", "Epithelial cell", "SPP1+CXCL8+CXCL1+ monocyte", "B cell",
                   "Smooth muscle cell", "Fibroblast", "Fibroblast", "Monocyte", "Epithelial cell",
                   "Plasma cell", "Fibroblast", "B cell", "SPP1+ monocyte", "Epithelial cell", "Mast cell", "Epithelial cell", "Neutrophil",
                   'Plasma cell') 
names(SCSA_clusters) <- levels(fov_stromal)
fov_stromal <- RenameIdents(fov_stromal, SCSA_clusters)

# Adding object metadata with
fov_stromal$SCSA_clusters <- factor(fov_stromal@active.ident)

dittoDimPlot(fov_stromal, var = "SCSA_clusters",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("SCSA epithelial clusters on UMAP")

multi_dittoDimPlot(fov_stromal, vars = c("SCATCH_clusters","SCSA_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
Idents(fov_stromal) <- "SCSA_clusters"
ImageDimPlot(fov_stromal, fov = "fov02", group.by = c("SCSA_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov23", group.by = c("SCSA_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_stromal, var = "SCSA_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "SCSA_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "SCSA_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_stromal) <- "Protein"
Idents(fov_stromal) <- "SCSA_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_stromal, features = features, sort = TRUE, ncol = 2)

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
library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=FALSE)
saveRDS(ref.data, "HumanPrimaryCellAtlasData.rds")
ref.data <- readRDS("./RDS_files/7_Reference_datasets/HumanPrimaryCellAtlasData/HumanPrimaryCellAtlasData.rds")

blueprint <- BlueprintEncodeData(ensembl=FALSE)
saveRDS(hemato, "BlueprintEncodeData.rds")
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

DefaultAssay(fov_stromal) <- "RNA" 
Idents(fov_stromal) <- 'st_harmony_clusters'
levels(fov_stromal)
sce <- as.SingleCellExperiment(fov_stromal, assay = "RNA")
saveRDS(sce, "fov_stromal_sce.rds")

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
saveRDS(pred, "./fov_stromal_SingleR.rds")
pred <- readRDS("./fov_stromal_SingleR.rds")

# 15.1.3- Using the Blueprintenconde ref

# main labels
pred3 <- SingleR(test = sce, ref = ref.data2, 
                 labels = ref.data2$label.main, assay.type.test=1)
colnames(pred3)
saveRDS(pred3, "./fov_stromal_SingleR_pred3.rds")
pred3 <- readRDS("./fov_stromal_SingleR_pred3.rds")


# 15.2- Cell Annotation with CIPR

# Reference: https://github.com/atakanekiz/CIPR-Shiny
library(CIPR)
DefaultAssay(fov_immune) <- "RNA" 
Idents(fov_immune) <- 'im_harmony_clusters'
levels(fov_immune)

# 15.2.1- Calculate average gene expression per cluster
# This is the step where we generate the input for CIPR's all-genes correlation methods
avgexp <- AverageExpression(fov_immune)
avgexp <- as.data.frame(avgexp$RNA)
avgexp$gene <- rownames(avgexp)
write.csv(avgexp, "avgexp_immune.csv")

# CIPR analysis

# 15.2.2- Standard logFC comparison method
# Immunological Genome Project (ImmGen) microarray data from sorted mouse immune cells. This dataset is prepared by using both V1 and V2 ImmGen releases and it contains 296 samples from 20 different cell types (253 subtypes).
# Mouse RNAseq data from sorted cells reported in Benayoun et al. (2019). This dataset contains 358 sorted immune and nonimmune samples from 18 different lineages (28 subtypes).
# Blueprint/Encode RNAseq data that contains 259 sorted human immune and nonimmune samples from 24 different lineages (43 subtypes).
# Human Primary Cell Atlas (hpca) that contains microarray data from 713 sorted immune and nonimmune cells (37 main cell types and 157 subtypes).
# DICE (Database for Immune Cell Expression(/eQTLs/Epigenomics) that contains 1561 human immune samples from 5 main cell types (15 subtypes). To reduce object sizes, mean TPM values per cell type is used.
# Human microarray data from sorted hematopoietic cells reported in Novershtern et al. (2011). This dataset contains data from 211 samples and 17 main cell types (38 subtypes)
# Human RNAseq (hsrnaseq) data from sorted cells reported in Monaco et al. (2019). This dataset contains 114 samples originating from 11 main cell types (29 subtypes)
# A custom reference dataset provided by the user. This dataset can be obtained from a number of high througput methods including microarray and bulk RNAseq. For details about how to prepare custom reference, please see the How-to tab on the Shiny website.
        
allmarkers <- read.csv("./Plots/7_Clustering/Option 2/Immune/Harmony/Markers_harmony_immune_clusters_RNA.csv")

CIPR(input_dat = allmarkers,
     comp_method = "logfc_dot_product", 
     reference = "hpca", 
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
#write.csv(CIPR_top_results, "CIPR_top_results_avgexp_correlation.csv")

# 15.2.3- Standard all-genes correlation method - SKIP THIS
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
Idents(fov_immune) <- "im_harmony_clusters"
levels(fov_immune)

fov_immune$CIPR_hpca_clusters <- NULL
CIPR_hpca_clusters <- c("Fibroblast", "Epithelial cell", "Epithelial cell", "SPP1+MMP12+ Macrophage", "Endothelial cell",
                        "Fibroblast", "Epithelial cell", "SPP1+CXCL8+CXCL1+ monocyte", "B cell",
                        "Endothelial cell", "Fibroblast", "Fibroblast", "Epithelial cell", "Epithelial cell",
                        "Plasma cell", "Fibroblast", "B cell", "SPP1+ DC", "Epithelial cell", "Mast cell", "Epithelial cell", "Neutrophil",
                        'Plasma cell') 
names(CIPR_hpca_clusters) <- levels(fov_immune)
fov_immune <- RenameIdents(fov_immune, CIPR_hpca_clusters)

# Adding object metadata with
fov_stromal$CIPR_hpca_clusters <- factor(fov_stromal@active.ident)

dittoDimPlot(fov_stromal, var = "CIPR_hpca_clusters",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

multi_dittoDimPlot(fov_stromal, vars = c("SCATCH_clusters","SCSA_clusters", "CIPR_hpca_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

# Visualise labelled clusters in spatial context
Idents(fov_stromal) <- "CIPR_hpca_clusters"
ImageDimPlot(fov_stromal, fov = "fov02", group.by = "CIPR_hpca_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov23", group.by = "CIPR_hpca_clusters", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_stromal, var = "CIPR_hpca_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "CIPR_hpca_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "CIPR_hpca_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_stromal) <- "Protein"
Idents(fov_stromal) <- "CIPR_hpca_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_stromal, features = features, sort = TRUE, ncol = 2)

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
gut_ref_immune <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/gut_ref_immune.csv")
colon_ref <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/colon_ref.csv")
rectum_ref <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/rectum_ref.csv")

merged <- merge(gut_ref_immune, colon_ref, by.x = "Gene", by.y = "Gene")
merged2 <- merge(merged, rectum_ref, by.x = "Gene", by.y = "Gene")# gave the best performance for the annotation
write.csv(merged2, "./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/colon_rectum_gut_ref_immune.csv")
colon_rectum_gut_ref_immune <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CIPR/colon_rectum_gut_ref_immune.csv")

# Annotation with CIPR
allmarkers <- read.csv("./Plots/7_Clustering/Option 2/Immune/Harmony/Markers_harmony_immune_clusters_RNA.csv")
allmarkers <- as.data.frame(allmarkers)
avgexp <- read.csv("./Plots/8_Annotation/Immune/avgexp_immune.csv", sep=",", row.names=1, header = TRUE)

# correlation method
CIPR(input_dat = avgexp,
     comp_method = "all_genes_spearman", 
     reference = "custom",
     custom_reference = gut_ref_immune,  
     plot_ind = F,
     plot_top = T, 
     global_results_obj = T, 
     global_plot_obj = F)

write.csv(CIPR_top_results, "CIPR_top_results_colon_rectum_gut_ref_immune.csv")

# 15.2.4- Add CIPR cluster lables
Idents(fov_immune) <- "im_harmony_clusters"
levels(fov_immune)

fov_immune.markers <- FindAllMarkers(fov_immune, only.pos = TRUE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_immune.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

macrophage.markers <- FindMarkers(fov_immune, ident.1 = "2", ident.2 = "19", only.pos = TRUE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(macrophage.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

#fov_stromal$CIPR_colon_clusters <- NULL
fov_immune$CIPR_gut_clusters <- NULL
#fov_stromal$CIPR_colon_gut_clusters <- NULL

CIPR_gut_clusters <- c("CD8+ T cell", "Tumor epithelial cell", "SPP1+ macrophage", "Fibroblast", "Fibroblast", "Pericyte", "B cell", "B cell",
                       "Plasma cell", "ChemokinesHigh monocyte", "Plasma cell", "Plasma cell", "SPP1+ macrophage",
                       "Goblet cell", "Plasma cell", "CD4+ T cell", "Plasma cell", "CD4+ T cell", "Goblet cell", "ChemokinesHigh macrophage", "Fibroblast", "SPP1+ monocyte")  
names(CIPR_gut_clusters) <- levels(fov_immune)
fov_immune <- RenameIdents(fov_immune, CIPR_gut_clusters)
# Adding object metadata with
fov_immune$CIPR_gut_clusters <- factor(fov_immune@active.ident)

dittoDimPlot(fov_immune, var = "CIPR_gut_clusters",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Immune clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
Idents(fov_immune) <- "CIPR_gut_clusters"
ImageDimPlot(fov_immune, fov = "fov02", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_immune, fov = "fov23", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_immune, var = "CIPR_gut_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_immune, var = "CIPR_gut_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_immune, var = "CIPR_gut_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_immune) <- "Protein"
Idents(fov_immune) <- "CIPR_gut_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_immune, features = features, sort = TRUE, ncol = 2)

DefaultAssay(fov_immune) <- "RNA" 
fov_immune.markers <- FindAllMarkers(fov_immune, only.pos = TRUE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_immune.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_immune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap((subset(fov_immune, downsample = 1000)), features = top10$gene, size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 8))

dittoHeatmap(fov_immune, top10$gene, 
             cluster_rows = FALSE, 
             annot.by = "CIPR_gut_clusters", 
             fontsize_row = 8)

avgexp <- AverageExpression(subset(fov_immune, features = top10$gene))
avgexp <- as.data.frame(avgexp$RNA)
write.csv(avgexp, "Average_expression_Top_expressed_genes_CIPR_immune_clusters_RNA.csv")
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

DefaultAssay(fov_stromal) <- "RNA" 
Idents(fov_stromal) <- 'st_harmony_clusters'
levels(fov_stromal)
sce <- as.SingleCellExperiment(fov_stromal, assay = "RNA")
saveRDS(sce, "fov_stromal_sce.rds")

# Customised reference data:
colon_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/colon_ref_chetah.rds")
gut_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/gut_ref_chetah.rds")
colon_gut_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/CHETAH/colon_gut_ref_chetah.rds")

# Running CHETAH
sce <- CHETAHclassifier(input = sce, ref_cells = gut_ref, thresh = 0.1) # it took 2h to run...

# output
# CHETAH returns the input object, but added: input$celltype_CHETAH - a named character vector that can directly be used in any other workflow/method.
# “hidden” int_colData and int_metadata, not meant for direct interaction, but which can all be viewed and interacted with using: PlotCHETAH and CHETAHshiny
PlotCHETAH(input = sce, redD = "UMAP.HARMONY")

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
# Customised colon reference data:

# Customised gut reference data:
gut_ref_immune <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_immune.csv", sep=",", row.names=1)
gut_ref_immune <- as.matrix(gut_ref_immune)
class(gut_ref_immune)
str(gut_ref_immune)
saveRDS(gut_ref_immune, "./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_immune_insitutype.rds")
gut_ref_immune <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_immune_insitutype.rds")

# Customised colon+rectum+gut immune:
ref1 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/gut_ref_immune.csv")
ref2 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/colon_ref.csv")
ref3 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/rectum_ref.csv")
merged <- merge(ref1, ref2, by.x = "Gene", by.y = "Gene")
merged2 <- merge(merged, ref3, by.x = "Gene", by.y = "Gene")# gave the best performance for the annotation
write.csv(x = merged2, file = "./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/colon_rectum_gut_ref_immune.csv", row.names = TRUE, quote = FALSE)
ref4 <- read.csv("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/colon_rectum_gut_ref_immune.csv", sep=",", row.names=1)
ref5 <- as.matrix(ref4)
str(ref5)
saveRDS(ref5, "./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/colon_rectum_gut_ref_immune_insitutype.rds")
colon_rectum_gut_ref <- readRDS("./RDS_files/7_Reference_datasets/CellProfileLibrary/Human/Insitutype/colon_rectum_gut_ref_immune_insitutype.rds")

# 16.2.2- Extract the count matrix with the negative probe expression and calculate cells background

fov_subset <- readRDS("./RDS_files/4_Subset/fov_subset_neg_probes.rds")
Idents(fov_subset) <- "IF_type"
levels(fov_subset)
fov_subset <- subset(fov_subset, idents = c("Stromal", "CD3 Lymphocyte"))

counts <- LayerData(fov_subset, assay = "RNA")
saveRDS(counts, "immune_all_counts.rds")
neg <- counts[(which(rownames(counts) %in% c('NegPrb3',	'NegPrb5', 'NegPrb6',	'NegPrb7', 'NegPrb8',	'NegPrb9', 'NegPrb10',	'NegPrb11',	'NegPrb12',	'NegPrb13',	'NegPrb14',	'NegPrb15',	'NegPrb16', 'NegPrb18',	'NegPrb19',	'NegPrb20',	'NegPrb21',	'NegPrb22',	'NegPrb23'))),]
saveRDS(neg, "immune_negative_counts.rds")

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
# 1. Identify “anchor cells”: cells that can confidently be assigned to a reference cell type. Anchor cells can be hand-selected based on marker genes or morphology. But for most cases, they will be selected using the data-driven approach implemented below. 
# 2. Estimate the mean profile of each cell type’s anchors. These mean profiles form the updated reference matrix.
# This step is slow for big datasets. It's recommended to run this once and save the results. Then you can iteratively adjust your choices when selecting anchor cells)

astats <- get_anchor_stats(counts = counts,
                           neg = negmean,
                           profiles = colon_ref)
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
                    reference_profiles = gut_ref_immune)   

Insitutype_gut <- sup$clust
write.csv(Insitutype_gut, "insituclusters_supervised_gut_ref_immune.csv")

# Run with ileum ref
sup2 <- insitutypeML(x = counts,
                    neg = negmean,
                    reference_profiles = colon_rectum_gut_ref)   

Insitutype_colon_rectum_gut <- sup2$clust
write.csv(Insitutype_colon_rectum_gut, "insituclusters_supervised_colon_rectum_gut_ref_immune.csv")

#fov_stromal$Insitutype_colon <- NULL
#fov_stromal$Insitutype_gut <- NULL
#fov_stromal$Insitutype_gut2 <- NULL
#fov_stromal$Insitutype_colon_gut <- NULL

DefaultAssay(fov_immune) <- "RNA" 
Idents(fov_immune) <- 'im_harmony_clusters'
levels(fov_immune)
fov_immune@meta.data <-cbind(fov_immune@meta.data, Insitutype_gut)
fov_immune@meta.data <-cbind(fov_immune@meta.data, Insitutype_colon_rectum_gut)

Idents(fov_immune) <- "Insitutype_colon_rectum_gut"
levels(fov_immune)

dittoDimPlot(fov_immune, var = "Insitutype_colon_rectum_gut",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Clusters on UMAP")

# Visualise SCSA labelled clusters in a spatial context
ImageDimPlot(fov_immune, fov = "fov02", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_immune, fov = "fov23", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

fp <- flightpath_plot(flightpath_result = NULL, insitutype_result = sup)
class(fp)
print(fp)

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_immune, var = "Insitutype_gut", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_immune, var = "Insitutype_gut", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_immune, var = "Insitutype_gut", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_immune) <- "Protein"
Idents(fov_immune) <- "Insitutype_gut"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

RidgePlot(fov_immune, features = features, sort = TRUE, ncol = 2)

multi_dittoDimPlot(fov_immune, vars = c("CIPR_gut_clusters", "Insitutype_gut"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")


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
p1 <- ImageDimPlot(fov_integrated, fov = "fov01", group.by = "cluster_ID", axes = TRUE, cols = cluster_ID_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("Epithelial cell", "PIGR+ Epithelial cell", "VEGFA+KRT15+ Epithelial cell", "ChemokineHighVEGFA+KRT15+ Epithelial cell", "Tuft cell")))
p2 <- ImageDimPlot(fov_integrated, fov = "fov02", group.by = "cluster_ID", axes = TRUE, cols = cluster_ID_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("Epithelial cell", "PIGR+ Epithelial cell", "VEGFA+KRT15+ Epithelial cell", "ChemokineHighVEGFA+KRT15+ Epithelial cell", "Tuft cell")))
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

# 18- Merging CIPR cluster annotated cell types

# reference: https://github.com/satijalab/seurat/issues/1748

fov_integrated <- readRDS("./RDS_files/5_Integration/All/fov_integrated.rds")
fov_epithelial <- readRDS("./RDS_files/6_Gated/fov_epithelial.rds")
fov_stromal <- readRDS("./RDS_files/6_Gated/fov_stromal.rds")
fov_immune <- readRDS("./RDS_files/6_Gated/fov_immune.rds")

#saveRDS(fov_integrated, "./fov_integrated.rds")
DefaultAssay(fov_integrated) <- "RNA"
Idents(fov_integrated) <- 'harmony_clusters'
levels(fov_integrated)

DefaultAssay(fov_epithelial) <- "RNA"
Idents(fov_epithelial) <- 'CIPR_rectum_gut_clusters'
levels(fov_epithelial)

DefaultAssay(fov_stromal) <- "RNA"
Idents(fov_stromal) <- 'CIPR_colon_gut_clusters'
levels(fov_stromal)

DefaultAssay(fov_immune) <- "RNA"
Idents(fov_immune) <- 'CIPR_gut_clusters'
levels(fov_immune)

# Generate a new column called sub_cluster in the metadata
fov_integrated$CIPR_cluster_annotation <- as.character(Idents(fov_integrated))

fov_integrated$CIPR_cluster_annotation[Cells(fov_epithelial)] <- paste(Idents(fov_epithelial))
fov_integrated$CIPR_cluster_annotation[Cells(fov_stromal)] <- paste(Idents(fov_stromal))
fov_integrated$CIPR_cluster_annotation[Cells(fov_immune)] <- paste(Idents(fov_immune))

#fov_integrated$CIPR_cluster_annotation <- NULL
#fov_integrated$CIPR_annotated_clusters <- NULL

Idents(fov_integrated) <- 'CIPR_cluster_annotation'
levels(fov_integrated)

CIPR_annotated_clusters <- c("Fibroblast", "Plasma cell", "Pericyte", "SPP1+ macrophage", "Treg", "CD8+ T cell", "Goblet cell", 
                             "CD4+ T cell", "B cell", "Endothelial cell", "Tumor epithelial cell", "Myofibroblast", "SPP1+ monocyte", "ChemokinesHigh monocyte",
                             "Tuft cell", "IGF2+AREG+ Tumor epithelial cell",  "IGF2+AREG+ stem cell", "Tumor epithelial cell", "ChemokinesHigh Tumor epithelial cell", 
                             "Stem cell", "Neutrophil", "ChemokinesHigh tuft cell", "ChemokinesHigh macrophage")

names(CIPR_annotated_clusters) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, CIPR_annotated_clusters)


fov_integrated$CIPR_annotated_clusters <- factor(fov_integrated@active.ident)

Idents(fov_integrated) <- 'CIPR_annotated_clusters'
levels(fov_integrated)

dittoDimPlot(fov_integrated, var = "CIPR_annotated_clusters",
             reduction.use = "umap.harmony", size = 0.5,
             do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Clusters on UMAP")

Idents(fov_integrated) <- 'CIPR_annotated_clusters'

# Visualise in a spatial context
ImageDimPlot(fov_integrated, fov = "fov02", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov23", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov10", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov19", axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# Plot frequencies of SCSA labelled clusters per sample
dittoBarPlot(fov_integrated, var = "CIPR_annotated_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_integrated, var = "CIPR_annotated_clusters", group.by = "CMS",  scale = c("percent")) 
dittoBarPlot(fov_integrated, var = "CIPR_annotated_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 

# Plot normalised protein expression per cluster
DefaultAssay(fov_integrated) <- "Protein"
Idents(fov_integrated) <- "CIPR_annotated_clusters"
features <- c("MeanPanCK","MaxPanCK")
features2 <- c("MeanCD45", "MaxCD45")
features3 <- c("MeanCD3", "MaxCD3")

RidgePlot(fov_integrated, features = features, sort = TRUE, ncol = 2)
RidgePlot(fov_integrated, features = features2, sort = TRUE, ncol = 2)
RidgePlot(fov_integrated, features = features3, sort = TRUE, ncol = 2)

DefaultAssay(fov_integrated) <- "RNA" 
fov_integrated.markers <- FindAllMarkers(fov_integrated, only.pos = TRUE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_integrated.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

dittoHeatmap(fov_integrated, top10$gene, 
             cluster_rows = FALSE, 
             annot.by = "CIPR_annotated_clusters", 
             fontsize_row = 8)

avgexp <- AverageExpression(subset(fov_integrated, features = top10$gene))
avgexp <- as.data.frame(avgexp$RNA)
write.csv(avgexp, "Average_expression_Top_expressed_genes_CIPR_all_clusters_RNA.csv")

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

ImageFeaturePlot(fov_integrated, fov = "fov02", features = c('CD3E',"TNFRSF18", "FOXP3"), max.cutoff = "q95", boundaries = "segmentation")
ImageFeaturePlot(fov_integrated, fov = "fov23", features = c('CD3E',"TNFRSF18", "FOXP3"), max.cutoff = "q95", boundaries = "segmentation")

##############################################################################################################################################################################################################################################


# 19- Aesthetics correction

# 19.1- Manually define color code for each cluster label
cell_type_color <- setNames(c("#023FA5", '#66D9FF', '#8595e1', 
                              '#783F04', "#BF9000", 
                              "#A381EF", "#4d7191",'#b5bbe3', 
                              "#CCCCCC", '#444444', "#34E334", "#006633", "#c0cf0f",
                              '#FF9900', '#FFD966', '#FFF2CC', 'yellow',
                              "#FF0000" ,"#FF00FF",  '#EF8ECC', "#BF0A3D",   
                              '#850D85'
                              ),
                            c("Tumor epithelial cell", "ChemokinesHigh Tumor epithelial cell", "IGF2+AREG+ Tumor epithelial cell", 
                              "Stem cell", "IGF2+AREG+ stem cell",
                              "Goblet cell", "Tuft cell", "ChemokinesHigh tuft cell",
                              "Plasma cell", 'B cell', "CD8+ T cell", "CD4+ T cell", "Treg",
                              'Fibroblast', 'Pericyte', "Endothelial cell", 'Myofibroblast',
                              'SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage',
                              'ChemokinesHigh monocyte'))
#"#006633", 

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



# 19.2- Reorder the order of clusters for the bar graphs
levels(fov_integrated)
fov_integrated$CIPR_ordered_clusters <- factor(fov_integrated$CIPR_annotated_clusters,levels=c("ChemokinesHigh macrophage", "SPP1+ macrophage", "ChemokinesHigh monocyte", "SPP1+ monocyte", "Neutrophil",
                                                                                               "Treg", "CD8+ T cell", "CD4+ T cell", "B cell", "Plasma cell",
                                                                                               "Fibroblast", "Pericyte", "Endothelial cell", "Myofibroblast", 
                                                                                               "IGF2+AREG+ stem cell", "Stem cell", "ChemokinesHigh tuft cell", "Tuft cell", "Goblet cell", 
                                                                                               "IGF2+AREG+ Tumor epithelial cell", "ChemokinesHigh Tumor epithelial cell", "Tumor epithelial cell"))
Idents(fov_integrated) <- 'CIPR_ordered_clusters'
levels(fov_integrated)


#If plotting metaclusters...

Idents(fov_integrated) <- 'metaclusters_ordered'
levels(fov_integrated)

metaclusters <- c("Myeloid", "T lymphoid", "B lymphoid", "Stromal", "Epithelial")

names(metaclusters) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, metaclusters)
fov_integrated$metaclusters <- factor(fov_integrated@active.ident)

Idents(fov_integrated) <- 'metaclusters'
levels(fov_integrated)

metacluster_order <- match(levels(fov_integrated@meta.data[['metaclusters_ordered']]), metaLevels('metaclusters_ordered', fov_integrated))
metacluster_color <- setNames(c("#00924C", '#FDC77F', '#013F89', "#7853FF", '#AD0791'), c("Epithelial", "Stromal", "Myeloid",  "T lymphoid",  "B lymphoid"))


# 19.3- Define order of clusters for ditto bar plots
cluster_order <- match(levels(fov_integrated@meta.data[['CIPR_ordered_clusters']]), metaLevels('CIPR_ordered_clusters', fov_integrated))

dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "sample_ID",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "sample_ID",  scale = c("count"), color.panel = cell_type_color, var.labels.reorder = cluster_order)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "CMS",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "MMR_Status_Jen", scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order)

dittoDimPlot(fov_integrated, var = "CIPR_ordered_clusters",
             reduction.use = "umap.harmony", size = 0.5, color.panel = cell_type_color,
             do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Clusters on UMAP")

dittoDimPlot(fov_integrated, var = "metaclusters",
             reduction.use = "umap.harmony", size = 0.5, color.panel = metacluster_color,
             do.label = TRUE, labels.size = 8, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Clusters on UMAP")

# 19.4- Visualise clusters in a spatial context
ImageDimPlot(fov_integrated, fov = "fov01", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov02", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov03", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov04", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov05", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov06", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov07", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov08", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov09", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov10", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov11", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov12", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov13", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov14", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov15", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov16", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov17", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov18", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov19", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov20", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov21", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov22", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov23", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov24", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)
ImageDimPlot(fov_integrated, fov = "fov25", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", dark.background = FALSE)

# 19.5- Plot specific cell types
# Tumour epithelial cells + T cells or TME + Myeloid cells
p1 <- ImageDimPlot(fov_integrated, fov = "fov01", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov01", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov02", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov02", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2 
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov03", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov03", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov04", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov04", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov05", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov05", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov06", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov06", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov07", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov07", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov08", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov08", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2 
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov09", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov09", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov10", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov10", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov11", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov11", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov12", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov12", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov13", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov13", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov14", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov14", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov15", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov15", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov16", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov16", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov17", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov17", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov18", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov18", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov19", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov19", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov20", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov20", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov21", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov21", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov22", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov22", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov23", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov23", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p3 <- ImageDimPlot(fov_integrated, fov = "fov24", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p4 <- ImageDimPlot(fov_integrated, fov = "fov24", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2
p3/p4

p1 <- ImageDimPlot(fov_integrated, fov = "fov25", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c("CD8+ T cell", "CD4+ T cell", "Treg")), dark.background = FALSE)
p2 <- ImageDimPlot(fov_integrated, fov = "fov25", group.by = "CIPR_ordered_clusters", axes = TRUE, cols = cell_type_color, size = 1, boundaries = "segmentation", cells = WhichCells(fov_integrated, idents = c('SPP1+ macrophage', 'SPP1+ monocyte', 'Neutrophil', 'ChemokinesHigh macrophage', 'ChemokinesHigh monocyte')), dark.background = FALSE)
p1/p2


table(fov_integrated@meta.data$CIPR_ordered_clusters)


################################################################################################################################################################################\

# 20- Creating AnnData from Seurat

library(reticulate)
use_condaenv("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate")
use_python("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/bin/python")

library(sceasy)

py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)
reticulate::py_install(py_pkgs)
reticulate::import("scanpy")
sc <- import("scanpy")

BM_query <- readRDS("./BM_merged_query.rds")
DefaultAssay(BM_query) <- "RNA"
Idents(BM_query) <- 'SingleR.labels'
levels(BM_query)

X   = t(GetAssayData(fov_integrated)) # log normalised RNA data slot
write.csv(X, "fov_integrated_exprMat_file.csv")

obs = fov_integrated[[]] # metadata file
write.csv(obs, "fov_integrated_metadata_file.csv")

var = GetAssay(fov_integrated)[['data']]

ad_fov_integrated <- sc$AnnData(
  X   = t(GetAssayData(BM_query)),
  obs = BM_query[[]]
)

ad_fov_integrated$obsm$update(umap = Embeddings(BM_query, "proj.umap"))

ad_fov_integrated

sc$AnnData$write_h5ad(ad_fov_integrated, 'BM_merged_query.h5ad')

# 21- Creating AnnData from SCE

sce <- readRDS("./RDS_files/5_Integration/All/fov_integrated_sce.rds")

adata_sce <- sc$AnnData(
  X   = t(logcounts(sce)),
  obs = as.data.frame(colData(sce)),
  var = as.data.frame(rowData(sce))
)
adata_sce$obsm$update(umap = reducedDim(sce, "UMAP.HARMONY"))

adata_sce

sc$AnnData$write_h5ad(ad_fov_integrated, 'adata_fov_integrated_sce.h5ad')


sceasy::convertFormat(fov_integrated, from="seurat", to="anndata",
                      outFile='adata_fov_integrated.h5ad')


# Extract log normalised expression data and metadata for each fov from integrated seurat object

# Subset seurat object

Idents(fov_integrated) <- "sample_ID"
levels(fov_integrated)

fov_ID <- c("ROI1", "ROI2", "ROI3", "ROI4", "ROI5", "ROI6", "ROI7", "ROI8", "ROI9", "ROI10", "ROI11", "ROI12", "ROI13", "ROI14",
            "ROI15", "ROI16", "ROI17", "ROI18", "ROI19", "ROI20", "ROI21", "ROI22", "ROI23", "ROI24", "ROI25")

names(fov_ID) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, fov_ID)

# 9.2- Adding object metadata with cluster names in order - to be used in the bar plots
fov_integrated$fov_ID <- factor(fov_integrated@active.ident)

Idents(fov_integrated) <- "fov_ID"
levels(fov_integrated)

fov01 <- subset(fov_integrated, idents = c("ROI1"))

X   = t(GetAssayData(fov_integrated)) # log normalised RNA data slot
write.csv(X, "fov_integrated_exprMat_file.csv")
obs = fov_integrated[[]] # metadata file
write.csv(obs, "fov_integrated_metadata_file.csv")

DefaultAssay(fov_integrated) <- "RNA" 
Idents(fov_integrated) <- 'CIPR_ordered_clusters'
levels(fov_integrated)
sce <- as.SingleCellExperiment(fov_integrated, assay = "RNA")
saveRDS(sce, "fov_integrated_sce.rds")

##############################################################################################################################################################################################################################################################################################################################################################


# 22- PLOTTING FREQUENCY BAR PLOTS IN SEURAT WITH GGPLOT

library(Seurat)
library(scRNAseq)
library(SingleCellExperiment)
library(ggplot2)

fov_integrated <- readRDS("./RDS_files/5_Integration/fov_integrated.rds")

# 22.1- Plotting cell proportions per fov

#lung.immune <-ad[,ad$pheno_cluster_new %in% c('B cell', 'CD11c+ cell', 'EM CD4 T cell', 'Proliferative CD4 T cell', 'CD4 Treg cell',
#                                           'CD8 T cell', 'EM CD8 T cell', 'Proliferative CD8 T cell' , 'CD3+ cell',  
#                                           'SARSCoV2+ ArginaseHighVISTAHigh Activated Neutrophil', 
#                                           'ArginaseHighVISTAHigh Activated Neutrophil', 'ArginaseLowVISTALow Activated Neutrophil', 
#                                          'ArginaseLowVISTALow Neutrophil', 
#                                        'Classical Monocyte', 'SARSCoV2+ Interstitial Macrophage', 'Proliferating Interstitial Macrophage', 'Interstitial Macrophage', 
#                                       'Apoptotic SARSCoV2+  Alveolar Macrophage', 'Apoptotic Alveolar Macrophage', 'Alveolar Macrophage')]

freq.table.group <- as.data.frame(prop.table(x = table(fov_integrated$CIPR_ordered_clusters, fov_integrated$sample_ID), margin = 2))

freq.table.group[freq.table.group==0] <- NA
freq.table.group<-freq.table.group[complete.cases(freq.table.group),]
write.csv(freq.table.group, "frequency_immune_groups.csv")

# to plot the reversed order of cells
freq.table.group$Var1 <- as.ordered(factor(freq.table.group$Var1,
                                           levels=c('B cell', 'CD11c+ cell', 'EM CD4 T cell', 'Proliferative CD4 T cell', 'CD4 Treg cell',
                                                    'CD8 T cell', 'EM CD8 T cell', 'Proliferative CD8 T cell' , 'CD3+ cell',  
                                                    'SARSCoV2+ ArginaseHighVISTAHigh Activated Neutrophil', 
                                                    'ArginaseHighVISTAHigh Activated Neutrophil', 'ArginaseLowVISTALow Activated Neutrophil', 
                                                    'ArginaseLowVISTALow Neutrophil', 
                                                    'Classical Monocyte', 'SARSCoV2+ Interstitial Macrophage', 'Proliferating Interstitial Macrophage', 'Interstitial Macrophage', 
                                                    'Apoptotic SARSCoV2+  Alveolar Macrophage', 'Apoptotic Alveolar Macrophage', 'Alveolar Macrophage')))
# to plot the normal order of cells
library(forcats)
#freq.table.group$Var1 <- fct_rev(freq.table.group$Var1)
freq.table.group$Var1 <- freq.table.group$Var1

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2))

#freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2, levels=rev(c("HIV_Neg", "HIV_Pos"))))


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


ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Sample", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = (levels(freq.table.group$Var2))) +
  scale_fill_manual(values=cell_type_color) + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(face = "plain", size=12, angle=45, vjust=1, hjust=0.5))

# # 22.1- Plotting cell proportions per recurrence

freq.table.group <- as.data.frame(prop.table(x = table(fov_integrated$CIPR_ordered_clusters, fov_integrated$Recurrence), margin = 2))

freq.table.group[freq.table.group==0] <- NA
freq.table.group<-freq.table.group[complete.cases(freq.table.group),]

# to plot the normal order of cells
library(forcats)
#freq.table.group$Var1 <- fct_rev(freq.table.group$Var1)
freq.table.group$Var1 <- freq.table.group$Var1

freq.table.group$Var2 <- as.ordered(factor(freq.table.group$Var2, levels=(c("Negative", "Local", "Liver", "Brain", "Multi_site"))))


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


ggplot(data=freq.table.group, aes(x=freq.table.group$Var2, y = freq.table.group$Freq, fill=freq.table.group$Var1)) + geom_bar(stat="identity", color="black") +
  labs(x="Recurrence", y="Proportion of cells", fill="Cell Type") + scale_x_discrete(limits = (levels(freq.table.group$Var2))) +
  scale_fill_manual(values=cell_type_color) + theme_minimal(base_size = 12) + theme(axis.text.x = element_text(face = "plain", size=12, angle=45, vjust=1, hjust=0.5))


##############################################################################################################################################################################################################################################################################################################################################################

# PLOTTING CIRCULAR PLOT FOR UMAP

library(plot1cell)
library(Seurat)

Idents(fov_integrated) <- "CIPR_ordered_clusters"

###Prepare data for ploting - merged atlas
circ_data <- prepare_circlize_data(fov_integrated, scale = 0.75)
set.seed(1234)


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

cluster_colors <- cell_type_color
disease_colors<- c("#F1CDB1", "#CBDFB8", "#8FA9DF")
names(disease_colors) <- c("1", "2", "3")
hiv_colors<-c("lightblue", "red")
names(hiv_colors) <- c("Low", "High")

###plot and save figures
pdf('circlize_plot.pdf', width = 6, height = 6)
plot_circlize(circ_data, do.label = FALSE, pt.size = 0.1, contour.levels = c(0, 0), col.use = cluster_colors, bg.color = 'white', kde2d.n = 200, repel = TRUE, label.cex = 0,
              contour.nlevels = 0)
add_track(circ_data, group = "Epithelial_cluster", colors = disease_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
add_track(circ_data, group = "KM_score", colors = hiv_colors, track_num = 3)
dev.off()

saveRDS(circ_data, "data_circular_umap.rds")





