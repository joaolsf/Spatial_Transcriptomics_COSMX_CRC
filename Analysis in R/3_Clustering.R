
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
library(garnett)
library(monocle3)

options(ggrepel.max.overlaps = Inf)

fov_integrated <- readRDS("./RDS_files/5_Integration/All/fov_integrated.rds")
saveRDS(fov_integrated, "./fov_integrated.rds")

#fov_integrated$IF_type.1 <- NULL
#fov_integrated$seurat_clusters <- NULL

# 8- Clustering analysis option 1 - run clustering based on harmony and fastMNN-integrated reduction based on RNA data only.

# 8.1- Clustering Based on Harmony
set.seed(1234)
DefaultAssay(fov_integrated) <- "SCT"
fov_integrated <- FindNeighbors(fov_integrated, reduction = "harmony") #when using harmony, the default assay has to be SCT for this function to work
fov_integrated <- FindClusters(fov_integrated, resolution = 1, cluster.name = "harmony_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_integrated), 5)

# UMAP of the identified clusters 
DimPlot(fov_integrated, group.by = c("IF_type","harmony_clusters"), reduction = "umap.harmony", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_integrated, vars = c("IF_type","harmony_clusters"),
             reduction.use = "umap.harmony", size = 1, ncol=2,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

# 8.1.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers that define clusters via differential expression. 
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. # FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.
# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, 
# and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. 
# You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. 
# As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. 
# While there is generally going to be a loss in power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_integrated) <- "SCT" 
fov_integrated.markers <- FindAllMarkers(fov_integrated, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_integrated.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap((subset(fov_epithelial, downsample = 1000)), features = top10$gene, size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 8))

# Quantifications
# Code to quantify the numbers of cells in each cluster and the proportion of cells 
table(fov_integrated$seurat_clusters)
prop.table(table(fov_integrated$seurat_clusters))
# median number of RNA molecules per cluster
tapply(fov_integrated$nCount_RNA, fov_integrated$seurat_clusters,  median)

saveRDS(fov_integrated.markers, "fov_integrated.markers_SCT_harmony.rds")

# 8.1.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   

DefaultAssay(fov_integrated) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_integrated) <- "harmony_clusters"
fov_integrated.marker2 <- FindAllMarkers(object = fov_integrated, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_integrated.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_integrated.marker2 %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_integrated, features = top10$gene, size = 4) + NoLegend() + theme(axis.text.y = element_text(size = 6))

# 8.1.3- Visualise clusters in a spatial context
Idents(fov_integrated) <- 'harmony_clusters'
levels(fov_integrated)
Images(fov_integrated)

ImageDimPlot(fov_integrated, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov02", group.by = c("IF_type","harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov23", group.by = c("IF_type","harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")


# 8.2- Clustering based on fastMNN
set.seed(1234)
DefaultAssay(fov_integrated) <- "SCT"
fov_integrated <- FindNeighbors(fov_integrated, reduction = "integrated.mnn")
fov_integrated <- FindClusters(fov_integrated, resolution = 1, cluster.name = "mnn_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_integrated), 5)

## UMAP of the identified clusters 
DimPlot(fov_integrated, group.by = c("IF_type","mnn_clusters"), reduction = "umap.mnn", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_integrated, vars = c("IF_type","mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

# 8.2.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers that define clusters via differential expression. 
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. 
# FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, 
# and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. 
# You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. 
# As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. 
# While there is generally going to be a loss in power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_integrated) <- "SCT" 
Idents(fov_integrated) <- "mnn_clusters"
fov_integrated.markers <- FindAllMarkers(fov_integrated, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_integrated.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste


# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_integrated, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

# Quantifications
# Code to quantify the numbers of cells in each cluster and the proportion of cells 
table(fov_integrated$seurat_clusters)
prop.table(table(fov_integrated$seurat_clusters))
# median number of RNA molecules per cluster
tapply(fov_integrated$nCount_RNA, fov_integrated$seurat_clusters,  median)

saveRDS(fov_integrated.markers, "fov_integrated.markers_fastMNN.rds")

# 8.2.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   

DefaultAssay(fov_integrated) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_integrated) <- "mnn_clusters"
fov_integrated.marker2 <- FindAllMarkers(object = fov_integrated, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_integrated.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_integrated.marker2 %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_integrated, features = top10$gene, size = 4) + NoLegend() + theme(axis.text.y = element_text(size = 6))

# 8.2.3- Visualise clusters in a spatial context
Idents(fov_integrated) <- 'mnn_clusters'
levels(fov_integrated)
Images(fov_integrated)

ImageDimPlot(fov_integrated, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov02", group.by = c("IF_type","mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov23", group.by = c("IF_type","mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")


# 8.3- Comparing uncorrected, Harmony and mnn clusters - plot in each UMAP

#DimPlot(fov_integrated, group.by = c("harmony_clusters","mnn_clusters"), reduction = "umap.harmony", label = FALSE) #try other ploting packages and explore the function
#DimPlot(fov_integrated, group.by = c("harmony_clusters","mnn_clusters"), reduction = "umap.mnn", label = FALSE)

multi_dittoDimPlot(fov_integrated, vars = c("harmony_clusters","mnn_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

multi_dittoDimPlot(fov_integrated, vars = c("harmony_clusters","mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

ImageDimPlot(fov_integrated, fov = "fov02", group.by = c("harmony_clusters","mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_integrated, fov = "fov23", group.by = c("harmony_clusters","mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")


set.seed(1234)
DefaultAssay(fov_integrated) <- "SCT"
fov_integrated <- FindNeighbors(fov_integrated, reduction = "pca") #when using harmony, the default assay has to be SCT for this function to work
fov_integrated <- FindClusters(fov_integrated, resolution = 1, cluster.name = "uncorrected_clusters")


multi_dittoDimPlot(fov_integrated, vars = c("IF_type","uncorrected_clusters"),
                   reduction.use = "umap", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Uncorrected clusters on UMAP")

multi_dittoDimPlot(fov_integrated, vars = c("uncorrected_clusters","harmony_clusters","mnn_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=3,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Comparing clusters on UMAP")

# 8.4- Plot frequencies of IF types and clusters per sample

dittoBarPlot(fov_integrated, var = "IF_type", group.by = "MMR_Status_Jen",  scale = c("percent")) 
dittoBarPlot(fov_integrated, var = "uncorrected_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_integrated, var = "harmony_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_integrated, var = "mnn_clusters", group.by = "sample_ID",  scale = c("percent"))
dittoBarPlot(fov_integrated, var = "harmony_clusters", group.by = "IF_type",  scale = c("percent"))
dittoBarPlot(fov_integrated, var = "IF_type", group.by = "harmony_clusters",  scale = c("percent"))
dittoBarPlot(fov_integrated, var = "IF_marker", group.by = "harmony_clusters",  scale = c("percent"))

# 8.5- Plot normalised protein expression per cluster

DefaultAssay(fov_integrated) <- "Protein"
Idents(fov_integrated) <- "harmony_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

VlnPlot(fov_integrated, features = features, layer = "counts", pt.size = 0.001, ncol = 2) + NoLegend()
VlnPlot(fov_integrated, features = features, layer = "data", pt.size = 0.001, ncol = 2) + NoLegend() #CLR normalised protein exp
RidgePlot(fov_integrated, features = features, ncol = 2)
RidgePlot(fov_integrated, features = features, layer = "counts", ncol = 2)


########################################################################################################################################################################################################################################################################################################################################################################

# 9- Clustering analysis option 2 - run clustering based on harmony and fastMNN-integrated reduction based on epithelial vs stromal/CD3/other gated data.

# 9.1- Gate (subset) seurat object
# reference: https://github.com/satijalab/seurat/issues/1748
# Subset object based on the IF type
DefaultAssay(fov_integrated) <- "SCT"

# Rename idents of If type
Idents(fov_integrated) <- "IF_type"
levels(fov_integrated)
IF_type2 <- c("Stromal", "Immune", "CD3 Lymphocyte", "Epithelial")
names(IF_type2) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, IF_type2)

# Adding object metadata with
fov_integrated$IF_type2 <- factor(fov_integrated@active.ident)

#fov_integrated <- RenameIdents(object = fov_integrated, `Stromal` = "Immune")
#fov_integrated <- RenameIdents(object = fov_integrated, `Other` = "Stromal")

# gate IF epithelial
fov_epithelial <- subset(fov_integrated, idents = c("Epithelial"))
saveRDS(fov_epithelial, "./fov_epithelial.rds")

dittoDimPlot(fov_epithelial, var = "IF_type",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_epithelial, var = "IF_type",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_epithelial, var = "IF_type", split.by = "sample_ID",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE)

# gate IF stromal cells
fov_stromal <- subset(fov_integrated, idents = c("Stromal"))
saveRDS(fov_stromal, "./fov_stromal.rds")
fov_stromal <- readRDS("./fov_stromal.rds")

Idents(fov_stromal) <- "IF_type"
fov_stromal <- RenameIdents(object = fov_stromal, `Other` = "Stromal")

dittoDimPlot(fov_stromal, var = "IF_type",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_stromal, var = "IF_type",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_stromal, var = "IF_type", split.by = "sample_ID",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE)

# gate IF immune and CD3 lymphocyte cells
fov_immune <- subset(fov_integrated, idents = c("CD3 Lymphocyte", "Immune"))
saveRDS(fov_immune, "./fov_immune.rds")
fov_immune <- readRDS("./fov_immune.rds")

dittoDimPlot(fov_immune, var = "IF_type",
             reduction.use = "umap.harmony", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_immune, var = "IF_type",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE, labels.size = 1, labels.highlight = TRUE)

dittoDimPlot(fov_immune, var = "IF_type2", split.by = "sample_ID",
             reduction.use = "umap.mnn", size = 0.75,
             do.label = FALSE)

########################################################################################################################################################################################################################################################################################################################################################################

# 10- Clustering option 2- clustering within IF epithelial compartment

fov_epithelial <- readRDS("./RDS_files/6_gated/fov_epithelial.rds")

# 10.1- Clustering Based on Harmony
set.seed(1234)
DefaultAssay(fov_epithelial) <- "SCT"
fov_epithelial <- FindNeighbors(fov_epithelial, reduction = "harmony") #when using harmony, the default assay has to be SCT for this function to work
fov_epithelial <- FindClusters(fov_epithelial, resolution = 1, cluster.name = "ep_harmony_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_epithelial), 5)

# UMAP of the identified clusters 
DimPlot(fov_epithelial, group.by = c("IF_type","ep_harmony_clusters"), reduction = "umap.harmony", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_epithelial, vars = c("harmony_clusters","ep_harmony_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

dittoDimPlot(fov_epithelial, var = "ep_harmony_clusters", split.by = "sample_ID",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

# 10.1.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers that define clusters via differential expression. 
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. # FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.
# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, 
# and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. 
# You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. 
# As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. 
# While there is generally going to be a loss in power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_epithelial) <- "SCT" 
fov_epithelial.markers <- FindAllMarkers(fov_epithelial, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_epithelial.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_epithelial.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_epithelial, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

saveRDS(fov_epithelial.markers, "fov_epithelial.markers_SCT_harmony.rds")

# 10.1.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   

DefaultAssay(fov_epithelial) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_epithelial) <- "ep_harmony_clusters"
fov_epithelial.marker2 <- FindAllMarkers(object = fov_epithelial, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_epithelial.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste


# 10.1.3- Visualise clusters in a spatial context
Idents(fov_epithelial) <- 'ep_harmony_clusters'
levels(fov_epithelial)
Images(fov_epithelial)

ImageDimPlot(fov_epithelial, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov02", group.by = c("IF_marker","ep_harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", group.by = c("IF_marker","ep_harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")


# 10.2- Clustering based on fastMNN
set.seed(1234)
DefaultAssay(fov_epithelial) <- "SCT"
fov_epithelial <- FindNeighbors(fov_epithelial, reduction = "integrated.mnn")
fov_epithelial <- FindClusters(fov_epithelial, resolution = 1, cluster.name = "ep_mnn_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_epithelial), 5)

## UMAP of the identified clusters 
DimPlot(fov_epithelial, group.by = c("IF_type","ep_mnn_clusters"), reduction = "umap.mnn", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_epithelial, vars = c("IF_type","ep_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

multi_dittoDimPlot(fov_epithelial, vars = c("mnn_clusters","ep_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

# 10.2.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_epithelial) <- "SCT" 
Idents(fov_epithelial) <- "ep_mnn_clusters"
fov_epithelial.markers <- FindAllMarkers(fov_epithelial, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_epithelial.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_epithelial.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_epithelial, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

saveRDS(fov_epithelial.markers, "fov_epithelial.markers_fastMNN.rds")

# 10.2.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   
DefaultAssay(fov_epithelial) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_epithelial) <- "ep_mnn_clusters"
fov_epithelial.marker2 <- FindAllMarkers(object = fov_epithelial, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_epithelial.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# 10.2.3- Visualise clusters in a spatial context
Idents(fov_epithelial) <- 'ep_mnn_clusters'
levels(fov_epithelial)
Images(fov_epithelial)

ImageDimPlot(fov_epithelial, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov02", group.by = c("IF_type","ep_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", group.by = c("IF_type","ep_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# 10.3- Comparing uncorrected, Harmony and mnn clusters - plot in each UMAP

multi_dittoDimPlot(fov_epithelial, vars = c("ep_harmony_clusters","ep_mnn_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

multi_dittoDimPlot(fov_epithelial, vars = c("ep_harmony_clusters","ep_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Epithelial clusters on UMAP")

ImageDimPlot(fov_epithelial, fov = "fov02", group.by = c("ep_harmony_clusters","ep_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_epithelial, fov = "fov23", group.by = c("ep_harmony_clusters","ep_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

set.seed(1234)
DefaultAssay(fov_epithelial) <- "SCT"
fov_epithelial <- FindNeighbors(fov_epithelial, reduction = "pca") #when using harmony, the default assay has to be SCT for this function to work
fov_epithelial <- FindClusters(fov_epithelial, resolution = 1, cluster.name = "ep_uncorrected_clusters")

multi_dittoDimPlot(fov_epithelial, vars = c("IF_type","ep_uncorrected_clusters"),
                   reduction.use = "umap", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Uncorrected clusters on UMAP")

multi_dittoDimPlot(fov_epithelial, vars = c("ep_uncorrected_clusters","ep_harmony_clusters","ep_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=3,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Comparing epithelial clusters on UMAP")

# 10.4- Plot frequencies of IF types and clusters per sample

dittoBarPlot(fov_epithelial, var = "ep_uncorrected_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "ep_harmony_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "ep_mnn_clusters", group.by = "sample_ID",  scale = c("percent"))
dittoBarPlot(fov_epithelial, var = "ep_harmony_clusters", group.by = "IF_type",  scale = c("percent"))
dittoBarPlot(fov_epithelial, var = "IF_type", group.by = "ep_harmony_clusters",  scale = c("percent"))
dittoBarPlot(fov_epithelial, var = "IF_marker", group.by = "ep_harmony_clusters",  scale = c("percent"))
dittoBarPlot(fov_epithelial, var = "ep_harmony_clusters", group.by = "MMR_Status_Jen",  scale = c("percent")) 
dittoBarPlot(fov_epithelial, var = "ep_harmony_clusters", group.by = "CMS",  scale = c("percent"))

# 10.5- Plot normalised protein expression per cluster

DefaultAssay(fov_epithelial) <- "Protein"
Idents(fov_epithelial) <- "ep_harmony_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

FeaturePlot(fov_epithelial, features = features, reduction = "umap.mnn", slot = "counts") + NoLegend()
VlnPlot(fov_epithelial, features = features, layer = "data", pt.size = 0.001, ncol = 2) + NoLegend() #CLR normalised protein exp
RidgePlot(fov_epithelial, features = features, ncol = 2)

dittoScatterPlot(
  object = fov_epithelial,
  x.var = "MeanPanCK", y.var = "MeanCD45",
  color.var = "ep_harmony_clusters")

dittoScatterPlot(
  object = fov_epithelial,
  x.var = "MeanPanCK", y.var = "MeanCD3",
  color.var = "ep_harmony_clusters")

########################################################################################################################################################################################################################################################################################################################################################################

# 11- Clustering option 2- clustering within IF stromal compartment

fov_stromal <- readRDS("./RDS_files/6_gated/fov_stromal.rds")

# 11.1- Clustering Based on Harmony
set.seed(1234)
DefaultAssay(fov_stromal) <- "SCT"
fov_stromal <- FindNeighbors(fov_stromal, reduction = "harmony") #when using harmony, the default assay has to be SCT for this function to work
fov_stromal <- FindClusters(fov_stromal, resolution = 1, cluster.name = "st_harmony_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_stromal), 5)

# UMAP of the identified clusters 
DimPlot(fov_stromal, group.by = c("IF_type2","st_harmony_clusters"), reduction = "umap.harmony", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_stromal, vars = c("IF_type2","st_harmony_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

multi_dittoDimPlot(fov_stromal, vars = c("harmony_clusters","st_harmony_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

# 11.1.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_stromal) <- "SCT" 
fov_stromal.markers <- FindAllMarkers(fov_stromal, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_stromal.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_stromal.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_stromal, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

saveRDS(fov_stromal.markers, "fov_stromal.markers_SCT_harmony.rds")

# 11.1.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   

DefaultAssay(fov_stromal) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_stromal) <- "st_harmony_clusters"
fov_stromal.marker2 <- FindAllMarkers(object = fov_stromal, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_stromal.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# 11.1.3- Visualise clusters in a spatial context
Idents(fov_stromal) <- 'st_harmony_clusters'
levels(fov_stromal)
Images(fov_stromal)

ImageDimPlot(fov_stromal, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov02", group.by = c("IF_marker","st_harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov23", group.by = c("IF_marker","st_harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# 11.2- Clustering based on fastMNN
set.seed(1234)
DefaultAssay(fov_stromal) <- "SCT"
fov_stromal <- FindNeighbors(fov_stromal, reduction = "integrated.mnn")
fov_stromal <- FindClusters(fov_stromal, resolution = 1, cluster.name = "st_mnn_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_stromal), 5)

## UMAP of the identified clusters 
DimPlot(fov_stromal, group.by = c("IF_type","st_mnn_clusters"), reduction = "umap.mnn", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_stromal, vars = c("IF_type","st_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

multi_dittoDimPlot(fov_stromal, vars = c("mnn_clusters","st_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

# 11.2.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_stromal) <- "SCT" 
Idents(fov_stromal) <- "st_mnn_clusters"
fov_stromal.markers <- FindAllMarkers(fov_stromal, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_stromal.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_stromal.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_stromal, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

saveRDS(fov_stromal.markers, "fov_stromal.markers_fastMNN.rds")

# 11.2.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   
DefaultAssay(fov_stromal) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_stromal) <- "st_mnn_clusters"
fov_stromal.marker2 <- FindAllMarkers(object = fov_stromal, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_stromal.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# 11.2.3- Visualise clusters in a spatial context
Idents(fov_stromal) <- 'st_mnn_clusters'
levels(fov_stromal)
Images(fov_stromal)

ImageDimPlot(fov_stromal, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov02", group.by = c("IF_type","st_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov23", group.by = c("IF_type","st_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# 11.3- Comparing uncorrected, Harmony and mnn clusters - plot in each UMAP

multi_dittoDimPlot(fov_stromal, vars = c("st_harmony_clusters","st_mnn_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Stromall clusters on UMAP")

multi_dittoDimPlot(fov_stromal, vars = c("st_harmony_clusters","st_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Stromal clusters on UMAP")

ImageDimPlot(fov_stromal, fov = "fov02", group.by = c("st_harmony_clusters","st_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_stromal, fov = "fov23", group.by = c("st_harmony_clusters","st_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

set.seed(1234)
DefaultAssay(fov_stromal) <- "SCT"
fov_stromal <- FindNeighbors(fov_stromal, reduction = "pca") #when using harmony, the default assay has to be SCT for this function to work
fov_stromal <- FindClusters(fov_stromal, resolution = 1, cluster.name = "st_uncorrected_clusters")

multi_dittoDimPlot(fov_stromal, vars = c("IF_type","st_uncorrected_clusters"),
                   reduction.use = "umap", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Uncorrected clusters on UMAP")

multi_dittoDimPlot(fov_stromal, vars = c("st_uncorrected_clusters","st_harmony_clusters","st_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=3,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Comparing stromal clusters on UMAP")

# 11.4- Plot frequencies of IF types and clusters per sample

dittoBarPlot(fov_stromal, var = "st_uncorrected_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "st_harmony_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_stromal, var = "st_mnn_clusters", group.by = "sample_ID",  scale = c("percent"))
dittoBarPlot(fov_stromal, var = "st_harmony_clusters", group.by = "IF_type",  scale = c("percent"))
dittoBarPlot(fov_stromal, var = "IF_type2", group.by = "st_harmony_clusters",  scale = c("percent"))
dittoBarPlot(fov_stromal, var = "IF_marker", group.by = "st_harmony_clusters",  scale = c("percent"))

# 11.5- Plot normalised protein expression per cluster

DefaultAssay(fov_stromal) <- "Protein"
Idents(fov_stromal) <- "st_harmony_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

FeaturePlot(fov_stromal, features = features, reduction = "umap.mnn", slot = "data") + NoLegend()
VlnPlot(fov_stromal, features = features, layer = "data", pt.size = 0.001, ncol = 2) + NoLegend() #CLR normalised protein exp
RidgePlot(fov_stromal, features = features, ncol = 2)

dittoScatterPlot(
  object = fov_stromal,
  x.var = "MeanCD45", y.var = "MeanPanCK",
  color.var = "st_harmony_clusters")

dittoScatterPlot(
  object = fov_stromal,
  x.var = "MeanCD3", y.var = "MeanPanCK",
  color.var = "st_harmony_clusters")

########################################################################################################################################################################################################################################################################################################################################################################

# 12- Clustering option 2- clustering within IF stromal compartment

fov_immune <- readRDS("./RDS_files/6_gated/fov_immune.rds")
fov_immune <- readRDS("./fov_immune.rds")
saveRDS(fov_immune, "./fov_immune.rds")

# 12.1- Clustering Based on Harmony
set.seed(1234)
DefaultAssay(fov_immune) <- "SCT"
fov_immune <- FindNeighbors(fov_immune, reduction = "harmony") #when using harmony, the default assay has to be SCT for this function to work
fov_immune <- FindClusters(fov_immune, resolution = 1, cluster.name = "im_harmony_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_immune), 5)

# UMAP of the identified clusters 
DimPlot(fov_immune, group.by = c("IF_type2","im_harmony_clusters"), reduction = "umap.harmony", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_immune, vars = c("IF_type2","im_harmony_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

multi_dittoDimPlot(fov_immune, vars = c("harmony_clusters","im_harmony_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

# 12.1.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_immune) <- "SCT" 
fov_immune.markers <- FindAllMarkers(fov_immune, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_immune.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_immune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_immune, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

saveRDS(fov_immune.markers, "fov_immune.markers_SCT_harmony.rds")

# 11.1.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   

DefaultAssay(fov_immune) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_immune) <- "im_harmony_clusters"
fov_immune.marker2 <- FindAllMarkers(object = fov_immune, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_immune.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# 11.1.3- Visualise clusters in a spatial context
Idents(fov_immune) <- 'im_harmony_clusters'
levels(fov_immune)
Images(fov_immune)

ImageDimPlot(fov_immune, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_immune, fov = "fov02", group.by = c("IF_marker","im_harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_immune, fov = "fov23", group.by = c("IF_marker","im_harmony_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# 11.2- Clustering based on fastMNN
set.seed(1234)
DefaultAssay(fov_immune) <- "SCT"
fov_immune <- FindNeighbors(fov_immune, reduction = "integrated.mnn")
fov_immune <- FindClusters(fov_immune, resolution = 1, cluster.name = "im_mnn_clusters")

# Look at cluster IDs of the first 5 cells
head(Idents(fov_immune), 5)

## UMAP of the identified clusters 
DimPlot(fov_immunel, group.by = c("IF_type2","im_mnn_clusters"), reduction = "umap.mnn", label = FALSE) #try other ploting packages and explore the function

multi_dittoDimPlot(fov_immune, vars = c("IF_type2","im_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

multi_dittoDimPlot(fov_immune, vars = c("mnn_clusters","im_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("MNN clusters on UMAP")

# 11.2.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_immune) <- "SCT" 
Idents(fov_immune) <- "im_mnn_clusters"
fov_immune.markers <- FindAllMarkers(fov_immune, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_immune.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_immune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_immune, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

saveRDS(fov_immune.markers, "fov_immune.markers_fastMNN.rds")

# 11.2.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   
DefaultAssay(fov_immune) <- "RNA" # much more marker genes than with SCT slot
Idents(fov_immune) <- "im_mnn_clusters"
fov_immune.marker2 <- FindAllMarkers(object = fov_immune, only.pos = FALSE, logfc.threshold = 0.25)
clipr::write_clip(fov_immune.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# 11.2.3- Visualise clusters in a spatial context
Idents(fov_immune) <- 'im_mnn_clusters'
levels(fov_immune)
Images(fov_immune)

ImageDimPlot(fov_immune, fov = c("fov02", "fov23"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_immune, fov = "fov02", group.by = c("IF_type2","im_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_immune, fov = "fov23", group.by = c("IF_type2","im_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

# 11.3- Comparing uncorrected, Harmony and mnn clusters - plot in each UMAP

multi_dittoDimPlot(fov_immune, vars = c("im_harmony_clusters","im_mnn_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Stromall clusters on UMAP")

multi_dittoDimPlot(fov_immune, vars = c("im_harmony_clusters","im_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Stromal clusters on UMAP")

ImageDimPlot(fov_immune, fov = "fov02", group.by = c("im_harmony_clusters","im_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")
ImageDimPlot(fov_immune, fov = "fov23", group.by = c("im_harmony_clusters","im_mnn_clusters"), axes = TRUE, cols = "polychrome", size = 1, boundaries = "segmentation")

set.seed(1234)
DefaultAssay(fov_immune) <- "SCT"
fov_immune <- FindNeighbors(fov_immune, reduction = "pca") #when using harmony, the default assay has to be SCT for this function to work
fov_immune<- FindClusters(fov_immune, resolution = 1, cluster.name = "im_uncorrected_clusters")

multi_dittoDimPlot(fov_immune, vars = c("IF_type2","im_uncorrected_clusters"),
                   reduction.use = "umap", size = 1, ncol=2,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Uncorrected clusters on UMAP")

multi_dittoDimPlot(fov_immune, vars = c("im_uncorrected_clusters","im_harmony_clusters","im_mnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=3,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Comparing stromal clusters on UMAP")

# 11.4- Plot frequencies of IF types and clusters per sample

dittoBarPlot(fov_immune, var = "im_uncorrected_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_immune, var = "im_harmony_clusters", group.by = "sample_ID",  scale = c("percent")) 
dittoBarPlot(fov_immune, var = "im_mnn_clusters", group.by = "sample_ID",  scale = c("percent"))
dittoBarPlot(fov_immune, var = "im_harmony_clusters", group.by = "IF_type",  scale = c("percent"))
dittoBarPlot(fov_immune, var = "IF_type2", group.by = "im_harmony_clusters",  scale = c("percent"))
dittoBarPlot(fov_immune, var = "IF_marker", group.by = "im_harmony_clusters",  scale = c("percent"))

dittoDimPlot(fov_immune, var = "im_harmony_clusters", split.by = "sample_ID",
             reduction.use = "umap.harmony", size = 1,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Harmony clusters on UMAP")

# 11.5- Plot normalised protein expression per cluster

DefaultAssay(fov_immune) <- "Protein"
Idents(fov_immune) <- "im_harmony_clusters"
features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")

FeaturePlot(fov_immune, features = features, reduction = "umap.mnn", slot = "data") + NoLegend()
VlnPlot(fov_immune, features = features, layer = "data", pt.size = 0.001, ncol = 2) + NoLegend() #CLR normalised protein exp
RidgePlot(fov_immune, features = features, ncol = 2)

dittoScatterPlot(
  object = fov_immune,
  x.var = "MeanCD45", y.var = "MeanPanCK",
  color.var = "im_harmony_clusters")

dittoScatterPlot(
  object = fov_immune,
  x.var = "MeanCD3", y.var = "MeanPanCK",
  color.var = "im_harmony_clusters")

########################################################################################################################################################################################################################################################################################################################################################################

# 13- Clustering option 3- WNN protein+RNA data

# Reference: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-cite-seq-rna-adt-1

DefaultAssay(fov_integrated) <- "Protein"

prots = rownames(fov_integrated@assays[["Protein"]])[1:8]
VariableFeatures(fov_integrated) = prots

fov_integrated = NormalizeData(fov_integrated, normalization.method = 'CLR', margin = 2)
fov_integrated = ScaleData(fov_integrated) %>% 
  RunPCA(reduction.name = 'apca', verbose = FALSE)

# For each cell, we calculate its closest neighbors in the dataset based on a weighted combination of RNA and protein similarities. 
# The cell-specific modality weights and multimodal neighbors are calculated in a single function, which takes ~2 minutes to run on this dataset. 
# We specify the dimensionality of each modality (similar to specifying the number of PCs to include in scRNA-seq clustering), but you can vary these settings to see that small changes have minimal effect on the overall results.

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using [['weighted.nn']]
# The WNN graph can be accessed at [["wknn"]], 
# and the SNN graph used for clustering at [["wsnn"]]

fov_integrated  = FindMultiModalNeighbors(
  fov_integrated , reduction.list = list("harmony", "apca"), #originally pca, but resulted in non-integrated output, so use an integrated red dim
  dims.list = list(1:50, 1:7), 
  modality.weight.name = c("SCT.weight", "Protein.weight"),
  verbose = FALSE
)

# The number of provided modality.weight.name is not equal to the number of modalities. 
# SCT.weight and Protein.weight are used to store the modality weights

# We can now use these results for downstream analysis, such as visualization and clustering. 
# For example, we can create a UMAP visualization of the data based on a weighted combination of RNA and protein data.
# We can also perform graph-based clustering and visualize these results on the UMAP, alongside a set of cell annotations.
fov_integrated  <- RunUMAP(fov_integrated , nn.name = "weighted.nn", 
                    reduction.name = "wnn.umap", 
                    reduction.key = "wnnUMAP_")

fov_integrated  <- FindClusters(fov_integrated , graph.name = "wsnn", 
                         resolution = 1, cluster.name = "wnn_clusters",
                         verbose = FALSE, 
                         random.seed = 1990)

# Visualise WNN clusters in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(fov_integrated, var = "wnn_clusters",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

dittoDimPlot(fov_integrated, var = "IF_type",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

dittoDimPlot(fov_integrated, var = "sample_ID",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

dittoDimPlot(fov_integrated, var = "sample_ID", split.by = "sample_ID",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

features <- c("MeanPanCK", "MaxPanCK","MeanCD45", "MaxCD45", "MeanCD3", "MaxCD3")
FeaturePlot(fov_integrated, features = features, reduction = "wnn.umap", slot = "counts") + NoLegend()
VlnPlot(fov_integrated, features = features, group.by = "wnn_clusters" ,layer = "data", pt.size = 0.001, ncol = 2) + NoLegend() #CLR normalised protein exp
RidgePlot(fov_integrated, features = features, ncol = 2)

VlnPlot(fov_integrated, features = features, group.by = "sample_ID" ,layer = "data", pt.size = 0.001, ncol = 2) + NoLegend() #CLR normalised protein exp
RidgePlot(fov_integrated, features = features, ncol = 2)



multi_dittoDimPlot(fov_integrated, vars = c("harmony_clusters","mnn_clusters", "wnn_clusters"),
                   reduction.use = "umap.mnn", size = 1, ncol=3,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Comparing clusters on UMAP")

multi_dittoDimPlot(fov_integrated, vars = c("harmony_clusters","mnn_clusters", "wnn_clusters"),
                   reduction.use = "umap.harmony", size = 1, ncol=3,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Comparing clusters on UMAP")
 
dittoBarPlot(fov_integrated, var = "wnn_clusters", group.by = "sample_ID",  scale = c("percent"))
dittoBarPlot(fov_integrated, var = "wnn_clusters", group.by = "IF_type",  scale = c("percent"))
dittoBarPlot(fov_integrated, var = "IF_type", group.by = "wnn_clusters",  scale = c("percent"))
dittoBarPlot(fov_integrated, var = "IF_marker", group.by = "wnn_clusters",  scale = c("percent"))

VlnPlot(fov_integrated, features = "SCT.weight", group.by = 'wnn_clusters', sort = TRUE, pt.size = 0.01) +
  NoLegend()

VlnPlot(fov_integrated, features = "Protein.weight", group.by = 'wnn_clusters', sort = TRUE, pt.size = 0.01) +
  NoLegend()

dittoScatterPlot(
  object = fov_integrated,
  x.var = "MeanPanCK", y.var = "MeanCD45",
  color.var = "wnn_clusters")

dittoScatterPlot(
  object = fov_integrated,
  x.var = "MeanPanCK", y.var = "MeanCD3",
  color.var = "wnn_clusters")

dittoScatterPlot(
  object = fov_integrated,
  x.var = "MeanCD45", y.var = "MeanCD3",
  color.var = "wnn_clusters")

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(fov_integrated) <- "SCT" 
Idents(fov_integrated) <- "wnn_clusters"
fov_integrated.markers <- FindAllMarkers(fov_integrated, only.pos = FALSE, logfc.threshold = 0.25) #min.pct = 0.25,
clipr::write_clip(fov_integrated.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
fov_integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(fov_integrated, features = top10$gene, size = 6) + NoLegend() + theme(axis.text.y = element_text(size = 8))

saveRDS(fov_integrated.markers, "fov_integrated.markers_WNN.rds")




