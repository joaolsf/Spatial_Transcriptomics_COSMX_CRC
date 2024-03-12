

# 10- Differential expression testing

# 10.1- Perform default differential expression tests - comparisons between clusters within the integrated/mapped Seurat object

# The bulk of Seurat’s differential expression features can be accessed through the FindMarkers() function. 
# As a default, Seurat performs differential expression based on the non-parametric Wilcoxon rank sum test. 
# This replaces the previous default test (‘bimod’). 
# To test for differential expression between two specific groups of cells, specify the ident.1 and ident.2 parameters. 

# which assay should be used for DGE comparisons? 
# Reference: Current best practices in single-cell RNA-seq analysis: a tutorial" described to use the raw pre-normalised measured values.

# list options for groups to perform differential expression on
DefaultAssay(BM_query) <- "SCT"
Idents(BM_query) <- 'predicted.celltype.l2'
levels(BM_query)

# 10.1.1- Prefilter features or cells to increase the speed of DE testing

# Reference: https://satijalab.org/seurat/articles/de_vignette.html
# To increase the speed of marker discovery, particularly for large datasets, Seurat allows for pre-filtering of features or cells. 
# For example, features that are very infrequently detected in either group of cells, or features that are expressed at similar average levels, are unlikely to be differentially expressed. 
# Example use cases of the min.pct, logfc.threshold, min.diff.pct, and max.cells.per.ident parameters are demonstrated below.

# Use min.pct to pre-filter features that are detected at <X% frequency;
# Use logfc.threshold to pre-filter features that have less than a X-fold change between the average expression;
# Use min.diff.pct to pre-filter features whose detection percentages across the two groups are similar;
# Use max.cells.per.ident to subsample each group to a maximum of 200 cells. Can be very useful for large clusters, or computationally-intensive DE tests;
# Increasing min.pct, logfc.threshold, and min.diff.pct, will increase the speed of DE testing, but could also miss features that are prefiltered.

# 10.1.2- Using the RNA normalised slot

BM_RNA_data_markers1 <- FindAllMarkers(BM_query, assay = "RNA", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM_RNA_data_markers1) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste
saveRDS(BM_RNA_data_markers1, "BM_RNA_data_markers.rds")

# top 5 upregulated genes per cluster
BM1_RNA_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  print(n=100)

# 10.1.3- Using the SCT normalised slot

BM_SCT_data_markers <- FindAllMarkers(BM_query, assay = "SCT", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM_SCT_data_markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste
saveRDS(BM_SCT_data_markers, "BM_SCT_data_markers.rds")
BM_SCT_data_markers <- readRDS("./Plots/12_DGE_after_annotation/BM_SCT_data_markers.rds")

BM_SCT_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  print(n=100)

BM_SCT_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

DoHeatmap(BM_query, features = top10$gene, size = 2) + NoLegend() + theme(axis.text.y = element_text(size = 3))
# Subset Seurat object top plot heatmap of specific clusters
DoHeatmap(subset(BM_query, downsample = 500, idents = c("Late Eryth", "Early Eryth")), features = top10$gene, size = 2, assay = "RNA", combine = TRUE) + theme(axis.text.y = element_text(size = 3))

# Heatmap visualization - DittoHeatmap
dittoHeatmap(BM_query, genes = top10$gene,
             assay = "SCT", order.by = c("predicted.celltype.l2"), cluster_rows = FALSE,
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = inferno(100), 
             annot.by = c("predicted.celltype.l2"), fontsize_row = 2)


#  Find differentially expressed features between Late Eryth and all other cells, only search for positive markers
Late_Eryth_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "Late Eryth", only.pos = TRUE, logfc.threshold = 0.1)
print((Late_Eryth_markers))
clipr::write_clip(Late_Eryth_markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# Find all markers distinguishing late vs early erythrocytes
Late_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "Late Eryth", ident.2 = "Early Eryth", min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.1)
print(Late_markers)
clipr::write_clip(Late_markers)

# Find all markers distinguishing early erythrocytes vs EMP
Early_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "Early Eryth", ident.2 = "EMP", min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.1)
print(Early_markers)
clipr::write_clip(Early_markers)

# Find all markers distinguishing CD8 Effector 1 erythrocytes vs CD8 Effector 2
CD8_effector_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "CD8 Effector_2", ident.2 = "CD8 Effector_1", min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.1)
print(CD8_effector_markers)
clipr::write_clip(CD8_effector_markers)


DGE1 <- read.csv("./Plots/12_DGE_after_annotation/DGE_CD8Effector2_CD8Effector1.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = DGE1$Features,
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = FALSE,
                title = "DGE CD8 Effector_2 vs CD8 Effector_1",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black',
) + theme_light()


# Highlight specific genes
celltype1 <- c('TFRC','GYPA', 'ITGA4', 'CD34')
celltype2 <- c('ITGA4')

EnhancedVolcano(toptable = DGE1,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = c(celltype1, celltype2),
                lab = DGE1$Features,
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                shape = 42,
                boxedLabels = FALSE,
                title = "DGE late eryth vs early eryth",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black',
                # encircle
                encircle = celltype1,
                encircleCol = 'skyblue',
                encircleSize = 2.5,
                encircleFill = 'pink',
                encircleAlpha = 1/2,
                # shade
                shade = celltype2,
                shadeAlpha = 1/2,
                shadeFill = 'skyblue',
                shadeSize = 1,
                shadeBins = 5
) + theme_bw()

# Exploring the expression of canonical marker genes
Idents(BM1_query) <- 'predicted.celltype.l2'
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "RNA", slot = "data", sort = TRUE, ncol=1) + NoLegend()
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "SCT", slot = "data", sort = TRUE, ncol=1) + NoLegend()

# 10.2- Subset to compare and plot specific clusters within the same condition

# Reference: https://satijalab.org/seurat/articles/de_vignette.html
# Reference: https://satijalab.org/seurat/articles/visualization_vignette.html

BM_query_eryth <- subset(BM_query, idents = c("Late Eryth", "Early Eryth"))
Late_markers <- FindAllMarkers(BM_query_eryth, assay = "RNA", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1)

Late_markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC) -> top50

DoHeatmap(subset(BM_query_eryth, downsample = 500), assay = "RNA", features = top50$gene, size = 4) + NoLegend() + theme(axis.text.y = element_text(size = 6))

dittoHeatmap(subset(BM_query_eryth, downsample = 500), genes = top50$gene,
             assay = "RNA", slot = "scale.data",order.by = c("predicted.celltype.l2"), cluster_rows = FALSE,
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = inferno(50), 
             annot.by = c("predicted.celltype.l2"), fontsize_row = 6)

# 10.3- DGE between same clusters different samples/conditions

# Reference: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions-1
# Using the normalized datasets with known celltype annotation, we can ask what genes change in different conditions for cells of the same type. 
# First, we create a column in the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column.
BM_query$celltype.sample <- paste(BM_query$predicted.celltype.l2, BM_query$sample_id,
                                  sep = "_")
Idents(BM_query) <- "celltype.sample"
levels(BM_query)

# Prior to performing differential expression, we first run PrepSCTFindMarkers, which ensures that the fixed value is set properly. 
# Then we use FindMarkers(assay="SCT") to find differentially expressed genes. 
BM_query <- PrepSCTFindMarkers(BM_query)

CD8_markers <- FindMarkers(BM_query, assay = "SCT", ident.1 = "CD8 Effector_2_HD_BM1.1", ident.2 = "CD8 Effector_2_HD_BM4.1", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1,
                           verbose = FALSE)
head(CD8_markers, n = 20)
clipr::write_clip(CD8_markers)
DGE1 <- read.csv("./Plots/12_DGE_after_annotation/DGE_CD8_BM1_BM4.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = DGE1$Features,
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = FALSE,
                title = "DGE CD8 Effector_2 vs CD8 Effector_1",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black',
) + theme_light()


# If running on a subset of the original object after running PrepSCTFindMarkers(), 
# FindMarkers() should be invoked with recorrect_umi = FALSE to use the existing corrected counts:
BM_query_CD8 <- subset(BM_query, idents = c("CD8 Effector_2_HD_BM1.1", "CD8 Effector_2_HD_BM4.1"))
CD8_markers <- FindAllMarkers(BM_query_CD8, assay = "SCT", ident.1 = "CD8 Effector_2_HD_BM1.1",
                              ident.2 = "CD8 Effector_2_HD_BM4.1", verbose = FALSE, recorrect_umi = FALSE)

# 10.4- Identify conserved cell type markers

# To identify canonical cell type marker genes that are conserved across conditions, we provide the FindConservedMarkers() function. 
# This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. 
# For example, we can identify genes that are conserved markers irrespective of stimulation condition in NK cells. 
# Note that the PrepSCTFindMarkers command does not to be rerun here.
Idents(BM_query) <- "predicted.celltype.l2"
LE.markers <- FindConservedMarkers(BM_query, assay = "SCT", ident.1 = "Late Eryth", grouping.var = "sample_id",
                                   verbose = FALSE)
head(LE.markers)

# 20- Subset Seurat object to cluster major cell types and re-add the info to the original seurat object

# reference: https://github.com/satijalab/seurat/issues/1748
# Subset BM1_query object based on the clusters and sub-cluster
BM1_eryth <- subset(BM1_query, idents = c("Late Eryth", "Early Eryth"))
BM1_eryth <- FindNeighbors(BM1_eryth, dims = 1:10, k.param = 5)
BM1_eryth <- FindClusters(BM1_eryth) 

# Generate a new column called sub_cluster in the metadata
BM1_query$sub_cluster <- as.character(Idents(BM1_query))

# Change the information of cells containing sub-cluster information
WhichCells(BM1_query, c("Late Eryth", "Early Eryth"))
rownames(BM1_eryth@meta.data) #to get cells names from BM1_eryth

BM1_query$sub_cluster[Cells(BM1_eryth)] <- paste("BM1_eryth",Idents(BM1_eryth))
DimPlot(BM1_query, group.by = "sub_cluster")
