
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

fov_integrated <- readRDS("./fov_integrated.rds")

saveRDS(fov_integrated, "./fov_integrated.rds")

# 1- Adding new metadata
# 1.1- Samples stratified in epithelial clusters 

DefaultAssay(fov_integrated) <- "RNA"
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

Epithelial_cluster <- c("2", "1", "3", "1", "1", "2", "1", "1", "2", "1", "2", "2", "1", "2", "3", "1", "2", "1", "1", "1", "3", "1", "2", "1", "1")

names(Epithelial_cluster) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, Epithelial_cluster)
fov_integrated$Epithelial_cluster <- factor(fov_integrated@active.ident)

# 1.2- Samples stratified in TME clusters 
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

TME_cluster <- c("2", "2",  "3", "1", "2", "2",  "2", "1", "3", "2", "2", "NA", "1", "2", "3", "2", "1", "2", "1", "2", "2", "3", "2", "1", "1")

names(TME_cluster) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, TME_cluster)
fov_integrated$TME_cluster <- factor(fov_integrated@active.ident)

# 1.3- Samples stratified in aSMA (stromal) clusters 
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

aSMA_cluster <- c("1", "2", "NA", "1", "2", "2", "2", "3", "2", "1", "1", "NA", "NA", "1", "2", "1", "NA", "NA", "2", "1", "NA", "2", "2", "NA", "NA")

names(aSMA_cluster) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, aSMA_cluster)
fov_integrated$aSMA_cluster <- factor(fov_integrated@active.ident)


# 1.4- Add Tumor Stroma percentage as TSP score
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

TSP_score <- c("<50%", ">50%", "<50%", ">50%", "<50%", "<50%", ">50%", "<50%", ">50%", "<50%", "<50%", "<50%", "<50%", "<50%", "<50%", "<50%", "<50%", "<50%", "<50%", "<50%", "<50%", ">50%", "<50%", "<50%", "<50%")

names(TSP_score) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, TSP_score)
fov_integrated$TSP_score <- factor(fov_integrated@active.ident)

# 1.5- Add new KM score 
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

KM_score <- c("High", "Low", "Low", "Low", "Low", "Low","Low", "Low", "Low", "Low", "Low","Low", "Low", "Low", "Low", "High", "High", "High", "Low", "High", "Low", "Low", "High", "High", "High")

names(KM_score) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, KM_score)
fov_integrated$KM_score <- factor(fov_integrated@active.ident)

# 1.6- Add new recurrence labeling
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)
fov_integrated$Recurrence<- NULL
Recurrence <- c("Negative", "Liver", "Negative", "Negative", "Negative", "Local", "Negative", "Brain", "Negative", "Negative", "Negative", "Negative", "Multi_site", "Negative", "Negative", "Negative", "Negative", "Multi_site", "Negative", "Negative", "Negative", "Negative", "Negative", "Negative", "Negative")

names(Recurrence) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, Recurrence)
fov_integrated$Recurrence <- factor(fov_integrated@active.ident)

# 1.7- T stage
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

T_stage <- c("1", "3", "1", "3", "4", "3", "3", "3", "3", "3", "3",  "3", "4", "3", "3", "3", "3", "2", "2", "3", "3", "3", "3",  "3", "1")

names(T_stage) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, T_stage)
fov_integrated$T_stage <- factor(fov_integrated@active.ident)

# 1.8- N stage
Idents(fov_integrated) <- 'sample_ID'
levels(fov_integrated)

N_stage <- c("1", "0", "0", "1", "2", "2", "1", "0", "2", "1", "2", "1", "0",  "0", "0", "0", "1", "0", "0", "0", "1", "0", "0", "0", "0")

names(N_stage) <- levels(fov_integrated)
fov_integrated <- RenameIdents(fov_integrated, N_stage)
fov_integrated$N_stage <- factor(fov_integrated@active.ident)

# Cleaning the seurat object by removing some metadata columns

fov_integrated$SCSA_clusters_st1 <- NULL
fov_integrated$CIPR_clusters_st1 <- NULL
fov_integrated$CIPR_clusters_st1_v2 <- NULL
fov_integrated$Insitutype_clusters_st1 <- NULL
fov_integrated$Insitutype_clusters_st1_v2 <- NULL
fov_integrated$IST_cluster_annotation <- NULL

# 2- Comparisons of cell composition (frequencies and counts - direct based analysis) based on clinical annotations.

# Define order of CIPR clusters for ditto bar plots
cluster_order <- match(levels(fov_integrated@meta.data[['CIPR_ordered_clusters']]), metaLevels('CIPR_ordered_clusters', fov_integrated))
# Manually define color code for each cluster label
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

# Define order of metaclusters for ditto bar plots
Idents(fov_integrated) <- "metaclusters"
levels(fov_integrated)
fov_integrated$metaclusters_ordered <- factor(fov_integrated$metaclusters,levels=c("Myeloid", "T cell", "B cell", "Stromal", "Epithelial"))

Idents(fov_integrated) <- 'metaclusters_ordered'
levels(fov_integrated)
metacluster_order <- match(levels(fov_integrated@meta.data[['metaclusters_ordered']]), metaLevels('metaclusters_ordered', fov_integrated))

metacluster_color <- setNames(c("#00924C", '#FDC77F', '#013F89', "#7853FF", '#AD0791'), c("Epithelial", "Stromal", "Myeloid", "T cell",  "B cell"))

# 2.1- Comparisons based on epithelial cluster
Idents(fov_integrated) <- 'Epithelial_cluster'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "Epithelial_cluster",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "Epithelial_cluster",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.2- Comparisons based on TME cluster
Idents(fov_integrated) <- 'TME_cluster'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "TME_cluster",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "TME_cluster",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.3- Comparisons based on aSMA cluster
Idents(fov_integrated) <- 'aSMA_cluster'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "aSMA_cluster",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "aSMA_cluster",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.4- Comparisons based on TSP score
Idents(fov_integrated) <- 'TSP_score'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "TSP_score",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "TSP_score",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.5- Comparisons based on KM score
Idents(fov_integrated) <- 'KM_score'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "KM_score",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "KM_score",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.6- Comparisons based on CMS
Idents(fov_integrated) <- 'CMS'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "CMS",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "CMS",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.7- Comparisons based on MMR status
Idents(fov_integrated) <- 'MMR_Status_Jen'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "MMR_Status_Jen",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "MMR_Status_Jen",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.8- Comparisons based on recurrence location
Idents(fov_integrated) <- 'Recurrence_location'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "Recurrence_location",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "Recurrence_location",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.9- Comparisons based on recurrence code
Idents(fov_integrated) <- 'Recurrence_code'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "Recurrence_code",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "Recurrence_code",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.10- Comparisons based on recurrence - reorder this
Idents(fov_integrated) <- 'Recurrence'
levels(fov_integrated)
fov_integrated$Recurrence_ordered <- factor(fov_integrated$Recurrence,levels=c("No recurrence", "Local", "Single-site: Liver", "Single-site: Brain", "Multi-site"))

#fov_integrated$Recurrence_ordered <- NULL

Idents(fov_integrated) <- 'Recurrence_ordered'
levels(fov_integrated)
recurrence_ordered <- match(levels(fov_integrated@meta.data[['Recurrence_ordered']]), metaLevels('Recurrence_ordered', fov_integrated))

dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "Recurrence_ordered",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "Recurrence",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.9- Comparisons based on T stage
Idents(fov_integrated) <- 'T_stage'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "T_stage",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "T_stage",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

# 2.9- Comparisons based on N stage
Idents(fov_integrated) <- 'N_stage'
levels(fov_integrated)
dittoBarPlot(fov_integrated, var = "CIPR_ordered_clusters", group.by = "N_stage",  scale = c("percent"), color.panel = cell_type_color, var.labels.reorder = cluster_order) 
dittoBarPlot(fov_integrated, var = "metaclusters_ordered", group.by = "N_stage",  scale = c("percent"), color.panel = metacluster_color, var.labels.reorder = metacluster_order)

############################################################################################################################################################################################################################################################################################################################################################################################

# 3- DGE Analysis - at single-cell level

# Reference: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions
# Reference: https://www.sc-best-practices.org/conditions/differential_gene_expression.html#motivation

# Generally, both, pseudobulk methods with sum aggregation such as edgeR, DESeq2, or Limma[Ritchie et al., 2015] and mixed models such as MAST with random effect setting were found to be superior compared to naive methods, 
# such as the popular Wilcoxon rank-sum test or Seurat’s [Hao et al., 2021] latent models, which do not account for them[Junttila et al., 2022].


DefaultAssay(fov_integrated) <- "RNA"
Idents(fov_integrated) <- 'Epithelial_cluster'
levels(fov_integrated)

# 3.1- DGE Epithelial cluster 3 vs 1

# 3.1.1- DGE between CD8 T cells
fov_integrated$DGE_epithelial_cluster <- paste(fov_integrated$CIPR_ordered_clusters, fov_integrated$Epithelial_cluster, sep = "_")
Idents(fov_integrated) <- "DGE_epithelial_cluster"
Idents(fov_integrated)
levels(fov_integrated)

# Because this is a imaging-based approach I will use the normalised RNA data slot for DGE instead the SCT data. Similar to original COSMX nat biotech paper
#fov_integrated <- PrepSCTFindMarkers(fov_integrated)
#DGE_CD8 <- FindMarkers(fov_integrated, assay = "SCT", ident.1 = "CD8+ T cell_3", ident.2 = "CD8+ T cell_1", verbose = FALSE)
#clipr::write_clip(DGE_CD8)

DGE_CD8_v2 <- FindMarkers(fov_integrated, assay = "RNA", slot = "counts", ident.1 = "CD8+ T cell_3", ident.2 = "CD8+ T cell_1", test.use = "MAST", verbose = FALSE)
clipr::write_clip(DGE_CD8_v2)

DGE1 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_1/DGE_CD8_Epithelial_cluster_3_vs_1.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE1$Features, pCutoff = 0.01,  FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 Epithelial cluster 3 vs Epithelial cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'), 
                legendPosition = 'right', legendLabSize = 8.0,  legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open",  colConnectors = 'white') + theme_light()



# 3.1.2- DGE between CD4 T cells
DGE_CD4 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD4+ T cell_3", ident.2 = "CD4+ T cell_1", verbose = FALSE)
clipr::write_clip(DGE_CD4)

DGE2 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_1/DGE_CD4_Epithelial_cluster_3_vs_1.csv")
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(toptable = DGE2, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE2$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 Epithelial cluster 3 vs Epithelial cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.1.3- DGE between Neutrophils - there is no neutrophils in epithelial cluster 3

# 3.1.4- DGE between SPP1+ macrophages
DGE_Macro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ macrophage_3", ident.2 = "SPP1+ macrophage_1", verbose = FALSE)
clipr::write_clip(DGE_Macro)

options(ggrepel.max.overlaps = Inf)
DGE3 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_1/DGE_Macrophage_Epithelial_cluster_3_vs_1.csv")
EnhancedVolcano(toptable = DGE3, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE3$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ macrophage Epithelial cluster 3 vs Epithelial cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.1.5- DGE between SPP1+ monocytes
DGE_Mono <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ monocyte_3", ident.2 = "SPP1+ monocyte_1", verbose = FALSE)
clipr::write_clip(DGE_Mono)

options(ggrepel.max.overlaps = Inf)
DGE4 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_1/DGE_Monocyte_Epithelial_cluster_3_vs_1.csv")
EnhancedVolcano(toptable = DGE4, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE4$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ macrophage Epithelial cluster 3 vs Epithelial cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.1.6- DGE between Tumour epithelial cell
DGE_TEC <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Tumor epithelial cell_3", ident.2 = "Tumor epithelial cell_1", verbose = FALSE)
clipr::write_clip(DGE_TEC)

options(ggrepel.max.overlaps = Inf)
DGE1 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_1/DGE_TEC_Epithelial_cluster_3_vs_1.csv")
EnhancedVolcano(toptable = DGE1, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE1$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumour eptihelial cell - Epithelial cluster 3 vs Epithelial cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()


# 3.2- DGE Epithelial cluster 3 vs 2

# 3.2.1- DGE between CD8 T cells
fov_integrated$DGE_epithelial_cluster <- paste(fov_integrated$CIPR_ordered_clusters, fov_integrated$Epithelial_cluster, sep = "_")
Idents(fov_integrated) <- "DGE_epithelial_cluster"
Idents(fov_integrated)
levels(fov_integrated)

DGE_CD8 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD8+ T cell_3", ident.2 = "CD8+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD8)

DGE1 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_2/DGE_CD8_Epithelial_cluster_3_vs_2.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE1$Features, pCutoff = 0.01,  FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 Epithelial cluster 3 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'), 
                legendPosition = 'right', legendLabSize = 8.0,  legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open",  colConnectors = 'white') + theme_light()

# 3.2.2- DGE between CD4 T cells
DGE_CD4 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD4+ T cell_3", ident.2 = "CD4+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD4)

options(ggrepel.max.overlaps = 50000)
DGE2 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_2/DGE_CD4_Epithelial_cluster_3_vs_2.csv")
EnhancedVolcano(toptable = DGE2, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE2$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 Epithelial cluster 3 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.2.3- DGE between Neutrophils - there is no neutrophils in epithelial cluster 3

# 3.2.4- DGE between SPP1+ macrophages
DGE_Macro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ macrophage_3", ident.2 = "SPP1+ macrophage_2", verbose = FALSE)
clipr::write_clip(DGE_Macro)

options(ggrepel.max.overlaps = Inf)
DGE3 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_2/DGE_Macro_Epithelial_cluster_3_vs_2.csv")
EnhancedVolcano(toptable = DGE3, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE3$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ macrophage Epithelial cluster 3 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.2.5- DGE between SPP1+ monocytes
DGE_Mono <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ monocyte_3", ident.2 = "SPP1+ monocyte_2", verbose = FALSE)
clipr::write_clip(DGE_Mono)

options(ggrepel.max.overlaps = Inf)
DGE4 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_2/DGE_Mono_Epithelial_cluster_3_vs_2.csv")
EnhancedVolcano(toptable = DGE4, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE4$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ monocyte Epithelial cluster 3 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.2.6- DGE between Tumour epithelial cell
DGE_TEC <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Tumor epithelial cell_3", ident.2 = "Tumor epithelial cell_2", verbose = FALSE)
clipr::write_clip(DGE_TEC)

options(ggrepel.max.overlaps = Inf)
DGE5 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/3_vs_2/DGE_TEC_Epithelial_cluster_3_vs_2.csv")
EnhancedVolcano(toptable = DGE5, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE5$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumour eptihelial cell - Epithelial cluster 3 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.3- DGE Epithelial cluster 1 vs 2

# 3.3.1- DGE between CD8 T cells
Idents(fov_integrated) <- "DGE_epithelial_cluster"
Idents(fov_integrated)
levels(fov_integrated)

DGE_CD8 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD8+ T cell_1", ident.2 = "CD8+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD8)

DGE1 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/1_vs_2/DGE_CD8_Epithelial_cluster_1_vs_2.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE1$Features, pCutoff = 0.01,  FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 Epithelial cluster 1 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'), 
                legendPosition = 'right', legendLabSize = 8.0,  legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open",  colConnectors = 'white') + theme_light()

# 3.3.2- DGE between CD4 T cells
DGE_CD4 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD4+ T cell_1", ident.2 = "CD4+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD4)

options(ggrepel.max.overlaps = Inf)
DGE2 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/1_vs_2/DGE_CD4_Epithelial_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE2, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE2$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 Epithelial cluster 3 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.3.3- DGE between Neutrophils - there is no neutrophils in epithelial cluster 3
DGE_Neutro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Neutrophil_1", ident.2 = "Neutrophil_2", verbose = FALSE)
clipr::write_clip(DGE_Neutro)

options(ggrepel.max.overlaps = Inf)
DGE6 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/1_vs_2/DGE_Neutrophil_Epithelial_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE6, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE6$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Neutrophil Epithelial cluster 1 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.3.4- DGE between SPP1+ macrophages
DGE_Macro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ macrophage_1", ident.2 = "SPP1+ macrophage_2", verbose = FALSE)
clipr::write_clip(DGE_Macro)

options(ggrepel.max.overlaps = Inf)
DGE3 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/1_vs_2/DGE_Macrophage_Epithelial_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE3, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE3$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ macrophage Epithelial cluster 1 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.3.5- DGE between SPP1+ monocytes
DGE_Mono <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ monocyte_1", ident.2 = "SPP1+ monocyte_2", verbose = FALSE)
clipr::write_clip(DGE_Mono)

options(ggrepel.max.overlaps = Inf)
DGE4 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/1_vs_2/DGE_Monocyte_Epithelial_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE4, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE4$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ monocyte Epithelial cluster 1 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.3.6- DGE between Tumour epithelial cell
DGE_TEC <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Tumor epithelial cell_1", ident.2 = "Tumor epithelial cell_2", verbose = FALSE)
clipr::write_clip(DGE_TEC)

options(ggrepel.max.overlaps = Inf)
DGE5 <- read.csv("./Plots/10_DGE/1_Epithelial cluster/1_vs_2/DGE_TEC_Epithelial_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE5, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE5$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumour eptihelial cell - Epithelial cluster 1 vs Epithelial cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()


############################################################################################################################################################################################################################################################################################################################################################################################

# 4- DGE Analysis TME clusters

DefaultAssay(fov_integrated) <- "RNA"
Idents(fov_integrated) <- 'TME_cluster'
levels(fov_integrated)

# 4.1- DGE TME cluster 1 vs 2

# 4.1.1- DGE between CD8 T cells
fov_integrated$DGE_TME_cluster <- paste(fov_integrated$CIPR_ordered_clusters, fov_integrated$TME_cluster, sep = "_")
Idents(fov_integrated) <- "DGE_TME_cluster"
Idents(fov_integrated)
levels(fov_integrated)

DGE_CD8 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD8+ T cell_1", ident.2 = "CD8+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD8)

DGE1 <- read.csv("./Plots/10_DGE/2_TME cluster/1_vs_2/DGE_CD8_TME_cluster_1_vs_2.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE1$Features, pCutoff = 0.01,  FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 TME cluster 1 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'), 
                legendPosition = 'right', legendLabSize = 8.0,  legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open",  colConnectors = 'white') + theme_light()

# 4.1.2- DGE between CD4 T cells
DGE_CD4 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD4+ T cell_1", ident.2 = "CD4+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD4)

DGE2 <- read.csv("./Plots/10_DGE/2_TME cluster/1_vs_2/DGE_CD4_TME_cluster_1_vs_2.csv")
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(toptable = DGE2, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE2$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 TME cluster 1 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.1.3- DGE between Neutrophils - there is no neutrophils in epithelial cluster 3
DGE_Neutro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Neutrophil_1", ident.2 = "Neutrophil_2", verbose = FALSE)
clipr::write_clip(DGE_Neutro)

options(ggrepel.max.overlaps = Inf)
DGE3 <- read.csv("./Plots/10_DGE/2_TME cluster/1_vs_2/DGE_Neutrophil_TME_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE3, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE3$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Neutrophil TME cluster 1 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.1.4- DGE between SPP1+ macrophages
DGE_Macro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ macrophage_1", ident.2 = "SPP1+ macrophage_2", verbose = FALSE)
clipr::write_clip(DGE_Macro)

options(ggrepel.max.overlaps = Inf)
DGE4 <- read.csv("./Plots/10_DGE/2_TME cluster/1_vs_2/DGE_Macrophage_TME_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE4, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE4$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ macrophage TME cluster 1 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.1.5- DGE between SPP1+ monocytes
DGE_Mono <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ monocyte_1", ident.2 = "SPP1+ monocyte_2", verbose = FALSE)
clipr::write_clip(DGE_Mono)

options(ggrepel.max.overlaps = Inf)
DGE5 <- read.csv("./Plots/10_DGE/2_TME cluster/1_vs_2/DGE_Monocyte_TME_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE5, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE5$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ monocyte TME cluster 1 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.1.6- DGE between Tumour epithelial cell
DGE_TEC <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Tumor epithelial cell_1", ident.2 = "Tumor epithelial cell_2", verbose = FALSE)
clipr::write_clip(DGE_TEC)

options(ggrepel.max.overlaps = Inf)
DGE6 <- read.csv("./Plots/10_DGE/2_TME cluster/1_vs_2/DGE_TEC_TME_cluster_1_vs_2.csv")
EnhancedVolcano(toptable = DGE6, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE6$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumour eptihelial cell - TME cluster 1 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()


# 4.2- DGE TME cluster 3 vs 1

# 4.2.1- DGE between CD8 T cells
Idents(fov_integrated) <- "DGE_TME_cluster"
Idents(fov_integrated)
levels(fov_integrated)

DGE_CD8 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD8+ T cell_3", ident.2 = "CD8+ T cell_1", verbose = FALSE)
clipr::write_clip(DGE_CD8)

DGE1 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_1/DGE_CD8_TME_cluster_3_vs_1.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE1$Features, pCutoff = 0.01,  FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 TME cluster 3 vs TME cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'), 
                legendPosition = 'right', legendLabSize = 8.0,  legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open",  colConnectors = 'white') + theme_light()

# 4.2.2- DGE between CD4 T cells
DGE_CD4 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD4+ T cell_3", ident.2 = "CD4+ T cell_1", verbose = FALSE)
clipr::write_clip(DGE_CD4)

DGE2 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_1/DGE_CD4_TME_cluster_3_vs_1.csv")
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(toptable = DGE2, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE2$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 TME cluster 3 vs TME cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.2.3- DGE between Neutrophils - there is no neutrophils in TME cluster 3

# 4.2.4- DGE between SPP1+ macrophages
DGE_Macro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ macrophage_3", ident.2 = "SPP1+ macrophage_1", verbose = FALSE)
clipr::write_clip(DGE_Macro)

options(ggrepel.max.overlaps = Inf)
DGE4 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_1/DGE_Macrophage_TME_cluster_3_vs_1.csv")
EnhancedVolcano(toptable = DGE4, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE4$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ macrophage TME cluster 3 vs TME cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.2.5- DGE between SPP1+ monocytes
DGE_Mono <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ monocyte_3", ident.2 = "SPP1+ monocyte_1", verbose = FALSE)
clipr::write_clip(DGE_Mono)

options(ggrepel.max.overlaps = Inf)
DGE5 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_1/DGE_Monocyte_TME_cluster_3_vs_1.csv")
EnhancedVolcano(toptable = DGE5, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE5$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ monocyte TME cluster 3 vs TME cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 3.1.6- DGE between Tumour epithelial cell
DGE_TEC <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Tumor epithelial cell_3", ident.2 = "Tumor epithelial cell_1", verbose = FALSE)
clipr::write_clip(DGE_TEC)

options(ggrepel.max.overlaps = Inf)
DGE6 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_1/DGE_TEC_TME_cluster_3_vs_1.csv")
EnhancedVolcano(toptable = DGE6, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE6$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumour eptihelial cell - TME cluster 3 vs TME cluster 1", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()


# 4.3- DGE Epithelial cluster 3 vs 2

# 4.3.1- DGE between CD8 T cells
Idents(fov_integrated) <- "DGE_TME_cluster"
Idents(fov_integrated)
levels(fov_integrated)

DGE_CD8 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD8+ T cell_3", ident.2 = "CD8+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD8)

DGE1 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_2/DGE_CD8_TME_cluster_3_vs_2.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE1$Features, pCutoff = 0.01,  FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD8 TME cluster 3 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'), 
                legendPosition = 'right', legendLabSize = 8.0,  legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open",  colConnectors = 'white') + theme_light()

# 4.3.2- DGE between CD4 T cells
DGE_CD4 <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "CD4+ T cell_3", ident.2 = "CD4+ T cell_2", verbose = FALSE)
clipr::write_clip(DGE_CD4)

DGE2 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_2/DGE_CD4_TME_cluster_3_vs_2.csv")
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(toptable = DGE2, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE2$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE CD4 TME cluster 3 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.3.3- DGE between Neutrophils - there is no neutrophils in epithelial cluster 3

# 4.3.4- DGE between SPP1+ macrophages
DGE_Macro <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ macrophage_3", ident.2 = "SPP1+ macrophage_2", verbose = FALSE)
clipr::write_clip(DGE_Macro)

options(ggrepel.max.overlaps = Inf)
DGE4 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_2/DGE_Macrophage_TME_cluster_3_vs_2.csv")
EnhancedVolcano(toptable = DGE4, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE4$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ macrophage TME cluster 3 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.3.5- DGE between SPP1+ monocytes
DGE_Mono <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "SPP1+ monocyte_3", ident.2 = "SPP1+ monocyte_2", verbose = FALSE)
clipr::write_clip(DGE_Mono)

options(ggrepel.max.overlaps = Inf)
DGE5 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_2/DGE_Monocyte_TME_cluster_3_vs_2.csv")
EnhancedVolcano(toptable = DGE5, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE5$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE SPP1+ monocyte TME cluster 3 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

# 4.3.6- DGE between Tumour epithelial cell
DGE_TEC <- FindMarkers(fov_integrated, assay = "RNA", ident.1 = "Tumor epithelial cell_3", ident.2 = "Tumor epithelial cell_2", verbose = FALSE)
clipr::write_clip(DGE_TEC)

options(ggrepel.max.overlaps = Inf)
DGE6 <- read.csv("./Plots/10_DGE/2_TME cluster/3_vs_2/DGE_TEC_TME_cluster_3_vs_2.csv")
EnhancedVolcano(toptable = DGE6, x = "avg_log2FC", y = "p_val_adj",
                lab = DGE6$Features, pCutoff = 0.01, FCcutoff = 1, pointSize = 3, labSize = 4.0, labCol = 'black', boxedLabels = FALSE, 
                title = "DGE Tumour eptihelial cell - TME cluster 3 vs TME cluster 2", legendLabels = c('Not significant', 'Fold change (but do not pass padj cutoff)', 'Pass padj cutoff', 'Pass both padj & fold change'),
                legendPosition = 'right', legendLabSize = 8.0, legendIconSize = 2.0, drawConnectors = TRUE, widthConnectors = 0.1, typeConnectors = "open", colConnectors = 'white') + theme_light()

############################################################################################################################################################################################################################################################################################################################################################################################
#5- Extract data to generate Anndata object

DefaultAssay(fov_integrated) <- "RNA"

X   = t(GetAssayData(fov_integrated)) # log normalised RNA data slot
write.csv(X, "fov_integrated_exprMat_file.csv")

obs = fov_integrated[[]] # metadata file
write.csv(obs, "fov_integrated_metadata_file.csv")

###########################################################################################################################################################################################################################################

# 6-Subset seurat object according to recurrences:

# https://github.com/satijalab/seurat/issues/6409#issuecomment-1568808941

library(magrittr)
library(Seurat)
source("./subset_obj_seurat.R") # load the intermediate solution function

DefaultAssay(fov_integrated) <- "RNA"
Idents(fov_integrated) <- "Recurrence"
levels(fov_integrated)

fov_rec1 <- subset_opt(fov_integrated, subset = Recurrence == "Negative")
fov_rec1

fov_local <- subset_opt(fov_integrated, subset = Recurrence == "Local")
fov_local

fov_liver <- subset_opt(fov_integrated, subset = Recurrence == "Liver")
fov_liver

fov_brain <- subset_opt(fov_integrated, subset = Recurrence == "Brain")
fov_brain

fov_multisite <- subset_opt(fov_integrated, subset = Recurrence == "Multi_site")
fov_multisite

# checks..
fov_integrated$sample_ID %>% table
Images(fov_integrated) # all FOVs of the merged object
Images(fov_rec1) # subsetted/kept FOVs
fov_rec1[[Images(fov_rec1)[i]]] %>% names # FOV classes
fov_rec1[[Images(fov_rec1)[i]]][["centroids"]] %>% str # cell centroids
fov_rec1[[Images(fov_rec1)[i]]][["centroids"]] %>% Cells %>% str # cell ids

# Convert seurat subsetted objected to Anndata

library(sceasy)
library(reticulate)
use_condaenv("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate")
use_python("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/bin/python")

sceasy::convertFormat(fov_integrated, from="seurat", to="anndata",
                      outFile='adata_fov_integrated.h5ad')


















