

setwd("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/COSMX_colon_dataset/COSMX_Colonrectal_cancer_project/TMA_original/imcrtools")


#3.1- Testing some plotting from Bodernmiller's tutorial using the SingleCellExperiment
library(imcRtools)
library(cytomapper)
library(openxlsx)
library(stringr)
library(dittoSeq)
library(RColorBrewer)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(bluster)
library(BiocParallel)
library(ggplot2)
library(scran)
library(CATALYST)
library(kohonen)
library(ConsensusClusterPlus)
library(patchwork)
library(pheatmap)
library(gridExtra)
library(SingleCellExperiment)
library(SpatialExperiment)
library(tidyverse)
library(ggridges)
library(scater)
library(cowplot)
library(viridis)
library(scMerge)

# 1- Read in the spatial experiment object
# 1.1- steinbock generated data
sce1 <- readRDS("./1_RDS_files/ad_fov01_sce.rds")
sce2 <- readRDS("./1_RDS_files/ad_fov02_sce.rds")
sce3 <- readRDS("./1_RDS_files/ad_fov03_sce.rds")
sce4 <- readRDS("./1_RDS_files/ad_fov04_sce.rds")
sce5 <- readRDS("./1_RDS_files/ad_fov05_sce.rds")
sce6 <- readRDS("./1_RDS_files/ad_fov06_sce.rds")
sce7 <- readRDS("./1_RDS_files/ad_fov07_sce.rds")
sce8 <- readRDS("./1_RDS_files/ad_fov08_sce.rds")
sce9 <- readRDS("./1_RDS_files/ad_fov09_sce.rds")
sce10 <- readRDS("./1_RDS_files/ad_fov10_sce.rds")
sce11 <- readRDS("./1_RDS_files/ad_fov11_sce.rds")
sce12 <- readRDS("./1_RDS_files/ad_fov12_sce.rds")
sce13 <- readRDS("./1_RDS_files/ad_fov13_sce.rds")
sce14 <- readRDS("./1_RDS_files/ad_fov14_sce.rds")
sce15 <- readRDS("./1_RDS_files/ad_fov15_sce.rds")
sce16 <- readRDS("./1_RDS_files/ad_fov16_sce.rds")
sce17 <- readRDS("./1_RDS_files/ad_fov17_sce.rds")
sce18 <- readRDS("./1_RDS_files/ad_fov18_sce.rds")
sce19 <- readRDS("./1_RDS_files/ad_fov19_sce.rds")
sce20 <- readRDS("./1_RDS_files/ad_fov20_sce.rds")
sce21 <- readRDS("./1_RDS_files/ad_fov21_sce.rds")
sce22 <- readRDS("./1_RDS_files/ad_fov22_sce.rds")
sce23 <- readRDS("./1_RDS_files/ad_fov23_sce.rds")
sce24 <- readRDS("./1_RDS_files/ad_fov24_sce.rds")
sce25 <- readRDS("./1_RDS_files/ad_fov25_sce.rds")


cellcluster <- setNames(c("#00924C","#00EA7B","#D5FFEB","#783D0D",
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


metadata(sce1)$color_vectors$cellcluster <- cellcluster
metadata(sce2)$color_vectors$cellcluster <- cellcluster
metadata(sce3)$color_vectors$cellcluster <- cellcluster
metadata(sce4)$color_vectors$cellcluster <- cellcluster
metadata(sce5)$color_vectors$cellcluster <- cellcluster
metadata(sce6)$color_vectors$cellcluster <- cellcluster
metadata(sce7)$color_vectors$cellcluster <- cellcluster
metadata(sce8)$color_vectors$cellcluster <- cellcluster
metadata(sce9)$color_vectors$cellcluster <- cellcluster
metadata(sce10)$color_vectors$cellcluster <- cellcluster
metadata(sce11)$color_vectors$cellcluster <- cellcluster
metadata(sce12)$color_vectors$cellcluster <- cellcluster
metadata(sce13)$color_vectors$cellcluster <- cellcluster
metadata(sce14)$color_vectors$cellcluster <- cellcluster
metadata(sce15)$color_vectors$cellcluster <- cellcluster
metadata(sce16)$color_vectors$cellcluster <- cellcluster
metadata(sce17)$color_vectors$cellcluster <- cellcluster
metadata(sce18)$color_vectors$cellcluster <- cellcluster
metadata(sce19)$color_vectors$cellcluster <- cellcluster
metadata(sce20)$color_vectors$cellcluster <- cellcluster
metadata(sce21)$color_vectors$cellcluster <- cellcluster
metadata(sce22)$color_vectors$cellcluster <- cellcluster
metadata(sce23)$color_vectors$cellcluster <- cellcluster
metadata(sce24)$color_vectors$cellcluster <- cellcluster
metadata(sce25)$color_vectors$cellcluster <- cellcluster

sce1 <- buildSpatialGraph(sce1, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce2 <- buildSpatialGraph(sce2, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce3 <- buildSpatialGraph(sce3, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce4 <- buildSpatialGraph(sce4, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce5 <- buildSpatialGraph(sce5, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce6 <- buildSpatialGraph(sce6, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce7 <- buildSpatialGraph(sce7, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce8 <- buildSpatialGraph(sce8, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce9 <- buildSpatialGraph(sce9, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce10 <- buildSpatialGraph(sce10, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce11 <- buildSpatialGraph(sce11, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce12 <- buildSpatialGraph(sce12, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce13 <- buildSpatialGraph(sce13, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce14 <- buildSpatialGraph(sce14, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce15 <- buildSpatialGraph(sce15, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce16 <- buildSpatialGraph(sce16, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce17 <- buildSpatialGraph(sce17, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce18 <- buildSpatialGraph(sce18, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce19 <- buildSpatialGraph(sce19, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce20 <- buildSpatialGraph(sce20, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce21 <- buildSpatialGraph(sce21, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce22 <- buildSpatialGraph(sce22, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce23 <- buildSpatialGraph(sce23, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce24 <- buildSpatialGraph(sce24, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))
sce25 <- buildSpatialGraph(sce25, img_id = "sample_ID", type = "expansion", threshold = 60, coords = c("CenterX_global_px", "CenterY_global_px"))

sce_list <- list(sce1, sce2, sce3, sce4, sce5, sce6, sce7, sce8, sce9, sce10, sce11, sce12, sce13, sce14, sce15, sce16, sce17, sce18,
             sce19, sce20, sce21, sce22, sce23, sce24, sce25)

sce_cbind(
  sce_list,
  method = "intersect",
  cut_off_batch = 0.01,
  cut_off_overall = 0.01,
  exprs = c("X"),
  colData_names = TRUE,
  batch_names = TRUE
)


# Patch analysis T cells
sce1 <- patchDetection(sce1, patch_cells = sce1$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                      coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce2 <- patchDetection(sce2, patch_cells = sce2$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce3 <- patchDetection(sce3, patch_cells = sce3$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce4 <- patchDetection(sce4, patch_cells = sce4$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce5 <- patchDetection(sce5, patch_cells = sce5$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                     coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce6 <- patchDetection(sce6, patch_cells = sce6$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce7 <- patchDetection(sce7, patch_cells = sce7$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce8 <- patchDetection(sce8, patch_cells = sce8$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce9 <- patchDetection(sce9, patch_cells = sce9$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce10 <- patchDetection(sce10, patch_cells = sce10$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce11 <- patchDetection(sce11, patch_cells = sce11$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce12 <- patchDetection(sce12, patch_cells = sce12$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce13 <- patchDetection(sce13, patch_cells = sce13$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce14 <- patchDetection(sce14, patch_cells = sce14$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce15 <- patchDetection(sce15, patch_cells = sce15$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce16 <- patchDetection(sce16, patch_cells = sce16$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce17 <- patchDetection(sce17, patch_cells = sce17$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce18 <- patchDetection(sce18, patch_cells = sce18$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce19 <- patchDetection(sce19, patch_cells = sce19$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce20 <- patchDetection(sce20, patch_cells = sce20$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce21 <- patchDetection(sce21, patch_cells = sce21$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce22 <- patchDetection(sce22, patch_cells = sce22$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce23 <- patchDetection(sce23, patch_cells = sce23$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce24 <- patchDetection(sce24, patch_cells = sce24$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")
sce25 <- patchDetection(sce25, patch_cells = sce25$metaclusters == "T cell", img_id = "sample_ID", expand_by = 60, min_patch_size = 3, #change this
                       coords = c("CenterX_global_px", "CenterY_global_px"), colPairName = "expansion_interaction_graph")

#We can now measure the size of each patch using the patchSize function and visualize tumor patch distribution per patient.
patch_size <- patchSize(sce1, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size <- merge(patch_size, colData(sce1)[match(patch_size$patch_id, sce1$patch_id),], by = "patch_id")

patch_size2 <- patchSize(sce2, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size2 <- merge(patch_size2, colData(sce2)[match(patch_size2$patch_id, sce2$patch_id),], by = "patch_id")

patch_size3 <- patchSize(sce3, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size3 <- merge(patch_size3, colData(sce3)[match(patch_size3$patch_id, sce3$patch_id),], by = "patch_id")

patch_size4 <- patchSize(sce4, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size4 <- merge(patch_size4, colData(sce4)[match(patch_size4$patch_id, sce4$patch_id),], by = "patch_id")

patch_size5 <- patchSize(sce5, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size5 <- merge(patch_size5, colData(sce5)[match(patch_size5$patch_id, sce5$patch_id),], by = "patch_id")

patch_size6 <- patchSize(sce6, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size6 <- merge(patch_size6, colData(sce6)[match(patch_size6$patch_id, sce6$patch_id),], by = "patch_id")

patch_size7 <- patchSize(sce7, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size7 <- merge(patch_size7, colData(sce7)[match(patch_size7$patch_id, sce7$patch_id),], by = "patch_id")

patch_size8 <- patchSize(sce8, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size8 <- merge(patch_size8, colData(sce8)[match(patch_size8$patch_id, sce8$patch_id),], by = "patch_id")

patch_size9 <- patchSize(sce9, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size9 <- merge(patch_size9, colData(sce9)[match(patch_size9$patch_id, sce9$patch_id),], by = "patch_id")

patch_size10 <- patchSize(sce10, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size10 <- merge(patch_size10, colData(sce10)[match(patch_size10$patch_id, sce10$patch_id),], by = "patch_id")

patch_size11 <- patchSize(sce11, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size11 <- merge(patch_size11, colData(sce11)[match(patch_size11$patch_id, sce11$patch_id),], by = "patch_id")

patch_size12 <- patchSize(sce12, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size12 <- merge(patch_size12, colData(sce12)[match(patch_size12$patch_id, sce12$patch_id),], by = "patch_id")

patch_size13 <- patchSize(sce13, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size13 <- merge(patch_size13, colData(sce13)[match(patch_size13$patch_id, sce13$patch_id),], by = "patch_id")

patch_size14 <- patchSize(sce14, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size14 <- merge(patch_size14, colData(sce14)[match(patch_size14$patch_id, sce14$patch_id),], by = "patch_id")

patch_size15 <- patchSize(sce15, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size15 <- merge(patch_size15, colData(sce15)[match(patch_size15$patch_id, sce15$patch_id),], by = "patch_id")

patch_size16 <- patchSize(sce16, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size16 <- merge(patch_size16, colData(sce16)[match(patch_size16$patch_id, sce16$patch_id),], by = "patch_id")

patch_size17 <- patchSize(sce17, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size17 <- merge(patch_size17, colData(sce17)[match(patch_size17$patch_id, sce17$patch_id),], by = "patch_id")

patch_size18 <- patchSize(sce18, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size18 <- merge(patch_size18, colData(sce18)[match(patch_size18$patch_id, sce18$patch_id),], by = "patch_id")

patch_size19 <- patchSize(sce19, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size19 <- merge(patch_size19, colData(sce19)[match(patch_size19$patch_id, sce19$patch_id),], by = "patch_id")

patch_size20 <- patchSize(sce20, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size20 <- merge(patch_size20, colData(sce20)[match(patch_size20$patch_id, sce20$patch_id),], by = "patch_id")

patch_size21 <- patchSize(sce21, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size21 <- merge(patch_size21, colData(sce21)[match(patch_size21$patch_id, sce21$patch_id),], by = "patch_id")

patch_size22 <- patchSize(sce22, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size22 <- merge(patch_size22, colData(sce22)[match(patch_size22$patch_id, sce22$patch_id),], by = "patch_id")

patch_size23 <- patchSize(sce23, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size23 <- merge(patch_size23, colData(sce23)[match(patch_size23$patch_id, sce23$patch_id),], by = "patch_id")

patch_size24 <- patchSize(sce24, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size24 <- merge(patch_size24, colData(sce24)[match(patch_size24$patch_id, sce24$patch_id),], by = "patch_id")

patch_size25 <- patchSize(sce25, "patch_id", coords = c("CenterX_global_px", "CenterY_global_px"))
patch_size25 <- merge(patch_size25, colData(sce25)[match(patch_size25$patch_id, sce25$patch_id),], by = "patch_id")


patch_size <- as.data.frame(patch_size)
patch_size2 <- as.data.frame(patch_size2)
patch_size3 <- as.data.frame(patch_size3)
patch_size4 <- as.data.frame(patch_size4)
patch_size5 <- as.data.frame(patch_size5)
patch_size6 <- as.data.frame(patch_size6)
patch_size7 <- as.data.frame(patch_size7)
patch_size8 <- as.data.frame(patch_size8)
patch_size9 <- as.data.frame(patch_size9)
patch_size10 <- as.data.frame(patch_size10)
patch_size11 <- as.data.frame(patch_size11)
patch_size12 <- as.data.frame(patch_size12)
patch_size13 <- as.data.frame(patch_size13)
patch_size14 <- as.data.frame(patch_size14)
patch_size15 <- as.data.frame(patch_size15)
patch_size16 <- as.data.frame(patch_size16)
patch_size17 <- as.data.frame(patch_size17)
patch_size18 <- as.data.frame(patch_size18)
patch_size19 <- as.data.frame(patch_size19)
patch_size20 <- as.data.frame(patch_size20)
patch_size21 <- as.data.frame(patch_size21)
patch_size22 <- as.data.frame(patch_size22)
patch_size23 <- as.data.frame(patch_size23)
patch_size24 <- as.data.frame(patch_size24)
patch_size25 <- as.data.frame(patch_size25)

patch_size_merge <- rbind(patch_size, patch_size2, patch_size3, patch_size4, patch_size, patch_size2, patch_size3, patch_size4,
                          patch_size5, patch_size6, patch_size7, patch_size8, patch_size9, patch_size10, patch_size11, patch_size12,
                          patch_size13, patch_size14, patch_size15, patch_size16, patch_size17, patch_size18, patch_size19, patch_size20,
                          patch_size21, patch_size22, patch_size23, patch_size24, patch_size25)

ggplot(patch_size_merge) + 
  geom_boxplot(aes(sample_ID, log10(size))) +
  geom_point(aes(sample_ID, log10(size))) + theme(axis.text.x = element_text(angle=90))

library(dplyr)

patch_size_merge <- patch_size_merge %>% select_if(~ !any(is.na(.)))
write_csv(patch_size_merge,"Tcell_patchsize.csv")



# Define SC color scheme
col_PD <- setNames(colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(sce1$patch_id))), 
                   sort(unique(sce1$spe$patch_id)))

plotSpatial(sce1, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
           theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce2, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce3, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce4, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce5, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce6, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce7, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce8, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce9, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce10, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce11, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce12, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce13, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce14, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce15, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce16, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)
plotSpatial(sce17, node_color_by = "patch_id", img_id = "sample_ID", 
            node_size_fix = 3, coords = c("CenterX_global_px", "CenterY_global_px")) +
  theme(legend.position = "none") + scale_color_manual(values = col_PD)






ggplot(as.data.frame(patch_size)) + 
  geom_boxplot(aes(Patient, log10(size))) +
  geom_point(aes(Patient, log10(size))) 

ggplot(as.data.frame(patch_size)) + 
  geom_boxplot(aes(Group, log10(size))) +
  geom_point(aes(Group, log10(size))) 


#We can now calculate the fraction of cells within each patch to roughly estimate cell infiltration.
library(tidyverse)
colData(sce1) %>% as_tibble() %>%
  group_by(patch_id, sample_ID) %>%
  summarize(cell_count = sum(CIPR_annotated_clusters == "CD8+ T cell" | CIPR_annotated_clusters == "CD4+ T cell"),
            patch_size = n(),
            cell_freq = cell_count / patch_size) %>%
  ggplot() +
  geom_point(aes(log10(patch_size), cell_freq, color = patch_id)) +
  theme_classic()

library(tidyverse)
colData(adNonCM) %>% as_tibble() %>%
  group_by(patch_id, Patient) %>%
  summarize(Tcell_count = sum(cell_type == "CD8 T cells" | cell_type == "CD4 T cells"),
            patch_size = n(),
            Tcell_freq = Tcell_count / patch_size) %>%
  ggplot() +
  geom_point(aes(log10(patch_size), Tcell_freq, color = Patient)) +
  theme_classic()

#We can now calculate the fraction of B cells within each lymphoid patch to roughly estimate B cell infiltration.
colData(adCM2) %>% as_tibble() %>%
  group_by(patch_id, Patient) %>%
  summarize(Bcell_count = sum(cell_type == "B cells"),
            patch_size = n(),
            Bcell_freq = Bcell_count / patch_size) %>%
  ggplot() +
  geom_point(aes(log10(patch_size), Bcell_freq, color = Patient)) +
  theme_classic()

colData(adNonCM) %>% as_tibble() %>%
  group_by(patch_id, Patient) %>%
  summarize(Bcell_count = sum(cell_type == "B cells"),
            patch_size = n(),
            Bcell_freq = Bcell_count / patch_size) %>%
  ggplot() +
  geom_point(aes(log10(patch_size), Bcell_freq, color = Patient)) +
  theme_classic()

#We can now calculate the fraction of Macrophages within each patch to roughly estimate infiltration.
colData(ad) %>% as_tibble() %>%
  group_by(patch_id, Patient) %>%
  summarize(DC_count = sum(cell_type == "Dendritic cells"),
            patch_size = n(),
            DC_freq = DC_count / patch_size) %>%
  ggplot() +
  geom_point(aes(log10(patch_size), DC_freq, color = Patient)) +
  theme_classic()










