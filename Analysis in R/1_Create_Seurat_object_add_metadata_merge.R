
setwd("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/COSMX_colon_dataset/COSMX_Colonrectal_cancer_project/TMA_original/Seurat")

library(Seurat)
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

saveRDS(fov_merged, "fov_merged.rds")
fov_merged <- readRDS("./fov_merged.rds")

saveRDS(fov_subset, "fov_subset.rds")
fov_subset <- readRDS("./fov_subset.rds")

# 1- Create Seurat object

# 1.1- with LoadNanostring

# 1.1.1- fov01
fov01 <- LoadNanostring(data.dir = "./fovs/fov01", fov = "fov01", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov01, "Fov01_SEURAT.rds")
fov01 <- readRDS("./RDS_files/1_Raw/Fov01_SEURAT.rds")

# The resulting object has 4,156 cells. But metadata contains 4,159 cells. Which cells are missing?
options(max.print=1000000)
metadata <- read.csv("./fov01/Fov01_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
polygons <- read.csv("./fov01/Fov01-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fov01/Fov01_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fov01/Fov01_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cells "446","1445","3005" are not in the levels of tx file
tx[tx$cell_ID == "446", ]  
tx[tx$cell_ID == "1445", ] 
tx[tx$cell_ID == "3005", ] 
colnames(fov01) #cells "446","1445","3005" are not in the fov01 Seurat object

# Add protein expression as another assay
# Remove cells "446","1445","3005" from the metadata file and save as "_edited"
fov01.ptn <- as.sparse(t(read.csv(file = "./fovs/fov01/Fov01_metadata_file_edited.csv", sep = ",",
                               header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov01.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov01), colnames(fov01.ptn))

ptn_assay <- CreateAssayObject(counts = fov01.ptn)
# add this assay to the previously created Seurat object
fov01[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov01)
DefaultAssay(fov01) <- "Protein"
ImageFeaturePlot(fov01, boundaries = "segmentation",fov = "fov01", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov01/Fov01_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov01@meta.data <-cbind(fov01@meta.data, IF_marker)
fov01@meta.data <-cbind(fov01@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                             c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov01, boundaries = "segmentation", fov = "fov01", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov01, "./RDS_files/2_Add_metadata/Fov01_SEURAT.rds")

# 1.1.2- fov02
fov02 <- LoadNanostring(data.dir = "./fovs/fov02", fov = "fov02", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov02, "Fov02_SEURAT.rds")
fov02 <- readRDS("./RDS_files/1_Raw/Fov02_SEURAT.rds")

# Add protein expression as another assay
fov02.ptn <- as.sparse(t(read.csv(file = "./fovs/fov02/Fov02_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov02.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov02), colnames(fov02.ptn))

ptn_assay <- CreateAssayObject(counts = fov02.ptn)
# add this assay to the previously created Seurat object
fov02[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov02)
DefaultAssay(fov02) <- "Protein"
ImageFeaturePlot(fov02, boundaries = "segmentation",fov = "fov02", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov02/Fov02_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov02@meta.data <-cbind(fov02@meta.data, IF_marker)
fov02@meta.data <-cbind(fov02@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov02, boundaries = "segmentation", fov = "fov02", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov02, "./RDS_files/2_Add_metadata/Fov02_SEURAT.rds")

# 1.1.3- fov03
fov03 <- LoadNanostring(data.dir = "./fovs/fov03", fov = "fov03", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov03, "./RDS_files/1_Raw/Fov03_SEURAT.rds")
fov03 <- readRDS("./RDS_files/1_Raw/Fov03_SEURAT.rds")

# Add protein expression as another assay
fov03.ptn <- as.sparse(t(read.csv(file = "./fovs/fov03/Fov03_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov03.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov03), colnames(fov03.ptn))

ptn_assay <- CreateAssayObject(counts = fov03.ptn)
# add this assay to the previously created Seurat object
fov03[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov03)
DefaultAssay(fov03) <- "Protein"
ImageFeaturePlot(fov03, boundaries = "segmentation",fov = "fov03", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov03/Fov03_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov03@meta.data <-cbind(fov03@meta.data, IF_marker)
fov03@meta.data <-cbind(fov03@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov03, boundaries = "segmentation", fov = "fov03", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov03, "./RDS_files/2_Add_metadata/Fov03_SEURAT.rds")

# 1.1.4- fov04
fov04 <- LoadNanostring(data.dir = "./fovs/fov04", fov = "fov04", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov04, "Fov04_SEURAT.rds")
fov04 <- readRDS("./RDS_files/1_Raw/Fov04_SEURAT.rds")

# Add protein expression as another assay
fov04.ptn <- as.sparse(t(read.csv(file = "./fovs/fov04/Fov04_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov04.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov04), colnames(fov04.ptn))

ptn_assay <- CreateAssayObject(counts = fov04.ptn)
# add this assay to the previously created Seurat object
fov04[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov04)
DefaultAssay(fov04) <- "Protein"
ImageFeaturePlot(fov04, boundaries = "segmentation",fov = "fov04", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov04/Fov04_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov04@meta.data <-cbind(fov04@meta.data, IF_marker)
fov04@meta.data <-cbind(fov04@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov04, boundaries = "segmentation", fov = "fov04", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov04, "./RDS_files/2_Add_metadata/Fov04_SEURAT.rds")

# 1.1.5- fov05
fov05 <- LoadNanostring(data.dir = "./fovs/fov05", fov = "fov05", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov05, "./RDS_files/1_Raw/Fov05_SEURAT.rds")
fov05 <- readRDS("./RDS_files/1_Raw/Fov05_SEURAT.rds")

# Add protein expression as another assay
fov05.ptn <- as.sparse(t(read.csv(file = "./fovs/fov05/Fov05_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov05.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov05), colnames(fov05.ptn))

ptn_assay <- CreateAssayObject(counts = fov05.ptn)
# add this assay to the previously created Seurat object
fov05[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov05)
DefaultAssay(fov05) <- "Protein"
ImageFeaturePlot(fov05, boundaries = "segmentation",fov = "fov05", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov05/Fov05_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov05@meta.data <-cbind(fov05@meta.data, IF_marker)
fov05@meta.data <-cbind(fov05@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov05, boundaries = "segmentation", fov = "fov05", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov05, "./RDS_files/2_Add_metadata/Fov05_SEURAT.rds")

# 1.1.6- fov06
fov06 <- LoadNanostring(data.dir = "./fovs/fov06", fov = "fov06", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov06, "./RDS_files/1_Raw/Fov06_SEURAT.rds")
fov06 <- readRDS("./RDS_files/1_Raw/Fov06_SEURAT.rds")

# Add protein expression as another assay
fov06.ptn <- as.sparse(t(read.csv(file = "./fovs/fov06/Fov06_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov06.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov06), colnames(fov06.ptn))

ptn_assay <- CreateAssayObject(counts = fov06.ptn)
# add this assay to the previously created Seurat object
fov06[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov06)
DefaultAssay(fov06) <- "Protein"
ImageFeaturePlot(fov06, boundaries = "segmentation",fov = "fov06", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov06/Fov06_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov06@meta.data <-cbind(fov06@meta.data, IF_marker)
fov06@meta.data <-cbind(fov06@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov06, boundaries = "segmentation", fov = "fov06", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov06, "./RDS_files/2_Add_metadata/Fov06_SEURAT.rds")

# 1.1.7- fov07
fov07 <- LoadNanostring(data.dir = "./fovs/fov07", fov = "fov07", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov07, "./RDS_files/1_Raw/Fov07_SEURAT.rds")
fov07 <- readRDS("./RDS_files/1_Raw/Fov07_SEURAT.rds")

# Add protein expression as another assay
metadata <- read.csv("./fov07/Fov07_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
polygons <- read.csv("./fov07/Fov07-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fov07/Fov07_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cells "367","1517","2037", "2534" are not in the levels of tx file
tx[tx$cell_ID == "367", ]  
tx[tx$cell_ID == "1517", ] 
tx[tx$cell_ID == "2037", ] 
colnames(fov07) # cells "367","1517","2037", "2534" are not in the fov01 Seurat object

fov07.ptn <- as.sparse(t(read.csv(file = "./fovs/fov07/Fov07_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov07.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov07), colnames(fov07.ptn))

ptn_assay <- CreateAssayObject(counts = fov07.ptn)
# add this assay to the previously created Seurat object
fov07[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov07)
DefaultAssay(fov07) <- "Protein"
ImageFeaturePlot(fov07, boundaries = "segmentation",fov = "fov07", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov07/Fov07_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov07@meta.data <-cbind(fov07@meta.data, IF_marker)
fov07@meta.data <-cbind(fov07@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov07, boundaries = "segmentation", fov = "fov07", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov07, "./RDS_files/2_Add_metadata/Fov07_SEURAT.rds")

# 1.1.8- fov08
fov08 <- LoadNanostring(data.dir = "./fovs/fov08", fov = "fov08", assay = "RNA")
saveRDS(fov08, "./RDS_files/1_Raw/Fov08_SEURAT.rds")
fov08 <- readRDS("./RDS_files/1_Raw/Fov08_SEURAT.rds")

# Add protein expression as another assay
fov08.ptn <- as.sparse(t(read.csv(file = "./fovs/fov08/Fov08_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov08.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov08), colnames(fov08.ptn))

ptn_assay <- CreateAssayObject(counts = fov08.ptn)
# add this assay to the previously created Seurat object
fov08[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov08)
DefaultAssay(fov08) <- "Protein"
ImageFeaturePlot(fov08, boundaries = "segmentation",fov = "fov08", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov08/Fov08_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov08@meta.data <-cbind(fov08@meta.data, IF_marker)
fov08@meta.data <-cbind(fov08@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov08, boundaries = "segmentation", fov = "fov08", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov08, "./RDS_files/2_Add_metadata/Fov08_SEURAT.rds")

# 1.1.9- fov09
fov09 <- LoadNanostring(data.dir = "./fovs/fov09", fov = "fov09", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov09, "./RDS_files/1_Raw/Fov09_SEURAT.rds")
fov09 <- readRDS("./RDS_files/1_Raw/Fov09_SEURAT.rds")

# Add protein expression as another assay
metadata <- read.csv("./fov09/Fov09_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
tx <- read.csv("./fov09/Fov09_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cell "1562" is not in the levels of tx file
tx[tx$cell_ID == "1562", ]  
colnames(fov09) # cell "1562"is not in the fov01 Seurat object

fov09.ptn <- as.sparse(t(read.csv(file = "./fovs/fov09/Fov09_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov09.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov09), colnames(fov09.ptn))

ptn_assay <- CreateAssayObject(counts = fov09.ptn)
# add this assay to the previously created Seurat object
fov09[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov09)
DefaultAssay(fov09) <- "Protein"
ImageFeaturePlot(fov09, boundaries = "segmentation",fov = "fov09", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov09/Fov09_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov09@meta.data <-cbind(fov09@meta.data, IF_marker)
fov09@meta.data <-cbind(fov09@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov09, boundaries = "segmentation", fov = "fov09", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov09, "./RDS_files/2_Add_metadata/Fov09_SEURAT.rds")

# 1.1.10- fov10
fov10 <- LoadNanostring(data.dir = "./fovs/fov10", fov = "fov10", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov10, "./RDS_files/1_Raw/Fov10_SEURAT.rds")
fov10 <- readRDS("./RDS_files/1_Raw/Fov10_SEURAT.rds")

# Add protein expression as another assay
fov10.ptn <- as.sparse(t(read.csv(file = "./fovs/fov10/Fov10_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov10.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov10), colnames(fov10.ptn))

ptn_assay <- CreateAssayObject(counts = fov10.ptn)
# add this assay to the previously created Seurat object
fov10[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov10)
DefaultAssay(fov10) <- "Protein"
ImageFeaturePlot(fov10, boundaries = "segmentation",fov = "fov10", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov10/Fov10_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov10@meta.data <-cbind(fov10@meta.data, IF_marker)
fov10@meta.data <-cbind(fov10@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov10, boundaries = "segmentation", fov = "fov10", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov10, "./RDS_files/2_Add_metadata/Fov10_SEURAT.rds")

# 1.1.11- fov11
fov11 <- LoadNanostring(data.dir = "./fovs/fov11", fov = "fov11", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov11, "./RDS_files/1_Raw/Fov11_SEURAT.rds")
fov11 <- readRDS("./RDS_files/1_Raw/Fov11_SEURAT.rds")

# Add protein expression as another assay
fov11.ptn <- as.sparse(t(read.csv(file = "./fovs/fov11/Fov11_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov11.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov11), colnames(fov11.ptn))

ptn_assay <- CreateAssayObject(counts = fov11.ptn)
# add this assay to the previously created Seurat object
fov11[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov11)
DefaultAssay(fov11) <- "Protein"
ImageFeaturePlot(fov11, boundaries = "segmentation",fov = "fov11", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov11/Fov11_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov11@meta.data <-cbind(fov11@meta.data, IF_marker)
fov11@meta.data <-cbind(fov11@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov11, boundaries = "segmentation", fov = "fov11", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov11, "./RDS_files/2_Add_metadata/Fov11_SEURAT.rds")

# 1.1.12- fov12
fov12 <- LoadNanostring(data.dir = "./fovs/fov12", fov = "fov12", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov12, "./RDS_files/1_Raw/Fov12_SEURAT.rds")
fov12 <- readRDS("./RDS_files/1_Raw/Fov12_SEURAT.rds")

# Add protein expression as another assay
fov12.ptn <- as.sparse(t(read.csv(file = "./fovs/fov12/Fov12_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov12.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov12), colnames(fov12.ptn))

ptn_assay <- CreateAssayObject(counts = fov12.ptn)
# add this assay to the previously created Seurat object
fov12[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov12)
DefaultAssay(fov12) <- "Protein"
ImageFeaturePlot(fov12, boundaries = "segmentation",fov = "fov12", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov12/Fov12_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov12@meta.data <-cbind(fov12@meta.data, IF_marker)
fov12@meta.data <-cbind(fov12@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov12, boundaries = "segmentation", fov = "fov12", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov12, "./RDS_files/2_Add_metadata/Fov12_SEURAT.rds")

# 1.1.13- fov13
fov13 <- LoadNanostring(data.dir = "./fovs/fov13", fov = "fov13", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov13, "./RDS_files/1_Raw/Fov13_SEURAT.rds")
fov13 <- readRDS("./RDS_files/1_Raw/Fov13_SEURAT.rds")

# Add protein expression as another assay
fov13.ptn <- as.sparse(t(read.csv(file = "./fovs/fov13/Fov13_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov13.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov13), colnames(fov13.ptn))

ptn_assay <- CreateAssayObject(counts = fov13.ptn)
# add this assay to the previously created Seurat object
fov13[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov13)
DefaultAssay(fov13) <- "Protein"
ImageFeaturePlot(fov13, boundaries = "segmentation",fov = "fov13", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov13/Fov13_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov13@meta.data <-cbind(fov13@meta.data, IF_marker)
fov13@meta.data <-cbind(fov13@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov13, boundaries = "segmentation", fov = "fov13", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov13, "./RDS_files/2_Add_metadata/Fov13_SEURAT.rds")

# 1.1.14- fov14
fov14 <- LoadNanostring(data.dir = "./fovs/fov14", fov = "fov14", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov14, "./RDS_files/1_Raw/Fov14_SEURAT.rds")
fov14 <- readRDS("./RDS_files/1_Raw/Fov14_SEURAT.rds")

# Add protein expression as another assay
metadata <- read.csv("./fov14/Fov14_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
tx <- read.csv("./fov14/Fov14_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cell "1001" is not in the levels of tx file
tx[tx$cell_ID == "1001", ]  
colnames(fov14) # cell "1001"is not in the fov01 Seurat object

fov14.ptn <- as.sparse(t(read.csv(file = "./fovs/fov14/Fov14_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov14.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov14), colnames(fov14.ptn))

ptn_assay <- CreateAssayObject(counts = fov14.ptn)
# add this assay to the previously created Seurat object
fov14[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov14)
DefaultAssay(fov14) <- "Protein"
ImageFeaturePlot(fov14, boundaries = "segmentation",fov = "fov14", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov14/Fov14_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov14@meta.data <-cbind(fov14@meta.data, IF_marker)
fov14@meta.data <-cbind(fov14@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov14, boundaries = "segmentation", fov = "fov14", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov14, "./RDS_files/2_Add_metadata/Fov14_SEURAT.rds")

# 1.1.15- fov15
fov15 <- LoadNanostring(data.dir = "./fovs/fov15", fov = "fov15", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov15, "./RDS_files/1_Raw/Fov15_SEURAT.rds")
fov15 <- readRDS("./RDS_files/1_Raw/Fov15_SEURAT.rds")

# Add protein expression as another assay
fov15.ptn <- as.sparse(t(read.csv(file = "./fovs/fov15/Fov15_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov15.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov15), colnames(fov15.ptn))

ptn_assay <- CreateAssayObject(counts = fov15.ptn)
# add this assay to the previously created Seurat object
fov15[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov15)
DefaultAssay(fov15) <- "Protein"
ImageFeaturePlot(fov15, boundaries = "segmentation",fov = "fov15", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov15/Fov15_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov15@meta.data <-cbind(fov15@meta.data, IF_marker)
fov15@meta.data <-cbind(fov15@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov15, boundaries = "segmentation", fov = "fov15", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov15, "./RDS_files/2_Add_metadata/Fov15_SEURAT.rds")

# 1.1.16- fov16
fov16 <- LoadNanostring(data.dir = "./fovs/fov16", fov = "fov16", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov16, "./RDS_files/1_Raw/Fov16_SEURAT.rds")
fov16 <- readRDS("./RDS_files/1_Raw/Fov16_SEURAT.rds")

# Add protein expression as another assay
metadata <- read.csv("./fovs/fov16/Fov16_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
tx <- read.csv("./fov16/Fov16_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cells "1488", "1489" not in the levels of tx file
tx[tx$cell_ID == "1488", ]  
colnames(fov16) # cells "1488", "1489" are not in the fov01 Seurat object

fov16.ptn <- as.sparse(t(read.csv(file = "./fovs/fov16/Fov16_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov16.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov16), colnames(fov16.ptn))

ptn_assay <- CreateAssayObject(counts = fov16.ptn)
# add this assay to the previously created Seurat object
fov16[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov16)
DefaultAssay(fov16) <- "Protein"
ImageFeaturePlot(fov16, boundaries = "segmentation",fov = "fov16", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov16/Fov16_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov16@meta.data <-cbind(fov16@meta.data, IF_marker)
fov16@meta.data <-cbind(fov16@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov16, boundaries = "segmentation", fov = "fov16", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov16, "./RDS_files/2_Add_metadata/Fov16_SEURAT.rds")

# 1.1.17- fov17
fov17 <- LoadNanostring(data.dir = "./fovs/fov17", fov = "fov17", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov17, "./RDS_files/1_Raw/Fov17_SEURAT.rds")
fov17 <- readRDS("./RDS_files/1_Raw/Fov17_SEURAT.rds")

# Add protein expression as another assay
fov17.ptn <- as.sparse(t(read.csv(file = "./fovs/fov17/Fov17_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov17.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov17), colnames(fov17.ptn))

ptn_assay <- CreateAssayObject(counts = fov17.ptn)
# add this assay to the previously created Seurat object
fov17[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov17)
DefaultAssay(fov17) <- "Protein"
ImageFeaturePlot(fov17, boundaries = "segmentation",fov = "fov17", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov17/Fov17_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov17@meta.data <-cbind(fov17@meta.data, IF_marker)
fov17@meta.data <-cbind(fov17@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov17, boundaries = "segmentation", fov = "fov17", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov17, "./RDS_files/2_Add_metadata/Fov17_SEURAT.rds")

# 1.1.18- fov18
fov18 <- LoadNanostring(data.dir = "./fovs/fov18", fov = "fov18", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov18, "./RDS_files/1_Raw/Fov18_SEURAT.rds")
fov18 <- readRDS("./RDS_files/1_Raw/Fov18_SEURAT.rds")

# Add protein expression as another assay
metadata <- read.csv("./fovs/fov18/Fov18_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
tx <- read.csv("./fov18/Fov18_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cells "2546", "2918" not in the levels of tx file
tx[tx$cell_ID == "2546", ]  
colnames(fov18) # cells "2546", "2918" not in the fov01 Seurat object

fov18.ptn <- as.sparse(t(read.csv(file = "./fovs/fov18/Fov18_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov18.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov18), colnames(fov18.ptn))

ptn_assay <- CreateAssayObject(counts = fov18.ptn)
# add this assay to the previously created Seurat object
fov18[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov18)
DefaultAssay(fov18) <- "Protein"
ImageFeaturePlot(fov18, boundaries = "segmentation",fov = "fov18", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov18/Fov18_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov18@meta.data <-cbind(fov18@meta.data, IF_marker)
fov18@meta.data <-cbind(fov18@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov18, boundaries = "segmentation", fov = "fov18", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov18, "./RDS_files/2_Add_metadata/Fov18_SEURAT.rds")

# 1.1.19- fov19
fov19 <- LoadNanostring(data.dir = "./fovs/fov19", fov = "fov19", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov19, "./RDS_files/1_Raw/Fov19_SEURAT.rds")
fov19 <- readRDS("./RDS_files/1_Raw/Fov19_SEURAT.rds")

# Add protein expression as another assay
fov19.ptn <- as.sparse(t(read.csv(file = "./fovs/fov19/Fov19_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov19.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov19), colnames(fov19.ptn))

ptn_assay <- CreateAssayObject(counts = fov19.ptn)
# add this assay to the previously created Seurat object
fov19[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov19)
DefaultAssay(fov19) <- "Protein"
ImageFeaturePlot(fov19, boundaries = "segmentation",fov = "fov19", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov19/Fov19_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov19@meta.data <-cbind(fov19@meta.data, IF_marker)
fov19@meta.data <-cbind(fov19@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov19, boundaries = "segmentation", fov = "fov19", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov19, "./RDS_files/2_Add_metadata/Fov19_SEURAT.rds")

# 1.1.20- fov20
fov20 <- LoadNanostring(data.dir = "./fovs/fov20", fov = "fov20", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov20, "./RDS_files/1_Raw/Fov20_SEURAT.rds")
fov20 <- readRDS("./RDS_files/1_Raw/Fov20_SEURAT.rds")

# Add protein expression as another assay
fov20.ptn <- as.sparse(t(read.csv(file = "./fovs/fov20/Fov20_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov20.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov20), colnames(fov20.ptn))

ptn_assay <- CreateAssayObject(counts = fov20.ptn)
# add this assay to the previously created Seurat object
fov20[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov20)
DefaultAssay(fov20) <- "Protein"
ImageFeaturePlot(fov20, boundaries = "segmentation",fov = "fov20", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov20/Fov20_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov20@meta.data <-cbind(fov20@meta.data, IF_marker)
fov20@meta.data <-cbind(fov20@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov20, boundaries = "segmentation", fov = "fov20", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov20, "./RDS_files/2_Add_metadata/Fov20_SEURAT.rds")

# 1.1.21- fov21
fov21 <- LoadNanostring(data.dir = "./fovs/fov21", fov = "fov21", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov21, "./RDS_files/1_Raw/Fov21_SEURAT.rds")
fov21 <- readRDS("./RDS_files/1_Raw/Fov21_SEURAT.rds")

# Add protein expression as another assay
metadata <- read.csv("./fovs/fov21/Fov21_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
tx <- read.csv("./fov21/Fov21_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cells "3681" not in the levels of tx file
tx[tx$cell_ID == "3681", ]  
colnames(fov21) # cell "3681"is not in the fov21 Seurat object

fov21.ptn <- as.sparse(t(read.csv(file = "./fovs/fov21/Fov21_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov21.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov21), colnames(fov21.ptn))

ptn_assay <- CreateAssayObject(counts = fov21.ptn)
# add this assay to the previously created Seurat object
fov21[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov21)
DefaultAssay(fov21) <- "Protein"
ImageFeaturePlot(fov21, boundaries = "segmentation",fov = "fov21", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov21/Fov21_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov21@meta.data <-cbind(fov21@meta.data, IF_marker)
fov21@meta.data <-cbind(fov21@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov21, boundaries = "segmentation", fov = "fov21", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov21, "./RDS_files/2_Add_metadata/Fov21_SEURAT.rds")

# 1.1.22- fov22
fov22 <- LoadNanostring(data.dir = "./fovs/fov22", fov = "fov22", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov22, "./RDS_files/1_Raw/FFov22_SEURAT.rds")
fov22 <- readRDS("./RDS_files/1_Raw/Fov22_SEURAT.rds")

# Add protein expression as another assay
fov22.ptn <- as.sparse(t(read.csv(file = "./fovs/fov22/Fov22_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov22.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov22), colnames(fov22.ptn))

ptn_assay <- CreateAssayObject(counts = fov22.ptn)
# add this assay to the previously created Seurat object
fov22[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov22)
DefaultAssay(fov22) <- "Protein"
ImageFeaturePlot(fov22, boundaries = "segmentation",fov = "fov22", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov22/Fov22_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov22@meta.data <-cbind(fov22@meta.data, IF_marker)
fov22@meta.data <-cbind(fov22@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov22, boundaries = "segmentation", fov = "fov22", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov22, "./RDS_files/2_Add_metadata/Fov22_SEURAT.rds")

# 1.1.23- fov23
fov23 <- LoadNanostring(data.dir = "./fovs/fov23", fov = "fov23", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov23, "./RDS_files/1_Raw/FFov23_SEURAT.rds")
fov23 <- readRDS("./RDS_files/1_Raw/Fov23_SEURAT.rds")

# Add protein expression as another assay
metadata <- read.csv("./fov23/Fov23_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
tx <- read.csv("./fov23/Fov23_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
setdiff(levels(metadata$cell_ID), levels(tx$cell_ID)) # cells "1551" not in the levels of tx file
tx[tx$cell_ID == "1551", ]  
colnames(fov23) # cell "1551"is not in the fov21 Seurat object

fov23.ptn <- as.sparse(t(read.csv(file = "./fovs/fov23/Fov23_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov23.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov23), colnames(fov23.ptn))

ptn_assay <- CreateAssayObject(counts = fov23.ptn)
# add this assay to the previously created Seurat object
fov23[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov23)
DefaultAssay(fov23) <- "Protein"
ImageFeaturePlot(fov23, boundaries = "segmentation",fov = "fov23", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov23/Fov23_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov23@meta.data <-cbind(fov23@meta.data, IF_marker)
fov23@meta.data <-cbind(fov23@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov23, boundaries = "segmentation", fov = "fov23", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov23, "./RDS_files/2_Add_metadata/Fov23_SEURAT.rds")

# 1.1.24- fov24
fov24 <- LoadNanostring(data.dir = "./fovs/fov24", fov = "fov24", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov24, "./RDS_files/1_Raw/Fov24_SEURAT.rds")
fov24 <- readRDS("./RDS_files/1_Raw/Fov24_SEURAT.rds")

# Add protein expression as another assay
fov24.ptn <- as.sparse(t(read.csv(file = "./fovs/fov24/Fov24_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov24.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov24), colnames(fov24.ptn))

ptn_assay <- CreateAssayObject(counts = fov24.ptn)
# add this assay to the previously created Seurat object
fov24[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov24)
DefaultAssay(fov24) <- "Protein"
ImageFeaturePlot(fov24, boundaries = "segmentation",fov = "fov24", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov24/Fov24_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov24@meta.data <-cbind(fov24@meta.data, IF_marker)
fov24@meta.data <-cbind(fov24@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov24, boundaries = "segmentation", fov = "fov24", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov24, "./RDS_files/2_Add_metadata/Fov24_SEURAT.rds")

# 1.1.25- fov25
fov25 <- LoadNanostring(data.dir = "./fovs/fov25", fov = "fov25", assay = "RNA") #when loading using "Primary_fov01", could not plot the molecule dots. Use "fov01" to make it work and load all files properly. Same wtih fov22!
saveRDS(fov25, "./RDS_files/1_Raw/Fov25_SEURAT.rds")
fov25 <- readRDS("./RDS_files/1_Raw/Fov25_SEURAT.rds")

# Add protein expression as another assay
fov25.ptn <- as.sparse(t(read.csv(file = "./fovs/fov25/Fov25_metadata_file_edited.csv", sep = ",",
                                  header = TRUE, row.names = 1))) #use transpose so cells will be the colnames and proteins the rownames
colnames(fov25.ptn)
# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(fov25), colnames(fov25.ptn))

ptn_assay <- CreateAssayObject(counts = fov25.ptn)
# add this assay to the previously created Seurat object
fov25[["Protein"]] <- ptn_assay
# Validate that the object now contains multiple assays
Assays(fov25)
DefaultAssay(fov25) <- "Protein"
ImageFeaturePlot(fov25, boundaries = "segmentation",fov = "fov25", features = c("MeanPanCK","MaxPanCK", "MeanCD45","MaxCD45", "MeanCD3", "MaxCD3"), max.cutoff = "q95")

# Add metadata regarding PanCK (tumor), CD45 (stroma) or other labels
metadata_edited <- read.csv("./fovs/fov25/Fov25_metadata_file_edited2.csv")
IF_marker <- metadata_edited$IF_marker
IF_type <- metadata_edited$IF_type
fov25@meta.data <-cbind(fov25@meta.data, IF_marker)
fov25@meta.data <-cbind(fov25@meta.data, IF_type)

IF_typecolor <- setNames(c("#00CC66", "#FF00FF","#FFFF00", "#606060"),
                         c("Epithelial", "Stromal", "CD3 Lymphocyte", "Other"))
ImageDimPlot(fov25, boundaries = "segmentation", fov = "fov25", group.by = "IF_type", cols = IF_typecolor, coord.fixed = FALSE, dark.background = TRUE)
saveRDS(fov25, "./RDS_files/2_Add_metadata/Fov25_SEURAT.rds")

#######################################################################################################################################################################################################################################################################################


# 2- Add metadata
# Based on the file "Annotations.csv"

# 2.1- Sample ID
fov01$sample_ID <- "fov01"
fov02$sample_ID <- "fov02"
fov03$sample_ID <- "fov03"
fov04$sample_ID <- "fov04"
fov05$sample_ID <- "fov05"
fov06$sample_ID <- "fov06"
fov07$sample_ID <- "fov07"
fov08$sample_ID <- "fov08"
fov09$sample_ID <- "fov09"
fov10$sample_ID <- "fov10"
fov11$sample_ID <- "fov11"
fov12$sample_ID <- "fov12"
fov13$sample_ID <- "fov13"
fov14$sample_ID <- "fov14"
fov15$sample_ID <- "fov15"
fov16$sample_ID <- "fov16"
fov17$sample_ID <- "fov17"
fov18$sample_ID <- "fov18"
fov19$sample_ID <- "fov19"
fov20$sample_ID <- "fov20"
fov21$sample_ID <- "fov21"
fov22$sample_ID <- "fov22"
fov23$sample_ID <- "fov23"
fov24$sample_ID <- "fov24"
fov25$sample_ID <- "fov25"

# 2.2- CMS (R CMS Classifier)
fov01$CMS <- "Non classified"
fov02$CMS <- "CMS3"
fov03$CMS <- "Non classified"
fov04$CMS <- "CMS4"
fov05$CMS <- "Non classified"
fov06$CMS <- "Non classified"
fov07$CMS <- "Non classified"
fov08$CMS <- "CMS2"
fov09$CMS <- "CMS4"
fov10$CMS <- "CMS1"
fov11$CMS <- "CMS4"
fov12$CMS <- "Non classified"
fov13$CMS <- "Non classified"
fov14$CMS <- "Non classified"
fov15$CMS <- "Non classified"
fov16$CMS <- "Non classified"
fov17$CMS <- "Non classified"
fov18$CMS <- "Non classified"
fov19$CMS <- "Non classified"
fov20$CMS <- "CMS4"
fov21$CMS <- "CMS1"
fov22$CMS <- "CMS2"
fov23$CMS <- "Non classified"
fov24$CMS <- "Non classified"
fov25$CMS <- "CMS2"

# 2.3- Tumour Stroma Percentage
fov01$Tumor_Stroma_Percent <- "0"
fov02$Tumor_Stroma_Percent <- "1"
fov03$Tumor_Stroma_Percent <- "0"
fov04$Tumor_Stroma_Percent <- "1"
fov05$Tumor_Stroma_Percent <- "0"
fov06$Tumor_Stroma_Percent <- "0"
fov07$Tumor_Stroma_Percent <- "1"
fov08$Tumor_Stroma_Percent <- "0"
fov09$Tumor_Stroma_Percent <- "1"
fov10$Tumor_Stroma_Percent <- "0"
fov11$Tumor_Stroma_Percent <- "0"
fov12$Tumor_Stroma_Percent <- "0"
fov13$Tumor_Stroma_Percent <- "0"
fov14$Tumor_Stroma_Percent <- "0"
fov15$Tumor_Stroma_Percent <- "0"
fov16$Tumor_Stroma_Percent <- "0"
fov17$Tumor_Stroma_Percent <- "0"
fov18$Tumor_Stroma_Percent <- "0"
fov19$Tumor_Stroma_Percent <- "0"
fov20$Tumor_Stroma_Percent <- "0"
fov21$Tumor_Stroma_Percent <- "0"
fov22$Tumor_Stroma_Percent <- "1"
fov23$Tumor_Stroma_Percent <- "0"
fov24$Tumor_Stroma_Percent <- "0"
fov25$Tumor_Stroma_Percent <- "0"
  
# 2.4- KM Score (immune)
fov01$KM_Score <- "1"
fov02$KM_Score <- "0"
fov03$KM_Score <- "0"
fov04$KM_Score <- "0"
fov05$KM_Score <- "0"
fov06$KM_Score <- "0"
fov07$KM_Score <- "0"
fov08$KM_Score <- "0"
fov09$KM_Score <- "0"
fov10$KM_Score <- "0"
fov11$KM_Score <- "0"
fov12$KM_Score <- "0"
fov13$KM_Score <- "0"
fov14$KM_Score <- "0"
fov15$KM_Score <- "0"
fov16$KM_Score <- "1"
fov17$KM_Score <- "1"
fov18$KM_Score <- "1"
fov19$KM_Score <- "0"
fov20$KM_Score <- "1"
fov21$KM_Score <- "0"
fov22$KM_Score <- "0"
fov23$KM_Score <- "1"
fov24$KM_Score <- "1"
fov25$KM_Score <- "1"

# 2.5- mGPS (systemic inflammation)
fov01$mGPS <- "0"
fov02$mGPS <- "0"
fov03$mGPS <- "0"
fov04$mGPS <- "0"
fov05$mGPS <- "0"
fov06$mGPS <- "0"
fov07$mGPS <- "0"
fov08$mGPS <- "1"
fov09$mGPS <- "0"
fov10$mGPS <- "0"
fov11$mGPS <- "0"
fov12$mGPS <- "0"
fov13$mGPS <- "0"
fov14$mGPS <- "1"
fov15$mGPS <- "0"
fov16$mGPS <- "2"
fov17$mGPS <- "0"
fov18$mGPS <- "0"
fov19$mGPS <- "0"
fov20$mGPS <- "0"
fov21$mGPS <- "0"
fov22$mGPS <- "0"
fov23$mGPS <- "0"
fov24$mGPS <- "1"
fov25$mGPS <- "0"

# 2.6- MMR_Status_Jen (1 dmmr, 2, pmmr, 3 partial)
fov01$MMR_Status_Jen <- "PMMR"
fov02$MMR_Status_Jen <- "DMMR"
fov03$MMR_Status_Jen <- "DMMR"
fov04$MMR_Status_Jen <- "PMMR"
fov05$MMR_Status_Jen <- "PMMR"
fov06$MMR_Status_Jen <- "PMMR"
fov07$MMR_Status_Jen <- "PMMR"
fov08$MMR_Status_Jen <- "PMMR"
fov09$MMR_Status_Jen <- "Partial"
fov10$MMR_Status_Jen <- "DMMR"
fov11$MMR_Status_Jen <- "PMMR"
fov12$MMR_Status_Jen <- "PMMR"
fov13$MMR_Status_Jen <- "PMMR"
fov14$MMR_Status_Jen <- "Partial"
fov15$MMR_Status_Jen <- "PMMR"
fov16$MMR_Status_Jen <- "DMMR"
fov17$MMR_Status_Jen <- "PMMR"
fov18$MMR_Status_Jen <- "PMMR"
fov19$MMR_Status_Jen <- "PMMR"
fov20$MMR_Status_Jen <- "PMMR"
fov21$MMR_Status_Jen <- "DMMR"
fov22$MMR_Status_Jen <- "PMMR"
fov23$MMR_Status_Jen <- "PMMR"
fov24$MMR_Status_Jen <- "DMMR"
fov25$MMR_Status_Jen <- "PMMR"

# 2.7- Recurrence location
fov01$Recurrence_location <- "NA"
fov02$Recurrence_location <- "Liver"
fov03$Recurrence_location <- "NA"
fov04$Recurrence_location <- "NA"
fov05$Recurrence_location <- "NA"
fov06$Recurrence_location <- "soft tissue mass left rectus lesion anterior to left EIAs"
fov07$Recurrence_location <- "NA"
fov08$Recurrence_location <- "Brain"
fov09$Recurrence_location <- "NA"
fov10$Recurrence_location <- "NA"
fov11$Recurrence_location <- "NA"
fov12$Recurrence_location <- "NA"
fov13$Recurrence_location <- "Liver lung"
fov14$Recurrence_location <- "NA"
fov15$Recurrence_location <- "NA"
fov16$Recurrence_location <- "NA"
fov17$Recurrence_location <- "NA"
fov18$Recurrence_location <- "liver spleen skin abdo wound intraperitoneal"
fov19$Recurrence_location <- "NA"
fov20$Recurrence_location <- "NA"
fov21$Recurrence_location <- "NA"
fov22$Recurrence_location <- "NA"
fov23$Recurrence_location <- "NA"
fov24$Recurrence_location <- "NA"
fov25$Recurrence_location <- "NA"

# 2.8- Recurrence code
fov01$Recurrence_code <- "0"
fov02$Recurrence_code <- "3"
fov03$Recurrence_code <- "0"
fov04$Recurrence_code <- "0"
fov05$Recurrence_code <- "0"
fov06$Recurrence_code <- "1"
fov07$Recurrence_code <- "0"
fov08$Recurrence_code <- "4"
fov09$Recurrence_code <- "0"
fov10$Recurrence_code <- "0"
fov11$Recurrence_code <- "0"
fov12$Recurrence_code <- "0"
fov13$Recurrence_code <- "5"
fov14$Recurrence_code <- "0"
fov15$Recurrence_code <- "0"
fov16$Recurrence_code <- "0"
fov17$Recurrence_code <- "0"
fov18$Recurrence_code <- "5"
fov19$Recurrence_code <- "0"
fov20$Recurrence_code <- "0"
fov21$Recurrence_code <- "0"
fov22$Recurrence_code <- "0"
fov23$Recurrence_code <- "0"
fov24$Recurrence_code <- "0"
fov25$Recurrence_code <- "0"

# 2.1- The [[ operator can add columns to object metadata. This is a great place to stash QC stats - calculate and add the percentage of MT genes
fov01[["percent.mt"]] <- PercentageFeatureSet(fov01, pattern = "^MT-") 
fov02[["percent.mt"]] <- PercentageFeatureSet(fov02, pattern = "^MT-") 
fov03[["percent.mt"]] <- PercentageFeatureSet(fov03, pattern = "^MT-") 
fov04[["percent.mt"]] <- PercentageFeatureSet(fov04, pattern = "^MT-") 
fov05[["percent.mt"]] <- PercentageFeatureSet(fov05, pattern = "^MT-") 
fov06[["percent.mt"]] <- PercentageFeatureSet(fov06, pattern = "^MT-") 
fov07[["percent.mt"]] <- PercentageFeatureSet(fov07, pattern = "^MT-") 
fov08[["percent.mt"]] <- PercentageFeatureSet(fov08, pattern = "^MT-") 
fov09[["percent.mt"]] <- PercentageFeatureSet(fov09, pattern = "^MT-") 
fov10[["percent.mt"]] <- PercentageFeatureSet(fov10, pattern = "^MT-") 
fov11[["percent.mt"]] <- PercentageFeatureSet(fov11, pattern = "^MT-") 
fov12[["percent.mt"]] <- PercentageFeatureSet(fov12, pattern = "^MT-") 
fov13[["percent.mt"]] <- PercentageFeatureSet(fov13, pattern = "^MT-") 
fov14[["percent.mt"]] <- PercentageFeatureSet(fov14, pattern = "^MT-") 
fov15[["percent.mt"]] <- PercentageFeatureSet(fov15, pattern = "^MT-") 
fov16[["percent.mt"]] <- PercentageFeatureSet(fov16, pattern = "^MT-") 
fov17[["percent.mt"]] <- PercentageFeatureSet(fov17, pattern = "^MT-") 
fov18[["percent.mt"]] <- PercentageFeatureSet(fov18, pattern = "^MT-") 
fov19[["percent.mt"]] <- PercentageFeatureSet(fov19, pattern = "^MT-") 
fov20[["percent.mt"]] <- PercentageFeatureSet(fov20, pattern = "^MT-") 
fov21[["percent.mt"]] <- PercentageFeatureSet(fov21, pattern = "^MT-") 
fov22[["percent.mt"]] <- PercentageFeatureSet(fov22, pattern = "^MT-") 
fov23[["percent.mt"]] <- PercentageFeatureSet(fov23, pattern = "^MT-") 
fov24[["percent.mt"]] <- PercentageFeatureSet(fov24, pattern = "^MT-") 
fov25[["percent.mt"]] <- PercentageFeatureSet(fov25, pattern = "^MT-") 

# 2.2- Calculation of the percentage and absolute numbers of human transcripts/UMIs
fov01[["RNA_counts"]] <- fov01$nCount_RNA
fov02[["RNA_counts"]] <- fov02$nCount_RNA
fov03[["RNA_counts"]] <- fov03$nCount_RNA
fov04[["RNA_counts"]] <- fov04$nCount_RNA
fov05[["RNA_counts"]] <- fov05$nCount_RNA
fov06[["RNA_counts"]] <- fov06$nCount_RNA
fov07[["RNA_counts"]] <- fov07$nCount_RNA
fov08[["RNA_counts"]] <- fov08$nCount_RNA
fov09[["RNA_counts"]] <- fov09$nCount_RNA
fov10[["RNA_counts"]] <- fov10$nCount_RNA
fov11[["RNA_counts"]] <- fov11$nCount_RNA
fov12[["RNA_counts"]] <- fov12$nCount_RNA
fov13[["RNA_counts"]] <- fov13$nCount_RNA
fov14[["RNA_counts"]] <- fov14$nCount_RNA
fov15[["RNA_counts"]] <- fov15$nCount_RNA
fov16[["RNA_counts"]] <- fov16$nCount_RNA
fov17[["RNA_counts"]] <- fov17$nCount_RNA
fov18[["RNA_counts"]] <- fov18$nCount_RNA
fov19[["RNA_counts"]] <- fov19$nCount_RNA
fov20[["RNA_counts"]] <- fov20$nCount_RNA
fov21[["RNA_counts"]] <- fov21$nCount_RNA
fov22[["RNA_counts"]] <- fov22$nCount_RNA
fov23[["RNA_counts"]] <- fov23$nCount_RNA
fov24[["RNA_counts"]] <- fov24$nCount_RNA
fov25[["RNA_counts"]] <- fov25$nCount_RNA

# 2.3- Calculation of the number of  genes detected per cell
fov01[["Genes"]] <- fov01$nFeature_RNA
fov02[["Genes"]] <- fov02$nFeature_RNA
fov03[["Genes"]] <- fov03$nFeature_RNA
fov04[["Genes"]] <- fov04$nFeature_RNA
fov05[["Genes"]] <- fov05$nFeature_RNA
fov06[["Genes"]] <- fov06$nFeature_RNA
fov07[["Genes"]] <- fov07$nFeature_RNA
fov08[["Genes"]] <- fov08$nFeature_RNA
fov09[["Genes"]] <- fov09$nFeature_RNA
fov10[["Genes"]] <- fov10$nFeature_RNA
fov11[["Genes"]] <- fov11$nFeature_RNA
fov12[["Genes"]] <- fov12$nFeature_RNA
fov13[["Genes"]] <- fov13$nFeature_RNA
fov14[["Genes"]] <- fov14$nFeature_RNA
fov15[["Genes"]] <- fov15$nFeature_RNA
fov16[["Genes"]] <- fov16$nFeature_RNA
fov17[["Genes"]] <- fov17$nFeature_RNA
fov18[["Genes"]] <- fov18$nFeature_RNA
fov19[["Genes"]] <- fov19$nFeature_RNA
fov20[["Genes"]] <- fov20$nFeature_RNA
fov21[["Genes"]] <- fov21$nFeature_RNA
fov22[["Genes"]] <- fov22$nFeature_RNA
fov23[["Genes"]] <- fov23$nFeature_RNA
fov24[["Genes"]] <- fov24$nFeature_RNA
fov25[["Genes"]] <- fov25$nFeature_RNA

saveRDS(fov01, "./RDS_files/2_Add_metadata/Fov01_SEURAT.rds")
saveRDS(fov02, "./RDS_files/2_Add_metadata/Fov02_SEURAT.rds")
saveRDS(fov03, "./RDS_files/2_Add_metadata/Fov03_SEURAT.rds")
saveRDS(fov04, "./RDS_files/2_Add_metadata/Fov04_SEURAT.rds")
saveRDS(fov05, "./RDS_files/2_Add_metadata/Fov05_SEURAT.rds")
saveRDS(fov06, "./RDS_files/2_Add_metadata/Fov06_SEURAT.rds")
saveRDS(fov07, "./RDS_files/2_Add_metadata/Fov07_SEURAT.rds")
saveRDS(fov08, "./RDS_files/2_Add_metadata/Fov08_SEURAT.rds")
saveRDS(fov09, "./RDS_files/2_Add_metadata/Fov09_SEURAT.rds")
saveRDS(fov10, "./RDS_files/2_Add_metadata/Fov10_SEURAT.rds")
saveRDS(fov11, "./RDS_files/2_Add_metadata/Fov11_SEURAT.rds")
saveRDS(fov12, "./RDS_files/2_Add_metadata/Fov12_SEURAT.rds")
saveRDS(fov13, "./RDS_files/2_Add_metadata/Fov13_SEURAT.rds")
saveRDS(fov14, "./RDS_files/2_Add_metadata/Fov14_SEURAT.rds")
saveRDS(fov15, "./RDS_files/2_Add_metadata/Fov15_SEURAT.rds")
saveRDS(fov16, "./RDS_files/2_Add_metadata/Fov16_SEURAT.rds")
saveRDS(fov17, "./RDS_files/2_Add_metadata/Fov17_SEURAT.rds")
saveRDS(fov18, "./RDS_files/2_Add_metadata/Fov18_SEURAT.rds")
saveRDS(fov19, "./RDS_files/2_Add_metadata/Fov19_SEURAT.rds")
saveRDS(fov20, "./RDS_files/2_Add_metadata/Fov20_SEURAT.rds")
saveRDS(fov21, "./RDS_files/2_Add_metadata/Fov21_SEURAT.rds")
saveRDS(fov22, "./RDS_files/2_Add_metadata/Fov22_SEURAT.rds")
saveRDS(fov23, "./RDS_files/2_Add_metadata/Fov23_SEURAT.rds")
saveRDS(fov24, "./RDS_files/2_Add_metadata/Fov24_SEURAT.rds")
saveRDS(fov25, "./RDS_files/2_Add_metadata/Fov25_SEURAT.rds")

# 3- Create one merged Seurat object before pre-processing and QC as we will use the same parameters across all samples

# Reference 1: https://satijalab.org/seurat/articles/merge_vignette.html
# Reference 2: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

# Add dataset labels as cell.ids just in case you have overlapping barcodes between the datasets. 
# Merge datasets into one single seurat object
fov01 <- readRDS("./RDS_files/2_Add_metadata/Fov01_SEURAT.rds")
fov02 <- readRDS("./RDS_files/2_Add_metadata/Fov02_SEURAT.rds")
fov03 <- readRDS("./RDS_files/2_Add_metadata/Fov03_SEURAT.rds")
fov04 <- readRDS("./RDS_files/2_Add_metadata/Fov04_SEURAT.rds")
fov05 <- readRDS("./RDS_files/2_Add_metadata/Fov05_SEURAT.rds")
fov06 <- readRDS("./RDS_files/2_Add_metadata/Fov06_SEURAT.rds")
fov07 <- readRDS("./RDS_files/2_Add_metadata/Fov07_SEURAT.rds")
fov08 <- readRDS("./RDS_files/2_Add_metadata/Fov08_SEURAT.rds")
fov09 <- readRDS("./RDS_files/2_Add_metadata/Fov09_SEURAT.rds")
fov10 <- readRDS("./RDS_files/2_Add_metadata/Fov10_SEURAT.rds")
fov11 <- readRDS("./RDS_files/2_Add_metadata/Fov11_SEURAT.rds")
fov12 <- readRDS("./RDS_files/2_Add_metadata/Fov12_SEURAT.rds")
fov13 <- readRDS("./RDS_files/2_Add_metadata/Fov13_SEURAT.rds")
fov14 <- readRDS("./RDS_files/2_Add_metadata/Fov14_SEURAT.rds")
fov15 <- readRDS("./RDS_files/2_Add_metadata/Fov15_SEURAT.rds")
fov16 <- readRDS("./RDS_files/2_Add_metadata/Fov16_SEURAT.rds")
fov17 <- readRDS("./RDS_files/2_Add_metadata/Fov17_SEURAT.rds")
fov18 <- readRDS("./RDS_files/2_Add_metadata/Fov18_SEURAT.rds")
fov19 <- readRDS("./RDS_files/2_Add_metadata/Fov19_SEURAT.rds")
fov20 <- readRDS("./RDS_files/2_Add_metadata/Fov20_SEURAT.rds")
fov21 <- readRDS("./RDS_files/2_Add_metadata/Fov21_SEURAT.rds")
fov22 <- readRDS("./RDS_files/2_Add_metadata/Fov22_SEURAT.rds")
fov23 <- readRDS("./RDS_files/2_Add_metadata/Fov23_SEURAT.rds")
fov24 <- readRDS("./RDS_files/2_Add_metadata/Fov24_SEURAT.rds")
fov25 <- readRDS("./RDS_files/2_Add_metadata/Fov25_SEURAT.rds")

DefaultAssay(fov01) <- "RNA"
DefaultAssay(fov02) <- "RNA"
DefaultAssay(fov03) <- "RNA"
DefaultAssay(fov04) <- "RNA"
DefaultAssay(fov05) <- "RNA"
DefaultAssay(fov06) <- "RNA"
DefaultAssay(fov07) <- "RNA"
DefaultAssay(fov08) <- "RNA"
DefaultAssay(fov09) <- "RNA"
DefaultAssay(fov10) <- "RNA"
DefaultAssay(fov11) <- "RNA"
DefaultAssay(fov12) <- "RNA"
DefaultAssay(fov13) <- "RNA"
DefaultAssay(fov14) <- "RNA"
DefaultAssay(fov15) <- "RNA"
DefaultAssay(fov16) <- "RNA"
DefaultAssay(fov17) <- "RNA"
DefaultAssay(fov18) <- "RNA"
DefaultAssay(fov19) <- "RNA"
DefaultAssay(fov20) <- "RNA"
DefaultAssay(fov21) <- "RNA"
DefaultAssay(fov22) <- "RNA"
DefaultAssay(fov23) <- "RNA"
DefaultAssay(fov24) <- "RNA"
DefaultAssay(fov25) <- "RNA"

fov_merged <- merge(fov01, c(fov02, fov03, fov04, fov05, fov06, fov07, fov08, fov09, fov10,
                             fov11, fov12, fov13, fov14, fov15, fov16, fov17, fov18, fov19,
                             fov20, fov21, fov22, fov23, fov24, fov25), 
                    add.cell.ids = c("fov01", "fov02", "fov03", "fov04", "fov05", "fov06", "fov07", "fov08", "fov09", "fov10",
                                     "fov11", "fov12", "fov13", "fov14", "fov15", "fov16", "fov17", "fov18", "fov19", "fov20",
                                     "fov21", "fov22", "fov23", "fov24", "fov25"), project = "Colonrectum_Cancer_COSMX")
saveRDS(fov_merged, "fov_merged.rds")













