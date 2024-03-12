

#1- Converting anndata to SingleCellExperiment object (SCE)
library(reticulate)
use_condaenv("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate")
use_python("/Users/joaoluizsfilho/Library/r-miniconda/envs/r-reticulate/bin/python")

library(rhdf5)
library(zellkonverter)
ad_ep1 <- readH5AD('./2_h5ad_files/ad_ep1.h5ad') 
ad_ep2 <- readH5AD('./2_h5ad_files/ad_ep2.h5ad')
ad_ep3 <- readH5AD('./2_h5ad_files/ad_ep3.h5ad')

ad_tme1 <- readH5AD('./2_h5ad_files/ad_tme1.h5ad') 
ad_tme2 <- readH5AD('./2_h5ad_files/ad_tme2.h5ad')
ad_tme3 <- readH5AD('./2_h5ad_files/ad_tme3.h5ad')

ad_sma1 <- readH5AD('./2_h5ad_files/ad_sma1.h5ad') 
ad_sma2 <- readH5AD('./2_h5ad_files/ad_sma2.h5ad')
ad_sma3 <- readH5AD('./2_h5ad_files/ad_sma3.h5ad')

ad_cms1 <- readH5AD('./2_h5ad_files/ad_cms1.h5ad') 
ad_cms2 <- readH5AD('./2_h5ad_files/ad_cms2.h5ad')
ad_cms3 <- readH5AD('./2_h5ad_files/ad_cms3.h5ad')
ad_cms4 <- readH5AD('./2_h5ad_files/ad_cms4.h5ad')
ad_cms5 <- readH5AD('./2_h5ad_files/ad_cms5.h5ad')

ad_dmmr <- readH5AD('./2_h5ad_files/ad_dmmr.h5ad')
ad_pmmr <- readH5AD('./2_h5ad_files/ad_pmmr.h5ad')
ad_partial <- readH5AD('./2_h5ad_files/ad_partial.h5ad')

ad_norecurrence <- readH5AD('./2_h5ad_files/ad_norecurrence_giotto.h5ad') 
ad_local <- readH5AD('./2_h5ad_files/ad_local_giotto.h5ad')
ad_liver <- readH5AD('./2_h5ad_files/ad_liver_giotto.h5ad')
ad_brain <- readH5AD('./2_h5ad_files/ad_brain_giotto.h5ad')
ad_multisite <- readH5AD('./2_h5ad_files/ad_multisite_giotto.h5ad')

#Save object
saveRDS(ad_ep1, "./1_RDS_files/ad_ep1_sce.rds")
saveRDS(ad_ep2, "./1_RDS_files/ad_ep2_sce.rds")
saveRDS(ad_ep3, "./1_RDS_files/ad_ep3_sce.rds")

saveRDS(ad_tme1, "./1_RDS_files/ad_tme1_sce.rds")
saveRDS(ad_tme2, "./1_RDS_files/ad_tme2_sce.rds")
saveRDS(ad_tme3, "./1_RDS_files/ad_tme3_sce.rds")

saveRDS(ad_sma1, "./1_RDS_files/ad_sma1_sce.rds")
saveRDS(ad_sma2, "./1_RDS_files/ad_sma2_sce.rds")
saveRDS(ad_sma3, "./1_RDS_files/ad_sma3_sce.rds")

saveRDS(ad_cms1, "./1_RDS_files/ad_cms1_sce.rds")
saveRDS(ad_cms2, "./1_RDS_files/ad_cms2_sce.rds")
saveRDS(ad_cms3, "./1_RDS_files/ad_cms3_sce.rds")
saveRDS(ad_cms4, "./1_RDS_files/ad_cms4_sce.rds")
saveRDS(ad_cms5, "./1_RDS_files/ad_cms5_sce.rds")

saveRDS(ad_dmmr, "./1_RDS_files/ad_dmmr_sce.rds")
saveRDS(ad_pmmr, "./1_RDS_files/ad_pmmr_sce.rds")
saveRDS(ad_partial, "./1_RDS_files/ad_[sma3]partial_sce.rds")

saveRDS(ad_norecurrence, "./1_RDS_files/ad_norecurrence_sce.rds")
saveRDS(ad_local, "./1_RDS_files/ad_local_sce.rds")
saveRDS(ad_liver, "./1_RDS_files/ad_liver_sce.rds")
saveRDS(ad_brain, "./1_RDS_files/ad_brain_sce.rds")
saveRDS(ad_multisite, "./1_RDS_files/ad_multisite_sce.rds")

ad <- readRDS("adata_fov_integrated_sce.rds")

#Trying to convert SCE to SpatialExperiment object (SPE)
library(spatialLIBD)
spe <- sce_to_spe(ad, imageData = NULL) # it did not work because it needs the imgData dataframe

#Trying to convert the SCE to SPE as suggested in https://github.com/theislab/zellkonverter/issues/61
library(SpatialExperiment) #explore more this function
coords <- as.matrix(reducedDim(ad, "spatial"))
colnames(coords) = c("CenterX_global_px","CenterY_global_px")
spe2 <- SpatialExperiment(
  assay = list(counts = assay(ad)), 
  colData = ad@colData, 
  spatialCoords = coords
)

saveRDS(spe2, "adata_fov_integrated_spe.rds")
spe <- readRDS("adata_fov_integrated_spe.rds")


#3- Analysis using imcRtools functions
#The Bodernmiller's pipeline works with both, SPE and SCE, just need to figure out the correct slots of SCE that should be used in the functions.
#Have a look at the help of each function to check which slot of the SCE should be used.

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



