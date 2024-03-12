

library(Giotto)

# Package Accessibility
default_instrs <- createGiottoInstructions()

# 1- Setup

# Custom color palettes from rcartocolor
# pal10 = rcartocolor::carto_pal(n = 10, name = 'Pastel')
pal10 = c("#66C5CC","#F6CF71","#F89C74","#DCB0F2","#FF0000",   #87C55F
          "#9EB9F3","#FE88B1","#C9DB74","#8BE0A4","#B3B3B3", "black")
# viv10 = rcartocolor::carto_pal(n = 10, name = 'Vivid')
viv10 = c("#E58606","#5D69B1","#52BCA3","#99C945","#CC61B0",
          "#24796C","#DAA51B","#2F8AC4","#764E9F","#A5AA99")

# set working directory
results_folder = "./Plots"

# Optional: Specify a path to a Python executable within a conda or miniconda environment. 
# If set to NULL (default), the Python executable within the previously installed Giotto environment will be used.
my_python_path = NULL # alternatively, "/local/python/path/python" if desired.

## Set object behavior
# by directly saving plots, but not rendering them you will save a lot of time
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  return_plot = FALSE,
                                  python_path = my_python_path)


# [Expected Directory] This function generates a giotto object when given a link to a cosmx output directory. 
# It expects the following items within the directory where the bolded portions are what this function matches against:
# CellComposite (folder of images)
# CellLabels (folder of images)
# CellOverlay (folder of images)
# CompartmentLabels (folder of images)
# experimentname_exprMat_file.csv (file)
# experimentname_fov_positions_file.csv (file)
# experimentname_metadata_file.csv (file)
# experimentname_tx_file.csv (file)

# ‘all’ - loads and requires subcellular information from tx_file and fov_positions_file and also the existing aggregated information (expression, spatial locations, and metadata) from exprMat_file and metadata_file.
# ‘subcellular’ - loads and requires subcellular information from tx_file and fov_positions_file only.
# ‘aggregate’ - loads and requires the existing aggregate information (expression, spatial locations, and metadata) from exprMat_file and metadata_file.

# [Images] Images in the default CellComposite, CellLabels, CompartmentLabels, and CellOverlay folders will be loaded as giotto largeImage objects in all workflows as long as they are available. 
# Additionally, CellComposite images will be converted to giotto image objects, making plotting with these image objects more responsive when accessing them from a server. 
# showGiottoImageNames can be used to see the available images.

showGiottoFeatInfo(fov_join)
showGiottoSpatialInfo(fov_join)


# 2- Create Giotto object - option 1
## provide path to nanostring folder
data_path = './Files/'

## create giotto cosmx object
fov_join = createGiottoCosMxObject(cosmx_dir = data_path,
                                   data_to_use = 'subcellular', # only subcellular
                                   FOVs = c(01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25), #looks like it can load several fovs at once
                                   instructions = instrs)

# 3- Saving and Loading the giotto object
saveGiotto(gobject = fov_join,
           dir = './COSMX_Giotto', overwrite = TRUE)

fov_join <- loadGiotto('./COSMX_Giotto/')

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 3- Create Giotto object - option 2

# 3.1- Subcellular detections (points info)

# tx_file.csv contains the subcellular detections information. It contains information on each of the individual feature detections within the sample.
# fov which FOV the detection happened in
# cell_ID the ID of the cell the detection happened in
# x_global_px the global spatial x location in pixels
# y_global_px the global spatial y location in pixels
# x_local_px the spatial x location in pixels within the FOV
# y_local_px the spatial y location in pixels within the FOV
# z the z plane the detection was called in (-1 to 16)
# target the feature the probe is targeted against
# CellComp Cellular compartment the detection happened in (0, Cytoplasm, Membrane, Nuclear)

## provide path to nanostring folder
data_path = './Files/'

# load transcript coordinates
library(data.table)
tx_coord_all = fread(paste0(data_path, 'S1_TMA_tx_file.csv'))

colnames(tx_coord_all)
cat('\n')
# z planes
tx_coord_all[, table(z)]
cat('\n')
# Cell compartment
tx_coord_all[, table(CellComp)]

# 3.2- Split detections by features vs negative probes

# tx_file.csv contains information on both actual features (960 targeted gene probes in this dataset) and negative probes (20) 
# that are targeted to alien sequences defined by the External RNA Controls Consortium (ERCC) that do not exist in human tissue.
# These two types of detections will be treated as separate feature types (feat_type) and placed in separate expression matrices.

all_IDs = tx_coord_all[, unique(target)]
# negative probe IDs
neg_IDs = all_IDs[grepl(pattern = 'NegPrb', all_IDs)]
cat('Negative Probe IDs\n')
neg_IDs
cat('\nFeature IDs\n')
feat_IDs = all_IDs[!all_IDs %in% neg_IDs]
length(feat_IDs)

# split detections
feat_coords_all = tx_coord_all[target %in% feat_IDs]
neg_coords_all = tx_coord_all[target %in% neg_IDs]

cat('\nFeatures: ', feat_coords_all[, .N], '\n',
    'NegProbes: ', neg_coords_all[, .N])

feat_IDs

# 3.2.1- Preview negative probes (optional)

# Previewing the probe information can be done by converting to giottoPoints and then using plot().
# Here we show a preview of the negative probes.
# Note: if previewing the rna expression information, it is highly recommended to set a subset of features using the feats param. 
# The default is to plot all points, which can be very slow for large data.

neg_points = createGiottoPoints(
  x = neg_coords_all[, .(target, x_global_px, y_global_px)]
)
plot(neg_points, point_size = 0.5, feats = neg_IDs) # plot weird if using fovs from different samples!


# 3.3 FOV shifts - IF ANALYSING MORE THAN 1 FOV AT A TIME

# fov_positions_file.csv contains information on the x and y shifts needed in order to put the FOVs tiles together into a cohesive whole. 
# This information is needed during the gobject join process.

#  load field of vision (fov) positions
fov_offset_file = fread(paste0(data_path, 'S1_TMA_fov_positions_file.csv'))

# 2.4 Choose field of view for analysis - IF ANALYSING MORE THAN 1 FOV AT A TIME

# CosMx data is large and Giotto loads in the subcellular information by FOV. 
# This dataset includes 28 FOVs which can be difficult for most computers to handle at once.
# This tutorial will use FOVs ‘02’, ‘03’, and ‘04’ which correspond to the 3 FOVs visible on the bottom right in the negative probe preview above.

gobjects_list = list()

id_set = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25')

# 3.4- Create a Giotto Object for each FOV - IF ANALYSING MORE THAN 1 FOV AT A TIME

for(fov_i in 1:length(id_set)) {
  
  fov_id = id_set[fov_i]
  
  # 1. original composite image as png
  original_composite_image = paste0(data_path, 'CellComposite/CellComposite_F0', fov_id,'.jpg')
  
  # 2. input cell segmentation as mask file
  segmentation_mask = paste0(data_path, 'CellLabels/CellLabels_F0', fov_id, '.tif')
  
  # 3. input features coordinates + offset
  tx_coord = tx_coord_all[fov == as.numeric(fov_id)]
  tx_coord = tx_coord[,.(x_local_px, y_local_px, z, target)]
  colnames(tx_coord) = c('x', 'y', 'z', 'gene_id')
  tx_coord = tx_coord[,.(x, y, gene_id)]
  
  fovsubset = createGiottoObjectSubcellular(
    gpoints = list('rna' = tx_coord),
    gpolygons = list('cell' = segmentation_mask),
    polygon_mask_list_params = list(
      mask_method = 'guess',
      flip_vertical = TRUE,
      flip_horizontal = FALSE,
      shift_horizontal_step = FALSE
    ),
    instructions = instrs
  )
  
  
  # cell centroids are now used to provide the spatial locations
  fovsubset = addSpatialCentroidLocations(fovsubset,
                                          poly_info = 'cell')
  
  # create and add Giotto images
  composite = createGiottoLargeImage(raster_object = original_composite_image,
                                     negative_y = FALSE,
                                     name = 'composite')
  
  fovsubset = addGiottoImage(gobject = fovsubset,
                             largeImages = list(composite))
  
  
  fovsubset = convertGiottoLargeImageToMG(giottoLargeImage = composite,
                                          #mg_name = 'composite',
                                          gobject = fovsubset,
                                          return_gobject = TRUE)
  
  gobjects_list[[fov_i]] = fovsubset
  
  
}

# 3.5- Join FOV Giotto Objects - IF ANALYSING MORE THAN 1 FOV AT A TIME
new_names = paste0("fov0", id_set)

id_match = match(as.numeric(id_set), fov_offset_file$fov)
x_shifts = fov_offset_file[id_match]$x_global_px
y_shifts = fov_offset_file[id_match]$y_global_px

# Create Giotto object that includes all selected FOVs
fov_join = joinGiottoObjects(gobject_list = gobjects_list,
                             gobject_names = new_names,
                             join_method = 'shift',
                             x_shift = x_shifts,
                             y_shift = y_shifts)

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 4- Aggregate subcellular features

# Giotto supports working directly with the subcellular features in order to generate cell by feature matrices. 
# The data generated this way is then given the spatial unit 'cell'. 
# This workflow is recommended over loading the provided cell by feature (aggregated expression) matrix
# and then including the subcellular information as secondary data.
# When both the raw subcellular information and the pre-made expression matrix are loaded in at the same time, 
# the subcellular data and all data generated from it should be given the spatial unit 'cell' 
# and the pre-generated aggregated information should be given a different spatial unit such as 'cell_agg' to differentiate between the two sources of information.

# In this step, we will be aggregating the feature points of 'rna' and 'neg_probe' into the 'cell' spatial unit.
# Find the feature points overlapped by polygons. This overlap information is then
# returned to the relevant giottoPolygon object's overlaps slot.
fov_join = calculateOverlapRaster(fov_join, feat_info = 'rna')
fov_join = calculateOverlapRaster(fov_join, feat_info = 'neg_probe')

# Convert the overlap information into a cell by feature expression matrix which
# is then stored in the Giotto object's expression slot
fov_join = overlapToMatrix(fov_join, feat_info = 'rna')
fov_join = overlapToMatrix(fov_join, feat_info = 'neg_probe')

p <- fov_join@cell_metadata[["cell"]][["rna"]]@metaDT

write.csv(p, "giotto_cell_metadata.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

# 5- Creating Giotto object from Anndata

# Reference: https://giottosuite.readthedocs.io/en/latest/subsections/datasets/interoperability_04122023.html
# To convert an AnnData Object back into a Giotto object, it must first be saved as a .h5ad file. 
# The name of said file may then be provided to anndataToGiotto() for conversion.
# If a nearest neighbor network or spatial netowkr was created using the key_added argument, they may be provided to arguments n_key_added and/or spatial_n_key_added, respectively.


# 5.1- adata_fov_integrated
fov_integrated_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/adata_fov_integrated_giotto.h5ad",
                                 python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_integrated_gobject,
           dir = './COSMX_Giotto/fov_integrated', overwrite = TRUE)

fov_integrated <- loadGiotto('./COSMX_Giotto/fov_integrated/saveGiottoDir')

# 5.2- No recurrence
fov_no_recurrence_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_norecurrence_giotto.h5ad",
                                          python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_no_recurrence_gobject,
           dir = './COSMX_Giotto/fov_no_recurrence', overwrite = TRUE)

fov_no_recurrence <- loadGiotto('./COSMX_Giotto/fov_no_recurrence/saveGiottoDir')

# 5.3- Local recurrence
fov_local_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_local_giotto.h5ad",
                                          python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_local_gobject,
           dir = './COSMX_Giotto/fov_local', overwrite = TRUE)

fov_local <- loadGiotto('./COSMX_Giotto/fov_local/saveGiottoDir')

# 5.4- Liver recurrence
fov_liver_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_liver_giotto.h5ad",
                                          python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_liver_gobject,
           dir = './COSMX_Giotto/fov_liver', overwrite = TRUE)

fov_liver <- loadGiotto('./COSMX_Giotto/fov_liver/saveGiottoDir')

# 5.5- Brain recurrence
fov_brain_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_brain_giotto.h5ad",
                                          python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_brain_gobject,
           dir = './COSMX_Giotto/fov_brain', overwrite = TRUE)

fov_brain <- loadGiotto('./COSMX_Giotto/fov_brain/saveGiottoDir')

# 5.6- Multi-site recurrence
fov_multisite_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_multisite_giotto.h5ad",
                                          python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_multisite_gobject,
           dir = './COSMX_Giotto/fov_multisite', overwrite = TRUE)

fov_multisite <- loadGiotto('./COSMX_Giotto/fov_integrated/saveGiottoDir')

# 5.7- Epithelial cluster 1
fov_ep1_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_ep1_giotto.h5ad",
                                         python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_ep1_gobject,
           dir = './COSMX_Giotto/fov_ep1', overwrite = TRUE)

fov_ep1 <- loadGiotto('./COSMX_Giotto/fov_ep1/saveGiottoDir')

# 5.8- Epithelial cluster 2
fov_ep2_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_ep2_giotto.h5ad",
                                   python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_ep2_gobject,
           dir = './COSMX_Giotto/fov_ep2', overwrite = TRUE)

fov_ep2 <- loadGiotto('./COSMX_Giotto/fov_ep2/saveGiottoDir')

# 5.9- Epithelial cluster 3
fov_ep3_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_ep3_giotto.h5ad",
                                   python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_ep3_gobject,
           dir = './COSMX_Giotto/fov_ep3', overwrite = TRUE)

fov_ep3 <- loadGiotto('./COSMX_Giotto/fov_ep3/saveGiottoDir')

# 5.10- Stromal cluster 1
fov_sma1_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_sma1_giotto.h5ad",
                                   python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_sma1_gobject,
           dir = './COSMX_Giotto/fov_sma1', overwrite = TRUE)

fov_sma1 <- loadGiotto('./COSMX_Giotto/fov_sma1/saveGiottoDir')

# 5.11- Stromal cluster 2
fov_sma2_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_sma2_giotto.h5ad",
                                   python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_sma2_gobject,
           dir = './COSMX_Giotto/fov_sma2', overwrite = TRUE)

fov_sma2 <- loadGiotto('./COSMX_Giotto/fov_sma2/saveGiottoDir')

# 5.12- Stromal cluster 3
fov_sma3_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/ad_sma3_giotto.h5ad",
                                   python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_sma3_gobject,
           dir = './COSMX_Giotto/fov_sma3', overwrite = TRUE)

fov_sma3 <- loadGiotto('./COSMX_Giotto/fov_sma3/saveGiottoDir')


# 5.13- adata_fov_merged
fov_merged_gobject <- anndataToGiotto(anndata_path = "./2_h5ad_files/adata_fov_merged.h5ad",
                                          python_path = my_python_path, deluanay_spat_net = FALSE) #it cannot contain prior dim reductions but keep the spatial locations >> remove them python

saveGiotto(gobject = fov_merged_gobject,
           dir = './COSMX_Giotto/fov_merged', overwrite = TRUE)

fov_integrated <- loadGiotto('./COSMX_Giotto/fov_merged/saveGiottoDir')




