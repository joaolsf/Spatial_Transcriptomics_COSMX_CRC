
options(max.print=1000000)
# 1- Confirm the metadata and the expression matrix files contains the same cells
mat <- read.csv("./fovs/fov_integrated_exprMat_file.csv")
mat$fov <- as.factor(mat$fov)
mat$fov
mat_fov01 <- mat[mat$fov %in% 1,]
mat_fov02 <- mat[mat$fov %in% 2,]
mat_fov03 <- mat[mat$fov %in% 3,]
mat_fov04 <- mat[mat$fov %in% 4,]
mat_fov05 <- mat[mat$fov %in% 5,]
mat_fov06 <- mat[mat$fov %in% 6,]
mat_fov07 <- mat[mat$fov %in% 7,]
mat_fov08 <- mat[mat$fov %in% 8,]
mat_fov09 <- mat[mat$fov %in% 9,]
mat_fov10 <- mat[mat$fov %in% 10,]
mat_fov11 <- mat[mat$fov %in% 11,]
mat_fov12 <- mat[mat$fov %in% 12,]
mat_fov13 <- mat[mat$fov %in% 13,]
mat_fov14 <- mat[mat$fov %in% 14,]
mat_fov15 <- mat[mat$fov %in% 15,]
mat_fov16 <- mat[mat$fov %in% 16,]
mat_fov17 <- mat[mat$fov %in% 17,]
mat_fov18 <- mat[mat$fov %in% 18,]
mat_fov19 <- mat[mat$fov %in% 19,]
mat_fov20 <- mat[mat$fov %in% 20,]
mat_fov21 <- mat[mat$fov %in% 21,]
mat_fov22 <- mat[mat$fov %in% 22,]
mat_fov23 <- mat[mat$fov %in% 23,]
mat_fov24 <- mat[mat$fov %in% 24,]
mat_fov25 <- mat[mat$fov %in% 25,]

write.csv(mat_fov01, "Fov01_exprMat_file.csv")
write.csv(mat_fov02, "Fov02_exprMat_file.csv")
write.csv(mat_fov03, "Fov03_exprMat_file.csv")
write.csv(mat_fov04, "Fov04_exprMat_file.csv")
write.csv(mat_fov05, "Fov05_exprMat_file.csv")
write.csv(mat_fov06, "Fov06_exprMat_file.csv")
write.csv(mat_fov07, "Fov07_exprMat_file.csv")
write.csv(mat_fov08, "Fov08_exprMat_file.csv")
write.csv(mat_fov09, "Fov09_exprMat_file.csv")
write.csv(mat_fov10, "Fov10_exprMat_file.csv")
write.csv(mat_fov11, "Fov11_exprMat_file.csv")
write.csv(mat_fov12, "Fov12_exprMat_file.csv")
write.csv(mat_fov13, "Fov13_exprMat_file.csv")
write.csv(mat_fov14, "Fov14_exprMat_file.csv")
write.csv(mat_fov15, "Fov15_exprMat_file.csv")
write.csv(mat_fov16, "Fov16_exprMat_file.csv")
write.csv(mat_fov17, "Fov17_exprMat_file.csv")
write.csv(mat_fov18, "Fov18_exprMat_file.csv")
write.csv(mat_fov19, "Fov19_exprMat_file.csv")
write.csv(mat_fov20, "Fov20_exprMat_file.csv")
write.csv(mat_fov21, "Fov21_exprMat_file.csv")
write.csv(mat_fov22, "Fov22_exprMat_file.csv")
write.csv(mat_fov23, "Fov23_exprMat_file.csv")
write.csv(mat_fov24, "Fov24_exprMat_file.csv")
write.csv(mat_fov25, "Fov25_exprMat_file.csv")

mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID


metadata <- read.csv("./fovs/fov_integrated_metadata_file.csv")
metadata$cell_ID <- as.factor(metadata$cell_ID)
metadata$cell_ID
levels(metadata$cell_ID)
setdiff(levels(metadata$cell_ID), levels(mat$cell_ID))

metadata$fov <- as.factor(metadata$fov)
metadata$fov

metadata_fov01 <- metadata[metadata$fov %in% 1,]
write.csv(metadata_fov01, "Fov01_metadata_file_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################


# 2- Check which cells should be removed from the polygons and tx files according to the cells in the metadata

#fov01
polygons <- read.csv("./fovs/fov01/Fov01-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov01/Fov01_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov01/Fov01_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01
[1] "0"    "190"  "227"  "378"  "397"  "450"  "517"  "591"  "647"  "655"  "689"  "698"  "708"  "732"  "745"  "751"  "797" 
[18] "879"  "1002" "1236" "1519" "1555" "1826" "1876" "1880" "1950" "2223" "2316" "2389" "2452" "2484" "2529" "2571" "2595"
[35] "2605" "2658" "2686" "2691" "2749" "2750" "2787" "2799" "2806" "2812" "2854" "2863" "2879" "2892" "2893" "2903" "2906"
[52] "2915" "2924" "2925" "2938" "2949" "2999" "3012" "3022" "3092" "3160" "3220" "3233" "3363" "3394" "3395" "3436" "3490"
[69] "3537" "3641" "3776" "3780" "3790" "3876" "3880" "3886" "4041" "4066" "4108" "4120"

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
tx2[tx2$cell_ID == "190", ]  
write.csv(tx2, "Fov01_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
polygons2[polygons2$cellID == "190", ]  
write.csv(polygons2, "Fov01-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov02
polygons <- read.csv("./fovs/fov02/Fov02-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov02/Fov02_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov02/Fov02_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01
[1] "0"    "1407" "2139" "2259" "2685" "2851" "3188" "3190"

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
tx2[tx2$cell_ID == "1407", ]  
write.csv(tx2, "Fov02_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
polygons2[polygons2$cellID == "1407", ]  
write.csv(polygons2, "Fov02-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov03
polygons <- read.csv("./fovs/fov03/Fov03-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov03/Fov03_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov03/Fov03_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov03_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov03-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov04
polygons <- read.csv("./fovs/fov04/Fov04-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov04/Fov04_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov04/Fov04_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov04_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov04-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov05
polygons <- read.csv("./fovs/fov05/Fov05-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov05/Fov05_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov05/Fov05_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov05_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov05-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov06
polygons <- read.csv("./fovs/fov06/Fov06-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov06/Fov06_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov06/Fov06_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov06_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov06-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov07
polygons <- read.csv("./fovs/fov07/Fov07-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov07/Fov07_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov07/Fov07_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov07_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov07-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov08
polygons <- read.csv("./fovs/fov08/Fov08-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov08/Fov08_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov08/Fov08_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov08_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov08-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov09
polygons <- read.csv("./fovs/fov09/Fov09-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov09/Fov09_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov09/Fov09_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov09_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov09-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov10
polygons <- read.csv("./fovs/fov10/Fov10-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov10/Fov10_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov10/Fov10_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov10_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov10-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov11
polygons <- read.csv("./fovs/fov11/Fov11-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov11/Fov11_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov11/Fov11_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov11_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov11-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov12
polygons <- read.csv("./fovs/fov12/Fov12-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov12/Fov12_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov12/Fov12_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov12_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov12-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov13
polygons <- read.csv("./fovs/fov13/Fov13-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov13/Fov13_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov13/Fov13_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov13_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov13-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov14
polygons <- read.csv("./fovs/fov14/Fov14-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov14/Fov14_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov14/Fov14_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
#setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov14_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov14-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov15
polygons <- read.csv("./fovs/fov15/Fov15-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov15/Fov15_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov15/Fov15_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
#setdiff(levels(polygon$cellID), levels(tx$cell_ID))

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov15_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov15-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov16
polygons <- read.csv("./fovs/fov16/Fov16-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov16/Fov16_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov16/Fov16_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov16_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov16-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov17
polygons <- read.csv("./fovs/fov17/Fov17-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov17/Fov17_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov17/Fov17_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov17_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov17-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov18
polygons <- read.csv("./fovs/fov18/Fov18-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov18/Fov18_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov18/Fov18_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov18_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov18-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov19
polygons <- read.csv("./fovs/fov19/Fov19-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov19/Fov19_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov19/Fov19_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov19_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov19-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov20
polygons <- read.csv("./fovs/fov20/Fov20-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov20/Fov20_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov20/Fov20_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov20_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov20-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov21
polygons <- read.csv("./fovs/fov21/Fov21-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov21/Fov21_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov21/Fov21_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov21_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov21-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov22
polygons <- read.csv("./fovs/fov22/Fov22-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov22/Fov22_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov22/Fov22_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov22_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov22-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov23
polygons <- read.csv("./fovs/fov23/Fov23-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov23/Fov23_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov23/Fov23_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov23_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov23-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov24
polygons <- read.csv("./fovs/fov24/Fov24-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov24/Fov24_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov24/Fov24_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov24_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov24-polygons_v2.csv")

################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

#fov25
polygons <- read.csv("./fovs/fov25/Fov25-polygons.csv")
polygons$cellID <- as.factor(polygons$cellID)
polygons$cellID
tx <- read.csv("./fovs/fov25/Fov25_tx_file.csv")
tx$cell_ID <- as.factor(tx$cell_ID)
tx$cell_ID
mat <- read.csv("./fovs/fov25/Fov25_exprMat_file.csv")
mat$cell_ID <- as.factor(mat$cell_ID)
mat$cell_ID

setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 
p <- setdiff(levels(tx$cell_ID), levels(mat$cell_ID)) 

# cells to exclude from tx file fov01

tx2 <- tx[!(tx$cell_ID %in% p), ]
tx2$cell_ID <- as.factor(tx2$cell_ID)
tx2$cell_ID
setdiff(levels(tx2$cell_ID), levels(mat$cell_ID)) 
write.csv(tx2, "Fov25_tx_file_v2.csv")

polygons2 <- polygons[!(polygons$cellID %in% p), ]
polygons2$cellID <- as.factor(polygons2$cellID)
polygons2$cellID
write.csv(polygons2, "Fov25-polygons_v2.csv")

# Combine all tx files

library(data.table)
df <- 
  list.files(path = "./", pattern = "*.csv") %>% 
  map_df(~fread(.))
df
fwrite(df, "S1_TMA_tx_file.csv")





