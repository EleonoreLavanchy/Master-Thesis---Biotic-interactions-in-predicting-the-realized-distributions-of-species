####################################################################################################
############################################## VISUALISATION #######################################
####################################################################################################
setwd("~/UNIL/Master/2th semester/Master Project/Data")
library("raster")
library("xlsx")
library("ggplot2")

############################################## MAP #######################################

coord <- read.csv("all_plots_coords.csv",sep = ";")
coord <- coord[coord$Plot != 874,]
coord_s_ID <- coord[,2:3]
DEM_raster <- raster("mnt25pe100.tif")

coordinates(coord_s_ID) <- ~X+Y
crs(coord_s_ID) <- DEM_raster@crs

pal=colorRampPalette(c("black", "white"))(100) 

plot(DEM_raster, main = "Map of sampling plots", xlab = "Latitude", ylab = "Longitude", col = pal)
points(coord_s_ID, pch = ".", col = "darkorange", cex = 2)

############################################## Barplot sp. occ #######################################

Plant_occ <- read.xlsx("Plants_Sp_tot.xlsx", sheetIndex = 1)

barplot(Plant_occ[,4:103] == "1", names.arg = seq(1,100), col = "black")
segments(x0 = 0, y0=70, x1=120 , y1=70, col = "red", lwd = 10)

######################################### Barplot sp. nb per plot #######################################

sum_pres_plot <- apply(Plant_occ, 1, function(x) length(Plant_occ[Plant_occ[,4:103] == "1"]))

