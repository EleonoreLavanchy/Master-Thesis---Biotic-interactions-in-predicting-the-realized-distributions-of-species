####################################################################################################
############################################## VISUALISATION #######################################
####################################################################################################
setwd("~/UNIL/Master/2th semester/Master Project/Data")
library("raster")
library("xlsx")
library("ggplot2")
library(tidyr)

############################################## MAP #######################################

coord <- read.csv("all_plots_coords.csv",sep = ";")
coord <- coord[coord$Plot != 874,]
coord_s_ID <- coord[,2:3]
DEM_raster <- raster("mnt25pe100.tif")

coordinates(coord_s_ID) <- ~X+Y
crs(coord_s_ID) <- DEM_raster@crs

pal=colorRampPalette(c("black", "white"))(100) 

plot(DEM_raster, main = "Map of the sampled sites", xlab = "Latitude", ylab = "Longitude", col = pal)
points(coord_s_ID, pch = ".", col = "darkorange", cex = 2.5)

############################################## Barplot sp. occ #######################################

Plant_occ <- read.xlsx("Plants_Sp_tot.xlsx", sheetIndex = 1)

occ <- apply(Plant_occ[4:103], MARGIN = 2, FUN = sum)
occ <- cbind(occ, seq(1:100))
colnames(occ) <- c("occurrences", "species")
occ <- as.data.frame(occ)
occ$species <- as.factor(occ$species)

ggplot(occ, aes(x = species, y = occurrences)) + geom_col(fill = "black") +
  labs(x = "Species", y = "Number of occurences") + geom_abline(slope = 0, intercept = 70, col = "red", size = 1) +
  scale_x_discrete(breaks = seq(0,100,5)) +
  theme(axis.line = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
        axis.title = element_text(size = 18), axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

######################################### Barplot sp. nb per plot #######################################

sum_pres_plot <- apply(Plant_occ, 1, function(x) length(Plant_occ[Plant_occ[,4:103] == "1"]))

