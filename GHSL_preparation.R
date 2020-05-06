library(raster)
library(sp)

rm(list=ls())

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

### load data
GHSL_ESP_2014 <- raster("original/GHSL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
GHSL_ESP_1990 <- raster("original/GHSL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
sevilla_boundaries <- shapefile("created/GADM/Sevilla_boundaries.shp")


reprojectAndCrop <- function(ghsl, boundary, epsg) {
  
  # crs in proj.4 format
  coordSys <- paste("+init=epsg:", epsg, sep = "")
  
  # reproject
  reprojectedRaster <- projectRaster(ghsl, res=250, crs = coordSys, method = 'ngb')
  reprojectedBoundaries <- spTransform(boundary, coordSys)
  
  # crop and mask
  ghsl_crop <- crop(reprojectedRaster,reprojectedBoundaries)
  ghsl_bound <- mask(ghsl_crop,reprojectedBoundaries)

  return(ghsl_bound)
}


makeBinary <- function(ghsl, threshold) {
  # determine binary threshold
  class.m <- c(0, threshold, 0, threshold, 100, 1)
  rcl.m <- matrix(class.m, ncol = 3, byrow = T)
  
  # reclassify
  ghsl_threshold <- reclassify(ghsl, rcl.m)
  return(ghsl_threshold)
}


getChange <- function(ghsl_early, ghsl_late, boundary, epsg, threshold) {
  # reproject, crop and use threshold on ghsl data
  ghsl_e_crop_bin <- makeBinary(reprojectAndCrop(ghsl_early, boundary, epsg), threshold)
  ghsl_l_crop_bin <- makeBinary(reprojectAndCrop(ghsl_late, boundary, epsg), threshold)
  
  # find built-up change
  builtUpChange <- (ghsl_l_crop_bin - ghsl_e_crop_bin)
  
  plot(builtUpChange)
  
  return(builtUpChange)
}


change_Sevilla <- getChange(GHSL_ESP_1990, GHSL_ESP_2014, sevilla_boundaries, 25830, 50)

writeRaster(change_Sevilla, "created/inR/sevilla_change.tif")


###############################################################################################
##########  old stuff   #######################################################################
###############################################################################################

plot(sevilla_boundaries)

## reproject
GHSL_2014_ETRS <- projectRaster(GHSL_ESP_2014, res=250, crs = "+init=epsg:25830", method = 'ngb')
Sevilla_boundaries_ETRS <- spTransform(sevilla_boundaries, "+init=epsg:25830")

ext <- extent(Sevilla_boundaries_ETRS)

## crop 
GHSL_crop <- crop(GHSL_2014_ETRS, c(ext[1]-500, ext[2]+500, ext[3]-500, ext[4]+500))
GHSL_2014_Sevilla <- mask(GHSL_crop, Sevilla_boundaries_ETRS)


class.m <- c(0, 50, 0, 50, 100, 1)
rcl.m <- matrix(class.m, ncol = 3, byrow = T)
rcl.m

GHSL_2014_Sevilla_50per <- reclassify(GHSL_2014_Sevilla, rcl.m)
plot(GHSL_2014_Sevilla_50per)

spplot(GHSL_2014_Sevilla)



crs(GHSL_2014_ETRS)
crs(Sevilla_boundaries_ETRS)
identicalCRS(GHSL_2014_ETRS, Sevilla_boundaries_ETRS)
rcl.m
