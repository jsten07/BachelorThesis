library(raster)
library(sp)

# setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")



reprojectAndCrop <- function(raster, boundary, epsg, resolution) {
  
  # crs in proj.4 format
  coordSys <- paste("+init=epsg:", epsg, sep = "")
  
  
  ## crop for 30 m resolution to reduce needed computation resources
  boundaries_reproj <- spTransform(boundary, crs(raster))
  ext <- extent(boundaries_reproj)
  raster_cropped <- crop(raster, c(ext[1]-500, ext[2]+500, ext[3]-500, ext[4]+500))
  #raster_cropped <- crop(raster, boundaries_reproj)
  
  
  # reproject
  reprojectedRaster <- projectRaster(raster_cropped, res=resolution, crs = coordSys, method = 'ngb')
  reprojectedBoundaries <- spTransform(boundary, coordSys)
  
  # crop and mask
  raster_crop <- crop(reprojectedRaster,reprojectedBoundaries)
  raster_bound <- mask(raster_crop,reprojectedBoundaries)

  return(raster_bound)
}


makeBinary <- function(ghsl, threshold) {
  # determine binary threshold
  class.m <- c(0, threshold, 0, threshold, 100, 1)
  rcl.m <- matrix(class.m, ncol = 3, byrow = T)
  
  # reclassify
  ghsl_threshold <- reclassify(ghsl, rcl.m)
  return(ghsl_threshold)
}


getChange <- function(ghsl_early, ghsl_late, boundary, epsg, resolution, threshold) {
  # reproject, crop and use threshold on ghsl data
  ghsl_e_crop_bin <- makeBinary(reprojectAndCrop(ghsl_early, boundary, epsg, resolution), threshold)
  ghsl_l_crop_bin <- makeBinary(reprojectAndCrop(ghsl_late, boundary, epsg, resolution), threshold)
  
  # find built-up change
  builtUpChange <- (ghsl_l_crop_bin - ghsl_e_crop_bin)
  
  plot(builtUpChange)
  
  return(builtUpChange)
}


getChangeFromMultitemp <- function(ghsl, boundary, epsg, resolution) {
  ghsl_crop <- reprojectAndCrop(ghsl, boundary, epsg, resolution)
  
  # 3-4: changed from 1990 to 2014 -> 1
  # 0-2 & 5-6 not changed from 1990 to 2014 -> 0
  class.m <- c(0, 2, 0,
               2, 4, 1, 
               4, 6, 0)
  rcl.m <- matrix(class.m, ncol = 3, byrow = T)
  
  # reclassify
  ghsl_changed <- reclassify(ghsl_crop, rcl.m, include.lowest = T)
  
  plot(ghsl_changed)
  
  return(ghsl_changed)
}


DNtoPercentage <- function(DN) {
  # convert pixel values to radians
  rad <- (acos(DN/250))
  # convert radians to percentage slope
  percentage <- (tan(rad)*100)
  
  return(percentage)
}


citySlopeAsPercentage <- function(slope_raster, boundary, epsg) {
  # clip and reproject slope dataset
  slope_reprojected <- reprojectAndCrop(slope_raster, boundary, epsg, 25)
  
  # convert to percentage
  slope_per <- DNtoPercentage(slope_reprojected)
  
  plot(slope_per)
  return(slope_per)
}


reclassify_landuse <- function(landuse_raster) {
  # 1: artificial
  # 2: crop
  # 3: pasture
  # 4: forest
  # 5: open / bare land
  # 6: water
  class.m <- c(111, 142, 1, 
               211, 223, 2, 
               241, 244, 2,
               231, 231, 3,
               311, 313, 4,
               323, 324, 4,
               321, 322, 5,
               331, 335, 5,
               411, 422, 5,
               423, 523, 6)
  rcl.m <- matrix(class.m, ncol = 3, byrow = T)
  # +1 in the 'to' column to include this value (interval will be open on the right, to be closed left)
  rcl.m[,2] <- rcl.m[,2]+1
  
  # reclassify
  landuse_reclassified <- reclassify(landuse_raster, rcl.m, include.lowest=T, right=FALSE)
  
  return(landuse_reclassified)
}


cropAndReclassify_landuse <- function(landuse, boundary, epsg, resolution) {
  landuse_crop <- reprojectAndCrop(landuse, boundary, epsg, resolution)
  landuse_recl <- reclassify_landuse(landuse_crop)
  
  plot(landuse_recl)
  return(landuse_recl)
}


### load data

# GHSL_ESP_2014 <- raster("original/GHSL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
# GHSL_ESP_1990 <- raster("original/GHSL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
# GHSL_ESP_30m <- raster("original/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_13_8_ESP/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_13_8.tif")
# GHSL_POL_2014 <- raster("original/GHSL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_19_3_POL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_19_3.tif")
# GHSL_POL_1990 <- raster("original/GHSL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_19_3_POL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_19_3.tif")
# GHSL_POL_30m <- raster("original/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_14_7_POL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_14_7.tif")
# sevilla_boundaries <- shapefile("created/GADM/Sevilla_boundaries.shp")
# krakow_boundaries <- shapefile("created/GADM/Krakow_boundaries.shp")
# dresden_boundaries <- shapefile("created/GADM/Dresden_boundaries.shp")
# 
# GHSL_pop_ESP <- raster("original/GHSL/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_17_4_ESP/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_17_4.tif")
# GHSL_pop_POL <- raster("original/GHSL/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_19_3_POL/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_19_3.tif")
# 
# slope_ESP <- raster("original/slope/EUD_CP-SLOP_2500015000-AA_ESP/EUD_CP-SLOP_2500015000-AA.tif")
# slope_GER <- raster("original/slope/EUD_CP-SLOP_4500035000-AA_GER/EUD_CP-SLOP_4500035000-AA.tif")
# slope_POL <- raster("original/slope/EUD_CP-SLOP_5500035000-AA_POL/EUD_CP-SLOP_5500035000-AA.tif")



### process data

# change_Sevilla <- getChange(GHSL_ESP_1990, GHSL_ESP_2014, sevilla_boundaries, epsg = 32630, resolution = 250, threshold = 50)
# writeRaster(change_Sevilla, "created/inR/sevilla_change.tif", overwrite=T)
# change_Sevilla_30m <- getChangeFromMultitemp(GHSL_ESP_30m, sevilla_boundaries, 32630, resolution = 30)
# writeRaster(change_Sevilla_30m, "created/inR/sevilla_change_30_m.tif", overwrite=T)
# change_Krakow <- getChange(GHSL_POL_1990, GHSL_POL_2014, krakow_boundaries, epsg = 32634, resolution = 250, threshold = 50)
# writeRaster(change_Krakow, "created/inR/krakow_change.tif", overwrite=T)
# change_Krakow_30m <- getChangeFromMultitemp(GHSL_POL_30m, krakow_boundaries, 32634, resolution = 30)
# writeRaster(change_Krakow_30m, "created/inR/krakow_change_30_m.tif", overwrite=T)
# change_Dresden <- getChange(GHSL_POL_1990, GHSL_POL_2014, dresden_boundaries, epsg = 32633, resolution = 250, threshold = 50)
# writeRaster(change_Dresden, "created/inR/dresden_change.tif", overwrite=T)
# change_Dresden_30m <- getChangeFromMultitemp(GHSL_POL_30m, dresden_boundaries, 32633, resolution = 30)
# writeRaster(change_Dresden_30m, "created/inR/dresden_change_30_m.tif", overwrite=T)
# 
# population_Sevilla <- reprojectAndCrop(GHSL_pop_ESP, sevilla_boundaries, 25830, 250)
# writeRaster(population_Sevilla, "created/inR/sevilla_popDens.tif", overwrite=T)
# population_Krakow <- reprojectAndCrop(GHSL_pop_POL, krakow_boundaries, 25834, 250)
# writeRaster(population_Krakow, "created/inR/krakow_popDens.tif", overwrite=T)
# population_Dresden <- reprojectAndCrop(GHSL_pop_POL, dresden_boundaries, 25833, 250)
# writeRaster(population_Dresden, "created/inR/dresden_popDens.tif", overwrite=T)
# 
# slope_Sevilla <- citySlopeAsPercentage(slope_ESP, sevilla_boundaries, 25830)
# writeRaster(slope_Sevilla, "created/slope/sevilla_slope.tif", overwrite=T)
# slope_Dresden <- citySlopeAsPercentage(slope_GER, dresden_boundaries, 25833)
# writeRaster(slope_Dresden, "created/slope/dresden_slope.tif", overwrite=T)
# slope_Krakow <- citySlopeAsPercentage(slope_POL, krakow_boundaries, 25834)
# writeRaster(slope_Krakow, "created/slope/krakow_slope.tif", overwrite=T)





###############################################################################################
##########  old stuff   #######################################################################
###############################################################################################

# plot(sevilla_boundaries)
# 
# ## reproject
# GHSL_2014_ETRS <- projectRaster(GHSL_ESP_2014, res=250, crs = "+init=epsg:25830", method = 'ngb')
# Sevilla_boundaries_ETRS <- spTransform(sevilla_boundaries, "+init=epsg:25830")
# 
# ext <- extent(Sevilla_boundaries_ETRS)
# 
# ## crop 
# GHSL_crop <- crop(GHSL_2014_ETRS, c(ext[1]-500, ext[2]+500, ext[3]-500, ext[4]+500))
# GHSL_2014_Sevilla <- mask(GHSL_crop, Sevilla_boundaries_ETRS)
# 
# 
# class.m <- c(0, 50, 0, 50, 100, 1)
# rcl.m <- matrix(class.m, ncol = 3, byrow = T)
# rcl.m
# 
# GHSL_2014_Sevilla_50per <- reclassify(GHSL_2014_Sevilla, rcl.m)
# plot(GHSL_2014_Sevilla_50per)
# 
# spplot(GHSL_2014_Sevilla)
# 
# 
# 
# crs(GHSL_2014_ETRS)
# crs(Sevilla_boundaries_ETRS)
# identicalCRS(GHSL_2014_ETRS, Sevilla_boundaries_ETRS)
# rcl.m
# 
# 
# boundaries_reproj <- spTransform(sevilla_boundaries, crs(GHSL_ESP_30m))
# ext <- extent(boundaries_reproj)
# ghsl_cropped <- crop(GHSL_ESP_30m, c(ext[1]-500, ext[2]+500, ext[3]-500, ext[4]+500))
# 
# # reproject
# reprojectedRaster <- projectRaster(ghsl_cropped, res=250, crs = "+init=epsg:25830", method = 'ngb')
# reprojectedBoundaries <- spTransform(sevilla_boundaries, "+init=epsg:25830")
# 
# # crop and mask
# ghsl_crop <- crop(reprojectedRaster,reprojectedBoundaries)
# ghsl_bound <- mask(ghsl_crop,reprojectedBoundaries)
# 
# plot(ghsl_bound)