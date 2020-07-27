library(raster)
library(sp)



reprojectAndCrop <- function(raster, boundary, epsg, resolution) {
  
  # crs in proj.4 format
  coordSys <- paste("+init=epsg:", epsg, sep = "")
  
  
  ## crop for 30 m resolution to reduce needed computation resources
  boundaries_reproj <- spTransform(boundary, crs(raster))
  ext <- extent(boundaries_reproj)
  raster_cropped <- crop(raster, c(ext[1]-500, ext[2]+500, ext[3]-500, ext[4]+500))

  
  # reproject
  reprojectedRaster <- projectRaster(raster_cropped, res=resolution, crs = coordSys, method = 'ngb')


  return(reprojectedRaster)
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
  # 2: not built up in any epoch -> 0
  # 0, 1, 5, 6: no data, water, built up before
  class.m <- c(0, 1, NA, # 0-1
               1, 2, 0,  # 2
               2, 4, 1,  # 3-4
               4, 6, NA) # 5-6
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


reclassify_landuse <- function(landuse_raster, year = 1990) {
  # 1: artificial
  # 2: crop
  # 3: pasture
  # 4: forest
  # 5: open / bare land
  # 6: water
  if(year == 2000) {
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
  } else {
    class.m <- c(1, 11, 1, 
                 12, 17, 2, 
                 19, 22, 2,
                 18, 18, 3,
                 23, 25, 4,
                 28, 29, 4,
                 26, 27, 5,
                 30, 34, 5,
                 35, 38, 5,
                 39, 44, 6)
  }
  
  
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


calc_dist_raster <- function(osm, boundaries, resolution, epsg) {
  # reproject vector data
  coordSys <- paste("+init=epsg:", epsg, sep = "")
  osm_reproj <- spTransform(osm, coordSys)
  boundaries_reproj <- spTransform(boundaries, coordSys)
  
  # get and increase  extent of osm (or other vector) data
  ext <- extent(c(extent(boundaries_reproj)[1]-5000, extent(boundaries_reproj)[2]+5000,extent(boundaries_reproj)[3]-5000,extent(boundaries_reproj)[4]+5000))
  
  # create raster template
  raster_template <- raster(ext, resolution = resolution, crs = coordSys)
  
  # rasterize vector data
  rasterized <- rasterize(osm_reproj, raster_template, field = 1)
  
  # calculate euclidean distances
  distances <- distance(rasterized)
  
  city_dist <- reprojectAndCrop(distances, boundaries, epsg, resolution)
  
  plot(city_dist)
  
  return(city_dist)
}


calc_builtup_density <- function(ghsl_30m, boundary, epsg, window_size) {
  
  ghsl_reprojected <- reprojectAndCrop(ghsl_30m, boundary, epsg, 30)
  
  ### reclassify to (not) built-up
  # 5-6: built before 1990
  # 0-4: not built-up berfore 1990
  class.m <- c(0, 4, 0,
               4, 6, 1)
  rcl.m <- matrix(class.m, ncol = 3, byrow = T)
  ghsl_changed <- reclassify(ghsl_reprojected, rcl.m, include.lowest = T)
  
  # count cells within 7 x 7 window
  builtupCells <- focal(ghsl_changed, w=matrix(1, nc=window_size, nr=window_size))
  
  builtupDensity <- (builtupCells/(window_size*window_size-1))*100
  
  plot(builtupDensity)

  return(builtupDensity)
}


create_stack <- function(ghsl_30m, epsg, boundary, ghsl_pop, builtupDens_windowSizw, slope, landuse, road, primary_road, river, train_stations, city_center, airport) {
  
  change <- getChangeFromMultitemp(ghsl_30m, boundary, epsg, 30)
  
  builtup_density <- calc_builtup_density(ghsl_30m, boundary, epsg, builtupDens_windowSizw)
  pop_density <- reprojectAndCrop(ghsl_pop, boundary, epsg, 250)
  slope <- citySlopeAsPercentage(slope, boundary, epsg)
  landuse <- cropAndReclassify_landuse(landuse, boundary, epsg, 100)
  dist_mRoad <- calc_dist_raster(road, boundary, 120, epsg)
  dist_pRoad <- calc_dist_raster(primary_road, boundary, 120, epsg)
  dist_river <- calc_dist_raster(river, boundary, 120, epsg)
  dist_train <- calc_dist_raster(train_stations, boundary, 120, epsg)
  dist_center <- calc_dist_raster(city_center, boundary, 120, epsg)
  dist_airport <- calc_dist_raster(airport, boundary, 120, epsg)
  
  builtup_density <- projectRaster(builtup_density, change)
  pop_density <- projectRaster(pop_density, change)
  slope <- projectRaster(slope, change)
  landuse <- projectRaster(landuse, change, method = 'ngb')
  dist_mRoad <- projectRaster(dist_mRoad, change)
  dist_pRoad <- projectRaster(dist_pRoad, change)
  dist_river <- projectRaster(dist_river, change)
  dist_train <- projectRaster(dist_train, change)
  dist_center <- projectRaster(dist_center, change)
  dist_airport <- projectRaster(dist_airport, change)
  
  change_stack <- stack(change, builtup_density, pop_density, slope, landuse, dist_mRoad, dist_pRoad, dist_river, dist_train, dist_center, dist_airport)
  
  # crop and mask
  coordSys <- paste("+init=epsg:", epsg, sep = "")
  boundaries_reproj <- spTransform(boundary, coordSys)
  change_stack <- crop(change_stack, boundaries_reproj)
  change_stack <- mask(change_stack, boundaries_reproj)
  
  names(change_stack) <- c("change", "built_dens", "pop_dens", "slope", "landuse", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist")
  
  return(change_stack)
}
