rm(list=ls())

library(raster)
library(sp)
library(sf)
library(rgdal)

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

###############################################
# to use train stations (OSM point data) it first has to be converted in singlepart feature 
# -> QGIS: multipart to singlepart
###############################################

dresden_boundaries <- shapefile("created/GADM/Dresden_boundaries.shp")
dresden_roads <- shapefile("created/OSM/Dresden_mRoads.shp")
dresden_river <- shapefile("created/OSM/Dresden_mRiver.shp")
dresden_tStations <- shapefile("created/OSM/Dresden_trainStation_singlepart.shp")
dresden_center <- shapefile("created/OSM/Dresden_townhall.shp")
dresden_airport <- shapefile("created/OSM/Dresden_airport.shp")


sevilla_boundaries <- shapefile("created/GADM/Sevilla_boundaries.shp")
sevilla_roads <- shapefile("created/OSM/Sevilla_mRoads.shp")
sevilla_river <- shapefile("created/OSM/Sevilla_mRiver.shp")
sevilla_tStations <- shapefile("created/OSM/Sevilla_trainStation_singlepart.shp")
sevilla_center <- shapefile("created/OSM/Sevilla_townhall.shp")
sevilla_airport <- shapefile("created/OSM/Sevilla_airport.shp")

krakow_boundaries <- shapefile("created/GADM/Krakow_boundaries.shp")
krakow_roads <- shapefile("created/OSM/Krakow_mRoads.shp")
krakow_river <- shapefile("created/OSM/Krakow_mRiver.shp")
krakow_tStations <- shapefile("created/OSM/Krakow_trainStation_singlepart.shp")
krakow_center <- shapefile("created/OSM/Krakow_townhall.shp")
krakow_airport <- shapefile("created/OSM/Krakow_airport.shp")



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
  
  # clip on city boundaries
  city_cropped <- crop(distances, boundaries_reproj)
  city_dist <- mask(city_cropped, boundaries_reproj)
  
  plot(city_dist)
  
  return(city_dist)
}



### calculate and save distance rasters

dresden_road_dist <- calc_dist_raster(dresden_roads, dresden_boundaries, 250, 32633)
writeRaster(dresden_road_dist, "created/distances/dresden_roads_250m.tif", overwrite=T)
dresden_river_dist <- calc_dist_raster(dresden_river, dresden_boundaries, 250, 32633)
writeRaster(dresden_river_dist, "created/distances/dresden_river_250m.tif", overwrite=T)
dresden_tStation_dist <- calc_dist_raster(dresden_tStations, dresden_boundaries, 250, 32633)
writeRaster(dresden_tStation_dist, "created/distances/dresden_tStation_250m.tif", overwrite=T)
dresden_airport_dist <- calc_dist_raster(dresden_airport, dresden_boundaries, 250, 32633)
writeRaster(dresden_airport_dist, "created/distances/dresden_airport_250m.tif", overwrite=T)
dresden_center_dist <- calc_dist_raster(dresden_center, dresden_boundaries, 250, 32633)
writeRaster(dresden_center_dist, "created/distances/dresden_center_250m.tif", overwrite=T)

sevilla_road_dist <- calc_dist_raster(sevilla_roads, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_road_dist, "created/distances/sevilla_roads_250m.tif", overwrite=T)
sevilla_river_dist <- calc_dist_raster(sevilla_river, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_river_dist, "created/distances/sevilla_river_250m.tif", overwrite=T)
sevilla_tStation_dist <- calc_dist_raster(sevilla_tStations, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_tStation_dist, "created/distances/sevilla_tStation_250m.tif", overwrite=T)
sevilla_airport_dist <- calc_dist_raster(sevilla_airport, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_airport_dist, "created/distances/sevilla_airport_250m.tif", overwrite=T)
sevilla_center_dist <- calc_dist_raster(sevilla_center, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_center_dist, "created/distances/sevilla_center_250m.tif", overwrite=T)

krakow_road_dist <- calc_dist_raster(krakow_roads, krakow_boundaries, 250, 32634)
writeRaster(krakow_road_dist, "created/distances/krakow_roads_250m.tif", overwrite=T)
krakow_river_dist <- calc_dist_raster(krakow_river, krakow_boundaries, 250, 32634)
writeRaster(krakow_river_dist, "created/distances/krakow_river_250m.tif", overwrite=T)
krakow_tStation_dist <- calc_dist_raster(krakow_tStations, krakow_boundaries, 250, 32634)
writeRaster(krakow_tStation_dist, "created/distances/krakow_tStation_250m.tif", overwrite=T)
krakow_airport_dist <- calc_dist_raster(krakow_airport, krakow_boundaries, 250, 32634)
writeRaster(krakow_airport_dist, "created/distances/krakow_airport_250m.tif", overwrite=T)
krakow_center_dist <- calc_dist_raster(krakow_center, krakow_boundaries, 250, 32634)
writeRaster(krakow_center_dist, "created/distances/krakow_center_250m.tif", overwrite=T)






#########################################################################
###     old   ###########################################################
#########################################################################

roads_reproj <- spTransform(dresden_roads, crs(dresden_change))

ext <- extent(c(extent(roads_reproj)[1]-500, extent(roads_reproj)[2]+500,extent(roads_reproj)[3]-500,extent(roads_reproj)[4]+500))

raster_template <- raster(ext, resolution = 250, crs = crs(roads_reproj))

rastered <- rasterize(roads_reproj, raster_template, field=1)

distances <- distance(rastered)
plot(distances)

plot(dresden_boundaries)
dresden_dist <- mask(distances, dresden_boundaries, filename = "created/distances/dresden_roads.tif", overwrite = T)

plot(dresden_dist)