library(raster)
library(sp)
library(sf)

rm(list=ls())

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

dresden_change <- raster("created/inR/dresden_change_30_m.tif")
dresden_roads <- shapefile("created/OSM/Dresden_mRoads.shp")
dresden_boundaries <- spTransform(shapefile("created/GADM/Dresden_boundaries.shp"), crs(dresden_change))

roads_reproj <- spTransform(dresden_roads, crs(dresden_change))

ext <- extent(c(extent(roads_reproj)[1]-500, extent(roads_reproj)[2]+500,extent(roads_reproj)[3]-500,extent(roads_reproj)[4]+500))

raster_template <- raster(ext, resolution = 250, crs = crs(roads_reproj))

rastered <- rasterize(roads_reproj, raster_template, field=1)

distances <- distance(rastered)
plot(distances)

plot(dresden_boundaries)
dresden_dist <- mask(distances, dresden_boundaries, filename = "created/distances/dresden_roads.tif", overwrite = T)

plot(dresden_dist)
