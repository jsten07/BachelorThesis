rm(list=ls())

### load functions
source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/preparation_functions.R")
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

########################################################################################
### load data

sevilla_boundaries <- shapefile("created/GADM/Sevilla_boundaries.shp")
krakow_boundaries <- shapefile("created/GADM/Krakow_boundaries.shp")
dresden_boundaries <- shapefile("created/GADM/Dresden_boundaries.shp")

GHSL_ESP_2014 <- raster("original/GHSL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
GHSL_ESP_1990 <- raster("original/GHSL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
GHSL_POL_2014 <- raster("original/GHSL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_19_3_POL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_19_3.tif")
GHSL_POL_1990 <- raster("original/GHSL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_19_3_POL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_19_3.tif")
GHSL_ESP_30m <- raster("original/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_13_8_ESP/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_13_8.tif")
GHSL_GER_30m <- raster("original/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_14_7_POL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_14_7.tif")
GHSL_POL_30m_wrong <- raster("original/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_15_7_POL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_15_7.tif")
GHSL_POL_ext <- c(2200000,2260000,6430000,6480000)
GHSL_POL_30m <- mosaic(crop(GHSL_GER_30m, GHSL_POL_ext), crop(GHSL_POL_30m_2, GHSL_POL_ext), fun = mean)


GHSL_pop_ESP <- raster("original/GHSL/GHS_POP_E1990_GLOBE_R2019A_54009_250_V1_0_17_4_ESP/GHS_POP_E1990_GLOBE_R2019A_54009_250_V1_0_17_4.tif")
GHSL_pop_POL <- raster("original/GHSL/GHS_POP_E1990_GLOBE_R2019A_54009_250_V1_0_19_3_POL/GHS_POP_E1990_GLOBE_R2019A_54009_250_V1_0_19_3.tif")

slope_ESP <- raster("original/slope/EUD_CP-SLOP_2500015000-AA_ESP/EUD_CP-SLOP_2500015000-AA.tif")
slope_GER <- raster("original/slope/EUD_CP-SLOP_4500035000-AA_GER/EUD_CP-SLOP_4500035000-AA.tif")
slope_POL <- raster("original/slope/EUD_CP-SLOP_5500035000-AA_POL/EUD_CP-SLOP_5500035000-AA.tif")

dresden_roads <- shapefile("created/OSM/Dresden_mRoads.shp")
dresden_primary_roads <- shapefile("created/OSM/Dresden_pRoads.shp")
dresden_river <- shapefile("created/OSM/Dresden_mRiver.shp")
dresden_tStations <- shapefile("created/OSM/Dresden_trainStation_singlepart.shp")
dresden_center <- shapefile("created/OSM/Dresden_townhall.shp")
dresden_airport <- shapefile("created/OSM/Dresden_airport.shp")

sevilla_roads <- shapefile("created/OSM/Sevilla_mRoads.shp")
sevilla_primary_roads <- shapefile("created/OSM/Sevilla_pRoads.shp")
sevilla_river <- shapefile("created/OSM/Sevilla_mRiver.shp")
sevilla_tStations <- shapefile("created/OSM/Sevilla_trainStation_singlepart.shp")
sevilla_center <- shapefile("created/OSM/Sevilla_townhall.shp")
sevilla_airport <- shapefile("created/OSM/Sevilla_airport.shp")

krakow_roads <- shapefile("created/OSM/Krakow_mRoads.shp")
krakow_primary_roads <- shapefile("created/OSM/Krakow_pRoads.shp")
krakow_river <- shapefile("created/OSM/Krakow_mRiver.shp")
krakow_tStations <- shapefile("created/OSM/Krakow_trainStation_singlepart.shp")
krakow_center <- shapefile("created/OSM/Krakow_townhall.shp")
krakow_airport <- shapefile("created/OSM/Krakow_airport.shp")

landuse_2000 <- raster("original/land cover/clc2006_clc2000_v2018_20_raster100m/CLC2006_CLC2000_V2018_20.tif")
landuse <- raster("original/land cover/u2000_clc1990_v2020_20u1_raster100m/DATA/U2000_CLC1990_V2020_20u1.tif")


########################################################################################
# create stack and save grid and csv

stack_sevilla <- create_stack(GHSL_ESP_30m, 32630, sevilla_boundaries, GHSL_pop_ESP, 5, slope_ESP, landuse, 
                              sevilla_roads, sevilla_primary_roads, sevilla_river, sevilla_tStations, sevilla_center, sevilla_airport)
stack_sevilla.df <- as.data.frame(getValues(stack_sevilla))
stack_sevilla.df_noNA <- stack_sevilla.df[which(stack_sevilla.df$change!="NA"), ]
write.csv(stack_sevilla.df_noNA, "created/stack/sevilla_lu90.csv")
writeRaster(stack_sevilla, "created/stack/sevilla_lu90.grd", overwrite = TRUE)

stack_dresden <- create_stack(GHSL_GER_30m, 32633, dresden_boundaries, GHSL_pop_POL, 5, slope_GER, landuse, 
                              dresden_roads, dresden_primary_roads, dresden_river, dresden_tStations, dresden_center, dresden_airport)
stack_dresden.df <- as.data.frame(getValues(stack_dresden))
stack_dresden.df_noNA <- stack_dresden.df[which(stack_dresden.df$change!="NA"), ]
write.csv(stack_dresden.df_noNA, "created/stack/dresden_li90.csv")
writeRaster(stack_dresden, "created/stack/dresden_lu90.grd", overwrite = TRUE)

stack_krakow <- create_stack(GHSL_POL_30m, 32634, krakow_boundaries, GHSL_pop_POL, 5, slope_POL, landuse, 
                              krakow_roads, krakow_primary_roads, krakow_river, krakow_tStations, krakow_center, krakow_airport)
stack_krakow.df <- as.data.frame(getValues(stack_krakow))
stack_krakow.df_noNA <- stack_krakow.df[which(stack_krakow.df$change!="NA"), ]
write.csv(stack_krakow.df_noNA, "created/stack/krakow_lu90.csv")
writeRaster(stack_krakow, "created/stack/krakow_lu90.grd", overwrite = TRUE)



########################################################################################
########################################################################################
### the following code was at the end not used for the thesis
########################################################################################
########################################################################################


########################################################################################
### process and save data

# 250 m change

change_Sevilla <- getChange(GHSL_ESP_1990, GHSL_ESP_2014, sevilla_boundaries, epsg = 32630, resolution = 250, threshold = 50)
writeRaster(change_Sevilla, "created/GHSL_R/sevilla_change.tif", overwrite=T)
change_Krakow <- getChange(GHSL_POL_1990, GHSL_POL_2014, krakow_boundaries, epsg = 32634, resolution = 250, threshold = 50)
writeRaster(change_Krakow, "created/GHSL_R/krakow_change.tif", overwrite=T)
change_Dresden <- getChange(GHSL_POL_1990, GHSL_POL_2014, dresden_boundaries, epsg = 32633, resolution = 250, threshold = 50)
writeRaster(change_Dresden, "created/GHSL_R/dresden_change.tif", overwrite=T)

# 30 m change 

change_Sevilla_30m <- getChangeFromMultitemp(GHSL_ESP_30m, sevilla_boundaries, 32630, resolution = 30)
writeRaster(change_Sevilla_30m, "created/GHSL_R/sevilla_change_30_m.tif", overwrite=T)
change_Krakow_30m <- getChangeFromMultitemp(GHSL_POL_30m, krakow_boundaries, 32634, resolution = 30)
writeRaster(change_Krakow_30m, "created/GHSL_R/krakow_change_30_m.tif", overwrite=T)
change_Dresden_30m <- getChangeFromMultitemp(GHSL_GER_30m, dresden_boundaries, 32633, resolution = 30)
writeRaster(change_Dresden_30m, "created/GHSL_R/dresden_change_30_m.tif", overwrite=T)

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/GHSL_R/")
writeRaster(mask(reprojectAndCrop(GHSL_ESP_30m, sevilla_boundaries, 32630, 30), spTransform(sevilla_boundaries, "+init=epsg:32630")), "sevilla_built.tif", overwrite = T)
writeRaster(mask(reprojectAndCrop(GHSL_GER_30m, dresden_boundaries, 32633, 30), spTransform(dresden_boundaries, "+init=epsg:32633")), "dresden_built.tif", overwrite = T)
writeRaster(mask(reprojectAndCrop(GHSL_POL_30m, krakow_boundaries, 32634, 30), spTransform(krakow_boundaries, "+init=epsg:32634")), "krakow_built.tif", overwrite = T)

# population density

population_Sevilla <- reprojectAndCrop(GHSL_pop_ESP, sevilla_boundaries, 32630, 250)
writeRaster(population_Sevilla, "created/GHSL_R/sevilla_popDens.tif", overwrite=T)
population_Krakow <- reprojectAndCrop(GHSL_pop_POL, krakow_boundaries, 32634, 250)
writeRaster(population_Krakow, "created/GHSL_R/krakow_popDens.tif", overwrite=T)
population_Dresden <- reprojectAndCrop(GHSL_pop_POL, dresden_boundaries, 32633, 250)
writeRaster(population_Dresden, "created/GHSL_R/dresden_popDens.tif", overwrite=T)

# slope

slope_Sevilla <- citySlopeAsPercentage(slope_ESP, sevilla_boundaries, 32630)
writeRaster(slope_Sevilla, "created/slope/sevilla_slope.tif", overwrite=T)
slope_Dresden <- citySlopeAsPercentage(slope_GER, dresden_boundaries, 32633)
writeRaster(slope_Dresden, "created/slope/dresden_slope.tif", overwrite=T)
slope_Krakow <- citySlopeAsPercentage(slope_POL, krakow_boundaries, 32634)
writeRaster(slope_Krakow, "created/slope/krakow_slope.tif", overwrite=T)

# landuse

landuse_Dresden <- cropAndReclassify_landuse(landuse, dresden_boundaries, 32633, 100)
writeRaster(landuse_Dresden, "created/landuse/dresden_landuse.tif", overwrite=T)
landuse_Sevilla <- cropAndReclassify_landuse(landuse, sevilla_boundaries, 32630, 100)
writeRaster(landuse_Sevilla, "created/landuse/sevilla_landuse.tif", overwrite=T)
landuse_Krakow <- cropAndReclassify_landuse(landuse, krakow_boundaries, 32634, 100)
writeRaster(landuse_Krakow, "created/landuse/krakow_landuse.tif", overwrite=T)

# builtup density

built_density_Dresden <- calc_builtup_density(GHSL_GER_30m, dresden_boundaries, 32633)
writeRaster(built_density_Dresden, "created/builtup_density/dresden_density.tif", overwrite=T)
built_density_Krakow <- calc_builtup_density(GHSL_POL_30m, krakow_boundaries, 32634)
writeRaster(built_density_Krakow, "created/builtup_density/krakow_density.tif", overwrite=T)
built_density_Sevilla <- calc_builtup_density(GHSL_ESP_30m, sevilla_boundaries, 32630)
writeRaster(built_density_Sevilla, "created/builtup_density/sevilla_density.tif", overwrite=T)

# distances dresden

dresden_road_dist <- calc_dist_raster(dresden_roads, dresden_boundaries, 250, 32633)
writeRaster(dresden_road_dist, "created/distances/dresden_roads_250m.tif", overwrite=T)
dresden_p_road_dist <- calc_dist_raster(dresden_primary_roads, dresden_boundaries, 250, 32633)
writeRaster(dresden_p_road_dist, "created/distances/dresden_primary_roads_250m.tif", overwrite=T)
dresden_river_dist <- calc_dist_raster(dresden_river, dresden_boundaries, 250, 32633)
writeRaster(dresden_river_dist, "created/distances/dresden_river_250m.tif", overwrite=T)
dresden_tStation_dist <- calc_dist_raster(dresden_tStations, dresden_boundaries, 250, 32633)
writeRaster(dresden_tStation_dist, "created/distances/dresden_tStation_250m.tif", overwrite=T)
dresden_airport_dist <- calc_dist_raster(dresden_airport, dresden_boundaries, 250, 32633)
writeRaster(dresden_airport_dist, "created/distances/dresden_airport_250m.tif", overwrite=T)
dresden_center_dist <- calc_dist_raster(dresden_center, dresden_boundaries, 250, 32633)
writeRaster(dresden_center_dist, "created/distances/dresden_center_250m.tif", overwrite=T)

# distances sevilla

sevilla_road_dist <- calc_dist_raster(sevilla_roads, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_road_dist, "created/distances/sevilla_roads_250m.tif", overwrite=T)
sevilla_p_road_dist <- calc_dist_raster(sevilla_primary_roads, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_p_road_dist, "created/distances/sevilla_primary_roads_250m.tif", overwrite=T)
sevilla_river_dist <- calc_dist_raster(sevilla_river, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_river_dist, "created/distances/sevilla_river_250m.tif", overwrite=T)
sevilla_tStation_dist <- calc_dist_raster(sevilla_tStations, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_tStation_dist, "created/distances/sevilla_tStation_250m.tif", overwrite=T)
sevilla_airport_dist <- calc_dist_raster(sevilla_airport, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_airport_dist, "created/distances/sevilla_airport_250m.tif", overwrite=T)
sevilla_center_dist <- calc_dist_raster(sevilla_center, sevilla_boundaries, 250, 32630)
writeRaster(sevilla_center_dist, "created/distances/sevilla_center_250m.tif", overwrite=T)

#distances krakow

krakow_road_dist <- calc_dist_raster(krakow_roads, krakow_boundaries, 250, 32634)
writeRaster(krakow_road_dist, "created/distances/krakow_roads_250m.tif", overwrite=T)
krakow_p_road_dist <- calc_dist_raster(krakow_primary_roads, krakow_boundaries, 250, 32634)
writeRaster(krakow_p_road_dist, "created/distances/krakow_primary_roads_250m.tif", overwrite=T)
krakow_river_dist <- calc_dist_raster(krakow_river, krakow_boundaries, 250, 32634)
writeRaster(krakow_river_dist, "created/distances/krakow_river_250m.tif", overwrite=T)
krakow_tStation_dist <- calc_dist_raster(krakow_tStations, krakow_boundaries, 250, 32634)
writeRaster(krakow_tStation_dist, "created/distances/krakow_tStation_250m.tif", overwrite=T)
krakow_airport_dist <- calc_dist_raster(krakow_airport, krakow_boundaries, 250, 32634)
writeRaster(krakow_airport_dist, "created/distances/krakow_airport_250m.tif", overwrite=T)
krakow_center_dist <- calc_dist_raster(krakow_center, krakow_boundaries, 250, 32634)
writeRaster(krakow_center_dist, "created/distances/krakow_center_250m.tif", overwrite=T)
