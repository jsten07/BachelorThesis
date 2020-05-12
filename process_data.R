rm(list=ls())

### load functions
source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/preparation_functions.R")
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

### load data

sevilla_boundaries <- shapefile("created/GADM/Sevilla_boundaries.shp")
krakow_boundaries <- shapefile("created/GADM/Krakow_boundaries.shp")
dresden_boundaries <- shapefile("created/GADM/Dresden_boundaries.shp")

GHSL_ESP_2014 <- raster("original/GHSL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
GHSL_ESP_1990 <- raster("original/GHSL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4_ESP/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_17_4.tif")
GHSL_POL_2014 <- raster("original/GHSL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_19_3_POL/GHS_BUILT_LDS2014_GLOBE_R2018A_54009_250_V2_0_19_3.tif")
GHSL_POL_1990 <- raster("original/GHSL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_19_3_POL/GHS_BUILT_LDS1990_GLOBE_R2018A_54009_250_V2_0_19_3.tif")
GHSL_ESP_30m <- raster("original/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_13_8_ESP/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_13_8.tif")
GHSL_POL_30m <- raster("original/GHSL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_14_7_POL/GHS_BUILT_LDSMT_GLOBE_R2018A_3857_30_V2_0_14_7.tif")

GHSL_pop_ESP <- raster("original/GHSL/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_17_4_ESP/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_17_4.tif")
GHSL_pop_POL <- raster("original/GHSL/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_19_3_POL/GHS_POP_E2000_GLOBE_R2019A_54009_250_V1_0_19_3.tif")

slope_ESP <- raster("original/slope/EUD_CP-SLOP_2500015000-AA_ESP/EUD_CP-SLOP_2500015000-AA.tif")
slope_GER <- raster("original/slope/EUD_CP-SLOP_4500035000-AA_GER/EUD_CP-SLOP_4500035000-AA.tif")
slope_POL <- raster("original/slope/EUD_CP-SLOP_5500035000-AA_POL/EUD_CP-SLOP_5500035000-AA.tif")

landuse <- raster("original/land cover/clc2006_clc2000_v2018_20_raster100m/CLC2006_CLC2000_V2018_20.tif")


### process and save data

change_Sevilla <- getChange(GHSL_ESP_1990, GHSL_ESP_2014, sevilla_boundaries, epsg = 32630, resolution = 250, threshold = 50)
writeRaster(change_Sevilla, "created/GHSL_R/sevilla_change.tif", overwrite=T)
change_Krakow <- getChange(GHSL_POL_1990, GHSL_POL_2014, krakow_boundaries, epsg = 32634, resolution = 250, threshold = 50)
writeRaster(change_Krakow, "created/GHSL_R/krakow_change.tif", overwrite=T)
change_Dresden <- getChange(GHSL_POL_1990, GHSL_POL_2014, dresden_boundaries, epsg = 32633, resolution = 250, threshold = 50)
writeRaster(change_Dresden, "created/GHSL_R/dresden_change.tif", overwrite=T)

change_Sevilla_30m <- getChangeFromMultitemp(GHSL_ESP_30m, sevilla_boundaries, 32630, resolution = 30)
writeRaster(change_Sevilla_30m, "created/GHSL_R/sevilla_change_30_m.tif", overwrite=T)
change_Krakow_30m <- getChangeFromMultitemp(GHSL_POL_30m, krakow_boundaries, 32634, resolution = 30)
writeRaster(change_Krakow_30m, "created/GHSL_R/krakow_change_30_m.tif", overwrite=T)
change_Dresden_30m <- getChangeFromMultitemp(GHSL_POL_30m, dresden_boundaries, 32633, resolution = 30)
writeRaster(change_Dresden_30m, "created/GHSL_R/dresden_change_30_m.tif", overwrite=T)

population_Sevilla <- reprojectAndCrop(GHSL_pop_ESP, sevilla_boundaries, 32630, 250)
writeRaster(population_Sevilla, "created/GHSL_R/sevilla_popDens.tif", overwrite=T)
population_Krakow <- reprojectAndCrop(GHSL_pop_POL, krakow_boundaries, 32634, 250)
writeRaster(population_Krakow, "created/GHSL_R/krakow_popDens.tif", overwrite=T)
population_Dresden <- reprojectAndCrop(GHSL_pop_POL, dresden_boundaries, 32633, 250)
writeRaster(population_Dresden, "created/GHSL_R/dresden_popDens.tif", overwrite=T)

slope_Sevilla <- citySlopeAsPercentage(slope_ESP, sevilla_boundaries, 32630)
writeRaster(slope_Sevilla, "created/slope/sevilla_slope.tif", overwrite=T)
slope_Dresden <- citySlopeAsPercentage(slope_GER, dresden_boundaries, 32633)
writeRaster(slope_Dresden, "created/slope/dresden_slope.tif", overwrite=T)
slope_Krakow <- citySlopeAsPercentage(slope_POL, krakow_boundaries, 32634)
writeRaster(slope_Krakow, "created/slope/krakow_slope.tif", overwrite=T)


landuse_Dresden <- cropAndReclassify_landuse(landuse, dresden_boundaries, 32633, 100)
writeRaster(landuse_Dresden, "created/landuse/dresden_landuse.tif", overwrite=T)
landuse_Sevilla <- cropAndReclassify_landuse(landuse, sevilla_boundaries, 32630, 100)
writeRaster(landuse_Sevilla, "created/landuse/sevilla_landuse.tif", overwrite=T)
landuse_Krakow <- cropAndReclassify_landuse(landuse, krakow_boundaries, 32634, 100)
writeRaster(landuse_Krakow, "created/landuse/krakow_landuse.tif", overwrite=T)

