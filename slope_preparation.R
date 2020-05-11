rm(list=ls())

library(raster)
library(sp)

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

slope_ESP <- raster("original/slope/EUD_CP-SLOP_2500015000-AA_ESP/EUD_CP-SLOP_2500015000-AA.tif")
slope_GER <- raster("original/slope/EUD_CP-SLOP_4500035000-AA_GER/EUD_CP-SLOP_4500035000-AA.tif")
slope_POL <- raster("original/slope/EUD_CP-SLOP_5500035000-AA_POL/EUD_CP-SLOP_5500035000-AA.tif")


DNtoPercentage <- function(DN) {
  rad <- (acos(DN/250))
  percentage <- (tan(rad)*100)
  
  return(percentage)
}

