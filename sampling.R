rm(list=ls())

library(raster)
library(sp)
library(dplyr)
library(spdep)
library(lctools)

set.seed(21)


# source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/preparation_functions.R")
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

# load grids
stack_dresden <- stack("created/stack/dresden.grd")
stack_krakow <- stack("created/stack/krakow.grd")
stack_sevilla <- stack("created/stack/sevilla.grd")


### TODO: save pixel number
choose_samples <- function(raster_stack, sample_rate) {
  # convert stack to data frame
  df <- as.data.frame(raster_stack, optional = T, xy = T)
  df_noNA <- df[which(df$change!="NA"), ]
  
  # count changed records
  change_count <- count(df_noNA[which(df_noNA$change=="1"),])
  sample_amount <- as.integer(change_count/sample_rate)
  
  # samplechange records
  df_change <- df_noNA[which(df_noNA$change=="1"),]
  change_samples <- df_change[sample(nrow(df_change), sample_amount),]
  
  
  # sample no change records
  df_nochange <- df_noNA[which(df_noNA$change=="0"),]
  nochange_samples <- df_nochange[sample(nrow(df_nochange), sample_amount),]

  
  # join data frames
  samples <- rbind(change_samples, nochange_samples)
  
  return(samples)
  
}


dresden_samples <- choose_samples(stack_dresden, 25)
write.csv(dresden_samples, "created/samples/dresden_samples.csv")
krakow_samples <- choose_samples(stack_krakow, 25)
write.csv(krakow_samples, "created/samples/krakow_samples.csv")
sevilla_samples <- choose_samples(stack_sevilla, 25)
write.csv(sevilla_samples, "created/samples/sevilla_samples.csv")


# sample records with spatial sampling first
choose_samples_spatial <- function(raster_stack, sample_rate) {
  # count not NA records
  notNA_count <- freq(raster_stack$change)[1,2]+freq(raster_stack$change)[2,2]
  sample_size <- as.integer(notNA_count/sample_rate)
  
  # spatial sampling
  samples <- as.data.frame(sampleRandom(raster_stack, size = sample_size, na.rm = TRUE, cells = TRUE, rowcol = TRUE, xy = TRUE))
  
  #
  change_count <- as.integer(count(samples[which(samples$change=="1"),]))
  
  # sample no change records
  noChange_samples <- sample_n(samples[which(samples$change=="0"),], size=change_count)
  
  # join data frames
  samples <- rbind(samples[which(samples$change=="1"),], noChange_samples)
  
  return(samples)
}


dresden_samples_sp <- choose_samples_spatial(stack_dresden, 100)
write.csv(dresden_samples_sp, "created/samples/dresden_samples_sp.csv")
krakow_samples_sp <- choose_samples_spatial(stack_krakow, 100)
write.csv(krakow_samples_sp, "created/samples/krakow_samples_sp.csv")
sevilla_samples_sp <- choose_samples_spatial(stack_sevilla, 100)
write.csv(sevilla_samples_sp, "created/samples/sevilla_samples_sp.csv")


########################################################################################################
# Moran's I
########################################################################################################


calc_moransI <- function(samples) {
  # get coordinates as matrix
  coords <- as.matrix(cbind(samples$x, samples$y))
  
  # identify neighbours
  k1 <- knn2nb(knearneigh(coords))
  k1dists <- unlist(nbdists(k1, coords))
  # summary(k1dists)
  nb_maxdist <- dnearneigh(coords, 0, max(k1dists))
  print(nb_maxdist)
  
  # get neighbour list
  lw <- nb2listw(nb_maxdist,zero.policy = T)
  
  # calculate Morans I
  morans <- moran.test(samples$change, lw, randomisation = F, alternative = "two.sided", zero.policy = T)
  
  return(morans)
}



seed <- 22
set.seed(seed)
samples_sp <- choose_samples_spatial(stack_sevilla, 25)
calc_moransI(samples_sp)
set.seed(seed)
samples <- choose_samples(stack_sevilla, 25)
calc_moransI(samples)


# without sampling
df <- as.data.frame(stack_krakow, optional = T, xy = T)
df_noNA <- df[which(df$change!="NA"), ]

calc_moransI(df_noNA)

# ######################################################################################################
# # test stuff
# ######################################################################################################



# samples <- krakow_samples_sp
# # samples <- 
# 
# # inverse distace matrix
# samples.dists <- as.matrix(dist(cbind(samples$x, samples$y)))
# samples.dists.inv <- 1/(samples.dists/1000)
# diag(samples.dists.inv) <- 0
# 
# samples.dists.inv[1:10,1:10]
# 
# Moran.I(samples$change, samples.dists.inv, alternative="two.sided")
# 
# moransI(cbind(samples$x/1000, samples$y/1000), 5, samples$change, WType = "Binary")
# 
# 
# coords <- as.matrix(cbind(samples$x, samples$y))
# 
# k1 <- knn2nb(knearneigh(coords))
# k1dists <- unlist(nbdists(k1, coords))
# summary(k1dists)
# 
# nb_maxdist <- dnearneigh(coords, 0, max(k1dists))
# print(nb_maxdist)
# 
# lw <- nb2listw(nb_maxdist)
# 
# moran.test(samples$landuse, lw, randomisation = F, alternative = "greater")


# noChange_samples <- sample_n(df_noNA[which(df_noNA$change=="0"),], size=sample_amount)
# above does not keep the original cell number as row number

# change_samples <- sample_n(df_noNA[which(df_noNA$change=="1"),], size=sample_amount)
# above does not keep the original cell number as row number

# 
# 
# ### load grids
# stack_dresden <- stack("created/stack/dresden.grd")
# dresden_change <- raster("created/stack/dresden.grd")
# 
# plot(dresden)
# spplot(dresden)
# spplot(dresden_change)
# 
# # samples <- spsample(dresden, n=1000, "regular")
# 
# samples <- sampleRandom(dresden, size=7000, 
#                         na.rm = TRUE, 
#                         cells=TRUE, xy=TRUE, asRaster=TRUE)
# samples
# freq(samples$change)
# 
# ### spatial sampling
# samples_change <- sampleRandom(dresden_change, size=12000,asRaster=TRUE)
# writeRaster(samples_change, "created/test/samples.grd", overwrite=T)
# freq(samples_change)
# 
# 
# 
# # idea: 1. sample, 2. make data.frame, 3. sample not built-up to same amount as built-up
# #       4. merge data.frame and create raster from it
# 
# 
# 
# stack_dresden.df <- as.data.frame(stack_dresden, optional = T, xy = T)
# stack_dresden.df_noNA <- stack_dresden.df[which(stack_dresden.df$change!="NA"), ]
# head(stack_dresden.df_noNA)
# 
# stack_dresden.df_noNA[2,0]
# 
# ### data.frame sampling
# stack_dresden.samples.change <- sample_n(stack_dresden.df_noNA[which(stack_dresden.df_noNA$change=="1"),], size =1000)
# stack_dresden.samples.noChange <- sample_n(stack_dresden.df_noNA[which(stack_dresden.df_noNA$change=="0"),], size =1000)
# stack_dresden.samples.noChange
# 
# 
# count(dresden_samples_sp[which(dresden_samples_sp$change=="1"),])
# freq(dresden_samples_sp$change)
# count(dresden_samples_sp[which(dresden_samples_sp$change=="1"),])
