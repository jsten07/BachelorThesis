rm(list=ls())

library(raster)
#library(sp)
#library(dplyr)
library(spdep) # morans
#library(lctools)
#library(spcosa)
library(splitstackshape) # stratified

set.seed(21)


# source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/preparation_functions.R")
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

# load grids
stack_dresden <- stack("created/stack/dresden.grd")
stack_krakow <- stack("created/stack/krakow.grd")
stack_sevilla <- stack("created/stack/sevilla.grd")



###################################################################################################################
# 820 - 878 samples
# Morans I: 0.25 - 0.29 / 0.26 - 0.40 (k = 8)
# AUC: 0.86 - 0.89
###################################################################################################################
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")
set.seed(13)
write.csv((samples_dresden <- stratified_sampling(stack_dresden, 5)), ("created/samples/stratified/dresden.csv"))
write.csv((samples_sevilla <- stratified_sampling(stack_sevilla, 6)), ("created/samples/stratified/sevilla.csv"))
write.csv((samples_krakow <- stratified_sampling(stack_krakow, 9)), ("created/samples/stratified/krakow.csv"))




###################################################################################################################
# 860 - 2728 samples
# Morans I: 0.28.- 0.37
# ROC: 0.84-0.87
###################################################################################################################
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")
set.seed(13)
write.csv((samples_dresden <- stratified_sampling(stack_dresden, 5)), ("created/samples/stratified/dresden.csv"))
write.csv((samples_sevilla <- stratified_sampling(stack_sevilla, 5)), ("created/samples/stratified/sevilla.csv"))
write.csv((samples_krakow <- stratified_sampling(stack_krakow, 5)), ("created/samples/stratified/krakow.csv"))



###################################################################################################################
# 464 - 494
# 0.24 - 0.37 (k = 8)
###################################################################################################################
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")
set.seed(13)
write.csv((samples_dresden <- stratified_sampling(stack_dresden, 7)), ("created/samples/stratified/dresden.csv"))
write.csv((samples_sevilla <- stratified_sampling(stack_sevilla, 8)), ("created/samples/stratified/sevilla.csv"))
write.csv((samples_krakow <- stratified_sampling(stack_krakow, 12)), ("created/samples/stratified/krakow.csv"))



###################################################################################################################
# stratified 
# with changed and not changed value for every strata
# 950 - 994 samples (window size: 24, 16, ?)
# MOrans I: 0.04 - -0.17
# AUC: 0.67 - 0.76
###################################################################################################################
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")
set.seed(13)
write.csv((samples_dresden <- strata_sampling(stack_dresden, 24)), ("created/samples/stratified/dresden_new.csv"))
write.csv((samples_sevilla <- strata_sampling(stack_sevilla, 16)), ("created/samples/stratified/sevilla_new.csv"))
write.csv((samples_krakow <- strata_sampling(stack_krakow, 28)), ("created/samples/stratified/krakow_new.csv"))


###################################################################################################################
# stratified 
# with changed and not changed value for every strata
# 2370 - 2540 samples (window size: 13, 9, 16)
# MOrans I: 0.11 - -0.06
# 0.65 - 0.72
###################################################################################################################
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")
set.seed(13)
write.csv((samples_dresden <- strata_sampling(stack_dresden, 13)), ("created/samples/stratified/dresden_new.csv"))
write.csv((samples_sevilla <- strata_sampling(stack_sevilla, 9)), ("created/samples/stratified/sevilla_new.csv"))
write.csv((samples_krakow <- strata_sampling(stack_krakow, 16)), ("created/samples/stratified/krakow_new.csv"))




str(samples_dresden)
str(samples_sevilla)
str(samples_krakow)

system.time(calc_moransI(samples_sevilla_all))

calc_moransI(samples_dresden, dist = 1115)
calc_moransI(samples_sevilla, dist = 1325)
calc_moransI(samples_krakow, dist = 1332)
calc_moransI(samples_dresden, k = 8)
calc_moransI(samples_sevilla, k = 8)
calc_moransI(samples_krakow, k = 8)

calc_moransI(data.d, k = 8)
calc_moransI(data.s, k = 8)
calc_moransI(data.k, k = 8)


calc_moransI(samples_dresden_all, k = 8)
calc_moransI(samples_sevilla_all, k = 8)
calc_moransI(samples_krakow_all, k = 8)

# process in parallel
library(doParallel) 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

# turn parallel processing off and run sequentially again:
registerDoSEQ()

calc_moransI(samples_krakow_all)

########################################################################################
######################### functions ####################################################
########################################################################################



### convert grid / stack to data.frame and remove NA values
create_df_without_NA <- function(stack) {
  df <- as.data.frame(stack, optional = T, xy = T)
  df_noNA <- df[which(df$change!="NA"), ]
  
  return(df_noNA)
}

write.csv(samples_dresden_all <- create_df_without_NA(stack_dresden), "created/samples/dresden_all.csv")
write.csv(samples_sevilla_all <- create_df_without_NA(stack_sevilla), "created/samples/sevilla_all.csv")
write.csv(samples_krakow_all <- create_df_without_NA(stack_krakow), "created/samples/krakow_all.csv")


########################################################################################
# random sampling
########################################################################################

### randomly choose 1/sample_rate of changed and not changed data as samples and join them to one dataset
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


# dresden_samples <- choose_samples(stack_dresden, 50)
# write.csv(dresden_samples, "created/samples/dresden_samples.csv")
# krakow_samples <- choose_samples(stack_krakow, 50)
# write.csv(krakow_samples, "created/samples/krakow_samples.csv")
# sevilla_samples <- choose_samples(stack_sevilla, 50)
# write.csv(sevilla_samples, "created/samples/sevilla_samples.csv")


### first sample all records with sampleRandom (spatially), then choose same amount of not changed as was generated for changed
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


# dresden_samples_sp <- choose_samples_spatial(stack_dresden, 100)
# write.csv(dresden_samples_sp, "created/samples/dresden_samples_sp.csv")
# krakow_samples_sp <- choose_samples_spatial(stack_krakow, 100)
# write.csv(krakow_samples_sp, "created/samples/krakow_samples_sp.csv")
# sevilla_samples_sp <- choose_samples_spatial(stack_sevilla, 100)
# write.csv(sevilla_samples_sp, "created/samples/sevilla_samples_sp.csv")



########################################################################################################
# stratified sampling
########################################################################################################

# write.csv(sampleStratified(stack_dresden$change, size = 500, xy = TRUE), "created/samples/test_stratified.csv")

# systematic sampling
# choose every n-th record 
sampleSystematic <- function(stack, window_size) {
  
  df <- as.data.frame(stack, optional = T, xy = T)
  
  # initialize sample data frame
  samples <- df[1,]
  
  # get every n*n-th record (n = window_size)
  row_iterations <- c(1:round(nrow(stack)/window_size))
  column_iterations <- c(1:round(ncol(stack)/window_size))
  
  for (i in row_iterations) {
    # start every row with n-th sample
    s <- sample(1:window_size, 1)
    for (j in column_iterations){
      # choose cell number
      cell <- i * (window_size) * ncol(stack) + j * window_size + s
      samples[nrow(samples)+1,] <- df[cell,]
    }
  }
  
  # remove NA records
  samples <- samples[which(samples$change!="NA"), ]
  
  # remove excessive not changed samples
  changed <- samples[which(samples$change=="1"), ]
  change_amount <- nrow(changed)
  not_changed <- samples[which(samples$change=="0"), ]
  not_changed <- not_changed[sample(nrow(not_changed), change_amount),]
  
  samples <- rbind(changed, not_changed)
  
  return(samples)
}

# sevilla_syst <- sampleSystematic(stack_sevilla, 7)
# write.csv(sevilla_syst, "created/samples/sevilla_syst.csv")
# 
# calc_moransI(sevilla_syst)



# just take every nth record
# choose sample way: rate(choose every n record) or amount(choose amount of samples)
sampleEasy <- function(stack, way = "rate", size = 25) {
  df <- as.data.frame(stack, xy = T)
  df_noNA <- df[which(df$change!="NA"), ]
  
  df.change <- df[which(df$change=="1"), ]
  df.nochange <- df[which(df$change=="0"), ]
  
  if(way == "amount") {
    size <- size/2
    size <- nrow(df.change)/size
  }
  
  s <- sample(1:size, 1)
  df.change.new <- df.change[seq(s, nrow(df.change), size), ]
  df.nochange.new <- df.nochange[seq(s, nrow(df.nochange), nrow(df.nochange)/nrow(df.change.new)), ]
  
  df.new <- rbind(df.change.new, df.nochange.new)
  
  return(df.new)
}
# 
# write.csv((dresden_samples <- sampleEasy(stack_dresden, way="amount", size = 1000)), "created/samples/dresden_easySamples.csv")
# write.csv((krakow_samples <- sampleEasy(stack_krakow, way="amount", size = 1000)), "created/samples/krakow_easySamples.csv")
# write.csv((sevilla_samples <- sampleEasy(stack_sevilla, way="amount", size = 500)), "created/samples/sevilla_easySamples.csv")
# calc_moransI(sevilla_samples)



##############################################################################
# create strata as square windows
# choose one record in every strata
# remove NAs and bring changed and not changed to same amount

stratified_sampling <- function(stack, window_size = 5) {
  # width of raster (amount of columns)
  raster_cols <- ncol(stack)
  
  df <- as.data.frame(stack, xy = T)
  df$landuse <- as.factor(df$landuse)
  cell_no <- as.integer(row.names(df))
  
  # calc strata numbers 
  strata <- calc_strata_no(cell_no, window_size, raster_cols)
  df["strata"] <- strata     # and add to data frame
  strata.cv <- calc_strata_no(cell_no, 200, raster_cols)
  df["cv_strata"] <- strata.cv
  
  # sample dataframe stratified
  df.strat <- stratified(df, c("strata"), size = 1, keep.rownames = TRUE)
  
  
  # remove NA values
  samples <- df.strat[which(df.strat$change!="NA"), ]
  sampleschange <- samples[which(samples$change=="1"), ]
  samplesNochange <- samples[which(samples$change=="0"), ]
  samplesNochange <- samplesNochange[sample(nrow(samplesNochange), nrow(sampleschange)),]

  samples <- rbind(sampleschange, samplesNochange)

  # !!! general problem: samples in strata with few records overrepresented
  
  ######################################################
  # sampleStratified() just works for raster data
  ######################################################
  # possibility to use stratified() method in fifer package
  ######################################################
  
  return(samples)
  
} 


##################################################################################
# create square strata in window_size
# choose one changed and one not changed record from every strata (if available)
# randomly choose not changed samples to get same amount as of changed

strata_sampling <- function(stack, window_size = 5) {
  # width of raster (amount of columns)
  raster_cols <- ncol(stack)
  
  df <- as.data.frame(stack, xy = T)
  df$landuse <- as.factor(df$landuse)
  cell_no <- as.integer(row.names(df))
  
  # calc strata numbers 
  strata <- calc_strata_no(cell_no, window_size, raster_cols)
  df["strata"] <- strata     # and add to data frame
  strata.cv <- calc_strata_no(cell_no, 200, raster_cols)
  df["cv_strata"] <- strata.cv
  
  # sample dataframe stratified
  df.strat.1 <- stratified(df, c("change", "strata"), size = 1, select = list(change = "1"), keep.rownames = TRUE)
  df.strat.0 <- stratified(df, c("change", "strata"), size = 1, select = list(change = "0"), keep.rownames = TRUE)
  df.strat.0 <- df.strat.0[sample(nrow(df.strat.0), nrow(df.strat.1)),]
  
  samples <- rbind(df.strat.1, df.strat.0)
  
  return(samples)
  
} 

# samples <- stratified_sampling(stack_dresden, window_size = 5)
# write.csv(samples, "created/samples/sevilla_strat_test.csv")

# cell_vector: vector with numbers, representing cell numbers of a raster
# window_size: divide raster in n x n windows; n = window_size
# raster_cols: width of the raster
# return: strata number for every raster cell
calc_strata_no <- function(cell_vector, window_size, raster_cols){
  cells <- (ceiling((cell_vector%%raster_cols)/window_size)
            + ceiling((cell_vector/raster_cols)/window_size) * ceiling(raster_cols/window_size))
  
  return(cells)
}




########################################################################################################
# Moran's I
########################################################################################################


calc_moransI <- function(samples, dist = NULL, k = NULL) {
  
  # get coordinates as matrix
  coords <- as.matrix(cbind(samples$x, samples$y))
  
  if(!is.null(dist)) {
    
    nb <- dnearneigh(coords, 0, dist)
    
  } else if(!is.null(k)) {
    
    nb <- knn2nb(knearneigh(coords, k = k))
    
  } else {
    
    # identify neighbours
    # within the minimum distance so every sample has at least one neighbour
    k1 <- knn2nb(knearneigh(coords))
    k1dists <- unlist(nbdists(k1, coords))
    # summary(k1dists)
    nb <- dnearneigh(coords, 0, max(k1dists))
    print("max dist:")
    print(max(k1dists))
    
  }
  
  print(nb)
  
  # get neighbour list
  lw <- nb2listw(nb,zero.policy = T)
  
  # lw_k1 <- nb2listw(k1, zero.policy = T)
  
  # calculate Morans I
  morans <- moran.test(samples$change, lw, randomisation = F, alternative = "two.sided", zero.policy = T)
  
  return(morans)
}


# seed <- 22
# set.seed(seed)
# samples_sp <- choose_samples_spatial(stack_sevilla, 25)
# calc_moransI(samples_sp)
# set.seed(seed)
# samples <- choose_samples(stack_sevilla, 25)
# calc_moransI(samples)
# 
# 
# # without sampling
# df <- as.data.frame(stack_krakow, optional = T, xy = T)
# df_noNA <- df[which(df$change!="NA"), ]
# 
# calc_moransI(df_noNA)

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
