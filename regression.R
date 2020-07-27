rm(list=ls())

library(pROC)
library(caret)
library(varhandle)
library(spdep) # morans


############################################################################################################################
###################### functions #######################################################################################
############################################################################################################################

create_model <- function(data, significance = 0.05, split = F, train_part = 0.5, select = T) {
  
  data$landuse <- as.factor(data$landuse)

  # split data in train and test part
  if(split == T)   {
    trainIds <- createDataPartition(data$X, p = train_part, list = FALSE)
    data.train <- data[trainIds,]
    data.test <- data[-trainIds,]
  } else {
    data.train <- data
  }
  
  
  # create model with all variables
  model <- glm(formula = change ~ landuse + built_dens + pop_dens + slope + 
                 mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist,
                 family = "binomial",
                 data = data.train)
  
  
  # 
  formula <- "change ~ "
  
  if(TRUE %in% (summary(model)$coefficients[1:5,4] < significance)) {
    formula <- paste(formula, "landuse")
  }
  
  # find signifcant variables and add them to the formula for the new model
  for (i in c(6:length(summary(model)$coefficients[,4]))) {
    if(summary(model)$coefficients[i,4] < significance) {
      formula <- paste(formula, "+", rownames(summary(model)$coefficients)[i])
    }
  }
  
  
  # model with significant variables
  model.s <- glm(formula = formula,
                 family = "binomial",
                 data = data.train)
  
  if(select == F) {
    model.s <- model
  }
  
  # calculate roc of new model
  model.s.roc <- calc_roc(model.s)
                          #, test_data = data.test)
  print(model.s.roc)
  
  return(model.s)
}



#############################################################
# forward feature selection
#############################################################

library(CAST)

ffs_model <- function(data, significance = 0.1) {

  # split landuse 
  data.dummy <- data
  data.dummy$landuse <- as.factor(data.dummy$landuse)

  # create folds for cross validation
  cv_strata_l <- length(unique((data.dummy[,"cv_strata"])))
  indices <- CreateSpacetimeFolds(data.dummy, spacevar = "cv_strata", k = cv_strata_l)
  ctrl <- trainControl(method = "cv", index = indices$index)
  
  # check for predictors (because landuse 3 is not existent in sevilla data)
  predictors <- c("landuse", "built_dens", "pop_dens", "slope", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist")
  
  
  # initialize trained object
  trained <- deparse(substitute(data))
  
  # glm with all factors
  trained$glm <- create_model(data.dummy, split = F, select = F)

  # glm with significant factors  
  trained$glm$sig <- create_model(data.dummy, split = F, significance = significance)
  
  trained$cv_strata <- cv_strata_l
  
  # use forward feature selection to select best combination of determinants
  trained$ffs <- ffs(data.dummy[,predictors], data.dummy[,"change"],
                   method = "glm",
                   family = "binomial", 
                   trControl = ctrl)
  
  # create a model with cv containing all determinants
  trained$train <- train(data.dummy[,predictors], data.dummy[,"change"],
                         method = "glm",
                         family = "binomial",
                         trControl = ctrl)
  
  # if a ffs model doesnt contain landuse, calculate a model with landuse
  ffs_predictors <- names(trained$ffs$finalModel$coefficients)
  if(!("landuse2" %in% ffs_predictors)) {
    ffs_predictors <- append(ffs_predictors, "landuse")
    ffs_predictors <- ffs_predictors[-1]
    trained$ffs_plus <- train(data.dummy[,ffs_predictors], data.dummy[,"change"],
                           method = "glm",
                           family = "binomial",
                           trControl = ctrl)
  } else {
    trained$ffs_plus <- trained$ffs
  }
   
  
  return(trained)
  
}


bss_model <- function (data, city_name) {
  data.dummy <- data
  data.dummy$landuse <- as.factor(data.dummy$landuse)
  
  # create folds for cross validation
  cv_strata_l <- length(unique((data.dummy[,"cv_strata"])))
  indices <- CreateSpacetimeFolds(data.dummy, spacevar = "cv_strata", k = cv_strata_l)
  ctrl <- trainControl(method = "cv", index = indices$index)
  
  # check for predictors (because landuse 3 is not existent in sevilla data)
  predictors <- c("landuse", "built_dens", "pop_dens", "slope", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist")
  
  
  
  trained_bss <- bss(data.dummy[,predictors], data.dummy[,"change"],
                              method = "glm",
                              family = "binomial", 
                              trControl = ctrl)
  
  setwd("C:/Users/janst/sciebo/Bachelor Thesis/results/models_RDS/")
  saveRDS(trained_bss, paste("trained_bss_", city_name, ".rds", sep = ""))

  bss_predictors <- names(trained_bss$finalModel$coefficients)
  if(!("landuse2" %in% bss_predictors)) {
    bss_predictors <- append(bss_predictors, "landuse")
    bss_predictors <- bss_predictors[-1]
    trained_bss_plus <- train(data.dummy[,bss_predictors], data.dummy[,"change"],
                              method = "glm",
                              family = "binomial",
                              trControl = ctrl)
  } else {
    trained_bss_plus <- trained_bss
  }
  
  saveRDS(trained_bss_plus, paste("trained_bss_plus", city_name, ".rds", sep = ""))
  
}


############################################################################################################################
###################### ROC functions #######################################################################################
############################################################################################################################

# see: https://stackoverflow.com/questions/18449013/r-logistic-regression-area-under-curve
calc_roc <- function(model, test_data = NULL, main = "") {
  if(is.null(test_data)){
    prob = predict(model, type = c("response"))
    model$data$prob = prob
    if (is.null(model$model$.outcome)) {
      g <- roc(model$model$change ~ prob, data = model$data)
    } else {
      g <- roc(model$model$.outcome ~ prob, data = model$data)
    }
    
    print("model data")
  } else {
    test_data$landuse <- as.factor(test_data$landuse)
    prob = predict(model, newdata = test_data, type = c("response"))
    # test_data$prob = prob
    g <- roc(test_data$change ~ prob, data = test_data)
    print("test data")
  }
  
  plot.roc(g, print.auc = TRUE, print.auc.x = 0.3, print.auc.y = 0, print.auc.cex = 1.5, main = main, cex.lab =1.5, cex.axis = 1.5)
  
  return(g$auc)
}




########################################################################################
########################################################################################
### the following code was at the end not used for the thesis
########################################################################################
########################################################################################


split_landuse <- function(data) {
  if("cv_strata" %in% colnames(data)) {
    predictors <- c("change", "built_dens", "pop_dens", "slope", "landuse", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist", "cv_strata")
  } else {
    predictors <- c("change", "built_dens", "pop_dens", "slope", "landuse", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist")
  }
  
  data.dummy <- cbind(data[,predictors], to.dummy(data$landuse, prefix = "landuse"))
  data.dummy["landuse"] <- NULL
  
  return(data.dummy)
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
