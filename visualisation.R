rm(list=ls())


#########################################################################################
# results
#########################################################################################

source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/regression.R")

library(formattable)
 
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/stratified/")

data.k <- read.csv("krakow.csv")
data.d <- read.csv("dresden.csv")
data.s <- read.csv("sevilla.csv")

data.k2 <- read.csv("krakow_lu90.csv")
data.d2 <- read.csv("dresden_lu90.csv")
data.s2 <- read.csv("sevilla_lu90.csv")

data.k.i <- read.csv("krakow_i900.csv")
data.d.i <- read.csv("dresden_i900.csv")
data.s.i <- read.csv("sevilla_i900.csv")

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/")
data.k.a <- read.csv("krakow_all.csv")
data.d.a <- read.csv("dresden_all.csv")
data.s.a <- read.csv("sevilla_all.csv")

trained.s <- ffs_model(data.s2)
model.s.ffs <- trained.s$ffs$finalModel
model.s.glm <- trained.s$glm$sig
model.s.glmA <- trained.s$glm
summary(model.s.ffs)
calc_roc(model.s.ffs)
logit.plot.quad(model.s.ffs)


trained.d <- ffs_model(data.d2)
model.d.ffs <- trained.d$ffs$finalModel
model.d.glm <- trained.d$glm$sig
model.d.glmA <- trained.d$glm
summary(model.d.ffs)
calc_roc(model.d.ffs)
logit.plot.quad(model.d.ffs)


trained.k <- ffs_model(data.k2)
model.k.ffs <- trained.k$ffs$finalModel
model.k.glm <- trained.k$glm$sig
model.k.glmA <- trained.k$glm
summary(model.k.ffs)$coefficients[,1]
calc_roc(model.k.ffs)
logit.plot.quad(model.k.ffs)


# create table with all coefficients and the cities
# for every city show the value of the predictor if its part of the model
compare_coefs <- function(model1, model2, model3) {
  
  
  predictors <- c("(Intercept)", "built_dens", "pop_dens", "slope", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist",
                  "landuse1", "landuse2", "landuse3", "landuse4", "landuse5", "landuse6", "aic", "roc", "Morans I")
  
  coefs <- matrix(, nrow = length(predictors), ncol = 6)
  
  rownames(coefs) <- predictors
  colnames(coefs) <- c(deparse(substitute(model1)), "p1", deparse(substitute(model2)), "p2", deparse(substitute(model3)), "p3")
  
  #coefs.sign <- coefs
  
  # models als liste übergeben um darüber zu iterieren?
  for (i in predictors) {
    if(i %in% names(model1$coefficients)) {
      coefs[i,1] <- model1$coefficients[i]
      coefs[i,2] <- summary(model1)$coefficients[i,4]
    }
    if(i %in% names(model2$coefficients)) {
      coefs[i,3] <- model2$coefficients[i]
      coefs[i,4] <- summary(model2)$coefficients[i,4]
    }
    if(i %in% names(model3$coefficients)) {
      coefs[i,5] <- model3$coefficients[i]
      coefs[i,6] <- summary(model3)$coefficients[i,4]
    }
  }
  
  coefs["aic", 1] <- model1$aic
  coefs["aic", 3] <- model2$aic
  coefs["aic", 5] <- model3$aic
  
  coefs["roc", 1] <- calc_roc(model1)
  coefs["roc", 3] <- calc_roc(model2)
  coefs["roc", 5] <- calc_roc(model3)
  
  message("Check data for Morans!")
  coefs["Morans I", 1] <- calc_moransI(data.d2)$estimate[1]
  coefs["Morans I", 3] <- calc_moransI(data.s2)$estimate[1]
  coefs["Morans I", 5] <- calc_moransI(data.k2)$estimate[1]
  
  options(scipen = 99, digits = 3)
  
  coefs <- as.data.frame(coefs)
  #coefs.sign <- as.data.frame(coefs.sign)
  
  print(formattable(coefs, format = "f", digits = 5))
  # print(formattable(coefs.sign, digits = 3, title = "p-value"))
  
  return(coefs)
}

setwd("C:/Users/janst/sciebo/Bachelor Thesis/results/models/")

ffs_models <- compare_coefs(model.d.ffs, model.s.ffs, model.k.ffs)
write.csv(ffs_models, "ffs_models_lu90.csv")
glm_models <- compare_coefs(model.d.glm, model.s.glm, model.k.glm)
write.csv(glm_models, "glm_models_lu90.csv")
glm_allV_models <- compare_coefs(trained.d$glm, trained.s$glm, trained.k$glm)
write.csv(glm_allV_models, "glm_allV_models_lu90.csv")
train_models <- compare_coefs(trained.d$train$finalModel, trained.s$train$finalModel, trained.k$train$finalModel)
write.csv(train_models, "train_models_lu90.csv")
train_models <- compare_coefs(trained.d$ffs_plus$finalModel, trained.s$ffs_plus$finalModel, trained.k$ffs_plus$finalModel)
write.csv(train_models, "ffs_plus_models_lu90.csv")






#########################################################################################
# data
#########################################################################################
# source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/preparation_functions.R")
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")

# load grids
stack_dresden <- stack("created/stack/dresden.grd")
stack_krakow <- stack("created/stack/krakow.grd")
stack_sevilla <- stack("created/stack/sevilla.grd")

# add units
# change and landuse as categorical

# dependet in seperate plot
# show original data to used data
plot(stack_dresden)
plot(stack_krakow)
plot(stack_sevilla)


############################################################################
# landuse
############################################################################

# show absolute and relative values for (no) change / landuse ratio



############################################################################
# report plot function

## plot logistic curve, threshold, T/F +/-, sensitivity, specificity
## arguments
## model a fitted glm
## threshold cutoff for sensitivity/specificity, default 0.5
## title (optional)
logit.plot.quad <- function(model, threshold=0.5, title="Model success") {
  sf<-sort(fitted(model), index=T)
  # leave extra space at bottom
  par(mar=c(6,4,4,2)+.1); par(xaxs="i", yaxs="r")
  plot(sf$x, ylim=c(0,1), type="l", col="blue", lwd=3, xlab="",
       ylab="probability of change")
  abline(h=c(0,1), lty=1)
  # show threshold and crossover point
  abline(h=threshold,lty=2); text(0,threshold+.02,
                                  paste("threshold =", threshold), pos=4)
  crossover <- sum(fitted(model) < threshold)
  abline(v=crossover,lty=2)
  text(crossover,.05,"crossover",pos=4)
  text(crossover, threshold-.03,
       "fitted probability of change",col="blue",pos=4)
  # name of the response field
  field.name <- attr(attr(terms(formula(model)), "factors"),
                     "dimnames")[[1]][1]
  # extract the T/F from it
  eval(parse(text=paste("tmp <- ",
                        ifelse(class(model$data) == "data.frame", "model$data$", ""),
                        field.name, sep="")))
  # show T/F as vertical bars at the index
  # colours differ with T/F predictions
  points(1:length(tmp),tmp[sf$ix],
         pch="|",cex=1,
         col=ifelse((tmp[sf$ix] == (sf$x>threshold)),"green4","red"))
  # compute proportions
  tn <- sum((!tmp[sf$ix]) & (sf$x < threshold))
  fn <- sum((!tmp[sf$ix]) & (sf$x >= threshold))
  tp <- sum(tmp[sf$ix] & (sf$x >= threshold))
  fp <- sum(tmp[sf$ix] & (sf$x < threshold))
  right <- length(sf$x)*.65
  text(0,.1,paste("True negatives:",tn), col="green4",pos=4)
  text(right,.1,paste("False positives:", fn), col="red",pos=4)
  text(right,.9,paste("True positives:", tp), col="green4",pos=4)
  text(0,.9,paste("False negatives:", fp), col="red",pos=4)
  title(main=title)
  title(sub=paste("Sensitivity:", round(tp/(tp+fp),4),
                  "; Specificity:", round(tn/(tn+fn),4)), line=4)
}





