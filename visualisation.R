# rm(list=ls())
# source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/regression.R")

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/stratified/")

data.k <- read.csv("krakow.csv")
data.d <- read.csv("dresden.csv")
data.s <- read.csv("sevilla.csv")

data.k2 <- read.csv("krakow_2.csv")
data.d2 <- read.csv("dresden_2.csv")
data.s2 <- read.csv("sevilla_2.csv")

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/")
data.k.a <- read.csv("krakow_all.csv")
data.d.a <- read.csv("dresden_all.csv")
data.s.a <- read.csv("sevilla_all.csv")

trained.s <- ffs_model(data.s)
model.s.ffs <- trained.s$ffs$finalModel
summary(model.s.ffs)
calc_roc(model.s.ffs, test_data = split_landuse(data.s.a))
logit.plot.quad(model.s.ffs)


trained.d <- ffs_model(data.d)
model.d.ffs <- trained.d$ffs$finalModel
summary(model.d.ffs)
calc_roc(model.d.ffs, test_data = split_landuse(data.d.a))
logit.plot.quad(model.d.ffs)


trained.k <- ffs_model(data.k)
model.k.ffs <- trained.k$ffs$finalModel
summary(model.k.ffs)$coefficients[,1]
calc_roc(model.k.ffs)
logit.plot.quad(model.k.ffs)


# create table with all coefficients and the cities
# for every city show the value of the predictor if its part of the model
compare_coefs <- function(model1, model2, model3) {
  
  
  predictors <- c("built_dens", "pop_dens", "slope", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist",
                  "landuse.1", "landuse.2", "landuse.3", "landuse.4", "landuse.6")
  
  coefs <- matrix(, nrow = length(predictors), ncol = 3)
  
  rownames(coefs) <- predictors
  
  # models als liste übergeben um darüber zu iterieren?
  for (i in models) {
    coefs <- 
  }
  
  as.table(coefs)
  
}


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
