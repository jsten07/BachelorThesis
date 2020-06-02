rm(list=ls())

# library(Deducer)
library(pROC)
library(caret)
library(varhandle)

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/stratified/")

# load data
data.k.n <- read.csv("krakow_new.csv")
trainIds <- createDataPartition(data.k.n$X, p = 0.5, list = FALSE)
data.k.train <- data.k.n[trainIds,]
data.k.test <- data.k.n[-trainIds,]

data.d.n <- read.csv("dresden_new.csv")
# trainIds <- createDataPartition(data.d$X, p = 0.5, list = FALSE)
# data.d.train <- data.d[trainIds,]
# data.d.test <- data.d[-trainIds,]

data.s.n <- read.csv("sevilla_new.csv")
# trainIds <- createDataPartition(data.d$X, p = 0.5, list = FALSE)
# data.s.train <- data.s[trainIds,]
# data.s.test <- data.s[-trainIds,]
# 
# data.k.a <- read.csv("krakow_all.csv")
# data.d.a <- read.csv("dresden_all.csv")
# data.s.a <- read.csv("sevilla_all.csv")

data.k <- read.csv("krakow.csv")
trainIds <- createDataPartition(data.k$X, p = 0.5, list = FALSE)
data.k.train <- data.k[trainIds,]
data.k.test <- data.k[-trainIds,]
calc_moransI(data.d)
data.d <- read.csv("dresden.csv")
data.s <- read.csv("sevilla.csv")

# data.k.test <- read.csv("krakow_test.csv")
##############################################################
# logistic regression
##############################################################

create_model(data.k)
create_model(data.d)
create_model(data.s)

create_model(data.k.n)
create_model(data.d.n)
create_model(data.s.n)




############################################################################################################################
###################### functions #######################################################################################
############################################################################################################################

create_model <- function(data, significance = 0.05, split = T, train_part = 0.5) {
  
  # split landuse in single binary variables
  # data.k.dummy <- cbind(data.k, to.dummy(data.k$landuse, prefix = "landuse"))
  
  # split data in train and test part
  if(split == T)   {
    trainIds <- createDataPartition(data$X, p = train_part, list = FALSE)
    data.train <- data[trainIds,]
    data.test <- data[-trainIds,]
  } else {
    data.train <- data
    data.test <- data
  }
  
  
  # create model with all variables
  model <- glm(formula = change ~ factor(landuse) + built_dens + pop_dens + slope + 
                 mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist,
               family = "binomial",
               data = data.train)
  
  # 
  formula <- "change ~ factor(landuse) "
  
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
  
  # calculate roc of new model
  model.s.roc <- calc_roc(model.s, test_data = data.test)
  print(model.s.roc)
  
  return(model.s)
}

############################################################################################################################
###################### ROC functions #######################################################################################
############################################################################################################################

# see: https://stackoverflow.com/questions/18449013/r-logistic-regression-area-under-curve
calc_roc <- function(model, test_data = NULL) {
  if(is.null(test_data)){
    prob = predict(model, type = c("response"))
    model$data$prob = prob
    g <- roc(model$model$.outcome ~ prob, data = model$data)
    print("model data")
  } else {
    prob = predict(model, newdata = test_data, type = c("response"))
    # test_data$prob = prob
    g <- roc(test_data$change ~ prob, data = test_data)
    print("test data")
  }
  
  plot.roc(g, print.auc = TRUE, print.auc.x = 0.3, print.auc.y = 0)
  
  return(g$auc)
}


##############################################################
##############################################################
# tests
##############################################################
##############################################################


#############################################################
# forward feature selection
library(CAST)

predictors <- c("change", "built_dens", "pop_dens", "slope", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist", "cv_strata")

# train model
data.dummy <- cbind(data.k.n[,predictors], to.dummy(data.k.n$landuse, prefix = "landuse"))
data.dummy["landuse"] <- NULL
trainIds <- createDataPartition(data.dummy$change, p = 0.5, list = FALSE)
data.k.train <- data.dummy[trainIds,]
data.k.test <- data.dummy[-trainIds,]
cv_strata_l <- length(unique((data.dummy[,"cv_strata"])))
indices <- CreateSpacetimeFolds(data.dummy, spacevar = "cv_strata", k = cv_strata_l)
ctrl <- trainControl(method = "cv", index = indices$index)

predictors <- c("built_dens", "pop_dens", "slope", "mRoads_dist", "pRoads_dist", "river_dist", "train_dist", "center_dist", "airport_dist",
                "landuse.1", "landuse.2", "landuse.3", "landuse.4", "landuse.6")

trained.k <- ffs(data.dummy[,predictors], data.dummy[,"change"],
                 method = "glm",
                 family = "binomial",
                 trControl = ctrl)
calc_roc(trained.k$finalModel, test_data = data.k.train)


trained.perfect.k <- trained.k

###########################################################
# manually 

glm.t.k <- glm(formula = change ~ factor(landuse) + built_dens + pop_dens + slope + 
               mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, 
               family = "binomial",
               data = data.k.train)
glm.t.d <- glm(formula = change ~ factor(landuse) + x+y+ built_dens + pop_dens + slope + 
               mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist,
               family = "binomial",
               data = data.d)
glm.t.s <- glm(formula = change ~ factor(landuse) + x+y+  built_dens + pop_dens + slope + 
               mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist,
               family = "binomial",
               data = data.s)
summary(glm.t.k)
summary(glm.t.d)
summary(glm.t.s)

calc_roc(glm.t.k, test_data = data.test)
calc_roc(glm.t.d)
calc_roc(glm.t.s)




##############################################################
# significant factors stratified samples
glm.t.k.s <- glm(formula = change ~ factor(landuse) +  built_dens + slope,# + slope + center_dist, 
                 family = "binomial",
                 data = data.k.train)
summary(glm.t.k.s)
calc_roc(glm.t.k.s)
calc_roc(glm.t.k.s, test_data = data.k.test)


glm.t.d.s <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + 
                 mRoads_dist + pRoads_dist+ center_dist,
                 family = "binomial",
                 data = data.d)
summary(glm.t.d.s)
calc_roc(glm.t.d.s)


glm.t.s.s <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope +
                 mRoads_dist + pRoads_dist,
                 family = "binomial",
                 data = data.s)
summary(glm.t.s.s)
calc_roc(glm.t.s.s)


##############################################################
# significant factors stratified samples new (one sample per strata and class)
data.k.dummy <- cbind(data.k, to.dummy(data.k$landuse, prefix = "landuse"))
glm.t.k.s <- glm(formula = change ~ factor(landuse) +  built_dens +  
                   mRoads_dist + train_dist, 
                 family = "binomial",
                 data = data.k)
summary(glm.t.k.s)
calc_roc(glm.t.k, test_data = data.k)
# calc_roc(glm.t.k.s, test_data = data.k.test)


glm.t.d.s <- glm(formula = change ~ factor(landuse) +  built_dens + 
                 pRoads_dist + center_dist,
                 family = "binomial",
                 data = data.d)
summary(glm.t.d.s)
calc_roc(glm.t.d.s)


glm.t.s.s <- glm(formula = change ~ factor(landuse) +  built_dens +
                   mRoads_dist + pRoads_dist,
                 family = "binomial",
                 data = data.s)
summary(glm.t.s.s)
calc_roc(glm.t.s.s)




# without sampling
glm.t.k.a <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope + 
                   mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.k.a)
glm.t.d.a <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope + 
                   mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.d.a)
glm.t.s.a <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope + 
                   mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.s.a)
summary(glm.t.k.a)
summary(glm.t.d.a)
summary(glm.t.s.a)

# systematic sampled
glm.t.dsys <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope + 
                    mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.d.sys)
summary(glm.t.dsys)
glm.t.dsys.s <- glm(formula = change ~ factor(landuse) +  pop_dens + pRoads_dist, data = data.d.sys)
summary(glm.t.dsys.s)








############################################################################################################################
# interaction between independent variables
############################################################################################################################

attach(data.k, name = "samples")
detach(samples)

##############
# pairwise correlations among the predictors 
cor(data.k.a[4:14])
cor(data.d.a[4:14])
cor(data.s.a[4:14])
(correlation.k <- (abs(cor(data.k[5:15])) > 0.7))
(correlation.d <- (abs(cor(data.d[5:15])) > 0.7))
(correlation.s <- (abs(cor(data.s[5:15])) > 0.7))

# continuous to categorical
summary(lm(mRoads_dist ~ factor(landuse)))
plot(mRoads_dist ~ factor(landuse))

# continuous to continuous
plot(data.k$pop_dens, data.k$built_dens)
abline(lm(pop_dens ~ built_dens, data = data.k))
plot(data.k$center_dist, data.k$airport_dist)
abline(lm(center_dist ~ airport_dist, data = data.k))
plot(data.d$pop_dens, data.d$built_dens)
plot(data.s$pop_dens, data.s$built_dens)
plot(data.s$center_dist, data.s$pRoads_dist)
plot(data.s$train_dist, data.s$pRoads_dist)
plot(data.s$center_dist, data.s$river_dist)
abline(lm(river_dist ~ center_dist, data = data.s))
plot(data.s$airport_dist, data.s$river_dist)
plot(data.s$center_dist, data.s$train_dist)

abline(h=mean(pRoads_dist),lty=2); abline(v=mean(built_dens),lty=2)
cor.test(pRoads_dist, built_dens)





############################################################################################################################
# from report


logit.roc <- function(model, steps=20) {
  # get the response field
  # from the model object
  field.name <- attr(attr(terms(formula(model)), "factors"),
                     "dimnames")[[1]][1]
  # and extract the T/F from it
  eval(parse(text=paste("tmp <- ",
                        ifelse(class(model$data) == "data.frame", "model$data$", ""),
                        field.name, sep="")))
  r <- data.frame(pts = seq(0, 1-(1/steps), by=1/steps),
                  sens = 0, spec=0);
  for (i in 0:steps) {
    thresh <- i/steps;
    r$sens[i] <- sum((fitted(model) >= thresh) & tmp)/sum(tmp);
    r$spec[i] <- sum((fitted(model) < thresh) & !tmp)/sum(!tmp)
  }
  return(r)
}

logit.roc.area <- function(r) {
  area <- 0;
  for (i in 1:(length(r$pts)-1))
    area <- area + ((1 - r$sens[i+1]) - (1 - r$sens[i])) *
      ((r$spec[i+1] + r$spec[i])/2);
  return(area)
}

logit.roc.plot <- function(r, title="ROC curve") {
  old.par <- par(no.readonly = TRUE); on.exit(par(old.par))
  par(xaxs="i", yaxs="i")
  plot(1 - r$spec, r$sens, xlim=c(0, 1), ylim=c(0,1), type="l",
       xlab="(1 - specificity): false positive rate",
       ylab="sensitivity: true positive rate",
       col="blue", lwd=2);
  points(1 - r$spec, r$sens, pch=20, cex=2, col="blue");
  abline(0, 1, lty=2);
  segments(1-r$spec, 1-r$spec, 1-r$spec, r$sens, lty=2)
  text(0, 0.9, paste("Area under ROC:",round(logit.roc.area(r),4)), pos=4)
  title(main = title)
}

r <- logit.roc(glm.t.k.s, steps=100)
logit.roc.area(r)
logit.roc.plot(r, "ROC for tenure, roads, settlements")





############################################################################################################################
####### descriptive analysis ###############################################################################################
############################################################################################################################


str(data)
data[1:5,]
summary(data)

attach(data)
search()

##############################################################
# cross classification
##############################################################
(ct <- table(change, landuse))
summary(ct)
(cs <- chisq.test(ct))
# expected change if random assignment
round(cs$expected)
ct - round(cs$expected)

# normalize to 1
(ct.p <- round(t(t(ct)/apply(ct,2,sum)),2))

# bar plot
# par(mfrow=c(2,2))
col.vec <- c("gray10", "gray60")
barplot(ct.p, col=col.vec, main="Proportion change",
        xlab="landuse")




