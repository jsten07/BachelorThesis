rm(list=ls())
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/")


# load data
data.k <- read.csv("krakow_samples.csv")
data.d <- read.csv("dresden_samples.csv")
data.s <- read.csv("sevilla_samples.csv")

data.k.a <- read.csv("krakow_all.csv")
data.d.a <- read.csv("dresden_all.csv")
data.s.a <- read.csv("sevilla_all.csv")

data.d.sys <- read.csv("dresden_syst.csv")

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



##############################################################
# logistic regression
##############################################################

glm.t.k <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope + 
               mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.k)
glm.t.d <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope + 
                 mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.d)
glm.t.s <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + slope + 
                 mRoads_dist + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.s)
summary(glm.t.k)
summary(glm.t.d)
summary(glm.t.s)

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

# significant factors
glm.t.k.s <- glm(formula = change ~ factor(landuse) +  built_dens + pop_dens + 
                 center_dist, data = data.k)
summary(glm.t.k.s)
glm.t.d.s <- glm(formula = change ~ factor(landuse) +  pop_dens + mRoads_dist + pRoads_dist, data = data.d)
summary(glm.t.d.s)
glm.t.s.s <- glm(formula = change ~ pop_dens + pRoads_dist + river_dist + train_dist + center_dist + airport_dist, data = data.s)
summary(glm.t.s.s)




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

r <- logit.roc(glm.t.dsys.s, steps=100)
logit.roc.area(r)
logit.roc.plot(r, "ROC for tenure, roads, settlements")
