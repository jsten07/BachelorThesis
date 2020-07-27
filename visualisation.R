rm(list=ls())


#########################################################################################
# results
#########################################################################################

source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/regression.R")

library(formattable)

#################################################################
# load data
#################################################################
 
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/stratified/")

data.k <- read.csv("krakow.csv")
data.d <- read.csv("dresden.csv")
data.s <- read.csv("sevilla.csv")

data.k2 <- read.csv("krakow_90.csv")
data.d2 <- read.csv("dresden_90.csv")
data.s2 <- read.csv("sevilla_90.csv")

str(data.k2)
str(data.d2)
str(data.s2)

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/samples/")
data.k.a <- read.csv("krakow_all.csv")
data.d.a <- read.csv("dresden_all.csv")
data.s.a <- read.csv("sevilla_all.csv")

str(data.d.a)
str(data.k.a)
str(data.s.a)


#################################################################
# calculate differnt models
#################################################################
trained.s <- ffs_model(data.s2)
model.s.ffs <- trained.s$ffs$finalModel
model.s.glm <- trained.s$glm$sig
model.s.glmA <- trained.s$glm
summary(model.s.ffs)
calc_roc(model.s.ffs)
logit.plot.quad(trained.s$ffs_plus$finalModel)
bss_model(data.s2, "s")


trained.d <- ffs_model(data.d2)
model.d.ffs <- trained.d$ffs$finalModel
model.d.glm <- trained.d$glm$sig
model.d.glmA <- trained.d$glm
summary(model.d.ffs)
calc_roc(model.d.ffs)
logit.plot.quad(trained.d$ffs_plus$finalModel)
bss_model(data.d2, "d")


trained.k <- ffs_model(data.k2)
model.k.ffs <- trained.k$ffs$finalModel
model.k.glm <- trained.k$glm$sig
model.k.glmA <- trained.k$glm
summary(model.k.ffs)
calc_roc(model.k.ffs)
logit.plot.quad(trained.k$ffs_plus$finalModel)
bss_model(data.k2, "k")


#################################################################
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
  
  coefs["roc", 1] <- calc_roc(model1, data.d.a)
  coefs["roc", 3] <- calc_roc(model2, data.s.a )
  coefs["roc", 5] <- calc_roc(model3, data.k.a)
  
  message("Check data for Morans!")
  coefs["Morans I", 1] <- calc_moransI(data.d2, k = 4)$estimate[1]
  coefs["Morans I", 3] <- calc_moransI(data.s2, k = 4)$estimate[1]
  coefs["Morans I", 5] <- calc_moransI(data.k2, k = 4)$estimate[1]
  
  options(scipen = 99, digits = 3)
  
  coefs <- as.data.frame(coefs)
  #coefs.sign <- as.data.frame(coefs.sign)
  
  print(formattable(coefs, format = "f", digits = 5))
  # print(formattable(coefs.sign, digits = 3, title = "p-value"))
  
  return(coefs)
}



### save tables

setwd("C:/Users/janst/sciebo/Bachelor Thesis/results/models/")

ffs_models <- compare_coefs(model.d.ffs, model.s.ffs, model.k.ffs)
write.csv(ffs_models, "ffs_models_90.csv")
glm_models <- compare_coefs(model.d.glm, model.s.glm, model.k.glm)
write.csv(glm_models, "glm_models_90.csv")
glm_allV_models <- compare_coefs(trained.d$glm, trained.s$glm, trained.k$glm)
write.csv(glm_allV_models, "glm_allV_models_90.csv")
train_models <- compare_coefs(trained.d$train$finalModel, trained.s$train$finalModel, trained.k$train$finalModel)
write.csv(train_models, "train_models_90.csv")
ffs_plus_models <- compare_coefs(trained.d$ffs_plus$finalModel, trained.s$ffs_plus$finalModel, trained.k$ffs_plus$finalModel)
write.csv(ffs_plus_models, "ffs_plus_models_90.csv")

setwd("C:/Users/janst/sciebo/Bachelor Thesis/results/models_RDS/")
trained_bss.k <- readRDS("trained_bss_k.rds")
trained_bss.d <- readRDS("trained_bss_d.rds")
trained_bss.s <- readRDS("trained_bss_s.rds")
trained_bss_plus.k <- readRDS("trained_bss_plus_k.rds")
trained_bss_plus.d <- readRDS("trained_bss_plus_d.rds")
trained_bss_plus.s <- readRDS("trained_bss_plus_s.rds")

setwd("C:/Users/janst/sciebo/Bachelor Thesis/results/models/")
bss_models <- compare_coefs(trained_bss.d$finalModel, trained_bss.s$finalModel, trained_bss.k$finalModel)
write.csv(bss_models, "bss_models_90.csv")
bss_plus_models <- compare_coefs(trained_bss_plus.d$finalModel, trained_bss_plus.s$finalModel, trained_bss_plus.k$finalModel)
write.csv(bss_plus_models, "bss_plus_models_90.csv")


### calculate model success
logit.plot.quad(trained_bss_plus.d$finalModel, title = "Model success Dresden")
logit.plot.quad(trained_bss_plus.k$finalModel, title = "Model success Krakow")
logit.plot.quad(trained_bss_plus.s$finalModel, title = "Model success Seville")

### calculate ROC
calc_roc(trained_bss_plus.d$finalModel, test_data = data.d.a, main = "ROC Dresden")
calc_roc(trained_bss_plus.s$finalModel, test_data = data.s.a, main = "ROC Seville")
calc_roc(trained_bss_plus.k$finalModel, test_data = data.k.a, main = "ROC Krakow")

calc_roc(trained.d$train$finalModel, test_data = data.d.a, main = "ROC Dresden train")
calc_roc(trained.s$train$finalModel, test_data = data.s.a, main = "ROC Seville train")
calc_roc(trained.k$train$finalModel, test_data = data.k.a, main = "ROC Krakow train")



#########################################################################################
# plots of prepared data
#########################################################################################
# source("C:/Users/janst/sciebo/Bachelor Thesis/R/BachelorThesis/preparation_functions.R")
setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/")
library(raster)
library(rasterVis)

# load grids
stack_dresden <- stack("created/stack/dresden.grd")
stack_krakow <- stack("created/stack/krakow.grd")
stack_sevilla <- stack("created/stack/sevilla.grd")

plot(stack_dresden)
plot(stack_krakow)
plot(stack_sevilla)

### plot non categorical data 
data_plots <- function(stack, city_name) {
  setwd("C:/Users/janst/sciebo/Bachelor Thesis/results/data/plots_font/")
  
  titles <- c("change", "built-up density (%)", "population density (people / cell)", 
                            "slope (%)", "landuse", "distance to major roads (m)", "distance to primary roads (m)", 
                             "distance to major river (m)", "distance to train stations (m)", 
                             "distance to city center (m)", "distance to airport (m)")

  
  for (i in c(1:nlayers(stack))) {
    plotname <- paste(paste(city_name, names(stack)[[i]], sep="_"), ".tiff", sep = "")
    
    if(names(stack[[i]]) == "change") {
      
     # plot(stack[[i]], main="Change to built up from 1990 to 2014")
      
    } else if(names(stack[[i]]) == "landuse") {
      
     # plot(stack[[i]], main = "landuse ...")
      
    } else {
      #tiff(plotname)
      par(cex = 1.5)
      plot(stack[[i]], main = titles[i], ylim = extent(stack)[3:4], xlim = extent(stack)[1:2])#, cex.main = 1.5, cex.axis = 1.5)
      #dev.off()
    }
  }
}

data_plots(stack_dresden, "Dresden")
data_plots(stack_krakow, "Krakow")
data_plots(stack_sevilla, "Sevilla")



### plot landuse and change data
# not pretty but levelplot does not work as expected within a function call but does when called manually !?!
# use extent and datasets for one of the cities and pass through the three plot functions

setwd("C:/Users/janst/sciebo/Bachelor Thesis/data/created/GHSL_R/")

y_l <- (y_d <- c(5645100, 5672900))
city_name <- "Dresden"
stack <- stack_dresden
built <- (dresden_built <- raster("dresden_built.tif"))

y_l <- (y_k <- c(5529000, 5560100))
city_name <- "Krakow"
stack <- stack_krakow
built <- (krakow_built <- raster("krakow_built.tif"))

y_l <- (y_s <- c(4131000, 4151000))
city_name <- "Sevilla"
stack <- stack_sevilla
built <- (sevilla_built <- raster("sevilla_built.tif"))


setwd("C:/Users/janst/sciebo/Bachelor Thesis/results/data/plots_font/")
plotname <- paste(city_name, "_change.tiff", sep = "")

change <- as.factor(stack$change)
lev <- levels(change)[[1]]
lev[["changed"]] <- c("not changed", "changed")
levels(change) <- lev
# tiff(plotname)
par(cex = 1.5)
levelplot(change, colorkey = list(labels = list(cex = 1.4)), col.regions=rev(terrain.colors(2)),  main=list(label = "change to built up from 1990 to 2014\n", cex = 1.7), ylim = y_l, scales=list(y=list(rot=90, cex = 1.5), x=list(cex = 1.54)))
# dev.off()


plotname <- paste(city_name, "_landuse.tiff", sep = "")

lu <- as.factor(stack$landuse)
lev <- levels(lu)[[1]]
if(city_name == "Sevilla") {
  lev[["landcover"]] <- c("artificial", "crop", "forest", "open", "water")
  lu_colors <- c("#F2F2F2FF","#EEB99FFF", "#00A600FF", "#E6E600FF", "blue")
} else {
  lev[["landcover"]] <- c("artificial", "crop", "pasture", "forest", "water")
  lu_colors <- c("#F2F2F2FF", "#EEB99FFF", "#EAB64EFF", "#00A600FF", "blue")
}
levels(lu) <- lev
#tiff(plotname)
levelplot(lu, colorkey = list(labels = list(cex = 1.4)), col.regions=lu_colors, main=list(label = "Landuse in 1990\n", cex = 1.7), ylim = y_l, scales=list(y=list(rot=90, cex = 1.5), x= list(cex = 1.5)))
#dev.off()


plotname <- paste(city_name, "_built.tiff", sep = "")

built <- as.factor(built)
lev <- levels(built)[[1]]
# 3-4: changed from 1990 to 2014 -> 1
# 2: not built up in any epoch -> 0
# 0, 1, 5, 6: no data, water, built up before
lev[["built"]] <- c("water", "not urban", "2014", "2000", "1990", "1975")
levels(built) <- lev
colors <- c("#66CCFF", "black", rev(heat.colors(4)))
# change_plot <- levelplot(change, col.regions=rev(terrain.colors(2)), main="change to built up from 1990 to 2014", ylim = y_d)
#tiff(plotname)
levelplot(built, colorkey = list(labels = list(cex = 1.4)), col.regions=colors, main=list(label = "built-up epochs\n", cex = 1.7), ylim = y_l, scales=list(y=list(rot=90, cex = 1.5), x= list(cex = 1.54)))
#dev.off()





############################################################################
# boxplots
############################################################################

# bind data for boxplots
bind.k <- data.k2
bind.k["city"] <- "Krakow"
bind.s <- data.s2
bind.s["city"] <- "Sevilla"
bind.d <- data.d2
bind.d["city"] <- "Dresden"
bind <- rbind(bind.d, bind.s, bind.k)

bind.k <- data.k.a
bind.k["city"] <- "Krakow"
bind.s <- data.s.a
bind.s["city"] <- "Sevilla"
bind.d <- data.d.a
bind.d["city"] <- "Dresden"
bind.a <- rbind(bind.d, bind.s, bind.k)



par(mar = c(6, 5.5, 1,1))

boxplot(bind$built_dens ~ bind$change:bind$city, 
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
 )      
title(ylab = "built-up density (%)", line = 4, cex.lab = 1.5)


(boxplot(bind$pop_dens ~ bind$change:bind$city, 
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
) )
title(ylab = "population density (people / cell)", line = 4, cex.lab = 1.5)


(boxplot(bind$slope ~ bind$change:bind$city, 
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),
        ylab = "",
        xlab = NULL,
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
) )
title(ylab = "slope (%)", line = 4, cex.lab = 1.5)


boxplot(bind$mRoads_dist ~ bind$change:bind$city, 
        #names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),  
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
) 
title(ylab = "distance to major roads (m)", line = 4, cex.lab = 1.5)



(boxplot(bind$pRoads_dist ~ bind$change:bind$city, 
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),  
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
) )
title(ylab = "distance to primary roads (m)", line = 4, cex.lab = 1.5)



boxplot(bind$river_dist ~ bind$change:bind$city, 
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
) 
title(ylab = "distance to major river (m)", line = 4, cex.lab = 1.5)




(boxplot(bind$train_dist ~ bind$change:bind$city, 
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
))
title(ylab = "distance to train stations (m)", line = 4, cex.lab = 1.5)



(boxplot(bind$center_dist ~ bind$change:bind$city, 
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
))
title(ylab = "distance to city center (m)", line = 4, cex.lab = 1.5)



(boxplot(bind$airport_dist ~ bind$change:bind$city, 
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),
        xlab = NULL,
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
))
title(ylab = "distance to airport (m)", line = 4, cex.lab = 1.5)



par(mar = c(6, 6, 1, 0.1))
levels <- c("artificial", "crop", "pasture", "forest", "open", "water")
lu_colors <- c("#F2F2F2FF", "#EEB99FFF", "#EAB64EFF", "#00A600FF","#E6E600FF", "blue")
lu_table <- table(paste(bind$city, bind$change), bind$landuse)
(lu_table_per <- round(((lu_table)/apply(lu_table,1,sum)),2))
barplot(t(lu_table_per), 
        xlim = c(0,9),
        col = lu_colors,
        # names = c("not changed\nDresden", "changed\nDresden", "not changed\nKrakow", "changed\nKrakow", "not changed\nSevilla", "changed\nSevilla"),
        ylab = "",
        names = c("not\nchanged\nDresden", "changed\nDresden", "not\nchanged\nKrakow", "changed\nKrakow", "not\nchanged\nSevilla", "changed\nSevilla"),
        las = 2,
        cex.axis=1.4,
        cex.names = 1.4,
        cex.lab = 1.5,
        boxwex = 0.5, 
        )
title(ylab = "sample rate per landuse class", line = 4, cex.lab = 1.5)
leg <- legend("right", levels, fill = lu_colors, cex = 1.3)



# all_table <- table(bind.a$city, bind.a$change)
# (all_table_per <- round(((all_table)/apply(all_table,1,sum)),2))





############################################################################
# report plot function
# taken from report: http://www.css.cornell.edu/faculty/dgr2/_static/files/R_PDF/lcc.pdf
# sizes modified

## plot logistic curve, threshold, T/F +/-, sensitivity, specificity
## arguments
## model a fitted glm
## threshold cutoff for sensitivity/specificity, default 0.5
## title (optional)
logit.plot.quad <- function(model, threshold=0.5, title="Model success") {
  sf<-sort(fitted(model), index=T)
  # leave extra space at bottom
  par(mar=c(5,5,2,1)+.1); par(xaxs="i", yaxs="r")
  plot(sf$x, ylim=c(0,1), type="l", col="blue", lwd=3, xlab="",
       ylab="probability of change", cex.axis = 1.5, cex.lab = 1.5)
  abline(h=c(0,1), lty=1)
  # show threshold and crossover point
  abline(h=threshold,lty=2); text(0,threshold+.02,
                                  paste("threshold =", threshold), pos=4, cex=1.5)
  crossover <- sum(fitted(model) < threshold)
  abline(v=crossover,lty=2)
  #text(crossover,.05,"crossover",pos=4)
  text(crossover, threshold-.03,
       "fitted probability \nof change",col="blue",pos=4, cex = 1.5)
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
  text(0,.1,paste("True negatives:",tn), col="green4",pos=4, cex = 1.5)
  text(right,.1,paste("False positives:\n", fn), col="red",pos=4, cex = 1.5)
  text(right,.9,paste("True positives:\n", tp), col="green4",pos=4, cex = 1.5)
  text(0,.9,paste("False negatives:", fp), col="red",pos=4, cex = 1.5)
  title(main=title, cex.main = 1.5)
  title(sub=paste("Sensitivity:", round(tp/(tp+fp),4),
                  "; Specificity:", round(tn/(tn+fn),4)), line=4, cex.sub = 1.5)
}





