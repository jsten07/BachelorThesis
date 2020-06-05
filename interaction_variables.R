############################################################################################################################
# interaction between independent variables
############################################################################################################################

# attach(data.k, name = "samples")
# detach(samples)

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

