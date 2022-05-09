##--------------------------------------------------------------------------------------------------------
## SCRIPT : Model with 4 best variables
##
## Authors : Auriane Virgili
## Last update : 2020-11-18
##
## R version 4.0.2 (2020-06-22) -- "Taking Off Again"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##
##--------------------------------------------------------------------------------------------------------

library(mgcv)
library(maptools)
library(raster)

#----------------------------------------------------------------------------------------------
# 1- Data
#----------------------------------------------------------------------------------------------
setwd("")
obs<-read.csv("Effort_5km_ATL_BW_Covariates.csv", h=T, sep=";", dec=".")

#----------------------------------------------------------------------------------------------
# 2- Data formatting
#----------------------------------------------------------------------------------------------

dat_full <- obs[,c(1:33,34:35,37:38,43:44,46:47,49:50,55:56,58:59,61:62,67:68,70:71,73:74,79:80,82:83,85:86,91:92)]

#----------------------------------------------------------------------------------------------
# 3- Remove outliers 
#----------------------------------------------------------------------------------------------

dat_full[which(dat_full$Seg_Length_km == 0),17] <- NA 
dat_full$SEASTATE<-as.numeric(as.character(dat_full$SEASTATE))
dat_full[which(dat_full$SEASTATE == 4),11] <- NA 

saturate <- function(x, threshold = NULL, lower = TRUE, alpha = 0.01, both = FALSE) {
  if(both) {
    if(is.null(threshold)) {
      threshold <- quantile(x, c(alpha/2, 1-alpha/2), na.rm = TRUE)
    }
    if(length(threshold) != 2){
      stop("must provide a lower and upper thresholds")
    }
    x <- ifelse(x < min(threshold), min(threshold), x)
    x <- ifelse(x > max(threshold), max(threshold), x)
  }
  else{
    if(is.null(threshold)) {
      threshold <- quantile(x, alpha, na.rm = TRUE)
    }
    if(lower) { x <- ifelse(x < threshold, threshold, x) }
    else {
      threshold <- quantile(x, 1-alpha, na.rm=TRUE)
      x <- ifelse(x > threshold, threshold, x) }
  }
  return(x)
}

saturation <- apply(dat_full[30:63],2, saturate, alpha = 0.01, both = TRUE, lower = FALSE)

dat_full <- cbind(dat_full[1:29], saturation)

dat<-dat_full[-which(is.na(dat_full[,c(11,17,30:63)]),arr.ind=T),] 

dat <- dat[,c(1:57)]

dat$Effort <- dat$Seg_Length_m/1000

setwd("")
best4var <- read.table("RankedVariablesByAkaikeWeight.txt",sep="\t",h=T)

formula <- paste("PODSIZE~ 1 +s(", best4var[1,1], ", k = 4, bs = 'tp') + s(", best4var[2,1], ", k = 4, bs = 'tp') + s(", 
                 best4var[3,1], ", k = 4, bs = 'tp') + s(", best4var[5,1], ", k = 4, bs = 'tp')", sep="") ## on choisit la 5eme car 3 et 4 sont correlees

modBest4var <- gam(as.formula(as.character(formula)), data = dat, 
                   offset = log(2*dat$esw*dat$Effort), 
                   family = tw(), method = "REML",
                   weights = ifelse(dat$PLATFORM == "Boat", 1, ifelse(dat$PODSIZE == 0, 1/5, 1))/mean(ifelse(dat$PLATFORM == "Boat", 1, ifelse(dat$PODSIZE == 0, 1/5, 1)))
)
summary(modBest4var)

### plot model ###

theta <- function(df, var_name, gam_model, unlog = TRUE) {
  # df is the original dataset used to calibrate the model
  # var_name is a vector of the variable name used in the model
  # gam_model is dsm model you want to use
  # set unlog to TRUE to obtain curves on the density/abundance scale
  lower <- function(x) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = 0.95)[1]) }
  upper <- function(x) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = 0.95)[2]) }
  nx <- 1e3
  n_sim <- 1e4
  X <- as.data.frame(sapply(var_name, function(id) {rep(mean(df[, id]), nx)}))
  Y <- NULL
  for(j in var_name) {
    Z <- X
    Z[, j] <- seq(min(df[, j]), max(df[, j]), length.out = nx)
    Z <- predict(gam_model, newdata = Z, off.set = 1, type = "lpmatrix")
    beta <- mvtnorm::rmvnorm(n_sim, mean = gam_model$coefficients, sigma = gam_model$Vc)
    linpred <- beta %*% t(Z); rm(Z, beta)
    if(unlog) { linpred <- exp(linpred) }
    Y <- rbind(Y,
               data.frame(x = seq(min(df[, j]), max(df[, j]), length.out = nx),
                          y = apply(linpred, 2, mean),
                          lower = apply(linpred, 2, lower),
                          upper = apply(linpred, 2, upper),
                          param = rep(j, nx)
               )
    )
  }
  return(Y)
}

pred <- (do.call('rbind', lapply(list(modBest4var), theta, df = dat, var_name = c("Depth", "Rough", "SstMSurf", "SstSd0_200"), unlog = TRUE)))

pred$response <- c("Densité relative")

pred <- pred
pred[which(pred$param=="SstSd0_200" & pred$x < 0.2),"x"] <-NA
pred[which(pred$param=="SstSd0_200" & pred$x > 0.5),"x"] <-NA
pred[which(pred$param=="SstMSurf" & pred$x > 16),"x"] <-NA
pred[which(pred$param=="SstMSurf" & pred$x < 13),"x"] <-NA

### plot 

library(ggplot2)
library(gridExtra)

theme_set(theme_bw(base_size = 18))
g1 <- ggplot(pred,
             aes(x = x, y = y), color = "midnightblue") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red") +
  geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "lightblue", color = "lightblue", alpha = 0.5) +
  geom_line(color = "midnightblue") +
  facet_grid(response~param, scales = "free") +
  ylab("") + xlab("") +
  coord_cartesian(ylim = c(0, 0.003)) +
  ggtitle("BW - ATLANTIQUE") +
  theme(legend.position = "top", 
        plot.title = element_text(lineheight = 0.8, face = "bold", size = 16), 
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
g1
ggsave( "/ModBest4variables_BWbob.png",g1,
        width = 11, height = 6, dpi = 300)


### model prediction ###

setwd("")
grille_var <- read.csv("PredictionGrid_WGS84_ATL_1-12.csv", h=T, sep=";", dec=".") 

best4variables <- all.vars(formula(modBest4var)[-1])

dat_pred_bestVar <- grille_var[,c("x","y",as.character(best4variables))]

dat_pred_bestVar[which(dat_pred_bestVar[,3]<min(obs[,best4variables[1]],na.rm=T)),3]<-NA
dat_pred_bestVar[which(dat_pred_bestVar[,3]>max(obs[,best4variables[1]],na.rm=T)),3]<-NA

dat_pred_bestVar[which(dat_pred_bestVar[,4]<min(obs[,best4variables[2]],na.rm=T)),4]<-NA
dat_pred_bestVar[which(dat_pred_bestVar[,4]>max(obs[,best4variables[2]],na.rm=T)),4]<-NA

dat_pred_bestVar[which(dat_pred_bestVar[,5]<min(obs[,best4variables[3]],na.rm=T)),5]<-NA
dat_pred_bestVar[which(dat_pred_bestVar[,5]>max(obs[,best4variables[3]],na.rm=T)),5]<-NA

dat_pred_bestVar[which(dat_pred_bestVar[,6]<min(obs[,best4variables[4]],na.rm=T)),6]<-NA
dat_pred_bestVar[which(dat_pred_bestVar[,6]>max(obs[,best4variables[4]],na.rm=T)),6]<-NA

pred_bestVar <- predict(modBest4var,newdata=dat_pred_bestVar,type="response",se.fit=T)

pred_bestVarCoord <- cbind(grille_var[,2:3],pred_bestVar$fit,pred_bestVar$se.fit)


setwd("")

raster_moy <- rasterFromXYZ(pred_bestVarCoord[,1:3], crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) 
writeRaster(raster_moy,"Prediction_ModBest4var_BW_bob.tif",formats=GTiff,overwrite=T)

raster_CV <- rasterFromXYZ(pred_bestVarCoord[,c(1,2,4)], crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) 
writeRaster(raster_CV,"Prediction_ModBest4var_BW_bob.tif",formats=GTiff,overwrite=T)
