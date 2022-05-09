##--------------------------------------------------------------------------------------------------------
## SCRIPT : Correlation matrix
##
## Authors : Auriane Virgili
## Last update : 2020-09-01
##
## R version 4.0.2 (2020-06-22) -- "Taking Off Again"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##
##--------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
#1- Data
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

dat_cor <- dat[,c(30:57)]

cor <- round(cor(dat_cor, method = c("pearson")),3)

library(corrplot)

names(dat_cor) <- c("Depth","Slope","Roughness","CanArea","mTsurf","sdTsurf","mGrTsurf","sdGrTsurf","mEKEsurf",     
                    "sdEKEsurf","mT0-200","sdT0-200","mGrT0-200","sdGrT0-200","mEKE0-200","sdEKE0-200","mT200-600","sdT200-600", 
                    "mGrT200-600","sdGrT200-600","mEKE200-600","sdEKE200-600","mT600-2000","sdT600-2000","mGrT600-2000","sdGrT600-2000","mEKE600-2000", 
                    "sdEKE600-2000") 

col_order <- c("Depth","Slope","Roughness","CanArea","mTsurf","sdTsurf","mGrTsurf","sdGrTsurf","mEKEsurf",     
               "sdEKEsurf","mT0-200","sdT0-200","mGrT0-200","sdGrT0-200","mEKE0-200","sdEKE0-200","mT200-600","sdT200-600", 
               "mGrT200-600","sdGrT200-600","mEKE200-600","sdEKE200-600","mT600-2000","sdT600-2000","mGrT600-2000","sdGrT600-2000","mEKE600-2000", 
               "sdEKE600-2000")

dat_cor <- dat_cor[, col_order]

mcor <- cor(dat_cor, method = c("pearson"))

mcor_ATL <- mcor

plot <- corrplot(mcor, type="upper", order="original", tl.col="black", tl.srt=45)

setwd("")
png("0_CorrelationMatrix_ATL.png",height=15,width=18,units="cm",res=600)
corrplot(mcor, type="upper", order="original", tl.col="black", tl.srt=45, tl.cex = 0.7, cl.cex = 0.8, mar=c(0,0,0,0))
dev.off()

