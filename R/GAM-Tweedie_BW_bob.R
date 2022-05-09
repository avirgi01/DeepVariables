##--------------------------------------------------------------------------------------------------------
## SCRIPT : Model selection
##
## Authors : Auriane Virgili
## Last update : 2020-09-01
##
## R version 4.0.2 (2020-06-22) -- "Taking Off Again"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##
##--------------------------------------------------------------------------------------------------------

library(mgcv)
library(maptools)

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

dat<-dat_full[-which(is.na(dat_full[,c(11,17,30:63)]),arr.ind=T),]  

dat <- dat[,c(1:63)]

dat_env<-dat[,c(30:63)]

#----------------------------------------------------------------------------------------------
# 4- Model selection
#----------------------------------------------------------------------------------------------

fit_all_gam <- function(envdata, outcome, predictors, 
                        family = "tweedie", esw = NULL, 
                        weight_by_g_naught = 1,
                        max_cor = 0.5, nb_max_pred = 4, complexity = 4
) {
  ## family must be one of "negative binomial", "poisson" or "tweedie"
  ## default is "negative binomial"
  
  rescale <- function(x) { (x - mean(x))/sd(x) } # Fonction qui normalise
  
  ## design matrix
  X <- envdata[, predictors]
  
  ## standardize
  envdata[, predictors] <- apply(X, 2, rescale) #Je normalise sur les colonnes
  
  ## prepare smooth terms
  writeLines("using gam() with thin-plate splines")
  smoothers <- paste("s(", predictors, ", k = ", complexity, ", bs = 'tp'", ")", sep = "")
  intercept <- "~ 1"
  
  ## all combinations among nb_max_pred
  all_x <- lapply(1:nb_max_pred, combn, x = length(predictors))
  
  ## check whether cross-correlation needs to be evaluated
  if(nb_max_pred == 1) {
    rm_combn <- c(rep(0, length(predictors) + 1))
  }
  else {
    ## identify which combination is to be removed
    rm_combn <- lapply(all_x[-1], function(mat) {
      sapply(1:ncol(mat), function(i) {
        rho <- cor(X[, mat[, i]], use="complete.obs") ; diag(rho) <- 0
        return(max(abs(as.numeric(rho))))
      })
    })
    rm_combn <- c(c(rep(0, length(predictors) + 1)), unlist(rm_combn))
  }
  
  ## Create list of models
  mlist <- function(n, y, predictors) {
    paste(y, 
          apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + "),
          sep = paste(intercept, "+", sep = " ")
    )
  }
  all_mods <- c(paste(outcome, intercept, sep = " "),
                unlist(lapply(1:nb_max_pred, mlist, y = outcome, predictors = smoothers))
  )
  
  ## remove combinations of variables that are too correlated
  all_mods <- all_mods[which(rm_combn < max_cor)]
  
  # suppress warnings
  options(warn=-1)
  
  ## fit the models
  if(is.null(esw)) { stop("Must Provide a value for esw") }
  else {
    if(length(esw) == 1 | length(esw) == nrow(envdata)) {
      if(any(names(envdata) == outcome) == FALSE) { stop("No response variable in envdata") }
      else {
        if(length(weight_by_g_naught) == 1 | length(weight_by_g_naught) == nrow(envdata)) {
          # rescale weights
          if(length(weight_by_g_naught) == 1) {
            w <- rep(1, nrow(envdata))
          }
          else{
            w <- weight_by_g_naught/mean(weight_by_g_naught)
          }
          my_dsm_fct <- function(x) {
            if(family == "negative binomial") {
              model <- gam(as.formula(x), data = envdata, offset = log(2*esw*envdata$Effort), family = nb(), method = "REML", weights = w)
            }
            if(family == "poisson") {
              model <- gam(as.formula(x), data = envdata, offset = log(2*esw*envdata$Effort), family = poisson(), method = "REML", weights = w)
            }
            if(family == "tweedie") {
              model <- gam(as.formula(x), data = envdata, offset = log(2*esw*envdata$Effort), family = tw(), method = "REML", weights = w)
            }
            ### store some results in a data frame
            data.frame(model = x,
                       Convergence = ifelse(model$converged, 1, 0),
                       AIC = model$aic,
                       GCV = model$gcv.ubre,
                       ResDev = model$deviance,
                       NulDev = model$null.deviance,
                       ExpDev = 100*round(1 - model$deviance/model$null.deviance, 3)
            )
          }
          all_fits <- lapply(all_mods, my_dsm_fct)
          ## Collapse to a data frame
          return(do.call(rbind, all_fits))
        }
        else { stop("Please check weight for g(0): provide either a single value or a vector with same length as envdata") }
      }
    }
    else {
      stop("Please check esw: provide either a single value or a vector with same length as envdata")
    }
  }
}

# change effort variable
dat$Effort <- dat$Seg_Length_m/1000

# model selection

Time_simul <- numeric(2)

Start <- Sys.time() 

dd <- fit_all_gam(envdata = dat, outcome = "PODSIZE",
                  predictors = names(dat_env),
                  esw = dat$esw,
                  nb_max_pred = 4,
                  weight_by_g_naught = ifelse(dat$PLATFORM == "Boat", 1, ifelse(dat$PODSIZE == 0, 1/5, 1)) ## !! valeur e changer en fonction de l'espece
)
Time_simul[1] <- difftime(Sys.time(), Start, units = "sec") # Time difference

# model results

dd.ord<-dd[order(dd$AIC),] 

library(qpcR)

aics<-akaike.weights(dd.ord[,3])

dd.ord$Delta.AIC<-aics$deltaAIC  

dd.ord$rel.likelihood<-aics$rel.LL 

dd.ord$Akaike.weight<-aics$weights 

setwd("")
write.table(dd.ord,"Results_fit_models_BWbob.txt",sep="\t",row.names=F)

Mod1 <- gam(as.formula(as.character(dd.ord$model[1])), data = dat, 
            offset = log(2*dat$esw*dat$Effort), 
            family = tw(), method = "REML",
            weights = ifelse(dat$PLATFORM == "Boat", 1, ifelse(dat$PODSIZE == 0, 1/5, 1))/mean(ifelse(dat$PLATFORM == "Boat", 1, ifelse(dat$PODSIZE == 0, 1/5, 1)))
)
Mod1
summary(Mod1)
plot(Mod1, shade=TRUE, pages=1, shade.col="lightblue", rug=T)
par(mfrow=c(2,2))
gam.check(Mod1)

#----------------------------------------------------------------------------------------------
# 5- Plot fonctions 
#----------------------------------------------------------------------------------------------
#### better plot
### spline curves
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

pred <- (do.call('rbind', lapply(list(Mod1), theta, df = dat, var_name = c("Slope", "SstMSurf", "SstSd0_200"), unlog = TRUE)))

pred$response <- c("Densité relative")

### plot

library(ggplot2)
library(gridExtra)

theme_set(theme_bw(base_size = 18))
g1 <- ggplot(pred,
             aes(x = x, y = y), color = "midnightblue") +
  geom_hline(yintercept = exp(Mod1$coefficients[1]), linetype = "dotted", color = "red") +
  geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "lightblue", color = "lightblue", alpha = 0.5) +
  geom_line(color = "midnightblue") +
  facet_grid(response~param, scales = "free") +
  ylab("") + xlab("") +
  coord_cartesian(ylim = c(0, 0.004)) +
  ggtitle("BW - surface") +
  theme(legend.position = "top", 
        plot.title = element_text(lineheight = 0.8, face = "bold", size = 16), 
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  )
g1
ggsave( "4_Mod1_BWbob_avec0_AllDepths_05.png",g1,
        width = 11, height = 6, dpi = 300)

# aic weight
library(data.table)

table <- data.frame(var = NA, count = NA, percent = NA, akaike.weight = NA)

var <- names(dat_env)

for (i in 1:length(var)) {
  select <- dd.ord[dd.ord$model %like% var[i], ]
  table[i,1] <- var[i]
  table[i,2] <- nrow(select) 
  table[i,3] <- nrow(select)/nrow(dd.ord)*100
  table[i,4] <- sum(select$Akaike.weight)*100
}
table <- table[order(-table$akaike.weight),]

