## Spatially-varying coefficients (SVC) model

* Here we use the R package *spBayes* to fit the SVC model to a subset of the BTH region in both seasons of 2015 for $PM_2.5$ concentration datasets; the SVC model is fitted to $PM_2.5$concentration on the square root scale, and the CMAQ covariate on the square root scale.
* The predictions are backtransformed to the original scale for predictive performance assessment.
* Parameters of the SVC model are estimated by running an MCMC algorithm for n.samples number of iterations. Posterior inference is based on the MCMC samples post burn-in (burn.in number of iterations).
* The SVC model employs an exponential covariance function to model the spatial dependence in observed PM$_{2.5}$ concentration measured at monitors with covariance parameters sigma.sq, phi and tau.sq.
* A flat priors is placed on the regression coefficients (e.g. intercept and coefficient of CMAQ), an Inverse Gamma prior is placed on sigma.sq and tau.sq, and a Uniform prior is specified for phi.
* The R workspace contains:

* Data:

- train_dat: square root PM$_2.5$concentration and covariates at 12 training cities from Jun 1, 2015 -Aug 31, 2015 and from Nov 1, 2015 -Jan 31, 2016
- test_dat:  square root PM$_2.5$ concentration and covariates at a testing cities from Jun 1, 2015 -Aug 31, 2015 and from Nov 1, 2015 -Jan 31, 2016


```r

rm(list=ls())
source("./LoadPackages/RDependPackages.R")
data("SiteData", package = "HDCM")
{
  
  PM25_2015w <- obs_PM25_2015w %>% setorderv(c("CITY", "ID","DATE_TIME"))
  setDF(PM25_2015w)
  ##################################################################
  ###################################################################
  #                           1. Data loading
  ###################################################################
  
  City.Name <- as.character(unique(PM25_2015w$CITY))
  setDF(PM25_2015w)
  DATE_TIME <- unique(PM25_2015w$DATE_TIME) %>% sort()
  Nt <- length(DATE_TIME)
  date.time <- data.frame(time.index = 1:Nt,
                          time.scale = seq(0, 1, , Nt),
                          DATE_TIME = DATE_TIME)
  PM25_2015w <- PM25_2015w  %>% left_join(date.time,  by = c("DATE_TIME"))
  PM25_2015w[, c("REAL_PM25")] <- sqrt(PM25_2015w[,  c("REAL_PM25")])
  PM25_2015w[, c("sim50_CMAQ_PM25")] <- sqrt(PM25_2015w[,  c("sim50_CMAQ_PM25")])
  
  Covariate <- c("sim50_CMAQ_PM25"
                 , "sim_TEMP"
                 , "sim_WIND_X"
                 , "sim_WIND_Y"
  );
  
  setDF(PM25_2015w)
  Cov.Index <- which(base::colnames(PM25_2015w) %in% Covariate)
  if(length(Cov.Index) > 1){
    mean_covariates <- apply(PM25_2015w[, Cov.Index], 2, mean)
    sd_covariates <- apply(PM25_2015w[, Cov.Index], 2, sd)
    PM25_2015w[, Cov.Index] <- scale(PM25_2015w[, Cov.Index],
                                     center = mean_covariates,
                                     scale = sd_covariates)
  }
}

# ###################################################################
# #                           2. Model
# ###################################################################
region <- unique(PM25_2015w$CITY)
###################################################################
#                       SVC model
###################################################################
{
  p <- length(Covariate) + 1
  # set parameters
  {
    n.samples <- 1e4
    
    starting <- list("phi" = 1e-5, "sigma.sq" = 1
                     , "tau.sq" = .1, "nu" = 0.5,
                     'beta.starting'= c(0, 0, rep(0, p - 2)))
    tun <- 0.1
    tuning <- list("phi" = 1e-8, "nu" = 1e-2, "sigma.sq" = tun
                   , "tau.sq" = tun, 'beta' = c(tun, tun))
    
    priors.1 <- list("beta.Norm" = list(rep(0, p), diag(1e5, p)),
                     "phi.Unif" = c(1/1E7, 1/1e4), #1/2e2, 1/1e1
                     "sigma.sq.IG" = c(2e0, 1e0),
                     "nu.Unif" = c(0.1, 3), #1/2e2, 1/1e1
                     "tau.sq.IG" = c(2, 1e0))
    
    cov.model <- "exponential"
    
    n.report <- 5000
    verbose <- F
  }
  
  
  region <- sort(as.character(unique(PM25_2015w$CITY)))
  region_num <- 1:length(region)
  setDT(PM25_2015w)
  
  Ens <- ceil(0.5*n.samples)
  
  colNames <- c("CITY","LON", "LAT", "DATE_TIME",
                "YEAR_MONTH", "YEAR","MONTH","DAY",
                "REAL_PM25")
  start.time <- proc.time()
  for(r in region_num)
  {
    
    year_range <- unique(PM25_2015w$YEAR)
    tem <- PM25_2015w[CITY %in% region[r], "ID"]
    
    for(Year in year_range)
    {
      Base_Table <- PM25_2015w[YEAR == Year,]
      month_range <- unique(Base_Table$MONTH)
      for(Month in month_range)
      {
        Base_Tab <- Base_Table[MONTH == Month,]
        day_range <- unique(Base_Tab$DAY)
        for(Day in day_range)
        {
          set.seed(1234)
          cat("\n\n   the ", r, "th region: ", region[r], "!!!\n\n")
          cat("   year: ", Year, "; month: ", Month, "; day: ", Day," \n\n")
          cat("...................SVC.................\n\n")
          # Database
          Da.mod <- Base_Tab[CITY %nin% region[r] & DAY == Day, ]
          Da.pre <- Base_Tab[CITY %in% region[r] & DAY == Day, ]
          
          setDF(Da.mod);setDF(Da.pre);
          Y = Da.mod[, "REAL_PM25"]
          X = as.matrix(cbind(Da.mod[, c(Covariate)]))
          # colnames(X) <- c("intercept", Covariate)
          colnames(X) <- c(Covariate)
          coords <- as.matrix(Da.mod[, c("LON_X", "LAT_Y")])
          
          
          m.1 <- spLM(Y ~ (X), coords = coords
                      , starting = starting, tuning = tuning
                      , priors = priors.1, cov.model = cov.model
                      , n.samples = n.samples, verbose = verbose
                      , n.report = n.report)
          # View parameter estimates
          # m.2 <- spRecover(m.1, start = Ens + 1, verbose = FALSE)
          # round(summary(m.2$p.theta.recover.samples
          #               )$quantiles, 2)
          X <-  as.matrix((Da.pre[, c(Covariate)]))
          X <- cbind(1, X)
          coords <- as.matrix(Da.pre[, c("LON_X", "LAT_Y")])
          y.pred <- spPredict(m.1, pred.covars = X,
                              pred.coords = coords,
                              start = Ens + 1,
                              verbose = F)
          testPred <- ifelse(y.pred$p.y.predictive.samples < 0, 0, 
                             y.pred$p.y.predictive.samples^2)
          y.pred <- apply(testPred, 1, quant)
          PM25.Pred.sd <- apply(testPred, 1, sd)
          
          temp.U95 <-  y.pred[3, ]
          temp.Pred <- y.pred[2, ]
          temp.L25 <-  y.pred[1, ]
          
          Da.pre$REAL_PM25 <- Da.pre$REAL_PM25^2
          spT <- spT_validation(Da.pre$REAL_PM25, temp.Pred,
                                sigma = NA,
                                zhat.Ens = NULL,
                                names = F, CC = F)[c(1, 4)]
          print(spT)
          
          if(r == region_num[1] & Year == year_range[1]  &
             Month == month_range[1] & Day == day_range[1])
          {
            SVC <- data.frame(Da.pre[, colNames],
                              PM25.L25 = temp.L25,
                              PM25.Pred = temp.Pred,
                              PM25.U95 = temp.U95,
                              Pred.sd = PM25.Pred.sd)
            
          }else{
            SVC <- rbind(SVC, data.frame(Da.pre[, colNames]
                                         , PM25.L25 = temp.L25
                                         , PM25.Pred = temp.Pred
                                         , PM25.U95 = temp.U95
                                         , Pred.sd = PM25.Pred.sd))
          }
          temp1 <- Validation.Group.Region(SVC, Sigma = NA,
                                           col = c("REAL_PM25", "PM25.Pred"),
                                           by = "CITY")
          cat("\n.............................\n")
          print(temp1)
        }
        cat("\n.............................\n")
      }
    }
    
    temp0 <- Validation.Group.Region(SVC, Sigma = NA,
                                     col = c("REAL_PM25", "PM25.Pred"),
                                     by = "CITY")
    cat("\n.............................\n")
    print(temp0)
  }
}
end.time <- proc.time()
SVC$run_time <- (end.time - start.time)[3]
# writexl::write_xlsx(temp0, path = "./Result/SVCx_w_cv.xlsx")
writexl::write_xlsx(SVC, path = "./Result/pred_SVCxz_w_cv.xlsx")
```
