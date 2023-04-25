## Authors: 
## Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca), 
## Arumugam Thiagarajan (arumugam.thiagarajan@ec.gc.ca)
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023

library(caret)
library(dplyr)
library(sensitivity)
library(boot)
library(doParallel)
library(parallel)
library(foreach)
library(purrr)
library(tidyverse)
library(ggplot2)

run.GSA.IPCCT2 <- function(main.dir,parameter_bounds,stocks_data,calibration_data,climate_data,initial_c,method,sample_size) {
  # read prior distribution from a csv file
  # (Required columns: Parameter, value, lower, upper)
  paramBounds <- parameter_bounds %>%
    arrange(Parameter)
  
  # names of parameters that are allowed to vary
  varSI       <- paramBounds$Parameter
  nParams     <- length(varSI)
  
  # sample size (10 used for illustration purposes)
  # (1024 used in Gurung et al., 2020)
  N <- sample_size
  
  # Sobols method required 2 random matrix
  m1 = matrix(runif(nParams*N), nrow=N);
  m2 = matrix(runif(nParams*N), nrow=N);
  M1 <- matrix(0, nrow = N, ncol = nParams)
  M2 <- matrix(0, nrow = N, ncol = nParams)
  
  # transform standard uniform to prior distribution 
  for(i in 1:nParams){
    pos <- which(paramBounds$Parameter == varSI[i])
    lower <- paramBounds[pos, "lower"]
    upper <- paramBounds[pos, "upper"]
    M1[, i] <- qunif(m1[, i], min = lower, max = upper)
    M2[, i] <- qunif(m2[, i], min = lower, max = upper)
  }
  X1 = data.frame(M1)
  X2 = data.frame(M2)
  names(X1) <- varSI
  names(X2) <- varSI
  
  if(method == "fast99") {
    fast99_qargs <- paramBounds %>%
      group_by(Parameter) %>%
      select(Parameter, min = lower, max = upper) %>%
      group_split %>%
      map(function(x) {y <- x %>%
        select(-Parameter) %>%
        as.list()
      })
  }
  
  if (method=="soboljansen"){
    si_obj2 <- sensitivity::soboljansen(model = NULL, X1 = X1, X2 = X2, nboot = 100)
  }
  if (method=="fast99"){
    si_obj2 <- sensitivity::fast99(model = NULL,
                                                   factors = paramBounds$Parameter,
                                                   n = N, M = 4,
                                                   q="qunif",
                                                   q.arg = fast99_qargs)
  }
  
  X <- si_obj2$X
  X <- cbind("SampleID" = 1:nrow(X), X)
  # NaN values can be in X if sample size was too small.
  if(any(is.nan(unlist(X)))) {warning("gsa: NaN values detected in sensitivity analysis. Try increasing sample_size")}
  
  
  params_list_sorted_names <- c("SampleID",varSI)
  params_list <- X %>%
    rowwise %>%
    group_split %>%
    map(function(x) {
      y <- x %>%
        pivot_longer(everything()) %>%
        deframe()
      z <- split(unname(y),names(y))
      return(z)
    }) %>%
    map(~ .[params_list_sorted_names])
  
  # Run the model and calculate log-likelihood
  # likelihoods were calculated assuming that the error (modeled - mesasured) are iid 
  
  Lkhood <- NULL
  Lkhood_list <- list()
  
  loglik=function(m,o){
    if(length(m)!=length(o)){
      print("Inequal number of modeled and observed values, cannot proceed")
      return()
    }
    
    res=log(m)-log(o)
    sigma=sqrt(mean(res^2))
    n=length(m)
    lk=-n*log(sigma)-(1/(2*sigma^2))*sum(res^2)
    return(lk)
    
  }
  
  site.list=unique(calibration_data$site)
  
  n_clust=detectCores()
  cl <- makeCluster(n_clust-4)
  registerDoParallel(cl)
  
  df_par <- foreach(site=site.list, .combine=rbind,  .inorder=TRUE) %dopar% {
  
    library(caret)
    library(dplyr)
    library(sensitivity)
    library(boot)
    library(doParallel)
    library(parallel)
    library(foreach)
    library(purrr)
    library(tidyverse)
    
    source(paste0(main.dir,"scripts//Recent//ipcct2_SOCTier2Model.r"))
    loglik=function(m,o){
      if(length(m)!=length(o)){
        print("Inequal number of modeled and observed values, cannot proceed")
        return()
      }
      
      res=log(m)-log(o)
      sigma=sqrt(mean(res^2))
      n=length(m)
      lk=-n*log(sigma)-(1/(2*sigma^2))*sum(res^2)
      return(lk)
      
    }
    
  for (i in 1:nrow(X)){
    site.sub=calibration_data[calibration_data$site==site,] %>%
      arrange(year) %>%
      distinct() %>%
      mutate(Diff = year -lag(year)) %>%
      replace_na(list(Diff=0))
    
    ## Address discontinuity in data ## Needs to be improved to interpolate gap in time series, now it just cut the time series when incomplete
    if (any(site.sub$Diff>1)) {
      problems=site.sub[site.sub$Diff>1,]
      site.sub <- site.sub %>% 
        slice(1:min(as.numeric(row.names(problems))-1))
    }
    
    parameters = X[i,]
    
    init.df <- initial_c[initial_c$site==site,]
    
    init_active = init.df[1,]$init_active
    init_slow = init.df[1,]$init_slow
    init_passive = init.df[1,]$init_passive
    
    climate_in <- climate_data[climate_data$POLYID==site.sub[1,]$POLYID,]
    
    climate_normal <- climate_in %>% group_by(month) %>% summarise_all(mean) %>% select(month,tavg,mappet,irrig)
    
    if (min(site.sub$year)<1981){
      climate_int <- data.frame(POLYID=site.sub[1,]$POLYID,
                                year=rep(min(site.sub$year):1980, each=12),
                                month= rep(1:12,1981-min(site.sub$year))) %>%
                                  merge(climate_normal,by=("month"))
      climate_in <- rbind(climate_int,climate_in)
    }
    
    
    modelled <- IPCCTier2SOMmodel(SiteData= site.sub,
                                    wth = climate_in,
                                    init.active = init_active,
                                    init.slow = init_slow,
                                    init.passive = init_passive,
                                  params=parameters)
    
    actuals <-  stocks_data %>%
      select(site, year,  actual = modelled_SOC)
    
    model_actual <- modelled %>%
      merge(actuals, by=c("site", "year")) %>%
      select(site, year, soc_total, actual) %>%
      filter(!is.na(actual))
    
    loglike <- loglik(model_actual$soc_total, model_actual$actual)
    
    out.df=data.frame(parameters,
                      site=site,
                      loglike=loglike)
    
    
    if (i==1){final.df <- out.df} else {final.df=rbind(final.df,out.df)}
    }
    return(final.df)
  }
  
  stopCluster(cl)
  
  Lkhood1 <- df_par %>%
    mutate(loglik = ifelse(loglike == -Inf, NA, loglike)) %>%
    group_by(id=SampleID) %>%
    summarise(loglik = mean(loglik, na.rm=T))
  
  si_obj2_llkhd <- sensitivity::tell(x = si_obj2, y = Lkhood1$loglik)
  
  
  if(method == "soboljansen"){
    # Calculate First-order and Total global sensitivity indices
    singleSI <- si_obj2_llkhd$S %>%
      select(SI = original, lci = `min. c.i.`, uci = `max. c.i.`) %>%
      rownames_to_column("params") %>%
      mutate(type="single")
    
    totalSI <- si_obj2_llkhd$T %>%
      select(SI = original, lci = `min. c.i.`, uci = `max. c.i.`) %>%
      rownames_to_column("params") %>%
      mutate(type="Total")
    
    combined_si <- rbind(singleSI,totalSI)
    
    return_si <- data.frame(combined_si)
  }
  
  if(method == "fast99"){
    return_si <- tibble(main = si_obj2_llkhd$D1 / si_obj2_llkhd$V,
                        interactions = 1 - si_obj2_llkhd$Dt / si_obj2_llkhd$V) %>%
      mutate(params =  varSI, .before=main)
  }
  
  return(return_si)
}



