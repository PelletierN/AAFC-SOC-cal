## Authors: 
## Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca), 
## Arumugam Thiagarajan (arumugam.thiagarajan@ec.gc.ca)
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023

#=================================================================================
# Sampling Importance Resampling algorithm with Latin hypercube sampling in R
#=================================================================================

run.SIR.IPCCT2 <- function(main.dir,
                           parameter.df,
                           stocks.df,
                           calibration.df,
                           climate_data,
                           initial_c,
                           sample_size,
                           resample_size) {
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

site.list=unique(calibration.df$site)

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
  
  source(paste0(main.dir,"scripts/Recent/ipcct2_SOCTier2Model.r"))
  
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
  
for (i in 1:nrow(parameter.df)){
    site.sub=calibration.df[calibration.df$site==site,] %>%
      arrange(year) %>%
      distinct() %>%
      mutate(Diff = year -lag(year)) %>%
      replace_na(list(Diff=0))
    
    ## Address discontinuity in data (mostly due to crop rotation) ##!! Something to improve here !!##
    if (any(site.sub$Diff>1)) {
      problems=site.sub[site.sub$Diff>1,]
      site.sub <- site.sub %>% 
        slice(1:min(as.numeric(row.names(problems))-1))
    }
    
    init_active = initial_c[initial_c$site==site,]$init_active
    init_slow = initial_c[initial_c$site==site,]$init_slow
    init_passive = initial_c[initial_c$site==site,]$init_passive
    
    parameters = parameter.df[i,]
    
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
    
    actuals <-  stocks.df %>%
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
  mutate(loglike = ifelse(loglike == -Inf, NA, loglike)) %>%
  group_by(SampleID) %>%
  summarise(loglike = mean(loglike, na.rm=T)) %>%
  mutate(weights = exp(loglike)/sum(exp(loglike)))

sampIndx <- sample(1:nrow(Lkhood1), 
                   size = SIR_resample_size,
                   replace = FALSE,
                   prob = Lkhood1$weights)

PostTheta <- as.data.frame(parameter.df[sampIndx,])
return(PostTheta)
}


