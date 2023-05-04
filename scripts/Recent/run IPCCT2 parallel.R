## Authors: 
## Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca), 
## Arumugam Thiagarajan (arumugam.thiagarajan@ec.gc.ca)
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023

#=================================================================================
# Sampling Importance Resampling algorithm with Latin hypercube sampling in R
#=================================================================================

run.IPCCT2.parralel <- function(main.dir,
                           parameter.df,
                           evaluation.df,
                           stocks.df,
                           RunMap,
                           climate_data) {

RunIndex.list=unique(evaluation.df$RunIndex)
parameters <- parameter.df[1,]

n_clust=detectCores()
cl <- makeCluster(n_clust-4)
registerDoParallel(cl)

df_par <- foreach(RunIndex=RunIndex.list, .combine=rbind,  .inorder=TRUE) %dopar% {
  
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
  
    RunIndex.sub=evaluation.df[evaluation.df$RunIndex==RunIndex,] %>%
      arrange(year) %>%
      distinct() %>%
      mutate(Diff = year -lag(year)) %>%
      replace_na(list(Diff=0))
    
    ## Address discontinuity in data (mostly due to crop rotation) ##!! Something to improve here !!##
    if (any(RunIndex.sub$Diff>1)) {
      problems=RunIndex.sub[RunIndex.sub$Diff>1,]
      RunIndex.sub <- RunIndex.sub %>% 
        slice(1:min(as.numeric(row.names(problems))-1))
    }
    
    init.df <- RunMap[RunMap$RunIndex==RunIndex,]
    
    init_active = init.df[1,]$init_active
    init_slow = init.df[1,]$init_slow
    init_passive = init.df[1,]$init_passive
    
    climate_in <- climate_data[climate_data$POLYID==init.df [1,]$POLYID,]
    
    climate_normal <- climate_in %>% group_by(month) %>% summarise_all(mean) %>% select(month,tavg,mappet,irrig)
    
    if (min(RunIndex.sub$year)<1981){
      climate_int <- data.frame(POLYID=RunIndex.sub[1,]$POLYID,
                                year=rep(min(RunIndex.sub$year):1980, each=12),
                                month= rep(1:12,1981-min(RunIndex.sub$year))) %>%
        merge(climate_normal,by=("month"))
      climate_in <- rbind(climate_int,climate_in)
    }
    
    
    modelled <- IPCCTier2SOMmodel(SiteData= RunIndex.sub,
                                  wth = climate_in,
                                  init.active = init_active,
                                  init.slow = init_slow,
                                  init.passive = init_passive,
                                  params=parameters)
    
    if (init.df [1,]$RunBy=="treatment"){
      actuals <-  stocks.df %>%
        select(TrtID_Final, year,  actual = modelled_SOC)
      
      model_actual <- modelled %>%
        rename(TrtID_Final=site) %>%
        merge(actuals, by=c("TrtID_Final", "year")) %>%
        select(TrtID_Final, year, soc_total, actual) %>%
        filter(!is.na(actual))
    } 
    if (init.df [1,]$RunBy=="stock"){
      actuals <-  stocks.df %>%
        select(site, year,  actual = modelled_SOC)
      
      model_actual <- modelled %>%
        merge(actuals, by=c("site", "year")) %>%
        select(site, year, soc_total, actual) %>%
        filter(!is.na(actual))
    }
    
    final.df=data.frame(model_actual)
    return(final.df)
}
return(df_par)
}



