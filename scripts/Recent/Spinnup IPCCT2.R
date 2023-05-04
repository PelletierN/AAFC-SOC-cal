## Authors: 
## Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca), 
## Arumugam Thiagarajan (arumugam.thiagarajan@ec.gc.ca)
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023
## modified may 2023 to run with the treatment run map

initiate.ipcc2 <- function (index.df,params,climate_data_ipcct2){
   
  poly.index=index.df
  
   for (i in 1:nrow(poly.index)){
    
  RunBy=as.character(poly.index[i,]$RunBy)
  initial_year=poly.index[i,]$first.year
  initial_soc=poly.index[i,]$initial.SOC.impute
  POLYID.in=as.character(poly.index[i,]$POLYID)
  
  if (RunBy=="treatment"){
    trt.on <- poly.index[i,]$TrtID_Final
    label=trt.on
    site.df=subset(site_data_ipcct2,site_data_ipcct2$TrtID_Final==trt.on) %>%
      group_by(year) %>%
      summarise(site=first(site),
                TrtID_Final=first(TrtID_Final),
                POLYID=first(POLYID),
                sand=first(sand),
                cinput=mean(cinput,na.rm=T),
                ligfrac=first(ligfrac),
                nfrac=first(nfrac),
                till=first(till),
                irrig=first(irrig))
  }
  if (RunBy=="stock"){
    site.on <- poly.index[i,]$UniqueStockID
    label=site.on
    site.df=subset(site_data_ipcct2,site_data_ipcct2$site==site.on) 
  }
  
  mean_c_input <- mean(site.df$cinput,na.rm=T)
  
  site.df <- site.df%>%
    arrange(year) %>%
    slice(1:10) %>%
    summarise(site=first(site),
              TrtID_Final=first(TrtID_Final),
              POLYID=first(POLYID),
              year=first(year),
              sand=first(sand),
              cinput=mean(cinput,na.rm=T),
              ligfrac=first(ligfrac),
              nfrac=first(nfrac),
              till=first(till),
              irrig=first(irrig))

  

    
  #Extract mean Cinput from first ten years of data

  if (initial_year>=1981+9){ # Process if we have 10 years of climate data to initiate the model
  climate.df=subset(climate_data_ipcct2,climate_data_ipcct2$POLYID==POLYID.in) %>%
    arrange(year,month) %>%
    filter(year<=initial_year) %>%
    filter(year>=initial_year-10) %>%
    group_by(month) %>%
    summarise(year=initial_year,
              tavg = mean(tavg),
              mappet = mean(mappet),
              irrig = mean(irrig))
  
  } else { # Process if experiment started too early and we don't have climate data (use climate normal)
    climate.df=subset(climate_data_ipcct2,climate_data_ipcct2$POLYID==POLYID.in) %>%
      arrange(year,month) %>%
      #filter(year<=1991) %>%
      group_by(month) %>%
      summarise(year=initial_year,
                tavg = mean(tavg),
                mappet = mean(mappet),
                irrig = mean(irrig))
  }
  
    
  # Calculate pools fraction with steady-state model by running the first 10 years of data

  ss <- IPCCTier2SOMmodel(SiteData = site.df,
                          wth = climate.df,
                          init.active = 0,
                          init.slow = 0,
                          init.passive = 0,
                          params)
  
    # get fraction in each pool
    ss.frac <- data.frame(init_active_frac=ss[1,]$soc_active/ss[1,]$soc_total,
                        init_slow_frac=ss[1,]$soc_slow/ss[1,]$soc_total,
                        init_passive_frac=ss[1,]$soc_passive/ss[1,]$soc_total)
    

  # calculate steady state pool size
    poly.index[i,]$init_active <- ss.frac[1,]$init_active_frac*initial_soc
    poly.index[i,]$init_slow <- ss.frac[1,]$init_slow_frac*initial_soc
    poly.index[i,]$init_passive <- ss.frac[1,]$init_passive_frac*initial_soc
    
    rm(climate.df)
    
    poly.index[i,]$mean_c_input=mean_c_input
    poly.index[i,]$RunIndex=paste0(RunBy,"_",label)

  }
   return(poly.index)
}

