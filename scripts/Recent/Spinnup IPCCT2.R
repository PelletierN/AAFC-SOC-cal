## Authors: 
## Francis Durnin-Vermette (francis.durnin-vermette@agr.gc.ca), 
## Arumugam Thiagarajan (arumugam.thiagarajan@ec.gc.ca)
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023

initiate.ipcc2<- function (index.df,params){
   poly.index=index.df
  for (i in 1:nrow(poly.index)){
  site.in=as.character(poly.index[i,]$site)
  POLYID.in=as.character(poly.index[i,]$POLYID)
  initial_year=poly.index[i,]$year
  initial_soc=poly.index[i,]$modelled_SOC
    
  #Extract mean Cinput from first ten years of data
  #!# probably need to change the first year of spin up to the first year with C stock measurement 
  site.df=subset(site_data_ipcct2,site_data_ipcct2$site==site.in) %>%
    arrange(year) %>%
    slice(1:10) %>%
    summarise(site=first(site),
              POLYID=first(POLYID),
              year=first(year),
              sand=first(sand),
              cinput=mean(cinput,na.rm=T),
              ligfrac=first(ligfrac),
              nfrac=first(nfrac),
              till=first(till),
              irrig=first(irrig))
  
  if (site.df[1,]$year>=1981+9){
  climate.df=subset(climate_data_ipcct2,climate_data_ipcct2$POLYID==POLYID.in) %>%
    arrange(year,month) %>%
    filter(year<=site.df[1,]$year) %>%
    filter(year>=site.df[1,]$year-10) %>%
    group_by(month) %>%
    summarise(year=site.df[1,]$year,
              tavg = mean(tavg),
              mappet = mean(mappet),
              irrig = mean(irrig))
  } else {
    climate.df=subset(climate_data_ipcct2,climate_data_ipcct2$POLYID==POLYID.in) %>%
      arrange(year,month) %>%
      filter(year<=1991) %>%
      group_by(month) %>%
      summarise(year=site.df[1,]$year,
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
    poly.index[i,]$init_active_frac <- ss.frac[[1]]*initial_soc
    poly.index[i,]$init_slow_frac <- ss.frac[[2]]*initial_soc
    poly.index[i,]$init_passive_frac <- ss.frac[[3]]*initial_soc

  }
   return(poly.index)
   }
