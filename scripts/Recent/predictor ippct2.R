predictor.ipcct2<- function(df){
for (i in 1:nrow(df)){
  RunBy=df[i,]$RunBy
  
    if (RunBy=="stock"){
      
      site.on <- df[i,]$UniqueStockID
      label <- site.on
      site.df=subset(site_data_ipcct2,site_data_ipcct2$site==site.on)
      full.list=seq(min(site.df$year),max(site.df$year))
      data.list=site.df$year
      interp=FALSE
      
      if (length(setdiff(full.list,data.list))>0){
        print(paste0("interpolation for stock ",site.on," (index ",i,"): missing ",setdiff(full.list,data.list)))
        interp=TRUE
        }
      if (any(duplicated(data.list))){
        print(paste0("error for stock ",site.on," (index ",i,"): double value for ",data.list[duplicated(data.list)]))
        interp=TRUE
        }
    }
  
    if (RunBy=="treatment"){
      trt.on <- df[i,]$TrtID_Final
      label <- trt.on
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
      full.list=seq(min(site.df$year),max(site.df$year))
      data.list=site.df$year
      interp=FALSE
      
      if (length(setdiff(full.list,data.list))>0){
        print(paste0("interpolation for treatment ",trt.on," (index ",i,"): missing ",setdiff(full.list,data.list)))
        interp=TRUE
        }
      if (any(duplicated(data.list))){
        print(paste0("error for treatment ",trt.on," (index ",i,"): double value for ",data.list[duplicated(data.list)]))
        interp=TRUE
        }
    }
  if (interp==TRUE){
    
    foo<- data.frame(year=seq(min(site.df$year),max(site.df$year),1))
    
    site.df <- foo %>% 
      left_join(site.df) %>%
      mutate(site=ifelse(is.na(site),"interim",site)) %>%
      mutate(cinput=na_interpolation(cinput)) %>%
      fill(c(TrtID_Final,POLYID,sand,ligfrac,nfrac,till,irrig) ,.direction="down")
  }
  
  site.df$RunIndex=paste0(RunBy,"_",label)
  
  if(i==1){return.df=site.df} else{return.df=rbind(return.df,site.df)}
  
  }
return(return.df)
}

