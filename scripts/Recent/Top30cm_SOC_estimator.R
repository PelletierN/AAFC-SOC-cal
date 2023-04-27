## Authors: 
## Nicolas Pelletier (nicolas.pelletier@agr.gc.ca)

## Date: April 2023

get.stocks <- function(site_data) {
  stocks.df <- site_data %>%
    mutate(
      site = UniqueStockID,
      year = year_name,
      depthID = paste0(soil_depth_min_cm,"-",soil_depth_max_cm)) %>%
    select(site,year,Exp_ID,TrtID,soc_tha,soil_depth_min_cm,soil_depth_max_cm,depthID) %>%
    distinct() %>%
    arrange(site,year) %>%
    mutate(soc_tha=as.numeric(soc_tha)) 
  
  missing.SOC.list <- stocks.df %>% 
    group_by(site) %>% 
    summarise(mean_soc=mean(soc_tha,na.rm=T)) %>% 
    filter(is.nan(mean_soc)) %>% 
    select(site)
  if(nrow(missing.SOC.list)>0){print(paste0("no C stock measurements available for StockID: ",missing.SOC.list$site))}
  
  # remove NAs
  stocks.df <- stocks.df %>%
    filter(!is.na(soc_tha))

  return (stocks.df)
}


obtain.dsm.correction <- function(stocks.df,dsm_data,merge.key) {

# Join the DSM data to the stocks
dsm=merge(stocks.df,dsm_data,by=merge.key,left=T) %>%
  select(site,Can_SOC30) %>% 
  group_by(site) %>%
  summarise(Can_SOC30=mean(Can_SOC30))

## Get correction factors for top 30 cm for each Exp_ID and create QA flag
correct.fac <- stocks.df %>%
  merge(dsm,by=c("site"),left=T) %>%
  mutate(dif=soc_tha-Can_SOC30,min=soil_depth_min_cm,max=soil_depth_max_cm) %>%
  group_by(site,depthID) %>%
  arrange(year) %>%
  slice(1) %>%
  group_by(depthID,Exp_ID,min,max) %>%
  summarise(SOC_correction_mean=mean(dif,na.rm=T)*-1,SOC_correction_median=median(dif,na.rm=T)*-1) %>%
  filter(!is.na(SOC_correction_mean)) %>%
  mutate(QA_flag1=case_when(min==0 & max<30 & SOC_correction_mean>0 ~ 1,
                            min==0 & max<30 & SOC_correction_mean<0 ~ 0,
                            min==0 & max>30 & SOC_correction_mean<0 ~ 1,
                            min==0 & max>30 & SOC_correction_mean>0 ~ 0)) 

return(correct.fac)
}


#test 12-12
# Check the counts of depth measurements
#stocks.df$depthID=paste0(stocks.df$soil_depth_min_cm,"-",stocks.df$soil_depth_max_cm)
#table(stocks.df$depthID)

# Check the counts of depth measurements by experiment
#table(stocks.df$depthID,stocks.df$Exp_ID)

# For debugging...
#df=stocks.df
#site=unique(df$site)[63]#[211] #
#year=1999
#depth_cm=30

# Function to estimate 30 cm SOC at all sites and all dates
top.SOC.estimate <- function(df,correct.fac,depth_cm){
  df <- df[!is.na(df$soc_tha),]
  df <- df[df$soc_tha!=0,]
  
  for (site in unique(df$site)){
    sub.df=df[df$site==site,]
    for (year in unique(sub.df$year)){
      sub.df2=sub.df[sub.df$year==year,] 
      type_string=unique(sub.df2$depthID) %>% sort %>% paste(collapse = " + ")
      TrtID=unique(sub.df2$TrtID)
      
      res.df=data.frame(depth.min=seq(0,depth_cm-0.5,0.5),
                        depth.max=seq(0.5,depth_cm,0.5),
                        measured_c=NA)
      
      
      if (sub.df2[1,]$soil_depth_min_cm==0 & sub.df2[1,]$soil_depth_max_cm>depth_cm | # Catch cases with measurements exceeding the target depth increment
          max(sub.df2$soil_depth_max_cm)<depth_cm) { # Catch cases with measurements short of the target depth increment
        
        cor.sub=correct.fac %>% filter(Exp_ID==sub.df2[1,]$Exp_ID, min==sub.df2[1,]$soil_depth_min_cm, max==sub.df2[1,]$soil_depth_max_cm)
        if (is.na(cor.sub[1,]$QA_flag1)) {
          print(paste0("missing DSM data for Exp_ID = ",sub.df[1,]$Exp_ID))
          modelled_SOC= sub.df2[1,]$soc_tha
        } else if (cor.sub[1,]$QA_flag1==1) {modelled_SOC= sub.df2[1,]$soc_tha + cor.sub[1,]$SOC_correction_mean} else {modelled_SOC= sub.df2[1,]$soc_tha}
        
      } else {
      for (row in 1:nrow(sub.df2)){
        data=sub.df2[row,]
        res.df[,row+3]=ifelse(res.df$depth.min>=data$soil_depth_min_cm & res.df$depth.max<=data$soil_depth_max_cm,
                              data$soc_tha/((data$soil_depth_max_cm-data$soil_depth_min_cm)*2),
                              NA)
      
    res.df <- replace(res.df, res.df==0, NA)
    
    res.df$measured_c=rowMeans(data.frame(res.df[,4:length(res.df)]),na.rm=T)
    
    # Filling missing values in top / bottom section(s)
    # Can we use the DSM estimation instead?
    res.df <- res.df %>%
      fill(measured_c, .direction = "down")%>%
      fill(measured_c, .direction = "up")
    
    modelled_SOC=sum(res.df$measured_c)
      }}
    out=data.frame(site,
                   experiment=sub("\\_.*", "", site),
                   TrtID,
                   year,
                   input_type=type_string,
                   modelled_SOC)
    
    if (site==unique(stocks.df$site)[1] & year==unique(sub.df$year)[1]){final=out } else {final=rbind(final,out)}
      
    }
  }
  return(final)
}


top.SOC.raw <- function(df){
  df <- df[!is.na(df$soc_tha),]
  df <- df[df$soc_tha!=0,]
  
  for (site in unique(df$site)){
    sub.df=df[df$site==site,]
    for (year in unique(sub.df$year)){
      sub.df2=sub.df[sub.df$year==year,] 
      type_string=unique(sub.df2$depthID) %>% sort %>% paste(collapse = " + ")
      TrtID=unique(sub.df2$TrtID)
      if (length(TrtID)>1){print(paste0("error: more than one treatment ID for stockID --> ",site))}
      d.max=max(sub.df2$soil_depth_max_cm)
      res.df=data.frame(depth.min=seq(0,d.max-0.5,0.5),
                        depth.max=seq(0.5,d.max,0.5),
                        measured_c=NA)
      
      for (row in 1:nrow(sub.df2)){
        data=sub.df2[row,]
        res.df[,row+3]=ifelse(res.df$depth.min>=data$soil_depth_min_cm & res.df$depth.max<=data$soil_depth_max_cm,
                              data$soc_tha/((data$soil_depth_max_cm-data$soil_depth_min_cm)*2),
                              NA)
        
        res.df <- replace(res.df, res.df==0, NA)
        
        res.df$measured_c=rowMeans(data.frame(res.df[,4:length(res.df)]),na.rm=T)
        
        # Filling missing values in top / bottom section(s)
        # Can we use the DSM estimation instead?
        #res.df <- res.df %>%
        #fill(measured_c, .direction = "down")%>%
        #fill(measured_c, .direction = "up")
        
        modelled_SOC=sum(res.df$measured_c)
      }
      out=data.frame(site,
                     experiment=sub("\\_.*", "", site),
                     TrtID,
                     year,
                     input_type=type_string,
                     modelled_SOC)
      
      if (site==unique(stocks.df$site)[1] & year==unique(sub.df$year)[1]){final=out } else {final=rbind(final,out)}
      
    }
  }
  return(final)
}

top.SOC.naive <- function(df,depth_cm){
  df <- df[!is.na(df$soc_tha),]
  df <- df[df$soc_tha!=0,]
  
  for (site in unique(df$site)){
    sub.df=df[df$site==site,]
    for (year in unique(sub.df$year)){
      sub.df2=sub.df[sub.df$year==year,] 
      type_string=unique(sub.df2$depthID) %>% sort %>% paste(collapse = " + ")
      TrtID=unique(sub.df2$TrtID)
      if (length(TrtID)>1){print(paste0("error: more than one treatment ID for stockID --> ",site))}
      d.max=max(sub.df2$soil_depth_max_cm)
      res.df=data.frame(depth.min=seq(0,depth_cm-0.5,0.5),
                        depth.max=seq(0.5,depth_cm,0.5),
                        measured_c=NA)
      
      for (row in 1:nrow(sub.df2)){
        data=sub.df2[row,]
        res.df[,row+3]=ifelse(res.df$depth.min>=data$soil_depth_min_cm & res.df$depth.max<=data$soil_depth_max_cm,
                              data$soc_tha/((data$soil_depth_max_cm-data$soil_depth_min_cm)*2),
                              NA)
        
        res.df <- replace(res.df, res.df==0, NA)
        
        res.df$measured_c=rowMeans(data.frame(res.df[,4:length(res.df)]),na.rm=T)
        
        # Filling missing values in top / bottom section(s)
        res.df <- res.df %>%
        fill(measured_c, .direction = "down")%>%
        fill(measured_c, .direction = "up")
        
        modelled_SOC=sum(res.df$measured_c)
      }
      out=data.frame(site,
                     experiment=sub("\\_.*", "", site),
                     TrtID,
                     year,
                     input_type=type_string,
                     modelled_SOC)
      
      if (site==unique(stocks.df$site)[1] & year==unique(sub.df$year)[1]){final=out } else {final=rbind(final,out)}
      
    }
  }
  return(final)
}
