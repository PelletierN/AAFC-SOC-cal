site_data_raw=site_data
cinput.df=cinput
yield.df=yield
poly.df=polyids_all
save.interim.df=FALSE
interim.df.file=NA

data.prep.ipcct2 <- function (site_data_raw,
                              cinput.df,
                              yield.df,
                              poly.df,
                              save.interim.df=FALSE,
                              interim.df.file=NA) {
  
site_data_ipcct2 <- site_data_raw %>%
  merge(poly.df, by=c("Exp_ID")) %>%
  mutate(tillage = ifelse(is.na(tillage), "CT", tillage)) %>%
  mutate(tillage = ifelse(tillage == "CT", "FT", tillage)) %>%
  mutate(
    site = UniqueStockID,
    POLYID = POLYID,
    year = year_name,
    sand = sand_px/100,
    till = tillage,
    irrig=0) %>%
  rowwise() %>% mutate(
    cinput = sum(as.numeric(as.character(crop_residue_kgha))*0.001, as.numeric(as.character(roots_residue_kgha))*0.001, as.numeric(as.character(above_ground_residue_kgha))*0.001,na.rm=T),
    yield = sum(as.numeric(as.character(grain_yield_kgha))*0.001, as.numeric(as.character(hay_yield_kgha))*0.001,na.rm=T)) %>%
  select(site,Exp_ID,POLYID,TrtID,year,sand,cinput,yield,till,irrig,AAFC_code) %>%
  group_by(site,POLYID,Exp_ID,TrtID,year,AAFC_code,till) %>% #Get mean values when multiple depth measurements are available
  summarise(sand=mean(sand,na.rm=T),
          cinput=mean(cinput,na.rm=T),
          yield=mean(yield,na.rm=T),
          irrig=mean(irrig,na.rm=T))

  foo <- yield.df %>%
  rename(POLYID=SL,
         year=YEAR) %>%
  select(-ECODISTRIC,-PROVINCE) %>%
  reshape2::melt(id.vars=c("POLYID","year")) %>%
  rename(AAFC_code=variable,
         pred_yield=value)

yield.df <- site_data_ipcct2 %>%
  left_join(foo,by=c("POLYID","year","AAFC_code")) %>% 
  mutate(yield=ifelse(yield==0,NA,yield),
         cinput=ifelse(cinput==0,NA,cinput))

all_sites <- unique(yield.df$site)

yield.df <- yield.df %>% 
  mutate(AAFC_code_adapted=AAFC_code) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="ALFALFA+BROME","ALFALFA",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="ALFALFA+CLOV","ALFALFA",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="ALFALFA+CWG+OHAYFD","ALFALFA",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="BARLEY+ALFALFA+BROME","BARLEY",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="BARLEY+CLOV","BARLEY",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="FALLOW","SUMMRF",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="SUMMRFG","SUMMRF",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="CSFG","SUMMRF", AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="GML","SUMMRF", AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="GRS","SUMMRF", AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="OAT+OHAYFD","OATS",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="OATS+PEA","OATS",AAFC_code_adapted)) %>%
  mutate(AAFC_code_adapted=ifelse(AAFC_code=="PEA","PEAS",AAFC_code_adapted))

test <- yield.df %>% left_join(cinput.df %>% rename(AAFC_code_adapted=crop_aafc),by="AAFC_code_adapted")

# Remove the treatments containing perennial crops
perenials=c("CLOVER","CLOV","CWG+HAY","ALFALFA","OHAYFD")

annual_crops <- test %>% group_by(TrtID) %>%
  mutate(contains_peren = ifelse(any(AAFC_code_adapted %in% perenials),TRUE,FALSE)) %>%
  ungroup %>%
  filter(contains_peren ==FALSE) %>%
  select(-contains_peren)

# Use available SLC yield table data to fill treatments with no yield data

annual_crops <-  annual_crops %>% group_by(TrtID) %>%
  mutate(some_yield = ifelse(any(!is.na(yield)),TRUE,FALSE)) %>%
  mutate(AdjProd_t=ifelse(is.na(yield) & some_yield==FALSE,pred_yield/1000,yield)) %>%
  ungroup()
  
## Calculate Harvest Index and Estimate C input
annual_crops <- annual_crops %>%
  mutate(HI=((AdjProd_t*(1-mc/100))*slope)+intercept) %>%
  mutate(AGR_t=ifelse(AdjProd_t>0.4, ((AdjProd_t*(1-mc/100)/HI)-(AdjProd_t*(1-mc/100)))*0.45, 0.15)) %>%
  mutate(BGR_t=((AdjProd_t*(1-mc/100))/HI)*0.45*rsr_wp) %>%
  mutate(c_input_wo_ex = AGR_t+BGR_t) %>%
  mutate(c_input_w_ex = AGR_t+BGR_t*1.65)

## Add constant value for summerfallow
annual_crops <- annual_crops %>%
  mutate(c_input_wo_ex=ifelse(is.na(c_input_wo_ex)&AAFC_code_adapted=="SUMMRF",0.25,c_input_wo_ex))

if (exists("interpolated.df")){rm(interpolated.df)}

for (i in unique(annual_crops$site)){
  sub.df=annual_crops[annual_crops$site==i,]
  int.df=data.frame(year=seq(min(sub.df$year),max(sub.df$year),1))
  int.df<- left_join(int.df,sub.df,by="year",multiple="first") 
  int.df2 <- int.df %>%
    mutate(site = first(site),
           POLYID = first(POLYID),
           Exp_ID = first(Exp_ID)) %>%
    fill(till,.direction="down") %>%
    fill(sand,.direction="down") %>%
    fill(TrtID,.direction="down") %>%
    mutate(irrig = na_interpolation(irrig)) 
  if (!exists("interpolated.df")) {interpolated.df=int.df2} else {interpolated.df=rbind(interpolated.df,int.df2)}
}

## Interpolation of missing data

interpolate.cinput <- function(x) {
  if (length(na.omit(x))>=10) {out=na_seasplit(x,find_frequency=TRUE)} 
  else if (length(na.omit(x))>=3) {out=na_kalman(x)} 
  else {out=na_replace(x,fill=mean(x,na.rm=T))}
  return(out)
}

mean.data <- interpolated.df %>%
  group_by(Exp_ID) %>%
  summarise(mean_sand=mean(sand,na.rm=T),
            mean_ligfrac=mean(ligfrac,na.rm=T),
            mean_nfrac=mean(nfrac,na.rm=T)) %>%
  summarise(mean_sand=mean(mean_sand,na.rm=T),
            mean_ligfrac=mean(mean_ligfrac,na.rm=T),
            mean_nfrac=mean(mean_nfrac,na.rm=T))

interpolated.df <- interpolated.df %>%
  group_by(site) %>%
  mutate(c_input_int= interpolate.cinput(c_input_wo_ex),
         ligfrac_int = interpolate.cinput(ligfrac),
         nfrac_int = interpolate.cinput(nfrac)) %>%
  ungroup()

interpolated.df <- interpolated.df %>%
  mutate(sand=ifelse(is.na(sand)|is.nan(sand),mean.data[1,]$mean_sand,sand),
         ligfrac=ifelse(is.na(ligfrac)|is.nan(ligfrac),mean.data[1,]$mean_ligfrac,ligfrac),
         nfrac=ifelse(is.na(nfrac)|is.nan(nfrac),mean.data[1,]$mean_nfrac,nfrac))

ipcct2.out=data.frame(
  site=interpolated.df$site,
  TrtID=interpolated.df$TrtID,
  POLYID=interpolated.df$POLYID,
  year=interpolated.df$year,
  sand=interpolated.df$sand,
  cinput=interpolated.df$c_input_int,
  ligfrac=interpolated.df$ligfrac_int,
  nfrac=interpolated.df$nfrac_int,
  till=interpolated.df$till,
  irrig=interpolated.df$irrig) %>%
  mutate(sand=ifelse(is.na(sand)|is.nan(sand),mean.data[1,]$mean_sand,sand),
         ligfrac=ifelse(is.na(ligfrac)|is.nan(ligfrac),mean.data[1,]$mean_ligfrac,ligfrac),
         nfrac=ifelse(is.na(nfrac)|is.nan(nfrac),mean.data[1,]$mean_nfrac,nfrac))

if (save.interim.df==TRUE) { write.csv(interpolated.df,interim.df.file,row.names=FALSE)}

return (ipcct2.out)
}

################# Climate data ###################

data.prep.climate <- function (dir) {
  
  climate.dir <- dir
  # Climate data
  ## open and combine files
  files <- list.files(climate.dir, full.names = TRUE, pattern = "csv$")
  
  climate_data <- files %>%
    lapply(read.csv) %>% 
    bind_rows %>%
    rename(Tavg = Tmean,
           PREC = Precip,
           JulianDay = Julian)
  
  ## Add potential evapotransipration
  holos_calculate_pet <- function(meanDailyTemperature,
                                  solarRadiation,
                                  relativeHumidity){
    # This is a vectorized implementation of the Holos reference PET calculator.
    # https://github.com/holos-aafc/Holos/blob/main/H.Core/Calculators/Climate/EvapotranspirationCalculator.cs
    
    term1 = 0.013
    term2 = meanDailyTemperature / (meanDailyTemperature + 15)
    term3 = (23.8856 * solarRadiation) + 50
    term4 = 1 + ((50 - relativeHumidity) / 70)
    
    result <- ifelse(relativeHumidity >= 50,
                     term1 * term2 * term3,
                     term1 * term2 * term3 * term4)
    result <- ifelse(result < 0, 0, result)
    result <- ifelse(meanDailyTemperature <= 0, 0, result)
    
    return(result)
  }
  
  climate_data <- climate_data %>%
    mutate(PET = holos_calculate_pet(meanDailyTemperature=Tavg,
                                     solarRadiation=Rad,
                                     relativeHumidity=RH)) 
  
  # Format climate data and add irrigation to this dataset (required for IPCCT2 wth file)
  climate_data_ipcct2 <- climate_data %>%
    group_by(POLYID, year = Year, month = Month) %>%
    summarise(tavg = mean(Tavg),
              mappet = ifelse(is.infinite(sum(PREC) / sum(PET)),0,sum(PREC) / sum(PET))
    ) %>%
    mutate(irrig=0) # at this point we don't have any irrigation data, so set it to zero
  
  return(climate_data_ipcct2)
}
