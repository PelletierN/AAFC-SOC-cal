library(imputeTS)

#key=read.csv("C://Users//pelletierni//OneDrive - AGR-AGR//Documents//SOC modelling//soc-modelling-framework-NP//data//input data//c_input//exp_location_SLC_linkage.csv")
yield=read.csv("C://Users//pelletierni//OneDrive - AGR-AGR//Documents//SOC modelling//soc-modelling-framework-NP//data//input data//c_input//slyield_merged_final.csv")

polyids_all<- read.csv("C://Users//pelletierni//OneDrive - AGR-AGR//Documents//SOC modelling//soc-modelling-framework-NP//data//input data//c_input//exp_location_SLC_linkage.csv") %>%
  select(Exp_ID,POLYID)

cinput.df <- read.csv("C://Users//pelletierni//OneDrive - AGR-AGR//Documents//SOC modelling//soc-modelling-framework-NP//data//input data//c_input//cninput.csv")

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
  int.df<- left_join(int.df,sub.df,by="year") 
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


test <- ipcct2.out[!complete.cases(ipcct2.out), ]

test2 <- test %>% filter(is.na(c_input_int))
  mutate(c_input_wo_ex_int = ifelse(any(!is.na(c_input_wo_ex))&n_obs>2,na_seasplit(find_frequency=TRUE,c_input_wo_ex),NA))

  test2 <- test %>% filter (site=="E40_CT")
  test2 <- annual_crops %>% filter (site=="E4_119")
#### To integrate
x=test.df$c_input_wo_ex

  
#all_sites <- unique(df$site)
#for (i in all_sites){
#sub<-df[df$site==i,]
    if (nrow(sub[!is.na(sub$yield),])>=10){sub$yield_int=na_seasplit(find_frequency=TRUE,sub$yield)} 
      else if (nrow(sub[!is.na(sub$pred_yield),])>=10) {sub$pred_yield_int=na_seasplit(find_frequency=TRUE,sub$pred_yield)}
      else if (nrow(sub[!is.na(sub$yield),])<10 & nrow(sub[!is.na(sub$yield),])>1) {sub$yield_int=na_interpolation(sub$yield)}
      else if (nrow(sub[!is.na(sub$pred_yield),])<10 & nrow(sub[!is.na(sub$pred_yield),])>1) {sub$pred_yield_int=na_interpolation(sub$pred_yield)}
        else {print(paste0("problem with site --> ",i))}
    if (i==all_sites[1]){out.df=sub} else {out.df=rbind(out.df,sub)}
  }
  return(out.df$pred_yield)
}

out.df2 <- out.df %>% 
  mutate(pred_yield=pred_yield/1000) %>%
  mutate(pred_yield_int=pred_yield_int/1000) %>%
  mutate(yield_final=ifelse(!is.na(yield),yield,NA)) %>%
  mutate(yield_final=ifelse(is.na(yield) & !is.na(yield_int),yield_int,yield_final)) %>%
  mutate(yield_final=ifelse(is.na(yield) & is.na(yield_int) & !is.na(pred_yield),pred_yield,yield_final)) %>%
  mutate(yield_final=ifelse(is.na(yield) & is.na(yield_int) & is.na(pred_yield) & !is.na(pred_yield_int),pred_yield_int,yield_final))
           
foo <- out.df2 %>% filter(is.na(yield_final))

ggplot(data=out.df2,aes(x=yield_final,y=cinput))+
  geom_point()+
  facet_wrap(~AAFC_code,scales="free")

out.df2$c_perc=out.df2$cinput/out.df2$yield_final

ggplot(data=out.df2[!is.na(out.df2$c_perc),],aes(y=c_perc,x=as.factor(AAFC_code)))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle=90))+
  ylab(expression(Yield~to~C_input~ratio))+
  xlab("AAFC code")

rm(out.df)
i="E4_102W"


## Adapt AAFC Code 
interpolated.df <- interpolated.df %>% 
  mutate(AAFC_code_adapted = ifelse(grepl("ALFALFA", AAFC_code),"ALFALFA",AAFC_code))








table(interpolated.df$AAFC_code_adapted)


unique(interpolated.df$Exp_ID)
dif.df=setdiff(out.df,out.df2)
table(dif.df$Exp_ID)
table(dif.df$site)

clipr::write_clip(table(site_data$TrtID,site_data$year_name))

########### Q/A for treatments imbalance

summary <- site_data %>%
  group_by(Exp_ID) %>%
  mutate(depth_ID=paste0(soil_depth_min_cm,"-",soil_depth_max_cm)) %>%
  summarise(length=max(year_name)-min(year_name)+1,
            n_treatments = length(unique(TrtID)),
            n_depths = length(unique(depth_ID)),
            n_replicates = length(unique(replication_number)),
            count=n()) %>%
  mutate(expected=length*n_treatments*n_depths*n_replicates) %>%
  filter(count!=expected)

clipr::write_clip(summary)

foo<- site_data %>% filter(Exp_ID=="E8")

table(foo$crop_rotation,foo$tillage)
table(foo$crop_rotation_code,foo$year_name,foo$tillage)
table(foo$crop_rotation,foo$replication_number)
table(foo$treatment_name,foo$replication_number)
table(foo$crop_rotation_code,foo$year_name,foo$treatment_name)


summary <- interpolated.df%>%
  group_by(Exp_ID) %>%
  summarise(length=max(year)-min(year)+1,
            n_treatments = length(unique(site)),
            count=n()) %>%
  mutate(expected=length*n_treatments) %>%
  filter(count!=expected)

clipr::write_clip(summary)

foo<- site_data %>% filter(Exp_ID=="E4")

table(foo$crop_rotation,foo$tillage)
table(foo$crop_rotation,foo$replication_number)
table(foo$treatment_name,foo$replication_number)
table(foo$treatment_name,foo$year_name)
########### Old stuff



not_interpolatable <- test %>%
  group_by(site) %>%
  mutate(flag=ifelse(all(is.na(pred_yield))&all(is.na(yield))&all(is.na(cinput)),1,0)) %>%
  filter(flag==1)

interpolatable <- test %>%
  group_by(site) %>%
  mutate(flag1=ifelse(all(is.na(pred_yield))&all(is.na(yield))&all(is.na(cinput)),1,0)) %>%
  filter(flag1==0)

big_problems <- problems %>% filter(site %in% unique(not_interpolatable$site))
table(big_problems$AAFC_code)


ggplot(data=test,aes(x=pred_yield/1000,y=yield))+
  geom_point()+
  facet_wrap(~AAFC_code,scales = "free")+ 
  geom_smooth(method="lm")+
  geom_abline(slope=1,intercept = 0,color="red")+
  ylab(expression(measured~yield~(T~ha^-1)))+
  xlab(expression(predicted~yield~(T~ha^-1)))+
  theme_bw()

ggplot(data=test[test$POLYID==425023,],aes(x=pred_yield/1000,y=cinput))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)+
  facet_wrap(~AAFC_code,scales = "free")+ 
  geom_smooth(method="lm")+
ylab(expression(measured~C~input~(T~ha^-1)))+
  xlab(expression(predicted~yield~(T~ha^-1)))+
  theme_bw()

i=unique(test$site)[1]
rm(out.df)


for (i in unique(test$site)){
  sub.df=test[test$site==i,]
  if (any(!is.na(sub.df$yield)) & any(is.na(sub.df$yield))){
    
    if (!exists("out.df")) {out.df=sub.df} else {out.df=rbind(out.df,sub.df)}
  }
}

library(imputeTS)
out.df <- out.df %>% group_by(site,Exp_ID) %>% mutate(yield2=na_interpolation(yield))

plot1<-ggplot(data=out.df,aes(x=year))+
  geom_point(aes(y=yield2),color="grey",size=2)+
  geom_point(aes(y=yield),color="black",size=1)+
  facet_wrap(~site,scales="free")+
  theme_minimal()

test2 <- test %>% group_by(site,Exp_ID) %>% mutate(yield2=ifelse(site %in% unique(out.df$site),na_interpolation(yield),yield))

plot2<- ggplot(data=test2[test2$site %in% unique(out.df$site),],aes(x=year))+
  geom_point(aes(y=yield2),color="grey",size=2)+
  geom_point(aes(y=yield),color="black",size=1)+
  facet_wrap(~site,scales="free")


### Interpolated dataframe

i=unique(test2$site)[1]
rm(out.df)

for (i in unique(test2$site)){
  sub.df=test[test$site==i,]
  int.df=data.frame(year=seq(min(sub.df$year),max(sub.df$year),1))
  int.df<- merge(int.df,sub.df,by="year")
  if (!exists("out.df")) {out.df=int.df} else {out.df=rbind(out.df,int.df)}
  
  # if (nrow(sub.df)!=nrow(int.df)){
  #       if (!exists("out.df")) {out.df=int.df} else {out.df=rbind(out.df,int.df)}
  # } 
}

unique(test2$Exp_ID)


ggplot(data=out.df,aes(x=year,y=imputed_yield))+geom_point()+geom_line()+facet_wrap(~site,scales="free")
ggplot(data=imputed.df,aes(x=year,y=imputed_yield))+geom_point()+geom_line()+facet_wrap(~site,scales="free")









out.df$na_flag=ifelse(is.na(out.df$yield),"na","value")
  
table(out.df$site,out.df$na_flag)

i=unique(out.df$site)[1]
for (i in unique(out.df$site)){
  sub.df=test[test$site==i,]
  if (any(!is.na(sub.df$cinput)) & any(is.na(sub.df$cinput))){
    if (!exists("out.df")) {out.df=sub.df} else {out.df=rbind(out.df,sub.df)}
  }
}
