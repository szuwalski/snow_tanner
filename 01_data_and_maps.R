#=============================
# is there a relationship between tanner and snow crab distribution?
# is there a relationship between large snow crab and small tanner crab?
# what is the appropriate lag to consider?
# can the impact of temperature be simultaneously tested for?
# where do they overlap?
# are there thermal niches shared for different sized crab? (are lg snow and small Tanner cooccuring?)
# are there relationships between centroids or time series of different sized crab?

library(maps)
library("rnaturalearth")
library(interp)
library(RColorBrewer)
library(reshape2) # for melt
library(mgcv)  
library(PBSmapping)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggnewscale)
#==============================
# snow crab data
#==============================
survDAT<-read.csv("data/snow/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)

drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
SurvYR<-unique(drvYear)
AllStation<-unique(survDAT$GIS_STATION)

#==plot GIS stations
AllStnLoc<-matrix(ncol=2,nrow=length(AllStation))
for(w in 1:length(AllStation))
{
  temp<-survDAT[survDAT$GIS_STATION==AllStation[w],]
  AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
  AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
}

#========================
#==find survey densities (total)
#========================
nmiSurv<-140350
DensityM_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityF_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_55_65<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_65_75<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_78_100<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityMge45<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

num_M_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
num_M_101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

StationYr<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
station_bot_temp<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
for(y in 1:length(SurvYR))
{
  yrDAT<-survDAT[drvYear==SurvYR[y],]
  fileyr<-SurvYR[y]
  stationsUNQ<-(unique(yrDAT$GIS_STATION))
  #==density at station
  for(j in 1:length(stationsUNQ))
  {
    stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
    StationYr[y,j]<-as.character(stationsUNQ[j])
    Hauls<-(unique(stationALL$HAUL))
    female_45_55<-stationALL[stationALL$SEX==2 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_45_55<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_55_65<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>55 & stationALL$WIDTH<66,]
    male_65_75<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>65 & stationALL$WIDTH<76,]
    male_45_85<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<86,]
    male_78_100<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>78 & stationALL$WIDTH<101,]
    Males101<-stationALL[stationALL$SEX==1& stationALL$WIDTH>101,]
    male_ge45<-stationALL[stationALL$SEX==1& stationALL$WIDTH>44,]
    #==densities across hauls in crabs per km^2
    tempDensM45<-NULL
    tempDensF45<-NULL
    tempDensM55<-NULL
    tempDensM65<-NULL
    tempDensM85<-NULL
    tempDensM75101<-NULL
    tempDensM101<-NULL
    tempDensF<-NULL
    temp_num_45<-0
    temp_num_101<-0
    tempDensMge45<-NULL
    for(k in 1:length(Hauls))
    {
      SampFactM<-female_45_55$SAMPLING_FACTOR[which(female_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_45_55$AREA_SWEPT[which(female_45_55$HAUL==Hauls[k])[1]]
      tempDensF45<-length(female_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_55$SAMPLING_FACTOR[which(male_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_55$AREA_SWEPT[which(male_45_55$HAUL==Hauls[k])[1]]
      tempDensM45<-length(male_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_55_65$SAMPLING_FACTOR[which(male_55_65$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_55_65$AREA_SWEPT[which(male_55_65$HAUL==Hauls[k])[1]]
      tempDensM55<-length(male_55_65$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_65_75$SAMPLING_FACTOR[which(male_65_75$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_65_75$AREA_SWEPT[which(male_65_75$HAUL==Hauls[k])[1]]
      tempDensM65<-length(male_65_75$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_85$SAMPLING_FACTOR[which(male_45_85$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_85$AREA_SWEPT[which(male_45_85$HAUL==Hauls[k])[1]]
      tempDensM85<-length(male_45_85$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_45<-temp_num_45+length(Males101$HAUL==Hauls[k])*SampFactM

      SampFactM<-male_78_100$SAMPLING_FACTOR[which(male_78_100$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_78_100$AREA_SWEPT[which(male_78_100$HAUL==Hauls[k])[1]]
      tempDensM78101<-length(male_78_100$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_ge45$SAMPLING_FACTOR[which(male_ge45$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_ge45$AREA_SWEPT[which(male_ge45$HAUL==Hauls[k])[1]]
      tempDensMge45<-length(male_ge45$HAUL==Hauls[k])*SampFactM/AreaSweptM

      SampFactM<-Males101$SAMPLING_FACTOR[which(Males101$HAUL==Hauls[k])[1]] 
      AreaSweptM<-Males101$AREA_SWEPT[which(Males101$HAUL==Hauls[k])[1]]
      tempDensM101<-length(Males101$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_101<-temp_num_101+length(Males101$HAUL==Hauls[k])*SampFactM
    }
    DensityF_45_55[y,j]<-mean(tempDensF45)
    DensityM_45_55[y,j]<-mean(tempDensM45)
    DensityM_55_65[y,j]<-mean(tempDensM55)
    DensityM_65_75[y,j]<-mean(tempDensM65)
    DensityM_45_85[y,j]<-mean(tempDensM85)
    DensityM_78_100[y,j]<-mean(tempDensM78101)
    DensityMge45[y,j]<-mean(tempDensMge45)
    DensityM101[y,j]<-mean(tempDensM101)
    num_M_45_85[y,j]<-temp_num_45    
    num_M_101[y,j]<-temp_num_101
    
    station_bot_temp[y,j]<-stationALL$GEAR_TEMPERATURE[1]
  }
}
plot(apply(num_M_101,1,sum,na.rm=T)[-1]~SurvYR[-1],type='l',ylab='total number of observed crab >101mm')

snow_dat<-data.frame(log_abund_101=log(c(t(DensityM101[-1,]))),
                    fish_recruits=log(c(t(DensityM_78_100[-1,]))),
                    bot_temp=c(t(station_bot_temp[-nrow(station_bot_temp),])),
                    station=c(t(StationYr[-1,])),
                    lat=rep(AllStnLoc[,1],(nrow(DensityM_78_100)-1)),
                    lon=rep(AllStnLoc[,2],(nrow(DensityM_78_100)-1)),
                    year=rep(SurvYR[-1],each=ncol(DensityM101)),
                    log_abund_101=log(c(t(DensityM101[-1,]))),
                    log_abund_45=log(c(t(DensityM_45_55[-1,]))),
                    log_abund_45_f=log(c(t(DensityF_45_55[-1,]))),
                    log_abund_55=log(c(t(DensityM_55_65[-1,]))),
                    log_abund_65=log(c(t(DensityM_65_75[-1,]))),
                    log_abund_85=log(c(t(DensityM_45_85[-1,]))),
                    log_abund_ge45=log(c(t(DensityMge45[-1,]))),
                    loc="EBS",
                    obs_num_101=c(t(num_M_101[-1,])),
                    obs_num_45_85=c(t(num_M_45_85[-1,])))

for(x in 1:length(AllStation))
{
  snow_dat$lat[snow_dat$station==AllStation[x]]<-AllStnLoc[x,1]
  snow_dat$lon[snow_dat$station==AllStation[x]]<-AllStnLoc[x,2]
}

#==============================
# Tanner crab data
#==============================
survDAT<-read.csv("data/tanner/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)
drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
SurvYR<-unique(drvYear)
AllStation<-unique(survDAT$GIS_STATION)

#==plot GIS stations
AllStnLoc<-matrix(ncol=2,nrow=length(AllStation))
for(w in 1:length(AllStation))
{
  temp<-survDAT[survDAT$GIS_STATION==AllStation[w],]
  AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
  AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
}

#========================
#==find survey densities (total)
#========================
nmiSurv<-140350
DensityM_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityF_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_55_65<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_65_75<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_78_100<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityMge45<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

num_M_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
num_M_101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

StationYr<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
station_bot_temp<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
for(y in 1:length(SurvYR))
{
  yrDAT<-survDAT[drvYear==SurvYR[y],]
  fileyr<-SurvYR[y]
  stationsUNQ<-(unique(yrDAT$GIS_STATION))
  #==density at station
  for(j in 1:length(stationsUNQ))
  {
    stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
    StationYr[y,j]<-as.character(stationsUNQ[j])
    Hauls<-(unique(stationALL$HAUL))
    female_45_55<-stationALL[stationALL$SEX==2 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_45_55<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_55_65<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>55 & stationALL$WIDTH<66,]
    male_65_75<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>65 & stationALL$WIDTH<76,]
    male_45_85<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<86,]
    male_78_100<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>78 & stationALL$WIDTH<101,]
    Males101<-stationALL[stationALL$SEX==1& stationALL$WIDTH>101,]
    male_ge45<-stationALL[stationALL$SEX==1& stationALL$WIDTH>44,]
    #==densities across hauls in crabs per km^2
    tempDensM45<-NULL
    tempDensF45<-NULL
    tempDensM55<-NULL
    tempDensM65<-NULL
    tempDensM85<-NULL
    tempDensM75101<-NULL
    tempDensM101<-NULL
    tempDensF<-NULL
    temp_num_45<-0
    temp_num_101<-0
    tempDensMge45<-NULL
    for(k in 1:length(Hauls))
    {
      SampFactM<-female_45_55$SAMPLING_FACTOR[which(female_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_45_55$AREA_SWEPT[which(female_45_55$HAUL==Hauls[k])[1]]
      tempDensF45<-length(female_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_55$SAMPLING_FACTOR[which(male_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_55$AREA_SWEPT[which(male_45_55$HAUL==Hauls[k])[1]]
      tempDensM45<-length(male_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_55_65$SAMPLING_FACTOR[which(male_55_65$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_55_65$AREA_SWEPT[which(male_55_65$HAUL==Hauls[k])[1]]
      tempDensM55<-length(male_55_65$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_65_75$SAMPLING_FACTOR[which(male_65_75$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_65_75$AREA_SWEPT[which(male_65_75$HAUL==Hauls[k])[1]]
      tempDensM65<-length(male_65_75$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_85$SAMPLING_FACTOR[which(male_45_85$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_85$AREA_SWEPT[which(male_45_85$HAUL==Hauls[k])[1]]
      tempDensM85<-length(male_45_85$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_45<-temp_num_45+length(Males101$HAUL==Hauls[k])*SampFactM
      
      SampFactM<-male_78_100$SAMPLING_FACTOR[which(male_78_100$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_78_100$AREA_SWEPT[which(male_78_100$HAUL==Hauls[k])[1]]
      tempDensM78101<-length(male_78_100$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_ge45$SAMPLING_FACTOR[which(male_ge45$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_ge45$AREA_SWEPT[which(male_ge45$HAUL==Hauls[k])[1]]
      tempDensMge45<-length(male_ge45$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-Males101$SAMPLING_FACTOR[which(Males101$HAUL==Hauls[k])[1]] 
      AreaSweptM<-Males101$AREA_SWEPT[which(Males101$HAUL==Hauls[k])[1]]
      tempDensM101<-length(Males101$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_101<-temp_num_101+length(Males101$HAUL==Hauls[k])*SampFactM
    }
    DensityF_45_55[y,j]<-mean(tempDensF45)
    DensityM_45_55[y,j]<-mean(tempDensM45)
    DensityM_55_65[y,j]<-mean(tempDensM55)
    DensityM_65_75[y,j]<-mean(tempDensM65)
    DensityM_45_85[y,j]<-mean(tempDensM85)
    DensityM_78_100[y,j]<-mean(tempDensM78101)
    DensityMge45[y,j]<-mean(tempDensMge45)
    DensityM101[y,j]<-mean(tempDensM101)
    num_M_45_85[y,j]<-temp_num_45    
    num_M_101[y,j]<-temp_num_101
    
    station_bot_temp[y,j]<-stationALL$GEAR_TEMPERATURE[1]
  }
}
plot(apply(num_M_101,1,sum,na.rm=T)[-1]~SurvYR[-1],type='l',ylab='total number of observed crab >101mm')

tan_dat<-data.frame(log_abund_101=log(c(t(DensityM101[-1,]))),
                    fish_recruits=log(c(t(DensityM_78_100[-1,]))),
                    bot_temp=c(t(station_bot_temp[-nrow(station_bot_temp),])),
                    station=c(t(StationYr[-1,])),
                    lat=rep(AllStnLoc[,1],(nrow(DensityM_78_100)-1)),
                    lon=rep(AllStnLoc[,2],(nrow(DensityM_78_100)-1)),
                    year=rep(SurvYR[-1],each=ncol(DensityM101)),
                    log_abund_101=log(c(t(DensityM101[-1,]))),
                    log_abund_45=log(c(t(DensityM_45_55[-1,]))),
                    log_abund_45_f=log(c(t(DensityF_45_55[-1,]))),
                    log_abund_55=log(c(t(DensityM_55_65[-1,]))),
                    log_abund_65=log(c(t(DensityM_65_75[-1,]))),
                    log_abund_85=log(c(t(DensityM_45_85[-1,]))),
                    log_abund_ge45=log(c(t(DensityMge45[-1,]))),
                    loc="EBS",
                    obs_num_101=c(t(num_M_101[-1,])),
                    obs_num_45_85=c(t(num_M_45_85[-1,])))

for(x in 1:length(AllStation))
{
  tan_dat$lat[tan_dat$station==AllStation[x]]<-AllStnLoc[x,1]
  tan_dat$lon[tan_dat$station==AllStation[x]]<-AllStnLoc[x,2]
}



#===================================================================
# plot one map with all locations observed over time for all crab
#==================================================================
t_temp<-tan_dat[,c(3,5,6,7,14)]
s_temp<-snow_dat[,c(3,5,6,7,14)]
t_temp<-t_temp[!is.na(t_temp[,5]),]
s_temp<-s_temp[!is.na(s_temp[,5]),]
s_temp$species<-"Snow"
t_temp$species<-"Tanner"

world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<-min(in_dat$lon,na.rm=T)
lon_2<-max(in_dat$lon,na.rm=T)*.99
lat_1<-min(in_dat$lat,na.rm=T)
lat_2<-max(in_dat$lat,na.rm=T)
in_fin<-melt(s_temp,id=c('lat','lon','year','species'))
in_fin<-in_fin[complete.cases(in_fin),]

p<-ggplot() + 
  geom_tile(data=filter(in_fin,variable=='log_abund_ge45'), aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  facet_wrap(~year) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/snow_dist.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

in_fin<-melt(t_temp,id=c('lat','lon','year','species'))
in_fin<-in_fin[complete.cases(in_fin),]
q<-ggplot() + 
  geom_tile(data=filter(in_fin,variable=='log_abund_ge45'), aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  facet_wrap(~year) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/tan_dist.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#============================
# small male maps
t_temp<-tan_dat[,c(3,5,6,7,9)]
s_temp<-snow_dat[,c(3,5,6,7,9)]
t_temp<-t_temp[!is.na(t_temp[,5]),]
s_temp<-s_temp[!is.na(s_temp[,5]),]
s_temp$species<-"Snow"
t_temp$species<-"Tanner"
in_fin<-melt(s_temp,id=c('lat','lon','year','species'))
in_fin<-in_fin[complete.cases(in_fin),]

p<-ggplot() + 
  geom_tile(data=filter(in_fin,variable=='log_abund_45'), aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  facet_wrap(~year) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/sm_snow_dist.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

in_fin<-melt(t_temp,id=c('lat','lon','year','species'))
in_fin<-in_fin[complete.cases(in_fin),]
q<-ggplot() + 
  geom_tile(data=filter(in_fin,variable=='log_abund_45'), aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  facet_wrap(~year) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/sm_tan_dist.png",height=8,width=8,res=400,units='in')
print(q)
dev.off()

#============================
# large male maps
t_temp<-tan_dat[,c(1,3,5,6,7)]
s_temp<-snow_dat[,c(1,3,5,6,7)]
t_temp<-t_temp[!is.na(t_temp[,1]),]
s_temp<-s_temp[!is.na(s_temp[,1]),]
s_temp$species<-"Snow"
t_temp$species<-"Tanner"
in_fin<-melt(s_temp,id=c('lat','lon','year','species'))
in_fin<-in_fin[complete.cases(in_fin),]

p<-ggplot() + 
  geom_tile(data=filter(in_fin,variable=='log_abund_101'), aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  facet_wrap(~year) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/lg_snow_dist.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

in_fin<-melt(t_temp,id=c('lat','lon','year','species'))
in_fin<-in_fin[complete.cases(in_fin),]
q<-ggplot() + 
  geom_tile(data=filter(in_fin,variable=='log_abund_101'), aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  facet_wrap(~year) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/lg_tan_dist.png",height=8,width=8,res=400,units='in')
print(q)
dev.off()

#=================================
# average the densities over time by species at station
# plot contours
#=====================================
t_temp<-tan_dat[,c(3,4,5,6,7,14)]
s_temp<-snow_dat[,c(3,4,5,6,7,14)]
#t_temp<-t_temp[!is.na(t_temp[,6]),]
#s_temp<-s_temp[!is.na(s_temp[,6]),]
s_temp$species<-"Snow"
t_temp$species<-"Tanner"
s_temp<-s_temp[-which(nchar(s_temp$station)>4),]
t_temp<-t_temp[-which(nchar(t_temp$station)>4),]

all_temp<-merge(s_temp,t_temp,by=c("year","station"))

plot(exp(all_temp$log_abund_ge45.x),exp(all_temp$log_abund_ge45.y))
cor.test(exp(all_temp$log_abund_ge45.x),exp(all_temp$log_abund_ge45.y))

mod<-gam((all_temp$log_abund_ge45.x)~s(all_temp$log_abund_ge45.y))
summary(mod)
plot(mod)

avg_dens_s<-s_temp%>%
  group_by(station,lat,lon)%>%
  summarise(avg_dens=mean((log_abund_ge45),na.rm=T))
 
avg_dens_t<-t_temp%>%
  group_by(station,lat,lon)%>%
  summarise(avg_dens=mean(log_abund_ge45,na.rm=T)) 

p<-ggplot() + 
  geom_tile(data=avg_dens_s, aes(x = lon, y = lat, fill = log(avg_dens)),width=.5,height=.25) +
  scale_fill_distiller(palette="Greens", na.value="grey",trans='reverse') 
  
q<-ggplot() + 
  geom_tile(data=avg_dens_t, aes(x = lon, y = lat, fill = (avg_dens)),width=.5,height=.25) +
  scale_fill_distiller(palette="Reds", na.value="grey",trans='reverse') 
  

p<-ggplot() + 
  geom_tile(data=avg_dens_s, aes(x = lon, y = lat, fill = (avg_dens)),width=.5,height=.25,alpha=.6) +
  scale_fill_distiller(palette="Greens", na.value="grey",trans='reverse') +
  new_scale_fill() +
  geom_tile(data=avg_dens_t, aes(x = lon, y = lat, fill = (avg_dens)),width=.5,height=.25,alpha=.6) +
  scale_fill_distiller(palette="Reds", na.value="grey",trans='reverse') +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")

print(p)



#=========================
# time series
#=======================
t_temp<-tan_dat[,c(1,3,4,5,6,7,9,14)]
s_temp<-snow_dat[,c(1,3,4,5,6,7,9,14)]
s_temp$species<-"Snow"
t_temp$species<-"Tanner"
s_temp<-s_temp[-which(nchar(s_temp$station)>4),]
t_temp<-t_temp[-which(nchar(t_temp$station)>4),]
s_temp[is.na(s_temp)]<-0
t_temp[is.na(t_temp)]<-0
#==check for zeroes to make densities right
#==or subsample for stations at which there are always samples by species

all_temp<-merge(s_temp,t_temp,by=c("year","station"))
all_temp2<-filter(all_temp,log_abund_101.x!=0 & all_temp$log_abund_45.y!=0)
plot((all_temp2$log_abund_101.x),(all_temp2$log_abund_45.y))
cor.test((all_temp2$log_abund_101.x),(all_temp2$log_abund_45.y))

mod<-gam(data=all_temp,(log_abund_45.y )~s(year)+s(log_abund_101.x)+s(lon.y,lat.y))
summary(mod)
plot(mod,pages=1,too.far=1,scheme=2)

avg_dens_s<-s_temp%>%
  group_by(year,species)%>%
  summarise(avg_dens_tot=mean(exp(log_abund_ge45),na.rm=T),
            avg_dens_sm=mean(exp(log_abund_45 ),na.rm=T),
            avg_dens_lg=mean(exp(log_abund_101 ),na.rm=T))

avg_dens_t<-t_temp%>%
  group_by(year,species)%>%
  summarise(avg_dens_tot=mean(exp(log_abund_ge45),na.rm=T),
            avg_dens_sm=mean(exp(log_abund_45 ),na.rm=T),
            avg_dens_lg=mean(exp(log_abund_101 ),na.rm=T))

#==add SDs later
plot_dat<-melt(rbind(avg_dens_t,avg_dens_s),id.vars=c("year","species"))
ggplot(plot_dat)+
  geom_line(aes(x=year,y=log(value),group=species,col=species),lwd=2)+
  facet_wrap(~variable,ncol=1)+theme_bw()

ggCcf(filter(plot_dat,species=="Tanner"&variable=='avg_dens_sm')$value,
      filter(plot_dat,species=="Snow"&variable=='avg_dens_lg')$value)

cor.test(filter(plot_dat,species=="Tanner"&variable=='avg_dens_sm')$value[7:47],
         filter(plot_dat,species=="Snow"&variable=='avg_dens_lg')$value[1:41])

#==cor matrices
avg_dens_s<-s_temp%>%
  group_by(year,species)%>%
  summarise(tot_snow=mean(exp(log_abund_ge45),na.rm=T),
            sm_snow=mean(exp(log_abund_45 ),na.rm=T),
            lg_snow=mean(exp(log_abund_101 ),na.rm=T))

avg_dens_t<-t_temp%>%
  group_by(year,species)%>%
  summarise(tot_tan=mean(exp(log_abund_ge45),na.rm=T),
            sm_tan=mean(exp(log_abund_45 ),na.rm=T),
            lg_tan=mean(exp(log_abund_101 ),na.rm=T))

library(GGally)
library(forecast)
ggpairs(log(cbind(avg_dens_t[,3:5],avg_dens_s[,3:5])))+theme_bw()

#==plot small tanner and lg snow by station
t_temp<-tan_dat[,c(1,3,4,5,6,7,9)]
s_temp<-snow_dat[,c(1,3,4,5,6,7,9)]
s_temp$species<-"Snow"
t_temp$species<-"Tanner"
s_temp<-s_temp[-which(nchar(s_temp$station)>4),]
t_temp<-t_temp[-which(nchar(t_temp$station)>4),]

t_temp$year<-t_temp$year-4
merged<-merge(s_temp,t_temp,by=c("station","year"))
plot(merged$log_abund_101.x,merged$log_abund_45.y)

ggplot(merged)+
  geom_point(aes(x=log_abund_101.x,y=log_abund_45.y))+
  facet_wrap(~year)+theme_bw()

#==check the average density over the previous five years of large snow crab and compare to sm tanner
#==checking total densities
# small to small
ggCcf(log(avg_dens_t[,4]),log(avg_dens_s[,4]))+theme_bw()
# lg Tan vs. sm snow
ggCcf(log(avg_dens_t[,5]),log(avg_dens_s[,4]))+theme_bw()
# sm Tan vs. lg snow
ggCcf(log(avg_dens_t[,4]),log(avg_dens_s[,5]))+theme_bw()
# lg Tan vs. lg snow
ggCcf(log(avg_dens_t[,5]),log(avg_dens_s[,5]))+theme_bw()

cor(avg_dens_s[,3:5],avg_dens_t[,3:5])


#============================================================================
# now let's look at relationships within the shared areas of each size class



#=======================================
# plot total numbers over time
#=======================================
in_dat<-rbind(ebs_dat,nbs_dat)
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<-min(in_dat$lon,na.rm=T)
lon_2<-max(in_dat$lon,na.rm=T)*.99
lat_1<-min(in_dat$lat,na.rm=T)
lat_2<-max(in_dat$lat,na.rm=T)
in_gam_dat<-in_dat
in_gam_dat_2<-in_gam_dat[,c(5,6,7,14)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]
nbs_bound<-data.frame(lat=c(62.1,62.1,61.1,61.1,60.5,60.5),
                      lon=c(-177,-172.5,-172.5,-171,-171,-165))
p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  facet_wrap(~year) +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
  strip.text.x = element_text(margin= margin(1,0,1,0)),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_total_crab.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#=======================================
# plot small male numbers over time
#=======================================
library(dplyr)
in_gam_dat_2<-in_gam_dat[,c(5,6,7,9)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]

p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  facet_wrap(~year) +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_small_crab.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#=======================================
# plot large male numbers over time
#=======================================
library(dplyr)
in_gam_dat_2<-in_gam_dat[,c(5,6,7,1)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]

p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  facet_wrap(~year) +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_large_crab.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#=======================================
# plot fishery recruits numbers over time
#=======================================
in_gam_dat_2<-in_gam_dat[,c(5,6,7,2)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]

p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  facet_wrap(~year) +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_fishery_recruits.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()


#=========================
# plot selected years
#=========================
#=======================================
# plot total numbers 
#=======================================
in_dat<-rbind(ebs_dat,nbs_dat)
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<-min(in_dat$lon,na.rm=T)
lon_2<-max(in_dat$lon,na.rm=T)*.99
lat_1<-min(in_dat$lat,na.rm=T)
lat_2<-max(in_dat$lat,na.rm=T)
in_gam_dat<-in_dat
in_gam_dat_2<-in_gam_dat[,c(5,6,7,14)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]
nbs_bound<-data.frame(lat=c(62.1,62.1,61.1,61.1,60.5,60.5),
                      lon=c(-177,-172.5,-172.5,-171,-171,-165))
in_fin<-filter(in_fin,year==2023)
p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_total_crab_recent.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#=======================================
# plot small male numbers 
#=======================================
library(dplyr)
in_gam_dat_2<-in_gam_dat[,c(5,6,7,9)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]
in_fin<-filter(in_fin,year==2023)
p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_small_crab_recent.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#=======================================
# plot large male numbers 
#=======================================
library(dplyr)
in_gam_dat_2<-in_gam_dat[,c(5,6,7,1)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]
in_fin<-filter(in_fin,year==2023)
p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_large_crab_recent.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#=======================================
# plot fishery recruits numbers 
#=======================================
in_gam_dat_2<-in_gam_dat[,c(5,6,7,2)]
in_fin<-melt(in_gam_dat_2,id=c('lat','lon','year'))
in_fin<-in_fin[complete.cases(in_fin),]
in_fin<-filter(in_fin,year==2023)
p<-ggplot() + 
  geom_tile(data=in_fin, aes(x = lon, y = lat, fill = value),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  geom_line(data=nbs_bound,aes(x=lon,y=lat),col='red')+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA))+labs(fill="Log density")
png("plots/ebs_nbs_fishery_recruits_recent.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

#########################################
#===calculate centroids of abundance
#===do these as heat maps or clusters instead
##########################################
Calc_centroids<-function(input)
{
  
  Centroids<-matrix(nrow=2,ncol=length(SurvYR))
  for(y in 1:length(SurvYR))
  { 
    trimStd<-which(nchar(as.character(StationYr[y,]))==4)
    inputDens<-input[y,trimStd]
    inputLat<-AllStnLoc[match(StationYr[y,trimStd],AllStation),1]
    inputLon<-AllStnLoc[match(StationYr[y,trimStd],AllStation),2]
    
    takers <-!is.na(inputDens)
    inputw <-inputDens[takers]
    inputXlon<-inputLon[takers]
    inputXlat<-inputLat[takers]
    
    Centroids[1,y]<- weighted.mean(inputXlon,w=inputw)
    Centroids[2,y]<- weighted.mean(inputXlat,w=inputw)
  }  
  
  return(Centroids)
}
useCol<-topo.colors(80,alpha=1)


#==large males
fmat<-Calc_centroids(input=DensityM101)
fmat<-fmat[,which(SurvYR>1988)]
in_dat<-as.data.frame(t(fmat))
colnames(in_dat)<-c("Lon","Lat")
in_dat$year<-c(seq(1989,2019),seq(2021,2023))
p<-ggplot() + 
  geom_text(data=in_dat, aes(x = Lon, y = Lat,label=year,color=year),cex=4) +
  scale_color_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),
        legend.position=c(.9,.7))+xlab('')+ylab('')

in_dat2<-melt(in_dat,id.var='year')
q<-ggplot(in_dat2)+
  geom_line(aes(x=year,y=value))+
  facet_wrap(~variable,scales='free_y')+
  theme_bw()+ylab('')

png("plots/lg_male_centroids.png",height=10,width=8,res=400,units='in')
(p/q)+plot_layout(heights=c(2,1))
dev.off()

#==mid size males
fmat<-Calc_centroids(input=DensityM_45_85)
fmat<-fmat[,which(SurvYR>1988)]
in_dat<-as.data.frame(t(fmat))
colnames(in_dat)<-c("Lon","Lat")
in_dat$year<-c(seq(1989,2019),seq(2021,2023))
p<-ggplot() + 
  geom_text(data=in_dat, aes(x = Lon, y = Lat,label=year,color=year),cex=4) +
  scale_color_distiller(palette="RdYlBu", na.value="grey") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),
        legend.position=c(.9,.7))+xlab('')+ylab('')

in_dat2<-melt(in_dat,id.var='year')
q<-ggplot(in_dat2)+
  geom_line(aes(x=year,y=value))+
  facet_wrap(~variable,scales='free_y')+
  theme_bw()+ylab('')

png("plots/mid_male_centroids.png",height=10,width=8,res=400,units='in')
(p/q)+plot_layout(heights=c(2,1))
dev.off()


#========================
#==plot clutch fullness, barren crab, and bitter crab over time
# do this above instead of right here
# color of the eggs (0 = no eggs, 2 =purple, 3 = brown, 4 = orange, 5 = purple-brown, 6 = pink), the condition of the eggs (0 = no 
# eggs, 1 = uneyed, 2 = eyed, 3 = dead, 4 = empty egg cases), and the size of the egg clutch (0 = 
# immature, 1 = mature female no eggs, 2 = trace to 1/8, 3 = 1/4, 4 = 1/2, 5 = 3/4, 6 = full).
#========================
library(ggplot2)
library(dplyr)
library(reshape)
library(ggridges)
library(gghighlight)
library(gridExtra)
clut_ful<-matrix(0,ncol=8,nrow=length(SurvYR))
disease_m_tot<-matrix(0,ncol=7,nrow=length(SurvYR))
disease_m_LE50<-matrix(0,ncol=7,nrow=length(SurvYR))
disease_m_GE50<-matrix(0,ncol=7,nrow=length(SurvYR))
disease_f_imm<-matrix(0,ncol=7,nrow=length(SurvYR))
disease_f_mat<-matrix(0,ncol=7,nrow=length(SurvYR))

BCS<-matrix(nrow=5,ncol=length(SurvYR))
calc_BCS<-function(input)
{
  denom<-    sum(input$SAMPLING_FACTOR,na.rm=T)
  num<-    sum(input$SAMPLING_FACTOR[input$DISEASE_CODE==2],na.rm=T)
  return(num/denom)
}
dis_type<-unique(survDAT$DISEASE_CODE)

for(y in 1:length(SurvYR))
{
  yrDAT<-survDAT[drvYear==SurvYR[y],]
  
  m_tot<-yrDAT[yrDAT$SEX==2,]
  temp<-table(m_tot$DISEASE_CODE)/nrow(m_tot)
  disease_m_tot[y,match(names(temp),dis_type)]<-temp
  BCS[1,y]<-calc_BCS(m_tot)
  
  m_smol<-yrDAT[yrDAT$SEX==2 & yrDAT$WIDTH <=50,]
  temp<-table(m_smol$DISEASE_CODE)/nrow(m_smol)
  disease_m_LE50[y,match(names(temp),dis_type)]<-temp
  BCS[2,y]<-calc_BCS(m_smol)
  
  m_big<-yrDAT[yrDAT$SEX==2 & yrDAT$WIDTH >50,]
  temp<-table(m_big$DISEASE_CODE)/nrow(m_big)
  disease_m_GE50[y,match(names(temp),dis_type)]<-temp
  BCS[3,y]<-calc_BCS(m_big)
  
  imm_f<-yrDAT[yrDAT$SEX==2 & yrDAT$EGG_CONDITION==0,]
  temp<-table(imm_f$DISEASE_CODE)/nrow(imm_f)
  disease_f_imm[y,match(names(temp),dis_type)]<-temp
  BCS[4,y]<-calc_BCS(imm_f)
  
  mat_f<-yrDAT[yrDAT$SEX==2 & yrDAT$EGG_CONDITION>0,]
  temp<-table(mat_f$DISEASE_CODE)/nrow(mat_f)
  disease_f_mat[y,match(names(temp),dis_type)]<-temp
  BCS[5,y]<-calc_BCS(mat_f)
  if(y<48)
     clut_ful[y,as.numeric(names(table(mat_f$CLUTCH_SIZE)))+1]<-table(mat_f$CLUTCH_SIZE)
  if(y==48)
  clut_ful[y,as.numeric(names(table(mat_f$CLUTCH_SIZE))[1:6])+1]<-table(mat_f$CLUTCH_SIZE)[1:6]
} 
png(paste("plots/bitter_crab.png",sep=''),height=8,width=8,res=350,units='in') 
par(mfrow=c(2,1),mar=c(.1,.1,.1,.1),oma=c(4,5,1,1))
plot(disease_m_LE50[,3]~SurvYR,type='b',las=1,ylab='',
     xaxt='n',ylim=c(0,0.025),xlim=c(1989,max(SurvYR)),pch=16)
lines(disease_m_tot[,3]~SurvYR,lty=2,col=2,type='b',pch=16)
lines(disease_m_GE50[,3]~SurvYR,lty=3,col=3,type='b',pch=16)

legend("topleft",bty='n',col=c(1,2,3),lty=c(1,2,3),
       legend=c("Less than 50mm carapace",
                "All crab",
                "Greater than 50mm carapace"))

plot(BCS[2,]~SurvYR,type='b',las=1,ylab='',xlab='Year',
     ylim=c(0,0.025),xlim=c(1989,max(SurvYR)),pch=16)
lines(BCS[1,]~SurvYR,lty=2,col=2,type='b',pch=16)
lines(BCS[3,]~SurvYR,lty=3,col=3,type='b',pch=16)

mtext(side=2,outer=T,"Prevalence of bitter crab syndrome",
      line=3.25)
dev.off()




clut_ful2<-sweep(clut_ful,1,apply(clut_ful,1,sum,na.rm=T),"/")
colnames(clut_ful2)<-seq(0,7) 
rownames(clut_ful2)<-c(seq(min(SurvYR),2019),seq(2021,SurvYR[length(SurvYR)]))
melted<-melt(clut_ful2)
colnames(melted)<-c("Year","Clutch","value")

avg_full<-rep(NA,nrow(clut_ful2))
for(x in 1:length(avg_full))
  avg_full[x]<-weighted.mean(x=seq(0,7),w=clut_ful2[x,])

png(paste("plots/clutch_fullness_full.png",sep=''),height=6,width=8,res=350,units='in') 
par(mfrow=c(2,1),mar=c(.1,4,.1,.1),oma=c(4,.1,1,1))
plot(avg_full~rownames(clut_ful2),
     xlim=c(1982,2023),type='b',las=1,ylab='Average clutch fullness score',
     xlab="Year",ylim=c(0,7),col='purple',pch=16,xaxt='n')
plot(apply(clut_ful2[,6:7],1,sum,na.rm=T)~rownames(clut_ful2),
     xlim=c(1982,2023),type='b',las=1,ylab='Proportion',
     xlab="Year",ylim=c(0,1),col='dark green',pch=16)
lines(clut_ful2[,2]~rownames(clut_ful2),type='b',col='blue',pch=16)
legend('center',bty='n',pch=16,col=c("dark green",'blue'),
       legend=c("Full clutches","Empty clutches"))
dev.off()


# color of the eggs (0 = no eggs, 2 =purple, 3 = brown, 4 = orange, 5 = purple-brown, 6 = pink), the condition of the eggs (0 = no 
# eggs, 1 = uneyed, 2 = eyed, 3 = dead, 4 = empty egg cases), and the size of the egg clutch (0 = 
# immature, 1 = mature female no eggs, 2 = trace to 1/8, 3 = 1/4, 4 = 1/2, 5 = 3/4, 6 = full).
##==ggridges
p <- ggplot(data=filter(melted,Year>1981) )
#p <- ggplot(data=melted)
p <- p + geom_density_ridges(aes(x=Clutch, y=Year, height = value, group = Year, 
                                 fill=stat(y),alpha=.9999),stat = "identity",scale=2,fill='purple') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90)) +
  labs(x="Clutch fullness score")
png(paste("plots/clutch_fullness_all.png",sep=''),height=8,width=6,res=350,units='in') 
print(p)
dev.off()
