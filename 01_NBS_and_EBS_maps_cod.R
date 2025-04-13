#==loking at NBS by size
#++combine the NBS and EBS data from 2017-2019
#==have a multipanel plot with columns of different slices of the population

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
#==============================
# EBS data
#==============================
survDAT<-read.csv("cod/race_cpue_by_haul.csv",header=T,skip=7)
len_dat<-read.csv("cod/race_length_by_haul.csv",header=T,skip=7)
#=======================================
# plot total numbers over time
#=======================================
in_dat<-survDAT
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<-min(in_dat$end_lon,na.rm=T)
lon_2<-max(in_dat$end_lon,na.rm=T)*.99
lat_1<-min(in_dat$end_lat,na.rm=T)
lat_2<-max(in_dat$end_lat,na.rm=T)
nbs_bound<-data.frame(lat=c(62.1,62.1,61.1,61.1,60.5,60.5),
                      lon=c(-177,-172.5,-172.5,-171,-171,-165))
p<-ggplot() + 
  geom_tile(data=in_dat, aes(x = end_lon, y = end_lat, fill = log(CPUE..number.km2.)),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="white") +
  facet_wrap(~Year) +
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
png("plots/ebs_nbs_cod.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()

library(dplyr)
use_len<-filter(len_dat,Length..mm.>500)

p<-ggplot() + 
  geom_tile(data=use_len, aes(x = Ending.Longitude..dd., y = Ending.Latitude..dd., fill = log(Frequency)),width=.5,height=.25) +
  scale_fill_distiller(palette="RdYlBu", na.value="white") +
  facet_wrap(~Year) +
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
png("plots/ebs_nbs_cod_len.png",height=8,width=8,res=400,units='in')
print(p)
dev.off()
