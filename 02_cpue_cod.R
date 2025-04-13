library(dplyr)
library(ggplot2)
library(reshape2)
library(ggridges)
len_dat<-read.csv("cod/race_length_by_haul.csv",header=T,skip=7)
surv_dat<-read.csv("cod/race_cpue_by_haul.csv",header=T,skip=7)

tmp<-filter(len_dat,Survey!="NBS"&Sex!="Unknown")%>%
  group_by(Year,Length..mm.)%>%
  summarize(tot_num=sum(Frequency,na.rm=T))
  
size_yr <- ggplot(dat=tmp) 
size_yr <- size_yr + geom_density_ridges(aes(x=Length..mm., y=Year, height = tot_num,
                                             group = Year, 
                                             alpha=.9999),stat = "identity",scale=2,fill='royalblue3') +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  labs(x="Carapace width (mm)") +
  xlim(25,900)+
  scale_y_continuous(name='Year',position='right')

tmp_sum<-tmp%>%
  group_by(Year)%>%
  summarize(tot_n=sum(tot_num))

ggplot(tmp_sum)+
  geom_line(aes(x=Year,y=tot_n))+theme_bw()

names(surv_dat)
hist(surv_dat$Number.of.Fish)
sum(surv_dat$Number.of.Fish==0)
unique(surv_dat$Survey)

tmp2<-filter(surv_dat,Survey!="NBS")%>%
  group_by(Year)%>%
  summarize(avg_cpue=mean(CPUE..number.km2.))  

ggplot(tmp2)+
  geom_line(aes(x=Year,y=avg_cpue))+
  theme_bw()

tmp2<-filter(surv_dat,Survey!="NBS")%>%
  group_by(Year)%>%
  summarize(avg_cpue=mean(Weight.CPUE..kg.km2.))  

ggplot(tmp2)+
  geom_line(aes(x=Year,y=avg_cpue))+
  theme_bw()
