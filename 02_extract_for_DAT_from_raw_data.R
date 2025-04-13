library(dplyr)
library(reshape2)
library(ggplot2)
library(png)
library(grid)
annotation_custom2 <-   function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf,ymax = Inf, data){ layer(data = data, 
                                                                                                      stat = StatIdentity, 
                                                                                                      position = PositionIdentity,
                                                                                                      geom = ggplot2:::GeomCustomAnn,
                                                                                                      inherit.aes = TRUE, 
                                                                                                      params = list(grob = grob,xmin = xmin, xmax = xmax,
                                                                                                                    ymin = ymin, ymax = ymax))}
in_png_opie<-readPNG('data/snowcrab.png')
EBSdata	<-read.csv("data/survey/EBSCrab_AB_Sizegroup.csv")
names(EBSdata)[1]<-"SURVEY_YEAR"

SurveyTotals<- filter(EBSdata,SIZE_CLASS_MM=='TOTAL')%>%
  group_by(SURVEY_YEAR,SEX) %>%
  summarise(Totals = sum(ABUNDANCE,na.rm=T))

#==make upper and lower CIs
EBSdata$upper_bio_ci<-exp(log(EBSdata$BIOMASS_MT)+1.96*(sqrt(log(1+EBSdata$BIOMASS_MT_CV^2))))
EBSdata$lower_bio_ci<-exp(log(EBSdata$BIOMASS_MT)-1.96*(sqrt(log(1+EBSdata$BIOMASS_MT_CV^2))))
EBSdata$upper_abnd_ci<-exp(log(EBSdata$ABUNDANCE)+1.96*(sqrt(log(1+EBSdata$ABUNDANCE_CV^2))))
EBSdata$lower_abnd_ci<-exp(log(EBSdata$ABUNDANCE)-1.96*(sqrt(log(1+EBSdata$ABUNDANCE_CV^2)))) 


plot_males<-c("MALE_GE25","MALE_GE78","MALE_GE95","MALE_GE102","MALE_TOTAL")

male_101<-filter(EBSdata,SIZE_GROUP=="MALE_GE102")
div_n<-1000
large_males<-ggplot(male_101)+
  geom_ribbon(aes(x=SURVEY_YEAR,ymin=lower_bio_ci/div_n,ymax=upper_bio_ci /div_n),fill='light grey')+
  geom_line(aes(x=SURVEY_YEAR,y=BIOMASS_MT/div_n))+
  geom_point(aes(x=SURVEY_YEAR,y=BIOMASS_MT/div_n))+
  theme_bw()+
  expand_limits(y=0)+
  xlab("")+
  scale_y_continuous(name="Biomass (1,000 t)",position='left')+
  annotation_custom2(rasterGrob(in_png_opie, interpolate=TRUE), 
                     xmin=2005, xmax=2020, ymin=250, ymax=400,data=male_101[1:1280,])
png("plots/obs_large_males.png",height=5,width=8,res=400,units='in')
print(large_males)
dev.off()


plot_males<-c("MALE_GE78","MALE_GE102","MALE_TOTAL")
male_com<-filter(EBSdata,SIZE_GROUP%in%plot_males)

div_n<-1000
exploit_males<-ggplot(male_com)+
  geom_ribbon(aes(x=SURVEY_YEAR,ymin=lower_bio_ci/div_n,ymax=upper_bio_ci /div_n),fill='light grey')+
  geom_line(aes(x=SURVEY_YEAR,y=BIOMASS_MT/div_n))+
  geom_point(aes(x=SURVEY_YEAR,y=BIOMASS_MT/div_n))+
  theme_bw()+
  expand_limits(y=0)+
  xlab("")+
  facet_wrap(~SIZE_GROUP)+
  scale_y_continuous(name="Biomass (1,000 t)",position='left')
png("plots/bio_time_series_males.png",height=8,width=8,res=400,units='in')
print(exploit_males)
dev.off() 

div_n<-1000
abund_males<-ggplot(male_com)+
  geom_ribbon(aes(x=SURVEY_YEAR,ymin=lower_abnd_ci/div_n,ymax=upper_abnd_ci /div_n),fill='light grey')+
  geom_line(aes(x=SURVEY_YEAR,y=ABUNDANCE/div_n))+
  geom_point(aes(x=SURVEY_YEAR,y=ABUNDANCE/div_n))+
  theme_bw()+
  expand_limits(y=0)+
  xlab("")+
  facet_wrap(~SIZE_GROUP,scales='free_y',ncol=1)+
  scale_y_continuous(name="Abundance",position='left')

abund_males2<-ggplot(male_com)+
  geom_ribbon(aes(x=SURVEY_YEAR,ymin=lower_abnd_ci/div_n,ymax=upper_abnd_ci /div_n),fill='light grey')+
  geom_line(aes(x=SURVEY_YEAR,y=ABUNDANCE/div_n))+
  geom_point(aes(x=SURVEY_YEAR,y=ABUNDANCE/div_n))+
  theme_bw()+
  expand_limits(y=0)+
  xlab("")+
  facet_wrap(~SIZE_GROUP,ncol=1)+
  scale_y_continuous(name="Abundance",position='right')

png("plots/abnd_time_series_males.png",height=8,width=8,res=400,units='in')
abund_males|abund_males2
dev.off() 


#=====================================
#==Female Survey Length compositions==
#=====================================
# this file is from EBC Crab, Large data Download, Abundance/Biomass, with MATURITY_SIZE_1MM_SHELLCON
EBSdata_in<-read.csv("data/survey/EBSCrab_Abundance_Biomass_female.csv",skip=7)
EBSdata<-EBSdata_in[EBSdata_in$SIZE_CLASS_MM>=25 & EBSdata_in$SEX=="FEMALE",]

female_total<- filter(EBSdata)%>%
  group_by(SURVEY_YEAR) %>%
  summarise(Totals = sum(ABUNDANCE))

female_s_totals<-filter(SurveyTotals,SEX=="FEMALE")
plot(female_s_totals$Totals~female_s_totals$SURVEY_YEAR,type='b')
lines(female_total$Totals~female_total$SURVEY_YEAR,type='b',col=2)

Years				<-sort(unique(EBSdata$SURVEY_YEAR))
LengthBins			<-seq(25,135,5)
FemaleMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
FemaleImmature		<-matrix(nrow=length(Years),ncol=length(LengthBins))

for(x in 1:length(Years))
{
 temp<-EBSdata[EBSdata$SURVEY_YEAR==Years[x],]
 for(y in 1:(length(LengthBins)-1))
 {
 FemaleMature[x,y]	<-sum(temp$ABUNDANCE[ temp$MATURITY=="MATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
 FemaleImmature[x,y]	<-sum(temp$ABUNDANCE[ temp$MATURITY=="IMMATURE" & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
 }
}

FemaleMature<-FemaleMature/1000000
FemaleImmature<-FemaleImmature/1000000

nal_fem_mat<-round(FemaleMature,4)[,1:(length(LengthBins)-1)]
sc_fem_mat<-sweep(nal_fem_mat,1,apply(nal_fem_mat,1,sum,na.rm=T),FUN="/")

apply(sc_fem_mat,1,sum)

nal_fem_imm<-round(FemaleImmature,4)[,1:(length(LengthBins)-1)]
sc_fem_imm<-sweep(nal_fem_imm,1,apply(nal_fem_imm,1,sum,na.rm=T),FUN="/")
apply(sc_fem_imm,1,sum)

write.table(round(sc_fem_mat,4), "data/derived/surv_len_comp_fem_mat.txt",row.names=FALSE,col.names=F)
write.table(round(sc_fem_imm,4), "data/derived/surv_len_comp_fem_imm.txt",row.names=FALSE,col.names=F)

#==========================
#==calc female biomasses===
#==========================
fem_bios<-EBSdata%>%
  group_by(SURVEY_YEAR,MATURITY) %>%
  summarise(biomass = round(sum(BIOMASS_MT,na.rm=T),2))

write.csv(dcast(fem_bios,SURVEY_YEAR~MATURITY)[-c(1:2),]/1000,"data/derived/index_female_biomass.csv",row.names=FALSE)

#===================================
#==Male survey Length compositions==
#===================================
# this file is from EBC Crab Abundance/Biomass, with SIZE_1MM_SHELLCON
EBSdata_in<-read.csv("data/survey/EBSCrab_Abundance_Biomass_male.csv",skip=7)
EBSdata<-EBSdata_in[EBSdata_in$SIZE_CLASS_MM>=25 & EBSdata_in$SEX=="MALE" &EBSdata_in$SURVEY_YEAR>1981,]

male_total<- filter(EBSdata)%>%
  group_by(SURVEY_YEAR) %>%
  summarise(Totals = sum(ABUNDANCE))

male_s_totals<-filter(SurveyTotals,SEX=="MALE")

plot(male_s_totals$Totals~male_s_totals$SURVEY_YEAR,type='b')
lines(male_total$Totals~male_total$SURVEY_YEAR,type='b',col=2)

Years				<-sort(unique(EBSdata$SURVEY_YEAR))
LengthBins	<-seq(25,135,5)
MaleNew			<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleOld			<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleNewMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleOldMature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleNewImmature		<-matrix(nrow=length(Years),ncol=length(LengthBins))
MaleOldImmature		<-matrix(0,nrow=length(Years),ncol=length(LengthBins))
NewShellIndex		<-3

for(x in 1:length(Years))
{
 temp<-EBSdata[EBSdata$SURVEY_YEAR==Years[x],]
 for(y in 1:(length(LengthBins)-1))
 {
 MaleNew[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION<NewShellIndex & temp$SEX=="MALE"  & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
 MaleOld[x,y]	<-sum(temp$ABUNDANCE[temp$SHELL_CONDITION>=NewShellIndex & temp$SEX=="MALE"  & temp$SIZE_CLASS_MM>=LengthBins[y] & temp$SIZE_CLASS_MM<LengthBins[y+1] ])
}
}

#==cap the sizes
TakeInd<-which(LengthBins==135)-1
MaleNew[,TakeInd]<-apply(MaleNew[,TakeInd:length(LengthBins)],1,sum,na.rm=T)
MaleOld[,TakeInd]<-apply(MaleOld[,TakeInd:length(LengthBins)],1,sum,na.rm=T)
MaleNew<-MaleNew[,1:TakeInd]/1000000	
MaleOld<-MaleOld[,1:TakeInd]/1000000	

#==make matrix of probability of maturity given newshell based on 
#==observed + average over years when observations are not available
male_mat<-read.csv("data/survey/EBS_male_opilio_mature_immature_ratio.csv")
male_mat$prop_mature<-male_mat$mat_abun /male_mat$tot
male_mat$prop_mature[male_mat$prop_mature=="NaN"]<-0
male_mat<-data.frame(male_mat)

library(ggplot2)
p<-ggplot(male_mat)+
  geom_line(aes(x=Size_1mm ,y=prop_mature,group=Year,col=Year))+
  xlim(35,135)
print(p)

library(mgcv)
all_yr<-seq(1982,max(male_mat$Year))
new_dat<-seq(27.5,132.5,5)
allmat<-matrix(nrow=length(all_yr),ncol=length(new_dat))
rownames(allmat)<-all_yr
colnames(allmat)<-new_dat
yrs<-unique(male_mat$Year)
for(x in 1:length(yrs))
{
temp<-filter(male_mat,Year==yrs[x])
temp$prop_mature[temp$Size_1mm>115]<-1
temp$prop_mature[temp$Size_1mm<40]<-0
temp<-data.frame(temp)
mod<-gam(data=temp,prop_mature~s(Size_1mm,k=20))
new_preds<-predict(mod,newdata=list(Size_1mm=new_dat))
allmat[match(yrs[x],all_yr),]<-new_preds
}

allmat[allmat<0]<-0
allmat[allmat>1]<-1
mean_mat<-apply(allmat,2,mean,na.rm=T)
for(x in which(is.na(allmat[,1])))
 allmat[x,]<-mean_mat

write.csv(allmat,"data/derived/prob_term_molt_males.csv")

male_mat_alt<-melt(allmat)
names(male_mat_alt)<-c("year","size","prop_mature")
png("plots/maturity_facet.png",height=8,width=8,res=400,units='in')
p<-ggplot(male_mat_alt)+
  geom_line(aes(x=size ,y=prop_mature,group=year,col=year),lwd=2)+
  scale_color_distiller(palette="RdYlBu", na.value="grey") +
  xlim(35,135)+theme_bw()+ylab("Proportion new shell mature")+
  xlab("Carapace width (mm)")+
  theme(legend.position=c(.9,.2))
print(p)
dev.off()

#==divide out to mature and immature
 in_mat<-allmat[-which(rownames(allmat)==2020),]
 MaleNewMature 	<-MaleNew*in_mat
 male_immature   <-MaleNew*(1-in_mat)
 male_mature<-MaleNewMature+MaleOld
 
 
 sc_male_mat<-sweep(male_mature,1,apply(male_mature,1,sum,na.rm=T),FUN="/")

 sc_male_imm<-sweep(male_immature,1,apply(male_immature,1,sum,na.rm=T),FUN="/")

write.table(round(sc_male_mat ,4), "data/derived/surv_len_comp_male_mat.txt",row.names=FALSE,col.names=F)
write.table(round(sc_male_imm,4), "data/derived/surv_len_comp_male_imm.txt",row.names=FALSE,col.names=F)

#==calculate biomass by different inputs
wt_at_size<-read.csv("data/survey/wt_at_size.csv")
male_mat_bio<-apply(sweep(male_mature,2,wt_at_size[,1],FUN="*"),1,sum)
male_imm_bio<-apply(sweep(male_immature,2,wt_at_size[,1],FUN="*"),1,sum)

plot(male_mat_bio,type='b',ylim=c(0,400))
lines(male_imm_bio,type='b',col=2)

write.table(male_mat_bio,"data/derived/index_mmb.txt",row.names=FALSE,col.names=F)
write.table(male_imm_bio,"data/derived/index_imm_male.txt",row.names=FALSE,col.names=F)

#==compare status quo data vs. recalculated
sq_dat<-read.csv("data/old_survey/sq_calc_mat_male_n_at_len.csv",header=T,check.names=FALSE)
sq_mmb<-apply(sweep(sq_dat[,-1],2,wt_at_size[,1],FUN="*"),1,sum)
plot(sq_mmb,type='b')
lines(male_mat_bio,type='b',col=2)




#============================================================
#  Retained catch length frequencies
# modified from "Item 1_20XX-20XX snow crab dockside
#============================================================
Data<-read.csv("data/Item 1_2022-23 snow crab dockside.csv")
binwidth<-2.5
LengthBins<-seq(27.5,162.5,5)
#LengthBins<-seq(25,165,5)
TotalNnew<-rep(0,length(LengthBins))
TotalNold<-rep(0,length(LengthBins))

upperBnd<-132.5
lowerBnd<-27.5

for(x in 1:length(LengthBins))
{
 tmp<-Data$New[Data$Size>(LengthBins[x]-binwidth) & Data$Size<=(LengthBins[x]+binwidth)]
 if(length(tmp)>0)
  TotalNnew[x]<-sum(tmp)
 tmp<-Data$Old[Data$Size>(LengthBins[x]-binwidth) & Data$Size<=(LengthBins[x]+binwidth)]
 if(length(tmp)>0)
  TotalNold[x]<-sum(tmp)
}

# TotalNnew[which(LengthBins==lowerBnd)]<-sum(TotalNnew[1:(which(LengthBins==lowerBnd))])
# TotalNold[which(LengthBins==lowerBnd)]<-sum(TotalNold[1:(which(LengthBins==lowerBnd))])

# TotalNnew[which(LengthBins==upperBnd)]<-sum(TotalNnew[(which(LengthBins==upperBnd)):length(LengthBins)])
# TotalNold[which(LengthBins==upperBnd)]<-sum(TotalNold[(which(LengthBins==upperBnd)):length(LengthBins)])
TotalNnew[is.na(TotalNnew)]<-0
TotalNold[is.na(TotalNold)]<-0

#==take the observed retained numbers, multiply to get the numbers at len by shell cond
#==observed retained numbers come from "Item 3 and 5...fish ticket snow crab.csv"
write.table(rbind(TotalNnew[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)],TotalNold[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)]),
 "data/derived/retained_catch_length_comps_new_then_old.txt",
 row.names=FALSE,col.names=F)
par(mfrow=c(2,1))
barplot(TotalNnew[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)])
barplot(TotalNold[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)])


ret_n_at_len_new<-TotalNnew
ret_n_at_len_old<-TotalNold
#ret_n_at_len_new<-TotalNnew[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)]
#ret_n_at_len_old<-TotalNold[which(LengthBins==lowerBnd):which(LengthBins==upperBnd)]
#================================================
# Total length frequencies
# modified from "item 2_20XX-20XX discard length frequency.xlx"
# take the male data for snow crab and save it in a .csv
#================================================
# #==males old and new shell total
#==from directed fishery
data<-read.csv("data/item_2_opilio_discard_len_m.csv")

library(dplyr)

#LengthBins			<-seq(25,165,5)
tot_len_comp_old		<-rep(0,length(LengthBins))
 for(y in 1:(length(LengthBins)-1))
 {
   tot_len_comp_old[y]	<-sum(data$Male.VeryOld[data$Size>=LengthBins[y] & data$Size<LengthBins[y+1]],
                         data$Male.Old[data$Size>=LengthBins[y] & data$Size<LengthBins[y+1]],na.rm=T)
 }

tot_len_comp_new		<-rep(0,length(LengthBins))
for(y in 1:(length(LengthBins)-1))
{
  tot_len_comp_new[y]	<-sum(data$Male.NewSoft[data$Size>=LengthBins[y] & data$Size<LengthBins[y+1]],
                            data$Male.Pliable[data$Size>=LengthBins[y] & data$Size<LengthBins[y+1]],
                            data$Male.New[data$Size>=LengthBins[y] & data$Size<LengthBins[y+1]],na.rm=T)
}

#=this number comes from "Item 4..."  number of observed Males caught in EBS snow fishery
total_crab_num<-60067664/1000
#=this number comes from item 3 and 5
ret_crab_num<-37492237/1000
disc_num<-total_crab_num-ret_crab_num

#==the raw numbers for both the observer and dockside length comps are not 'right', so we have to make sense of this
#==this requires getting them in a common currency
#==this will also be not necessary with GMACS

#==compositionify totals
tot_temp<-sum(tot_len_comp_new,tot_len_comp_old)
tot_comp_new<-tot_len_comp_new/tot_temp
tot_comp_old<-tot_len_comp_old/tot_temp

#==compositionify retained
ret_tot<-sum(ret_n_at_len_new,ret_n_at_len_old)
ret_comp_new<-ret_n_at_len_new/ret_tot
ret_comp_old<-ret_n_at_len_old/ret_tot

#==scale to the tot leng comps to total number of crab caught 
use_tot_new_n_len<-total_crab_num*tot_comp_new
use_tot_old_n_len<-total_crab_num*tot_comp_old

#==scale ret len comps to retained num caught
use_ret_new_n_len<-ret_crab_num*ret_comp_new
use_ret_old_n_len<-ret_crab_num*ret_comp_old

#==subtract retained from total
use_disc_new_len<-use_tot_new_n_len-use_ret_new_n_len
use_disc_old_len<-use_tot_old_n_len-use_ret_old_n_len

par(mfrow=c(2,1))
barplot(use_disc_new_len)
barplot(use_disc_old_len)

#==D'oh...there are negatives...'fix' them
use_disc_new_len[use_disc_new_len<0]<-0
use_disc_old_len[use_disc_old_len<0]<-0

plot(use_disc_new_len~LengthBins)
par(mfrow=c(3,1),mar=c(1,3,1,1))
barplot(use_tot_new_n_len,ylim=c(0,16000))
barplot(use_ret_new_n_len,ylim=c(0,16000))
barplot(use_disc_new_len,ylim=c(0,16000))

par(mfrow=c(4,1),mar=c(1,3,1,1))
barplot(use_tot_old_n_len,ylim=c(-1000,1600))
barplot(use_ret_old_n_len,ylim=c(-1000,1600))
barplot(use_disc_old_len,ylim=c(-1000,1600))
barplot(use_disc_old_len+use_ret_old_n_len,ylim=c(-1000,1600))

write.table(round(rbind(use_disc_new_len[1:22],use_disc_old_len[1:22]),1),"data/derived/male_disc_comps_new_then_old.txt",row.names=FALSE,col.names=F)

# #==males old and new shell total
#==from directed fishery
library(dplyr)
data<-read.csv("data/item_2_opilio_discard_len_f.csv")

total_disc_num_f<-10966/1000 # from item 4
LengthBins			<-seq(25,165,5)
tot_len_comp_fem		<-rep(0,length(LengthBins))
for(y in 1:(length(LengthBins)-1))
{
  tot_len_comp_fem[y]	<-sum(data$Fem.Total[data$Size>=LengthBins[y] & data$Size<LengthBins[y+1]],na.rm=T)
}

tot_f<-sum(tot_len_comp_fem)
tot_comp_f<-tot_len_comp_fem/tot_f
use_disc_n_len_f<-tot_comp_f*total_disc_num_f
plot(use_disc_n_len_f~LengthBins)
write.table(c(use_disc_n_len_f[1:22]),"data/derived/female_disc_comps.txt",row.names=FALSE,col.names=F)


#=======================================================================
# bycatch length composition
# these data come from the "Observer data" tab on AKfin
# Enter 'NORPAC Length Report - Haul & Length"
# Download the data for snow crab
# Time period should be July 1 (last year) to June 30 (this year)
#=======================================================================

LenDatBig<-read.csv("data/norpac_length_report/norpac_length_report.csv",skip=6)
LenDatBig$Haul.Offload.Date<-strptime(LenDatBig$Haul.Offload.Date,format="%d-%b-%y")
range(LenDatBig$Haul.Offload.Date)

#==CHECK DATES, EH!
LenDat<-LenDatBig[LenDatBig$Haul.Offload.Date >= "2022-07-01" & LenDatBig$Haul.Offload.Date <= "2023-06-30" & LenDatBig$Species.Name=="OPILIO TANNER CRAB",]

#==males old shell discard
LengthBins			<-seq(25,135,5)
BycatchFem			<-rep(0,length(LengthBins))
BycatchMale			<-rep(0,length(LengthBins))

 for(y in 1:(length(LengthBins)-1))
 {
  BycatchFem[y]	<-sum(LenDat$Frequency[LenDat$Length..cm.>=LengthBins[y] & LenDat$Length..cm.<LengthBins[y+1] & LenDat$Sex=="F"])
  BycatchMale[y]	<-sum(LenDat$Frequency[LenDat$Length..cm.>=LengthBins[y] & LenDat$Length..cm.<LengthBins[y+1] & LenDat$Sex=="M"])
 }

upperBnd							<-130
BycatchFem[which(LengthBins==upperBnd)]	<-sum(BycatchFem[(which(LengthBins==upperBnd)):length(LengthBins)])
BycatchFem						<-BycatchFem[1:which(LengthBins==upperBnd)]
BycatchMale[which(LengthBins==upperBnd)]	<-sum(BycatchMale[(which(LengthBins==upperBnd)):length(LengthBins)])
BycatchMale						<-BycatchMale[1:which(LengthBins==upperBnd)]

par(mfrow=c(1,2))
barplot(BycatchFem)
barplot(BycatchMale)

write.table(rbind(BycatchFem,BycatchMale),"data/derived/bycatch_len_comps_f_then_m.txt",
 row.names=FALSE,col.names=F)

#=======================================================================
# bycatch numbers
# these data come from the "Observer data" tab on AKfin
# Enter 'NORPAC catch Report"
# Download the data for snow crab in BS of BSAI
# Time period should be July 1 (last year) to June 30 (this year)
#=======================================================================


#==WHAT IS THE POINT OF THIS CHUNK OF CODE?
bycatch_dat<-read.csv("data/norpac_catch_report/norpac_catch_report.csv",skip=6)
temp<-strptime(bycatch_dat$Haul.Date,format="%d-%b-%y")
bycatch_dat$Haul.Date<-substr(temp,start=1,stop=10)

#==bycatch numbers total
  bycatchDat<-bycatch_dat[bycatch_dat$Haul.Date >= "2022-07-01" & bycatch_dat$Haul.Date <= "2023-06-30"& 
                            bycatch_dat$Species.Name=="OPILIO TANNER CRAB" & bycatch_dat$Gear.Description!="POT OR TRAP" ,]
  bycatch_num_tot<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
  bycatch_wt_tot<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)
  
  #==bycatch numbers total
  bycatchDat<-bycatch_dat[bycatch_dat$Haul.Date >= "2020-07-01" & bycatch_dat$Haul.Date <= "2021-06-30"& 
                            bycatch_dat$Species.Name=="OPILIO TANNER CRAB" & bycatch_dat$Gear.Description!="POT OR TRAP" ,]
  bycatch_num_tot<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
  
  bycatch_wt_tot<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)
#=========================================
# all years bycatch
#==================================

bycatch_dat_big<-bycatch_dat
#temp<-strptime(bycatch_dat_big$Haul.Date,format="%y-%b-%d")
#bycatch_dat_big$Haul.Date<-substr(temp,start=1,stop=10)
bycatch_dat_big$crab.year<-bycatch_dat_big$Year

#==bycatch numbers total
bycatch_year    <-sort(unique(bycatch_dat_big$Year))
bycatch_num_tot		<-rep(0,length(bycatch_year))
bycatch_wt_tot		<-rep(0,length(bycatch_year))

for(y in 1:(length(bycatch_year)-1))
{  
  bycatchDat<-bycatch_dat_big[bycatch_dat_big$Haul.Date >= paste(bycatch_year[y],"-07-01",sep="") & bycatch_dat_big$Haul.Date <= paste(bycatch_year[y]+1,"-06-30",sep="") & 
                                bycatch_dat_big$Species.Name=="OPILIO TANNER CRAB" ,]
  bycatch_num_tot[y]<-sum(bycatchDat$Extrapolated.Number,na.rm=T)
  bycatch_wt_tot[y]<-sum(bycatchDat$Extrapolated.Weight..kg.,na.rm=T)  
}

write.table(cbind(bycatch_year,bycatch_num_tot),"data/derived/bycatch_numbers_total.txt")
cbind(bycatch_year,bycatch_wt_tot/1000000)

plot(bycatch_wt_tot/1000000~bycatch_year,type='l')
write.table(cbind(bycatch_year,bycatch_wt_tot/1000000),"data/derived/bycatch_wt_total.txt")

#==bycatch by gear type
library(dplyr)
library(ggplot2)
for(y in 1:(length(bycatch_year)-1))
  bycatch_dat_big$crab.year[bycatch_dat_big$Haul.Date >= paste(bycatch_year[y],"-07-01",sep="") & bycatch_dat_big$Haul.Date <= paste(bycatch_year[y]+1,"-06-30",sep="")]<-bycatch_year[y]

in_dat<-bycatch_dat_big[,-24]
temp<- in_dat %>%
  group_by(crab.year,Gear.Description) %>%
  summarise(Bycatch = sum(Extrapolated.Number))

temp$Year<-as.numeric(temp$Year)

p<-qplot(x=crab.year,y=Bycatch,col=Gear.Description,data=temp)
png("plots/bycatch.png",height=8,width=8,res=400,units='in')
p + geom_line() + theme_bw()+theme(legend.position=c(.8,.8))
dev.off()


#==bycathc in other fisheries
crab_bycatch<-read.csv("data/crab_bycatch.csv")
plot_bycatch<-melt(crab_bycatch[,1:4],id.vars="Year")
colnames(plot_bycatch)<-c("Year","Fishery","Bycatch")
ggplot(plot_bycatch)+
  geom_line(aes(x=Year,y=Bycatch/1000000,color=Fishery),lwd=2)+
  theme_bw()+
  theme(legend.position=c(.5,.8))+
  ylab("Bycatch (kt)")
