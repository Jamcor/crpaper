#This code was written for the analysis of cCFU-F cell lifetime data as described in CR paper:
#'Quantifying intrinsic and extrinsic control of cell fates in cancer and stem/progenitor cell pedigrees with competing risks analysis'
#Code was written by R.E. Nordon and J.A. Cornwell 
#For questions please email: r.nordon@unsw.edu.au or james_cornwell@outlook.com

#The data is drawn from experiments: 'cCFUdata.txt'

#dependent libraries
library(R.matlab) #R.matlab version 3.1.1 used (newer or older versions may not be compatible)
library(timereg) #timereg version 1.7.0 used (newer or older versions may not be compatible)
library(mets) #mets version 0.1-13 used (newer or older versions may not be compatible)
library(cmprsk) #cmprisk version 2.2-7 used (newer or older versions may not be compatible)
library(survival) #survival version 2.37-7 used (newer or older versions may not be compatible)
library(ICC) #ICC version 2.2.1 used (newer or older versions may not be compatible)
library(coin) #coin version 1.0-24 use (newer or older versions may not be compatible)
library(GenABEL) #GenABEL version 1.8-0 used (newer or older versions may not be compatible)
library(psych) #psych version 1.5.1 used (newer or older versions may not be compatible)

#dependant files, 
#note: copy all .R files below into your working directory and substitute "L:/James/R/CR paper repository/" for the path of your working directory
source("L:/James/R/CR paper repository/ProbandwisePlot.R")
source("L:/James/R/CR paper repository/addSisterClusterID.R")
source("L:/James/R/CR paper repository/PlotPairedData.R")
source("L:/James/R/CR paper repository/ClusterMothersAndDaughters.R")
source("L:/James/R/CR paper repository/getMothersGFPStatus.R")
source("L:/James/R/CR paper repository/getMothers.R")
source("L:/James/R/CR paper repository/getCompletePairs.R")
source("L:/James/R/CR paper repository/SampleClones.R")
source('L:/James/R/CR paper repository/getpairs.R')
source('L:/James/R/CR paper repository/getMotherAndDaughterpairs.R')
source('L:/James/R/CR paper repository/ResetClustering.R')


#load data,  
#note: copy data file 'cCFU data.txt' into your working directory and substitute "L:/James/R/CR paper repository/" for the path of your working directory
P0<-read.delim("L:/James/R/CR paper repository/cCFU data.txt")


#create treatment groups and map fate outcomes 
P0$PDGF<-(P0$Group=='Condition 2')|(P0$Group=='Condition 5')|(P0$Group=='Condition 6')|(P0$Group=='Condition 8')
P0$FGF<-(P0$Group=='Condition 4')|(P0$Group=='Condition 6')|(P0$Group=='Condition 7')|(P0$Group=='Condition 8')
P0$TGF<-(P0$Group=='Condition 3')|(P0$Group=='Condition 5')|(P0$Group=='Condition 7')|(P0$Group=='Condition 8')
P0$CauseID<-"Censored"
P0$CauseID[P0$Cause==1]<-"Division" #CauseID = 1  is division
P0$CauseID[P0$Cause==2]<-"Apoptosis" #CauseID = 2 is death
P0$Cause[P0$CauseID=="Censored"]<-0; #CauseID = 0 is censored (lost of not complete)
# Group IDs
P0$GroupID[P0$Group=="Condition 1"]<-"_" #Serum free with no factors
P0$GroupID[P0$Group=="Condition 2"]<-"P" #PDGF only
P0$GroupID[P0$Group=="Condition 3"]<-"T" #TGF only
P0$GroupID[P0$Group=="Condition 4"]<-"F" #FGF only
P0$GroupID[P0$Group=="Condition 5"]<-"PT" #PDGF:TGF
P0$GroupID[P0$Group=="Condition 6"]<-"PF" #PDGF:FGF
P0$GroupID[P0$Group=="Condition 7"]<-"TF" #TGF:FGF
P0$GroupID[P0$Group=="Condition 8"]<-"PTF" #PDGF:TGF:FGF

#code for the analysis of cCFU-F data as presented in CR paper is shown below, 
#note: each section may be run independently  

# Generation 0 division --------------------------------------------------
ndx<-P0$Generation==0
tempdata<-P0[ndx,]
cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + TGF+PDGF+FGF+TGF:PDGF+TGF:FGF+PDGF:FGF+TGF:FGF:PDGF,
                                     data=tempdata,cause=tempdata$Cause,
                                     causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)

# predict
ksample<-cuminc(ftime=tempdata$Event_Time,fstatus=tempdata$CauseID,group=tempdata$GroupID,cencode="Censored")
cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF)+const(PDGF)+const(FGF)+const(TGF*PDGF)+const(TGF*FGF)
                              +const(PDGF*FGF)+const(TGF*FGF*PDGF),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
par(mfrow=c(1,1))
v<-qqnorm(cif.model.division$gamma,datax=TRUE,main="Generation 0",ylab="Relative risk of division",ylim=c(-0.5,3.0))
text(v$x,v$y,labels=c("T*","P**","F**","TP","TF","PF*","TFP***"),pos=4)
text(2,-0.6,pos=4,"*     p<0.05");text(2,-0.8,pos=4,"**   p<0.01");text(2,-1,pos=4,"*** p<0.001")


# probablity plot of covariates


# predict
division_G0<-predict(cif.model.division,X=1,Z=c(0,0,0,0,0,0,0),uniform=TRUE)
division_T_G0<-predict(cif.model.division,X=1,Z=c(1,0,0,0,0,0,0),uniform=TRUE)
division_P_G0<-predict(cif.model.division,X=1,Z=c(0,1,0,0,0,0,0),uniform=TRUE)
division_F_G0<-predict(cif.model.division,X=1,Z=c(0,0,1,0,0,0,0),uniform=TRUE)
division_PT_G0<-predict(cif.model.division,X=1,Z=c(1,1,0,1,0,0,0),uniform=TRUE)
division_FT_G0<-predict(cif.model.division,X=1,Z=c(1,0,1,0,1,0,0),uniform=TRUE)
division_PF_G0<-predict(cif.model.division,X=1,Z=c(0,1,1,0,0,1,0),uniform=TRUE)
division_PFT_G0<-predict(cif.model.division,X=1,Z=c(1,1,1,1,1,1,1),uniform=TRUE)
# plot division
attributes(ksample)
plot(ksample[["T Division"]]$time,ksample[["T Division"]]$est,
      type="s",,xlab="Time to Division (Days)",ylab="Cumulative Incidence",
      ylim=c(0,1.1),xlim=c(-0.5,4),lty=2,lwd=2,col="green",main="Generation 0")
lines(ksample[["_ Division"]]$time,ksample[["_ Division"]]$est,
      type="s",lty=2,lwd=2,col="grey")
lines(ksample[["P Division"]]$time,ksample[["P Division"]]$est,
      type="s",lty=2,lwd=2,col="red")
lines(ksample[["F Division"]]$time,ksample[["F Division"]]$est,
      type="s",lty=2,lwd=2,col="blue")
lines(ksample[["PTF Division"]]$time,ksample[["PTF Division"]]$est,
      type="s",lty=2,lwd=2,col="black")
legend(-0.65,1.2,legend=c("TGF+PDGF+bFGF","bFGF","PDGF","TGF","control"),lty=c(2,2,2,2,2),lwd=c(2,2,2,2,2,2),
       col=c("black","blue","red","green","grey"),bty="n")

plot(division_G0$time,division_G0$P1,type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1.1),xlim=c(-0.4,3.5),lty=1,lwd=2,col="grey",main="Generation 0")
lines(division_G0$time,division_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(division_G0$time,division_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(division_T_G0$time,division_T_G0$P1,type="s",lty=1,lwd=2,col="green")
lines(division_T_G0$time,division_T_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="green")
lines(division_T_G0$time,division_T_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="green")
lines(division_P_G0$time,division_P_G0$P1,type="s",lty=1,lwd=2,col="red")
lines(division_P_G0$time,division_P_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="red")
lines(division_P_G0$time,division_P_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="red")
lines(division_F_G0$time,division_F_G0$P1,type="s",lty=1,lwd=2,col="blue")
lines(division_F_G0$time,division_F_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="blue")
lines(division_F_G0$time,division_F_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="blue")
lines(division_PFT_G0$time,division_PFT_G0$P1,type="s",lty=1,lwd=2,col="black")
lines(division_PFT_G0$time,division_PFT_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="black")
lines(division_PFT_G0$time,division_PFT_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="black")
legend(-0.65,1.2,legend=c("TGF+PDGF+bFGF","SE","bFGF","PDGF","TGF","control"),lty=c(1,2,1,1,1,1),lwd=c(2,1,2,2,2,2),
       col=c("black","black","blue","red","green","grey"),bty="n")
# Generation >0 division --------------------------------------------------
ndx<-P0$Generation>0
tempdata<-P0[ndx,]
times=seq(from=0.05,by=0.05,to=3.5)
cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + (FGF+PDGF+TGF),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division)
ksample<-cuminc(ftime=tempdata$Event_Time,fstatus=tempdata$CauseID,group=tempdata$GroupID,cencode="Censored")
attributes(ksample)

cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF)+const(PDGF)+const(FGF),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division)
# predict
division_G0<-predict(cif.model.division,X=1,Z=c(0,0,0),uniform=TRUE)
division_T_G0<-predict(cif.model.division,X=1,Z=c(1,0,0),uniform=TRUE)
division_P_G0<-predict(cif.model.division,X=1,Z=c(0,1,0),uniform=TRUE)
division_F_G0<-predict(cif.model.division,X=1,Z=c(0,0,1),uniform=TRUE)
division_PT_G0<-predict(cif.model.division,X=1,Z=c(1,1,0),uniform=TRUE)
division_FT_G0<-predict(cif.model.division,X=1,Z=c(1,0,1),uniform=TRUE)
division_PF_G0<-predict(cif.model.division,X=1,Z=c(0,1,1),uniform=TRUE)
division_PFT_G0<-predict(cif.model.division,X=1,Z=c(1,1,1),uniform=TRUE)

# plot division
attributes(ksample)
plot(ksample[["T Division"]]$time,ksample[["T Division"]]$est,
     type="s",,xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=2,lwd=2,col="green",main="Generation>0")
lines(ksample[["P Division"]]$time,ksample[["P Division"]]$est,
      type="s",lty=2,lwd=2,col="red")
lines(ksample[["F Division"]]$time,ksample[["F Division"]]$est,
      type="s",lty=2,lwd=2,col="blue")
lines(ksample[["PTF Division"]]$time,ksample[["PTF Division"]]$est,
      type="s",lty=2,lwd=2,col="black")

plot(division_G0$time,division_G0$P1,type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,3.5),lty=1,lwd=2,col="grey",main="Generation > 0")
lines(division_G0$time,division_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(division_G0$time,division_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(division_P_G0$time,division_P_G0$P1,type="s",lty=1,lwd=2,col="red")
lines(division_P_G0$time,division_P_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="red")
lines(division_P_G0$time,division_P_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="red")
lines(division_F_G0$time,division_F_G0$P1,type="s",lty=1,lwd=2,col="blue")
lines(division_F_G0$time,division_F_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="blue")
lines(division_F_G0$time,division_F_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="blue")
lines(division_PFT_G0$time,division_PFT_G0$P1,type="s",lty=1,lwd=2,col="black")
lines(division_PFT_G0$time,division_PFT_G0$P1+division_G0$se.P1,type="s",lty=2,lwd=1,col="black")
lines(division_PFT_G0$time,division_PFT_G0$P1-division_G0$se.P1,type="s",lty=2,lwd=1,col="black")

legend(0,1,legend=c("TGF+PDGF+bFGF","SE","bFGF","PDGF","control"),lty=c(1,2,1,1,1),lwd=c(2,1,2,2,2),
       col=c("black","black","blue","red","grey"),bty="n")
# Generation 0 death ------------------------------------------------------

ndx<-P0$Generation==0
times<-seq(0.01,4,by=0.02)
tempdata<-P0[ndx,]
cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + TGF+PDGF+FGF+TGF:PDGF+TGF:FGF+PDGF:FGF+TGF:FGF:PDGF,
                              data=tempdata,cause=tempdata$Cause,times=times,
                              causeS=2,resample.iid=1,model="additive")
summary(cif.model.death)
cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF)+const(PDGF)+const(FGF)+const(TGF*PDGF)+const(TGF*FGF)
                              +const(PDGF*FGF)+const(TGF*FGF*PDGF),
                              data=tempdata,cause=tempdata$Cause,times=times,
                              causeS=2,resample.iid=1,model="additive")
cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF),
                           data=tempdata,cause=tempdata$Cause,times=times,
                           causeS=2,resample.iid=1,model="additive")
summary(cif.model.death)

ksample<-cuminc(ftime=tempdata$Event_Time,fstatus=tempdata$CauseID,group=tempdata$GroupID,cencode="Censored")
attributes(ksample)

# predict
death_G0<-predict(cif.model.death,X=1,Z=0,uniform=TRUE)
death_T_G0<-predict(cif.model.death,X=1,Z=1,uniform=TRUE)


# plot death
plot(ksample[["T Apoptosis"]]$time,ksample[["T Apoptosis"]]$est,
     type="s",,xlab="Time to Death (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=2,lwd=2,col="green")
lines(ksample[["P Apoptosis"]]$time,ksample[["P Apoptosis"]]$est,
      type="s",lty=2,lwd=2,col="red")
lines(ksample[["F Apoptosis"]]$time,ksample[["F Apoptosis"]]$est,
      type="s",lty=2,lwd=2,col="blue")
lines(ksample[["PTF Apoptosis"]]$time,ksample[["PTF Apoptosis"]]$est,
      type="s",lty=2,lwd=2,col="black")
# predict
death_G0<-predict(cif.model.death,X=1,Z=0,uniform=TRUE)
death_T_G0<-predict(cif.model.death,X=1,Z=1,uniform=TRUE)
plot(death_G0$time,death_G0$P1,type="s",xlab="Time to Death (Days)",ylab="Cumulative Incidence",
     ylim=c(0,0.2),xlim=c(0,3.5),lty=1,lwd=2,col="grey",main="Generation 0")
lines(death_G0$time,death_G0$P1+death_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(death_G0$time,death_G0$P1-death_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(death_T_G0$time,death_T_G0$P1,type="s",lty=1,lwd=2,col="green")
lines(death_T_G0$time,death_T_G0$P1+death_G0$se.P1,type="s",lty=2,lwd=1,col="green")
lines(death_T_G0$time,death_T_G0$P1-death_G0$se.P1,type="s",lty=2,lwd=1,col="green")

legend(0,0.2,legend=c("TGF","SE","Control"),lty=c(1,2,1),lwd=c(2,1,2),
       col=c("green","green","grey"),bty="n")
# Generatation >0 death ---------------------------------------------------
ndx<-P0$Generation>0
tempdata<-P0[ndx,]
times=seq(0.01,2.9,by=0.02)
cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + TGF+PDGF+FGF+TGF:PDGF+TGF:FGF+PDGF:FGF+TGF:FGF:PDGF,
                           data=tempdata,cause=tempdata$Cause,times=times,
                           causeS=2,resample.iid=1,model="additive")
cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + TGF+PDGF+FGF,
                           data=tempdata,cause=tempdata$Cause,times=times,
                           causeS=2,resample.iid=1,model="additive")
summary(cif.model.death)
cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF),times=times,
                           data=tempdata,cause=tempdata$Cause,
                           causeS=2,resample.iid=1,model="additive")
summary(cif.model.death)
# predict
death_G0<-predict(cif.model.death,X=1,Z=0,uniform=TRUE)
death_T_G0<-predict(cif.model.death,X=1,Z=1,uniform=TRUE)
ksample<-cuminc(ftime=tempdata$Event_Time,fstatus=tempdata$CauseID,group=tempdata$TGF,cencode="Censored")
attributes(ksample)
# plot death
plot(ksample[["TRUE Apoptosis"]]$time,ksample[["TRUE Apoptosis"]]$est,
     type="s",,xlab="Time to Death (Days)",ylab="Cumulative Incidence",
     ylim=c(0,0.3),xlim=c(0,4),lty=2,lwd=2,col="red")
lines(ksample[["FALSE Apoptosis"]]$time,ksample[["FALSE Apoptosis"]]$est,
      type="s",lty=2,lwd=2,col="grey")
legend(0,0.3,legend=c("TGF","Control"),lty=c(2,2),lwd=c(2,2),
       col=c("red","grey"),bty="n")
cif.model.division<-comp.risk(Surv(Event_Time,Cause>0) ~ 1 +(GFP+FGF)^2,
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="fg")
summary(cif.model.division)


plot(death_G0$time,death_G0$P1,type="s",xlab="Time to Death (Days)",ylab="Cumulative Incidence",
     ylim=c(0,0.3),xlim=c(0,3),lty=1,lwd=2,col="grey",main="Generation>0")
lines(death_G0$time,death_G0$P1+death_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(death_G0$time,death_G0$P1-death_G0$se.P1,type="s",lty=2,lwd=1,col="grey")
lines(death_T_G0$time,death_T_G0$P1,type="s",lty=1,lwd=2,col="green")
lines(death_T_G0$time,death_T_G0$P1+death_G0$se.P1,type="s",lty=2,lwd=1,col="green")
lines(death_T_G0$time,death_T_G0$P1-death_G0$se.P1,type="s",lty=2,lwd=1,col="green")

legend(0,0.3,legend=c("TGF","SE","Control"),lty=c(1,2,1),lwd=c(2,1,2),
       col=c("green","green","grey"),bty="n")
# Effect of GFP level on division -----------------------------------------
ndx<-(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")&(P0$Generation>0)
tempdata<-P0[ndx,]
cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(GFPAtDeath),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF)+const(PDGF)+const(FGF)+const(GFPAtDeath)+
                                const(GFPAtDeath*PDGF)+const(GFPAtDeath*TGF)+const(GFPAtDeath*FGF),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF)+const(PDGF)+const(FGF)+const(GFPAtBirth)+
                                const(GFPAtBirth*PDGF)+const(GFPAtBirth*TGF)+const(GFPAtBirth*FGF),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
# Effect of GFP level on death -----------------------------------------------
cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1 + const(TGF)+const(PDGF)+const(FGF)+const(GFPAtDeath)+
                                const(GFPAtDeath*PDGF)+const(GFPAtDeath*TGF)+const(GFPAtDeath*FGF),
                              data=tempdata,cause=tempdata$Cause,
                             causeS=2,resample.iid=1,model="additive")
summary(cif.model.death)
# Threshold of GFP for action of PDGF ------------------------------

P0$GFP<-P0$GFPAtDeath>100
P0$GFPTypeCause<-0
P0$GFPTypeCause[P0$GFP&(P0$Cause==1)]<-1
P0$GFPTypeCause[(!P0$GFP)&(P0$Cause==1)]<-2
P0$GFPTypeCause[P0$GFP&(P0$Cause==2)]<-3
P0$GFPTypeCause[(!P0$GFP)&(P0$Cause==2)]<-4
P0$GFPTypeCauseID[P0$GFPTypeCause==0]<-"Censored"
P0$GFPTypeCauseID[P0$GFPTypeCause==1]<-"Division (GFP+)"
P0$GFPTypeCauseID[P0$GFPTypeCause==2]<-"Division (GFP-)"
P0$GFPTypeCauseID[P0$GFPTypeCause==3]<-"Death (GFP+)"
P0$GFPTypeCauseID[P0$GFPTypeCause==4]<-"Death (GFP-)"
#Generation>0 division
ndx<-(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")&(P0$Generation>0)
tempdata<-P0[ndx,]
#get threshold for PDGFa effect
cif.model.division<-comp.risk(Surv(Event_Time,Cause>0) ~ 1 +(GFP+TGF)^2,
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="fg")
summary(cif.model.division)

# only GFP and PDGF interact 100%
times=seq(from=0.05,to=3.4,by=0.05)
cif.model.division.np<-comp.risk(Surv(Event_Time,Cause>0) ~ 1 +(GFP+PDGF)^2,
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division.np)
cif.model.division.p<-comp.risk(Surv(Event_Time,Cause>0) ~ 1 +(const(GFP)+const(PDGF))^2,
                                 data=tempdata,cause=tempdata$Cause,
                                 causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division.p)
# plot effect of PDGFRa and PDGF

# control 
marg.division<-predict(cif.model.division.p,X=1,Z=c(0,0,0))
plot(marg.division$time,marg.division$P1,type="s",,xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,3),lty=1,lwd=2,col="black",main="Generation>0")
lines(marg.division$time,marg.division$P1+marg.division$se.P1,type="s",lty=2,lwd=1,col="black")
lines(marg.division$time,marg.division$P1-marg.division$se.P1,type="s",lty=2,lwd=1,col="black")
# PDGFRa alone (+FGF)
marg.division<-predict(cif.model.division.p,X=1,Z=c(1,0,0))
lines(marg.division$time,marg.division$P1,type="s",lty=1,lwd=2,col="blue")
lines(marg.division$time,marg.division$P1+marg.division$se.P1,type="s",lty=2,lwd=1,col="blue")
lines(marg.division$time,marg.division$P1-marg.division$se.P1,type="s",lty=2,lwd=1,col="blue")
# PDGF (+FGF) alone
marg.division<-predict(cif.model.division.p,X=1,Z=c(0,1,0))
lines(marg.division$time,marg.division$P1,type="s",lty=1,lwd=2,col="red")
lines(marg.division$time,marg.division$P1+marg.division$se.P1,type="s",lty=2,lwd=1,col="red")
lines(marg.division$time,marg.division$P1-marg.division$se.P1,type="s",lty=2,lwd=1,col="red")
# PDGF and PDGFRa (+FGF)
marg.division<-predict(cif.model.division.p,X=1,Z=c(1,1,1))
lines(marg.division$time,marg.division$P1,type="s",lty=1,lwd=2,col="green")
lines(marg.division$time,marg.division$P1+marg.division$se.P1,type="s",lty=2,lwd=1,col="green")
lines(marg.division$time,marg.division$P1-marg.division$se.P1,type="s",lty=2,lwd=1,col="green")
legend(0,1,legend=c("PDGFRa-,PDGF-","PDGFRa+,PDGF-","PDGFRa-,PDGF+","PDGFRa+,PDGF+"),lty=c(1,1,1,1),lwd=c(2,2,2,2),
       col=c("black","blue","red","green"),bty="n")
# PDGFRa renewal ----------------------------------------------------------
#Generation 0 division
ndx<-(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")&(P0$Generation==0)
tempdata<-P0[ndx,]
times=seq(from=0.05,to=3.5,by=0.05)
cif.model.division<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(FGF+PDGF),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,times=times,resample.iid=1,model="additive")
summary(cif.model.division)
marg.division<-predict(cif.model.division,X=c(1,1,1))
plot(marg.division$time,marg.division$P1,type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,3.5),lty=2,lwd=2,col="black",main="GFP+ divisions")
cif.model.division<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(const(FGF)+const(PDGF)),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,times=times,resample.iid=1,model="additive")
summary(cif.model.division)
marg.division<-predict(cif.model.division,X=1,Z=c(1,1))
lines(marg.division$time,marg.division$P1,type="s",lty=1,lwd=2,col="black")
#Generation 1 division
ndx<-(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")&(P0$Generation==1)
tempdata<-P0[ndx,]
cif.model.division<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(FGF+PDGF),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,times=times,resample.iid=1,model="additive")
summary(cif.model.division)
marg.division<-predict(cif.model.division,X=c(1,1,1))
lines(marg.division$time,marg.division$P1,type="s",lty=2,lwd=2,col="blue")
cif.model.division<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(const(FGF)+const(PDGF)),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
marg.division<-predict(cif.model.division,X=1,Z=c(1,1))
lines(marg.division$time,marg.division$P1,type="s",lty=1,lwd=2,col="blue")
#Generation 2 division
ndx<-(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")&(P0$Generation==2)
tempdata<-P0[ndx,]
cif.model.division<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(FGF+PDGF),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,times=times,resample.iid=1,model="additive")
summary(cif.model.division)
marg.division<-predict(cif.model.division,X=c(1,1,1))
lines(marg.division$time,marg.division$P1,type="s",lty=2,lwd=2,col="red")
cif.model.division<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(const(FGF)+const(PDGF)),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
marg.division<-predict(cif.model.division,X=1,Z=c(1,1))
lines(marg.division$time,marg.division$P1,type="s",lty=1,lwd=2,col="red")
lines(c(0,3.5),c(0.5,0.5),type="l",lty="dotdash",lwd=1,col="black")

legend(0,1,legend=c("Gen0","Gen1","Gen2"),lty=c(1,1,1),lwd=c(2,2,2),
       col=c("black","blue","red"),bty="n")
# GFP+ vers GFP- division
ndx<-(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")&(P0$Generation>0)
tempdata<-P0[ndx,]
cif.model.division.pos<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(FGF+PDGF+TGF),
                                  data=tempdata,cause=tempdata$GFPTypeCause,times=times,
                                  causeS=1,resample.iid=1,model="fg")
summary(cif.model.division.pos)
cif.model.division.pos<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(const(FGF)+const(PDGF)),
                                  data=tempdata,cause=tempdata$GFPTypeCause,times=times,
                                  causeS=1,resample.iid=1,model="fg")
summary(cif.model.division.pos)
#GFP-
cif.model.division.neg<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(FGF+PDGF+TGF),
                                  data=tempdata,cause=tempdata$GFPTypeCause,times=times,
                                  causeS=2,resample.iid=1,model="fg")
summary(cif.model.division.neg)
cif.model.division.neg<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +const(FGF),
                                  data=tempdata,cause=tempdata$GFPTypeCause,times=times,
                                  causeS=2,resample.iid=1,model="fg")
summary(cif.model.division.neg)
marg.division.pos<-predict(cif.model.division.pos,X=1,Z=c(1,1))
marg.division.neg<-predict(cif.model.division.neg,X=1,Z=1)
plot(marg.division.pos$time,marg.division.pos$P1,type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,3.5),lty=1,lwd=2,col="green",main="PDGFRa+ versus PDGFRa- (Generation >0)")
lines(marg.division.pos$time,qnorm(0.025,marg.division.pos$P1,marg.division.pos$se.P1),type="s",lty=2,lwd=1,col="green")
lines(marg.division.pos$time,qnorm(0.975,marg.division.pos$P1,marg.division.pos$se.P1),type="s",lty=2,lwd=1,col="green")
lines(marg.division.neg$time,(marg.division.neg$P1),type="s",lty=1,lwd=2,col="black")
lines(marg.division.neg$time,qnorm(0.025,marg.division.neg$P1,marg.division.neg$se.P1),type="s",lty=2,lwd=1,col="black")
lines(marg.division.neg$time,qnorm(0.975,marg.division.neg$P1,marg.division.neg$se.P1),type="s",lty=2,lwd=1,col="black")
lines(c(0,3.5),c(0.5,0.5),type="l",lty="dotdash",lwd=1,col="black")
legend(0,1,legend=c("PDGFRa+","PDGFRa+ (95% CI)","PDGFRa-","PDGFRa- (95% CI)"),lty=c(1,2,1,2),lwd=c(2,1,2,1),
       col=c("green","green","black","black"),bty="n")
# PDGFRa survival (death) -------------------------------------------------
ndx<-(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")&(P0$Generation>0)
tempdata<-P0[ndx,]
times=seq(from=0.05,to=3.5,by=0.05)
# GFP+ death, generation >0
cif.model.death<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(FGF+PDGF+TGF),
                           data=tempdata,cause=tempdata$GFPTypeCause,
                           causeS=3,times=times,resample.iid=1,model="additive")
summary(cif.model.death)
cif.model.death<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +const(FGF)+const(PDGF)+const(TGF),
                           data=tempdata,cause=tempdata$GFPTypeCause,
                           causeS=3,times=times,resample.iid=1,model="additive")
summary(cif.model.death)
marg.death<-predict(cif.model.death,X=1,Z=c(1,1,1))
plot(marg.death$time,marg.death$P1,type="s",xlab="Time to Death (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,3.5),lty=1,lwd=2,col="green",main="GFP+ death")
lines(marg.death$time,qnorm(0.975,marg.death$P1,marg.death$se.P1),type="s",lty=2,lwd=1,col="green")
lines(marg.death$time,qnorm(0.025,marg.death$P1,marg.death$se.P1),type="s",lty=2,lwd=1,col="green")
# GFP- death, generation>0
cif.model.death<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +(FGF+PDGF+TGF),
                           data=tempdata,cause=tempdata$GFPTypeCause,
                           causeS=4,times=times,resample.iid=1,model="additive")
summary(cif.model.death)
cif.model.death<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 +const(FGF)+const(PDGF)+const(TGF),
                           data=tempdata,cause=tempdata$GFPTypeCause,
                           causeS=4,times=times,resample.iid=1,model="additive")
summary(cif.model.death)
marg.death<-predict(cif.model.death,X=1,Z=c(1,1,1))
lines(marg.death$time,marg.death$P1,type="s",lty=1,lwd=2,col="black")
lines(marg.death$time,qnorm(0.975,marg.death$P1,marg.death$se.P1),type="s",lty=2,lwd=1,col="black")
lines(marg.death$time,qnorm(0.025,marg.death$P1,marg.death$se.P1),type="s",lty=2,lwd=1,col="black")
# Concordance analysis of sisters -----------------------------------------
MothersAndDaughters<-ClusterMothersAndDaughters(P0,TRUE) #Get MotherAndDaughter clusters
P0<-ClusterMothersAndDaughters(P0,FALSE) #Get sibling clusters
ndx<-(P0$Generation>0)&(!is.na(P0$SisterClusterID))
tempdata<-P0[ndx,]
cif.model.division<-comp.risk(Surv(Event_Time,Cause>0) ~ 1+FGF+PDGF+cluster(SisterClusterID),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
division.cor<-cor.cif(cif.model.division,data=tempdata,cause1=1,cause2=1,sym=1)
summary(division.cor)

par(mfrow=c(1,1))
marg.division<-predict(cif.model.division,X=c(1,1,1))
ProbandwisePlot(division.cor,marg.cif1=marg.division,title="Division concordance",pos=c(1.3,0.2),
                eventtxt1="Sibling 1",eventtxt2="Sibling 2")

# # more concordance for shorter division times!
# cif.model.division<-comp.risk(Surv(Event_Time,Cause>0) ~ 1+cluster(SisterClusterID),
#                               data=tempdata,cause=tempdata$Cause,
#                               causeS=1,resample.iid=1,model="additive")
# summary(cif.model.division)
# t<-(tempdata$Event_Time-0.4)^3
# theta.des<-model.matrix(~1+t)
# division.cor<-cor.cif(cif.model.division,data=tempdata,theta.des=theta.des,cause1=1,cause2=1,sym=1)
# summary(division.cor)
# randt<-sample(tempdata$Event_Time,length(tempdata$Event_Time),replace=TRUE) # randomize data to check
# theta.des<-model.matrix(~1+randt)
# division.cor<-cor.cif(cif.model.division,data=tempdata,theta.des=theta.des,cause1=1,cause2=1,sym=1)
# summary(division.cor)

cif.model.death<-comp.risk(Surv(Event_Time,Cause>1) ~ 1+TGF+cluster(SisterClusterID),
                           data=tempdata,cause=tempdata$Cause,
                           causeS=2,resample.iid=1,model="additive")

summary(cif.model.death)
death.cor<-cor.cif(cif.model.death,data=tempdata,cause1=2,cause2=2,sym=1)
# create marginal distribution for death
marg.death<-predict(cif.model.death,X=c(1,1))
summary(death.cor)
ProbandwisePlot(death.cor,marg.cif1=marg.death,title="Death concordance",pos=c(0,1),
                eventtxt1="Sibling 1",eventtxt2="Sibling 2")

lifevsdeath.cor<-cor.cif(cif=cif.model.division,data=tempdata,sym=0,cause1=1,cause2=2,cif2=cif.model.death)
summary(lifevsdeath.cor)
# Concordance analysis of mothers and daughters ---------------------------
ndx<-MothersAndDaughters$Keep
ndx[is.na(ndx)]<-FALSE
tempdata<-MothersAndDaughters[ndx,]
cif.model.division<-comp.risk(Surv(Event_Time,Cause>1) ~ 1+PDGF+FGF+TGF+cluster(ClusterID),
                              data=tempdata,cause=tempdata$Cause,
                              causeS=1,resample.iid=1,model="fg")
summary(cif.model.division)
division.cor<-cor.cif(cif.model.division,data=tempdata,cause1=1,cause2=1,sym=1)
summary(division.cor)
marg.division<-predict(cif.model.division,X=c(1,1,1,1))
ProbandwisePlot(division.cor,marg.cif1=marg.division,title="Division concordance",pos=c(0,1),
                eventtxt1="Mother",eventtxt2="Daughter")

# Using separate marginals for mother and daughter 
tempdata$MotherDaughterCauseID<-0 
ndx<-(tempdata$IsMother)&(tempdata$Cause==1)
tempdata$MotherDaughterCauseID[ndx]<-1 # mother division
ndx<-(!tempdata$IsMother)&(tempdata$Cause==1)
tempdata$MotherDaughterCauseID[ndx]<-2 # daughter division
ndx<-(tempdata$Cause==2)
tempdata$MotherDaughterCauseID[ndx]<-3 # death
times=seq(from=0.1, to=3.5, by=0.05 )
cif.model.division.mother<-comp.risk(Surv(Event_Time,MotherDaughterCauseID>0) ~ 1+TGF+FGF+PDGF+cluster(ClusterID),
                                     data=tempdata,cause=tempdata$MotherDaughterCauseID,
                                     causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division.mother)
cif.model.division.daughter<-comp.risk(Surv(Event_Time,MotherDaughterCauseID>0) ~ 1+TGF+FGF+PDGF+cluster(ClusterID),
                                       data=tempdata,cause=tempdata$MotherDaughterCauseID,
                                       causeS=2,times=times,resample.iid=1,model="fg")
summary(cif.model.division.daughter)

division.cor<-cor.cif(cif.model.division.mother,
                      data=tempdata,cause1=1,cause2=1,cif2=cif.model.division.daughter,sym=0)
summary(division.cor)
marg.daughter<-predict(cif.model.division.daughter,X=c(1,1,1,1))
marg.mother<-predict(cif.model.division.mother,X=c(1,1,1,1))
ProbandwisePlot(division.cor,marg.cif1=marg.daughter,marg.cif2=marg.mother,title="Division concordance",pos=c(0,1),
                eventtxt1="Daughter",eventtxt2="Mother")



# bivariates plots (excluding censoring)
stats<-PlotPairedData(P0,MothersAndDaughters) 

# Check effect of generation on division times
ndx<-((P0$Generation>0)&(P0$Generation<4))
tempdata<-P0[ndx,]
times<-seq(from=0.1, to=3,by=0.05 )
cif.model.division<-comp.risk(Surv(Event_Time,Cause>0) ~ FGF+PDGF+TGF+Generation,
                              data=tempdata,cause=tempdata$Cause,times=times,
                              causeS=1,resample.iid=1,model="fg")
summary(cif.model.division)
cif.model.death<-comp.risk(Surv(Event_Time,Cause>0) ~ 1+FGF+PDGF+TGF+Generation,
                           data=tempdata,cause=tempdata$Cause,times=times,
                           causeS=2,resample.iid=1,model="additive")
summary(cif.model.death)
ndata<-data.frame(FGF=c(1,1,1),PDGF=c(1,1,1),TGF=c(1,1,1),Generation=c(1,2,3))

generations<-predict(cif.model.division,newdata=ndata,uniform=1,n.sim=1000)
attributes(generations)
plot(generations$time,generations$P1[1,],,xlab="Time to division (days)",ylab="Cumulative incidence",
     type="s",lty=1,lwd=2,col="black",main="Effect of generation number",ylim=c(0,1),xlim=c(0,3))
lines(generations$time,generations$P1[2,],type="s",lty=1,lwd=2,col="blue")
lines(generations$time,generations$P1[3,],type="s",lty=1,lwd=2,col="red")
# plot death
generations<-predict(cif.model.death,newdata=ndata,uniform=1,n.sim=1000)
lines(generations$time,generations$P1[1,],type="s",lty=2,lwd=2,col="black")
lines(generations$time,generations$P1[2,],type="s",lty=2,lwd=2,col="blue")
lines(generations$time,generations$P1[2,],type="s",lty=2,lwd=3,col="red")
legend(0.1,0.9,legend=c("G1 (division)","G2 (divisions)","G3 (division)","G1 (death)","G2 (death)","G3 (death)"),
       lty=c(1,1,1,2,2,2),lwd=c(2,2,2,2,2,2),
       col=c("black","blue","red","black","blue","red"),bty="n")
ksample=cuminc(ftime=tempdata$Event_Time,fstatus=tempdata$Cause,group=tempdata$Generation,cencode=0)
print(ksample)
attributes(ksample)

plot(ksample[["1 1"]]$time,ksample[["1 1"]]$est,xlab="Time to division (days)",ylab="Cumulative incidence",
     type="s",lty=1,lwd=2,col="black",main="Effect of generation number",ylim=c(0,1),xlim=c(0,3))
lines(ksample[["2 1"]]$time,ksample[["2 1"]]$est,type="s",lty=1,lwd=2,col="blue")
lines(ksample[["3 1"]]$time,ksample[["3 1"]]$est,type="s",lty=1,lwd=2,col="red")
lines(ksample[["1 2"]]$time,ksample[["1 2"]]$est,type="s",lty=2,lwd=2,col="black")
lines(ksample[["2 2"]]$time,ksample[["2 2"]]$est,type="s",lty=2,lwd=2,col="blue")
lines(ksample[["3 2"]]$time,ksample[["3 2"]]$est,type="s",lty=2,lwd=2,col="red")
legend(0.1,0.9,legend=c("G1 (division)","G2 (divisions)","G3 (division)","G1 (death)","G2 (death)","G3 (death)"),
       lty=c(1,1,1,2,2,2),lwd=c(2,2,2,2,2,2),
       col=c("black","blue","red","black","blue","red"),bty="n")
# Concordance analysis of GFP+ or GFP- sisters ------------------------------------
ndx<-(P0$Generation>0)&(!is.na(P0$SisterClusterID))
tempdata<-P0[ndx,]
cif.model.division<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1+FGF+PDGF+cluster(SisterClusterID),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,resample.iid=1,model="additive")
summary(cif.model.division)
division.cor<-cor.cif(cif.model.division,data=tempdata,cause1=1,cause2=1,sym=1)
summary(division.cor)

par(mfrow=c(1,1))
marg.division<-predict(cif.model.division,X=c(1,1,1))
ProbandwisePlot(division.cor,marg.cif1=marg.division,title="PDGFRa+ Division Concordance",pos=c(0,0.3),
                eventtxt1="Sibling 1",eventtxt2="Sibling 2")
# asymetric divisions
times=seq(from=0.05,to=3.5.by=0.05)
cif.model.division.pos<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1+FGF+PDGF+cluster(SisterClusterID),
                                  data=tempdata,cause=tempdata$GFPTypeCause,
                                  causeS=1,times=times,resample.iid=1,model="additive")
summary(cif.model.division.pos)
marg.division.pos<-predict(cif.model.division.pos,X=c(1,1,1))
cif.model.division.neg<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1+FGF+PDGF+cluster(SisterClusterID),
                                  data=tempdata,cause=tempdata$GFPTypeCause,
                                  causeS=2,times=times,resample.iid=1,model="additive")
summary(cif.model.division.neg)
mar.division.neg<-predict(cif.model.division.neg,X=c(1,1,1))
division.cor<-cor.cif(cif.model.division.pos,data=tempdata,cif2=cif.model.division.neg,cause1=1,cause2=2,sym=1)
summary(division.cor)
# no asymetric divisions
# symmetric negative divisions
cif.model.division.neg<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1+FGF+PDGF+cluster(SisterClusterID),
                                  data=tempdata,cause=tempdata$GFPTypeCause,
                                  causeS=2,times=times,resample.iid=1,model="additive")
summary(cif.model.division.neg)
division.cor<-cor.cif(cif.model.division.neg,data=tempdata,cause1=2,cause2=2,sym=1)
summary(division.cor)
marg.division.neg<-predict(cif.model.division.neg,X=c(1,1,1))
ProbandwisePlot(division.cor,marg.cif1=marg.division.neg,title="GFP- Division concordance",pos=c(0,1),
                eventtxt1="Sibling 1",eventtxt2="Sibling 2")
# Inheritance of GFP ------------------------------------------------------
P0<-getMothersGFPStatus(P0)
# GFP positive daughters from GFP positive mothers
times=seq(from=0.05, to= 3.5, by=0.05)
ndx<-(!is.na(P0$IsMotherGFPPos))&(P0$IsMotherGFPPos)&(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")
tempdata<-P0[ndx,]
cif.model.division.np<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (PDGF+FGF),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division.np)
cif.model.division.sp<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (const(PDGF)+const(FGF)),
                              data=tempdata,cause=tempdata$GFPTypeCause,
                              causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division.sp)
# plot
marg.division.np<-predict(cif.model.division.np,X=c(1,1,1))
marg.division.sp<-predict(cif.model.division.sp,X=1,Z=c(1,1))
par(mfcol=c(1,2))
plot(marg.division.sp$time,marg.division.sp$P1,type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=1,lwd=2,col="green",main="PDGFRa positive mother")
lines(marg.division.sp$time,qnorm(0.025,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="green")
lines(marg.division.sp$time,qnorm(0.975,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="green")
# GFP negative daughters from GFP positive mothers
times=seq(from=0.05, to= 3.5, by=0.05)
ndx<-(!is.na(P0$IsMotherGFPPos))&(P0$IsMotherGFPPos)&(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")
tempdata<-P0[ndx,]
cif.model.division.np<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (PDGF+FGF),
                                 data=tempdata,cause=tempdata$GFPTypeCause,
                                 causeS=2,times=times,resample.iid=1,model="fg")
summary(cif.model.division.np)
cif.model.division.sp<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (const(PDGF)+const(FGF)),
                                 data=tempdata,cause=tempdata$GFPTypeCause,
                                 causeS=2,times=times,resample.iid=1,model="fg")
summary(cif.model.division.sp)
# plot
marg.division.np<-predict(cif.model.division.np,X=c(1,1,1))
marg.division.sp<-predict(cif.model.division.sp,X=1,Z=c(1,1))
lines(marg.division.sp$time,marg.division.np$P1,type="s",lty=1,lwd=2,col="black")
lines(marg.division.sp$time,qnorm(0.025,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="black")
lines(marg.division.sp$time,qnorm(0.975,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="black")
legend(0,1.05,legend=c("PDGFRa+ daughter","95% CI","PDGFRa- daughter","95% CI"),lty=c(1,2,1,2),
       lwd=c(2,1,2,1),col=c("green","green","black","black"),bty="n")
# GFP positive daughters from GFP negative mothers
times=seq(from=0.05, to= 3.5, by=0.05)
ndx<-(!is.na(P0$IsMotherGFPPos))&(!P0$IsMotherGFPPos)&(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")
tempdata<-P0[ndx,]
cif.model.division.np<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (PDGF+FGF),
                                 data=tempdata,cause=tempdata$GFPTypeCause,
                                 causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division.np)
cif.model.division.sp<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (const(PDGF)+const(FGF)),
                                 data=tempdata,cause=tempdata$GFPTypeCause,
                                 causeS=1,times=times,resample.iid=1,model="fg")
summary(cif.model.division.sp)
# plot
marg.division.np<-predict(cif.model.division.np,X=c(1,1,1))
marg.division.sp<-predict(cif.model.division.sp,X=1,Z=c(1,1))
plot(marg.division.sp$time,marg.division.sp$P1,type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=1,lwd=2,col="green",main="PDGFRa negative mother")
lines(marg.division.sp$time,qnorm(0.025,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="green")
lines(marg.division.sp$time,qnorm(0.975,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="green")
# GFP negative daughters from GFP negative mothers
times=seq(from=0.05, to= 3.5, by=0.05)
ndx<-(!is.na(P0$IsMotherGFPPos))&(!P0$IsMotherGFPPos)&(P0$GFPAtBirth!="NaN")&(P0$GFPAtDeath!="NaN")
tempdata<-P0[ndx,]
cif.model.division.np<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (PDGF+FGF),
                                 data=tempdata,cause=tempdata$GFPTypeCause,
                                 causeS=2,times=times,resample.iid=1,model="fg")
summary(cif.model.division.np)
cif.model.division.sp<-comp.risk(Surv(Event_Time,GFPTypeCause>0) ~ 1 + (const(PDGF)+const(FGF)),
                                 data=tempdata,cause=tempdata$GFPTypeCause,
                                 causeS=2,times=times,resample.iid=1,model="fg")
summary(cif.model.division.sp)
# plot
marg.division.np<-predict(cif.model.division.np,X=c(1,1,1))
marg.division.sp<-predict(cif.model.division.sp,X=1,Z=c(1,1))
lines(marg.division.sp$time,marg.division.np$P1,type="s",lty=1,lwd=2,col="black")
lines(marg.division.sp$time,qnorm(0.025,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="black")
lines(marg.division.sp$time,qnorm(0.975,marg.division.sp$P1,marg.division.sp$se.P1),
      type="s",lty=2,lwd=1,col="black")
legend(0,1.05,legend=c("PDGFRa+ daughter","95% CI","PDGFRa- daughter","95% CI"),lty=c(1,2,1,2),
       lwd=c(2,1,2,1),col=c("green","green","black","black"),bty="n")
# Pie charts for fate outcomes, analysis of GFP+ and GFP- cells, histograms for GFP expression, pie charts GFP inheritance -----------------------------------------------
hist(P0$GFPAtDeath[P0$PDGF])
hist(P0$GFPAtDeath[P0$TGF])
hist(P0$GFPAtDeath[P0$FGF])

Threshold<-100
pdgf_gfp_pos<-sum(P0$GFPAtDeath[P0$PDGF]>=Threshold)
pdgf_gfp_neg<-sum(P0$GFPAtDeath[P0$PDGF]<Threshold)
fgf_gfp_pos<-sum(P0$GFPAtDeath[P0$FGF]>=Threshold)
fgf_gfp_neg<-sum(P0$GFPAtDeath[P0$FGF]<Threshold)
tgf_gfp_pos<-sum(P0$GFPAtDeath[P0$TGF]>=Threshold)
tgf_gfp_neg<-sum(P0$GFPAtDeath[P0$TGF]<Threshold)

pdgf.prop<-c(pdgf_gfp_pos, pdgf_gfp_neg)
fgf.prop<-c(fgf_gfp_pos, fgf_gfp_neg)
tgf.prop<-c(tgf_gfp_pos, tgf_gfp_neg)
pie(pdgf.prop, names(pdgf.prop)<-c("GFP positive" ,"GFP negative"),edges=5000)
pie(fgf.prop, names(fgf.prop)<-c("GFP positive" ,"GFP negative"),edges=5000)
pie(tgf.prop, names(tgf.prop)<-c("GFP positive" ,"GFP negative"),edges=5000)

total_gfp_pos<-sum(P0$GFPAtDeath>=Threshold)
total_gfp_neg<-sum(P0$GFPAtDeath<Threshold)

#pie charts for each generation
gen0_gfp_pos<-sum(P0$GFPAtDeath[P0$Generation==0]>=Threshold)
gen0_gfp_neg<-sum(P0$GFPAtDeath[P0$Generation==0]<Threshold)
gen1_gfp_neg<-sum(P0$GFPAtDeath[P0$Generation==1]<Threshold)
gen1_gfp_pos<-sum(P0$GFPAtDeath[P0$Generation==1]>=Threshold)
gen2_gfp_neg<-sum(P0$GFPAtDeath[P0$Generation==2]<Threshold)
gen2_gfp_pos<-sum(P0$GFPAtDeath[P0$Generation==2]>=Threshold)
gen3_gfp_neg<-sum(P0$GFPAtDeath[P0$Generation==3]<Threshold)
gen3_gfp_pos<-sum(P0$GFPAtDeath[P0$Generation==3]>=Threshold)
total.prop<-c(total_gfp_pos,total_gfp_neg)
gen0.prop<-c(gen0_gfp_pos,gen0_gfp_neg)
gen1.prop<-c(gen1_gfp_pos,gen1_gfp_neg)
gen2.prop<-c(gen2_gfp_pos,gen2_gfp_neg)
gen3.prop<-c(gen3_gfp_pos,gen3_gfp_neg)
percentotal<-(total_gfp_pos/(total_gfp_pos+total_gfp_neg))
percentG0<-(gen0_gfp_pos/(gen0_gfp_pos+gen0_gfp_neg))*100
percentG1<-(gen1_gfp_pos/(gen1_gfp_pos+gen1_gfp_neg))*100
percentG2<-(gen2_gfp_pos/(gen2_gfp_pos+gen2_gfp_neg))*100
percentG3<-(gen3_gfp_pos/(gen3_gfp_pos+gen3_gfp_neg))*100

pie(total.prop, names(total.prop)<-c("GFP positive" ,"GFP negative"),edges=5000,main="All generations")
pie(gen0.prop, names(gen0.prop)<-c("GFP positive" ,"GFP negative"),edges=5000,main="Gen0")
pie(gen1.prop, names(gen1.prop)<-c("GFP positive" ,"GFP negative"),edges=5000,main="Gen1")
pie(gen2.prop, names(gen2.prop)<-c("GFP positive" ,"GFP negative"),edges=5000,main="Gen2")
pie(gen3.prop, names(gen3.prop)<-c("GFP positive" ,"GFP negative"),edges=5000,main="Gen3")
P0$GFPAtDeath[P0$Generation==0]

##find number of GFP negative and positive cells that do and don't divide
gfp_pos<-P0[P0$GFPAtDeath>=Threshold,]
gfp_pos_total<-length(gfp_pos$Cause)
gfp_pos_div<-sum(gfp_pos$Cause==1)
gfp_pos_nc<-sum(gfp_pos$Cause==0)
# gfp_pos_lost<-sum(gfp_pos$Cause==3)
gfp_pos_death<-sum(gfp_pos$Cause==2)
percentdiv<-(gfp_pos_div/(gfp_pos_total))*100
percentnc<-(gfp_pos_nc/(gfp_pos_total))*100
percentlost<-(gfp_pos_lost/(gfp_pos_total))*100
percentdeath<-(gfp_pos_death/(gfp_pos_total))*100
gfp_pos.prop<-c(gfp_pos_div,gfp_pos_nc,gfp_pos_death)
pie(gfp_pos.prop, names(gfp_pos.prop)<-c("Division" ,"Right censored","Death"),edges=5000,main="GFP+")


gfp_neg<-P0[P0$GFPAtDeath<Threshold,]
gfp_neg_total<-length(gfp_neg$Cause)
gfp_neg_div<-sum(gfp_neg$Cause==1)
gfp_neg_nc<-sum(gfp_neg$Cause==0)
# gfp_neg_lost<-sum(gfp_neg$Cause==3)
gfp_neg_death<-sum(gfp_neg$Cause==2)
percentdiv<-(gfp_neg_div/(gfp_neg_total))*100
percentnc<-(gfp_neg_nc/(gfp_neg_total))*100
percentlost<-(gfp_neg_lost/(gfp_neg_total))*100
percentdeath<-(gfp_neg_death/(gfp_neg_total))*100
gfp_neg.prop<-c(gfp_neg_div,gfp_neg_nc,gfp_neg_death)
pie(gfp_neg.prop, names(gfp_neg.prop)<-c("Division" ,"Right censored","Death"),edges=5000,main="GFP-")

####number of transitioning cells 
gfp_pos_to_neg<-sum(P0$GFPAtBirth>=Threshold&P0$GFPAtDeath<Threshold)
gfp_neg_to_pos<-sum(P0$GFPAtBirth<Threshold&P0$GFPAtDeath>=Threshold)
gfp_neg_to_neg<-sum(P0$GFPAtBirth<Threshold&P0$GFPAtDeath<Threshold)
gfp_pos_to_pos<-sum(P0$GFPAtBirth>=Threshold&P0$GFPAtDeath>=Threshold)
total<-gfp_pos_to_neg+gfp_pos_to_pos+gfp_neg_to_neg+gfp_neg_to_pos
# 
# Comparative statistics (Pearson's correlation, intraclass correlation, bivariate dotplot, frequency of observed fate outcomes-----

MothersAndDaughters<-ClusterMothersAndDaughters(P0,TRUE) #Get MotherAndDaughter clusters
P0<-ClusterMothersAndDaughters(P0,FALSE) #Get sibling clusters
# ndx<-(P0$Generation>0)&(!is.na(P0$SisterClusterID)) 
ndx<-(P0$Generation>0)&(!is.na(P0$SisterClusterID)&!is.na(P0$GFPAtDeath)&(P0$GFPAtDeath>=100))
tempdata<-P0[ndx,]
sib1_time<-0
sib2_time<-0
ndx<-(tempdata$Progeny!=1)
ndx<-ndx&(tempdata$Cause==1)
tempdata<-tempdata[ndx,]

##bivariate plot of sibling cell cycle times (change the above ndx to gate for all cells or for GFP+ cells only)
i<-0
j<-0
for (i in 1:max(tempdata$SisterClusterID)) {
  sibndx<-(tempdata$SisterClusterID==i)
  if (sum(sibndx)==2){
    sibdata<-tempdata[sibndx,]
    sib1_time[j]<-sibdata[1,]$Event_Time
    sib2_time[j]<-sibdata[2,]$Event_Time
    j<-j+1
  }
}
plot(sib1_time,sib2_time,main="Concordance in Sibling Cycle Time",
     xlab=("Sibling 1 Cycle Time (Days)"),ylab=("Sibling 2 Cycle Time (Days)"),
     cex.main=0.75,cex.lab=0.75,col="blue",pch=19)

#Pearson's correlation coefficient (PCC)
r<-cor(sib1_time,sib2_time)
pearson<-cor.test(sib1_time,sib2_time,method="pearson")
print(pearson)

#Intraclass correlation coefficient (ICC)
ndx<-!is.na(P0$SisterClusterID)&!is.na(P0$Cause)&!is.na(P0$GFPAtDeath)&(P0$GFPAtDeath>=100) #change ndx to gate for GFP+ cells or not
sisters<-P0[ndx,]
icc_test<-ICCest(SisterClusterID,Event_Time,data=sisters,alpha=.05)
print(icc_test)

#get all unique sibling pairs 
# ndx<-!is.na(P0$SisterClusterID)&!is.na(P0$Cause)&!is.na(P0$GFPAtDeath)&(P0$GFPAtDeath>=100)
ndx<-!is.na(P0$SisterClusterID)&!is.na(P0$Cause)&!is.na(P0$GFPAtDeath)
sisters<-P0[ndx,]
uniquepairs<-unique(sisters$SisterClusterID)

#create groups for possible sibling pair fate outcomes 
#can run code for GFP+ cells and GFP- cells separately or all cells together including 'P0$GFPAtDeath>=100' in the above ndx

div_div<-0 #both siblings divide
death_death<-0 #both siblings dead
fate_censored<-0 #at least one sibling censored
div_death<-0 #one sibling divided, one sibling dead
gfp_nogfp<-0 #one sibling GFP+, one sibling GFP-
gfp_gfp<-0 #both siblings GFP+
nogfp_nogfp<-0 #both siblings GFP-

#loop through all sibling pairs to assign paired fate outcomes for use in binomial test and Yule's Q 
for(i in 1:length(uniquepairs)){
  pairdata<-sisters[sisters$SisterClusterID==uniquepairs[i],]
  
  if (length(pairdata$Cause)==2){
    
    if (pairdata$Cause[1]==1 & pairdata$Cause[2]==1) {
      div_div<-div_div+1
      
    }
    else if (pairdata$Cause[1]==2 & pairdata$Cause[2]==2){
      death_death<-death_death+1
    }
    else if (pairdata$Cause[1]==0 | pairdata$Cause[2]==0){
      fate_censored<-fate_censored+1
    }
    else if (pairdata$Cause[1]!=pairdata$Cause[2] & pairdata$Cause[1]!=0 & pairdata$Cause[2]!=0){
      div_death<-div_death+1
    }
    
    if (pairdata$GFPAtDeath[1]>=100 & pairdata$GFPAtDeath[2]<100){
      gfp_nogfp<-gfp_nogfp+1
      
    }
    else if (pairdata$GFPAtDeath[1]>=100& pairdata$GFPAtDeath[2]>=100){
      gfp_gfp<-gfp_gfp+1
      
    }
    else if (pairdata$GFPAtDeath[1]<100 & pairdata$GFPAtDeath[2]<100){
      nogfp_nogfp<-nogfp_nogfp+1
      
    }
    
  }
}


#caclulate the proportion of observed gfp sibling outcomes
total<-gfp_nogfp+gfp_gfp+nogfp_nogfp
prop.gfp_gfp<-(gfp_gfp/total)*100 #both daughters GFP+
prop.nogfp_nogfp<-(nogfp_nogfp/total)*100 #both daughters GFP-
prop.gfp_nogfp<-(gfp_nogfp/total)*100 #one daughter GFP+, one daughter GFP-


#binomial test (either including or excluding censored data)
binomial_div_cens<-binom.test(div_div,length(uniquepairs))
binomial_death_cens<-binom.test(death_death,length(uniquepairs))

binomial_div<-binom.test(div_div,(div_div+death_death))
binomial_death<-binom.test(death_death,(div_div+death_death))

#Yule's Q test for association in sibling cell fate 
x<-c(div_div,death_death,div_death,div_death)
x<-matrix(0,2,2)
x[1,1]<-div_div
x[1,2]<-death_death
x[2,1]<-div_death
x[2,2]<-div_death
yules<-Yule(x)
print(yules)
#Permutation test and Monte Carlo Simulation----

#initialise variables
pair_diff<-0
samp<-0
ndx<-tempdata$Cause==1  #change ndx to include or exclude GFP+ cells by adding or removing '&tempdata$GFPAtDeath>=100'
ndx<-tempdata$Cause==1&tempdata$GFPAtDeath>=100
sib_div<-tempdata[ndx,]
sib_div$SisterClusterID<-as.character(sib_div$SisterClusterID)
uniquepairs<-unique(sib_div$SisterClusterID) #get unique sibling pairs where both cells divide


#loop to remove any sibling pairs with only one sibling
for(i in 1:length(uniquepairs)){
  pairdata<-sib_div[sib_div$SisterClusterID==uniquepairs[i],]
  if(length(pairdata$Cause)==1){
    str=paste('\\b',uniquepairs[i],sep="")
    str=paste(str,'\\b',sep="")
    remove<-grep(str,sib_div$SisterClusterID)
    sib_div<-sib_div[-c(remove),]
  }
}
sib_div$SisterClusterID<-as.integer(sib_div$SisterClusterID)
uniquepairs<-unique(sib_div$SisterClusterID)
sibtimes<-matrix(0,2,(length(sib_div$X)/2))

#loop through all unique sibling pairs and calculate the difference in sibling pair cycle times
for (i in 1:length(uniquepairs)){
  pairdata<-sib_div[sib_div$SisterClusterID==uniquepairs[i],]
  sibtimes[1,i]<-(pairdata$Event_Time[1]*24)
  sibtimes[2,i]<-(pairdata$Event_Time[2]*24)
  pair_diff[i]<-((pairdata$Event_Time[1]*24)-(pairdata$Event_Time[2]*24))
#   print(pair_diff[i])
}


mean_sib_diff<-mean(abs(pair_diff)) #mean difference in sibling cycle time
sd_sib_diff<-sd(abs(pair_diff)) #standard deviation in sibling cycle time

#one way K-sample permutation test with 100,000 permutations
y<-factor(sib_div$SisterClusterID) #
pt<-oneway_test(sib_div$Event_Time~y,distribution=approximate(B=1000))
print(pt)

#plot sibling cycle times
plot((sibtimes[1,])/24,sibtimes[2,]/24,xlab=("Sibling 1 cycle time (days)"),
     ylab="(Sister 2 cycle time (days)",xlim=c(0,3.5),ylim=c(0,3.5),col="blue")

#use this to overlay GFP+ cells over the top of All cells, after changing the ndx and running above code again
points((sibtimes[1,])/24,sibtimes[2,]/24,col="green")

#pearsons correlation coefficient for all cells and for GFP+ cells
pearson<-cor.test(sibtimes[1,],sibtimes[2,],method="pearson")
print(pearson)

#plot density and distribution of the standardized test statistic 
layout(matrix(1:2, nrow = 2))
s <- support(pt)
d <- sapply(s, function(x) dperm(pt, x))
p <- sapply(s, function(x) pperm(pt, x))
plot(s, d, type = "S", xlab = "Time (days)", ylab = "Density")
plot(s, p, type = "S", xlab = "Time (days)", ylab = "Cumm. Probability")

#Monte Carlo simulation (bootstrapped) to calculate the mean difference in 10,000 random permutations of cell pairs
#i.e. the mean difference in sibling cell cycle times is calculated over 10,000 randomly generated data sets (by label shuffling)
#initialise variables
permutations<-matrix(0,3,length(sib_div$X)/2)
iterations<-10000 #set number of iterations
mean_diff<-matrix(0,1,iterations)
k<-seq(from=1,to=(length(sib_div$X)-1),by=2) 
for(i in 1:iterations){
  samp<-sample(sib_div$X,length(sib_div$X),replace=FALSE) #random sample from sib-div
  cell1_ndx<-samp[k] 
  cell2_ndx<-samp[k+1] 
#loop through the number of unique sibling pairs to create a data set the same as the original size  
  for(j in 1:length(sib_div$X)/2){
    cell1_div<-sib_div$Event_Time[sib_div$X==cell1_ndx[j]]*24 
    cell2_div<-sib_div$Event_Time[sib_div$X==cell2_ndx[j]]*24
    diff1<-abs(cell1_div-cell2_div) #take absolute value difference in cycle time
    diff2<-(cell1_div-cell2_div) #take difference in sibling cell cycle time
    permutations[1,j]<-cell1_div #random sib1 cycle time
    permutations[2,j]<-cell2_div #random sib2 cycle time
    permutations[3,j]<-diff1 #random difference in sibling cell cycle time
    # print(diff1)
    # print(diff2)
  }
  mean_diff[i]<-mean(permutations[3,])
  #   print(mean(permutations[3,]))
}
par(mfcol=c(1,1))
mean_pair_diff<-mean(abs(mean_diff)) #mean of mean pairwise differences over 10000 permutations 
sd_pair_diff<-sd(abs(mean_diff)) #standard deviation of the mean of mean pairwise differences
mean_sib_diff<-mean(abs(pair_diff))  #mean of sibling pair differences 
sd_sib_diff<-sd(abs(pair_diff)) #standard deviation of sibling pair differences 

#histograms for simulated data and empirical data
# hist(mean_diff,breaks=50,main="Mean pairwise difference distribution (random)",freq=TRUE,xlim=range(5,25),ylab="Cell number",xlab="Mean difference in cycle time from random permutations (hours)")
# hist(mean_diff,breaks=50,main="Mean pairwise difference distribution (random)",freq=TRUE,ylim=c(0,600),ylab="Cell number",xlab="Mean difference in cycle time from random permutations (hours)")
# lines(x=h$mids, y=h$density, type="l")
# hist(pair_diff,breaks=50,main="Sibling pairwise difference")
# hist(permutations[3,],breaks=50)

#z transform
z_distribution<-ztransform(mean_diff[,])
z_score<-(mean_sib_diff-mean_pair_diff)/(sd(mean_diff))
hist(z_distribution,main="Z-transformed distribution",xlim=range(-7,7))
p_value<-pnorm(z_score)

# hist(test)
# test2<-rntransform(pair_diff)
# hist(test2)
# icc_test<-ICCest(SisterClusterID,Event_Time,data=sib_div)
# mcfunction<-mcstoc(rempiricalD,)




