#This code was written for the analysis of breast cancer cell lifetime data as described in CR paper:
#'Quantifying intrinsic and extrinsic control of cell fates in cancer and stem/progenitor cell pedigrees with competing risks analysis'
#Code was written by J.A. Cornwell and R.E. Nordon
#For questions please email: r.nordon@unsw.edu.au or james_cornwell@outlook.com

#The data used in this code: 'BCC data.txt'

#dependant libraries
library(R.matlab) #R.matlab version 3.1.1 used (newer or older versions may not be compatible)
library(timereg) #timereg version 1.7.0 used (newer or older versions may not be compatible)
library(mets) #mets version 0.1-13 used (newer or older versions may not be compatible)
library(cmprsk) #cmprisk version 2.2-7 used (newer or older versions may not be compatible)

#load data,  
#note: copy data file 'BCC data.txt' into your working directory and substitute "L:/James/R/CR paper repository/" for the path of your working directory
P0<-read.delim("L:/James/R/CR paper repository/BCC data.txt")

#Map fate outcomes and create treatment groups
Data$CauseNum[Data$StopReason=='Division']<-1
Data$CauseNum[Data$StopReason=='Apoptosis']<-2
Data$CauseNum[Data$StopReason=='Lost cell']<-0
Data$CauseNum[Data$StopReason=='Not complete']<-0
Data$CauseNum[Data$StopReason=='Endomitosis']<-0
Data$CauseNum[Data$StopReason=='Nil']<-0
Data$CauseID[Data$CauseNum==1]<-"Division" #division is 1 
Data$CauseID[Data$CauseNum==2]<-"Death" #death is 2
Data$CauseID[Data$CauseNum==0]<-"RC"#right censored is 3


# Data=Data[Data$Generation==0,] # ndx to only look at generation 0
Data=Data[(Data$Treatment=="control")|(Data$Treatment=="dox")|(Data$Treatment=="nut"),] # only look at dox and nut versus control

#Compare wild type cell lines, no treatment----
par(mfcol=c(2,2))
ndx<-(Data$TP53=="wt")&(Data$Treatment=="control")
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$Cell_Line,cencode="RC")
print(ksample)
plot(ksample,main="Wild type p53, Control",xlab="days",ylim=c(0,1))

ndx<-(Data$TP53=="wt")&(Data$Treatment=="dox")
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$Cell_Line,cencode="RC")
print(ksample)
plot(ksample,main="Wild type p53, Doxorubicin",xlab="days",ylim=c(0,1))
## Compare mutant cell lines, no treatment----
ndx<-(Data$TP53=="mut")&(Data$Treatment=="control")
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$Cell_Line,cencode="RC")
print(ksample)
plot(ksample,main="Mutant p53, Control",xlab="days",ylim=c(0,1))

ndx<-(Data$TP53=="mut")&(Data$Treatment=="dox")
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$Cell_Line,cencode="RC")
print(ksample)
plot(ksample,main="Mutant p53, Doxorubicin",xlab="days",ylim=c(0,1))
#Compare mutant with wild type cell lines, no treatment----
ndx<-(Data$Cell_Line=="BT474")|(Data$Cell_Line=="mda-mb231")
ndx<-ndx&Data$Treatment=="control"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="BT474(wt) versus mda-mb231 (mut), Control",xlab="days",ylim=c(0,1))

ndx<-(Data$Cell_Line=="BT474")|(Data$Cell_Line=="T47D",xlab="days",ylim=c(0,1))
ndx<-ndx&Data$Treatment=="control"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="BT474(wt) versus T47D(mut), Control",xlab="days",ylim=c(0,1))

ndx<-(Data$Cell_Line=="MCF7")|(Data$Cell_Line=="mda-mb231")
ndx<-ndx&Data$Treatment=="control"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="MCF7(wt) versus mda-mb231 (mut), Control",xlab="days",ylim=c(0,1))


ndx<-(Data$Cell_Line=="MCF7")|(Data$Cell_Line=="T47D")
ndx<-ndx&Data$Treatment=="control"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="MCF7(wt) versus T47D (mut), Control",xlab="days",ylim=c(0,1))
#Compare mutant with wild type cell lines, Doxo----
ndx<-(Data$Cell_Line=="BT474")|(Data$Cell_Line=="mda-mb231")
ndx<-ndx&Data$Treatment=="dox"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="BT474(wt) versus mda-mb231 (mut), Doxorubicin",xlab="days",ylim=c(0,1))

ndx<-(Data$Cell_Line=="BT474")|(Data$Cell_Line=="T47D",xlab="days",ylim=c(0,1))
ndx<-ndx&Data$Treatment=="dox"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="BT474(wt) versus T47D(mut), Doxorubicin",xlab="days",ylim=c(0,1))

ndx<-(Data$Cell_Line=="MCF7")|(Data$Cell_Line=="mda-mb231")
ndx<-ndx&Data$Treatment=="dox"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="MCF7(wt) versus mda-mb231 (mut), Doxorubicin",xlab="days",ylim=c(0,1))


ndx<-(Data$Cell_Line=="MCF7")|(Data$Cell_Line=="T47D")
ndx<-ndx&Data$Treatment=="dox"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="MCF7(wt) versus T47D (mut), Doxorubicin",xlab="days",ylim=c(0,1))
#Compare mutant with wild type cell lines, nut----
ndx<-(Data$Cell_Line=="BT474")|(Data$Cell_Line=="mda-mb231")
ndx<-ndx&Data$Treatment=="nut"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="BT474(wt) versus mda-mb231 (mut), Nutlin3a",xlab="days",ylim=c(0,1))

ndx<-(Data$Cell_Line=="BT474")|(Data$Cell_Line=="T47D",xlab="days",ylim=c(0,1))
ndx<-ndx&Data$Treatment=="nut"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="BT474(wt) versus T47D(mut), Nutlin3a",xlab="days",ylim=c(0,1))

ndx<-(Data$Cell_Line=="MCF7")|(Data$Cell_Line=="mda-mb231")
ndx<-ndx&Data$Treatment=="nut"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="MCF7(wt) versus mda-mb231 (mut), Nutlin3a",xlab="days",ylim=c(0,1))


ndx<-(Data$Cell_Line=="MCF7")|(Data$Cell_Line=="T47D")
ndx<-ndx&Data$Treatment=="nut"
tempdata<-Data[ndx,]
ksample<-cuminc(ftime=tempdata$Age,fstatus=tempdata$CauseID,group=tempdata$TP53,cencode="RC")
print(ksample)
plot(ksample,main="MCF7(wt) versus T47D (mut), Nutlin3a",xlab="days",ylim=c(0,1))
#Regression----
ndx<-(Data$Cell_Line=="MCF7"|Data$Cell_Line=="mda-mb231")
tempdata<-Data[ndx,]
tempdata$Nut=0;
tempdata$Nut[tempdata$Treatment=="nut"]=1;
tempdata$Dox=0;
tempdata$Dox[tempdata$Treatment=="dox"]=1;
tempdata$wt=0;
tempdata$wt[tempdata$TP53=="wt"]=1;
Division<-comp.risk(Surv(Age,CauseNum>0) ~ 1+wt+Dox+Nut+wt*Nut+wt*Dox,data=tempdata,
                    cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division)
Division_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(wt)+const(Dox)+const(Nut)+const(wt*Nut)+const(wt*Dox),data=tempdata,
                       cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division_SP)

Death<-comp.risk(Surv(Age,CauseNum>0) ~ 1+wt+Dox+Nut+wt*Nut+wt*Dox,data=tempdata,
                 cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death)

Death_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(wt)+const(Dox)+const(Nut)+const(wt*Nut)+const(wt*Dox),data=tempdata,
                    cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death_SP)

# plot

division_control_NP<-predict(Division,X=c(1,0,0,0,0,0),uniform=TRUE)
division_dox_NP<-predict(Division,X=c(1,0,1,0,0,0),uniform=TRUE)
division_wt_nut_NP<-predict(Division,X=c(1,1,0,1,1,0),uniform=TRUE)
death_control_NP<-predict(Death,X=c(1,0,0,0,0,0),uniform=TRUE)
death_dox_NP<-predict(Death,X=c(1,0,1,0,0,0),uniform=TRUE)
death_wt_dox_NP<-predict(Death,X=c(1,1,1,0,0,1),uniform=TRUE)

division_control<-predict(Division_SP,X=1,Z=c(0,0,0,0,0),uniform=TRUE)
division_dox<-predict(Division_SP,X=1,Z=c(0,1,0,0,0),uniform=TRUE)
division_wt_nut<-predict(Division_SP,X=1, Z=c(1,0,1,1,0),uniform=TRUE)
death_control<-predict(Death_SP,X=1, Z=c(0,0,0,0,0),uniform=TRUE)
death_dox<-predict(Death_SP,X=1, Z=c(0,1,0,0,0),uniform=TRUE)
death_wt_dox<-predict(Death_SP,X=1,Z=c(1,1,0,0,1),uniform=TRUE)
par(mfcol=c(1,2))
plot(division_control_NP$time,division_control_NP$P1,
     type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=1,lwd=2,col="green")
lines(division_dox_NP$time,division_dox_NP$P1,
      type="s",lty=1,lwd=2,col="red")
lines(division_wt_nut_NP$time,division_wt_nut_NP$P1,
      type="s",lty=1,lwd=2,col="blue")
legend(1,0.5,legend=c("P53mut","P53mut+dox","P53wt+nut"),lty=c(1,1,1),lwd=c(2,2,2),col=c("green","red","blue"),bty="n")

plot(death_control_NP$time,death_control_NP$P1,
     type="s",xlab="Time to Death (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=1,lwd=2,col="green")
lines(death_dox_NP$time,death_dox_NP$P1,
      type="s",lty=1,lwd=2,col="red")
lines(death_wt_dox_NP$time,death_wt_dox_NP$P1,
      type="s",lty=1,lwd=2,col="blue")
legend(1,0.5,legend=c("P53mut","P53mut+dox","P53wt+dox"),lty=c(1,1,1),lwd=c(2,2,2),col=c("green","red","blue"),bty="n")

# refine regression model
summary(Division)
# make dox and intercept non-parametric
Division_Ver1<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(wt)+Dox+const(Nut)+const(wt*Nut)+const(wt*Dox),data=tempdata,
                         cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division_Ver1)
# goodness of fit?
division_control_ver1<-predict(Division_Ver1,X=c(1,0),Z=c(0,0,0,0),uniform=TRUE)
division_dox_ver1<-predict(Division_Ver1,X=c(1,1),Z=c(0,0,0,0),uniform=TRUE)
division_wt_nut_ver1<-predict(Division_Ver1,X=c(1,0),Z=c(1,1,1,0),uniform=TRUE)

par(mfcol=c(1,1))
plot(division_control_NP$time,division_control_NP$P1,
     type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=1,lwd=2,col="green")
lines(division_dox_NP$time,division_dox_NP$P1,
      type="s",lty=1,lwd=2,col="red")
lines(division_wt_nut_NP$time,division_wt_nut_NP$P1,
      type="s",lty=1,lwd=2,col="blue")
lines(division_control_ver1$time,division_control_ver1$P1,
      type="s",lty=2,lwd=2,col="green")
lines(division_dox_ver1$time,division_dox_ver1$P1,
      type="s",lty=2,lwd=2,col="red")
lines(division_wt_nut_ver1$time,division_wt_nut_ver1$P1,
      type="s",lty=2,lwd=2,col="blue")
#compare wild types using regression model----
ndx<-(Data$Cell_Line=="MCF7"|Data$Cell_Line=="BT474")
tempdata<-Data[ndx,]
tempdata$ID=0;
tempdata$ID[tempdata$Cell_Line=="BT474"]=1;
tempdata$Nut=0;
tempdata$Nut[tempdata$Treatment=="nut"]=1;
tempdata$Dox=0;
tempdata$Dox[tempdata$Treatment=="dox"]=1;
Division<-comp.risk(Surv(Age,CauseNum>0) ~ 1+ID+Dox+Nut+ID*Nut+ID*Dox,data=tempdata,
                    cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division)
Division_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(ID)+const(Dox)+const(Nut)+const(ID*Nut)+const(ID*Dox),data=tempdata,
                       cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division_SP)

Death<-comp.risk(Surv(Age,CauseNum>0) ~ 1+ID+Dox+Nut+ID*Nut+ID*Dox,data=tempdata,
                 cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death)

Death_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(ID)+const(Dox)+const(Nut)+const(ID*Nut)+const(ID*Dox),data=tempdata,
                    cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death_SP)
#compare P53 mut types using regression model----
ndx<-(Data$Cell_Line=="mda-mb231"|Data$Cell_Line=="T47D")
tempdata<-Data[ndx,]
tempdata$ID=0;
tempdata$ID[tempdata$Cell_Line=="T47D"]=1;
tempdata$Nut=0;
tempdata$Nut[tempdata$Treatment=="nut"]=1;
tempdata$Dox=0;
tempdata$Dox[tempdata$Treatment=="dox"]=1;
Division<-comp.risk(Surv(Age,CauseNum>0) ~ 1+ID+Dox+Nut+ID*Nut+ID*Dox,data=tempdata,
                    cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division)
Division_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+ID+Dox+const(Nut)+const(ID*Nut)+ID*Dox,data=tempdata,
                       cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division_SP)

tempdata<-tempdata[tempdata$Age<3.6,]
Death<-comp.risk(Surv(Age,CauseNum>0) ~ 1+ID+Dox+Nut+ID*Nut+ID*Dox,data=tempdata,
                 cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death)

Death_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(ID)+const(Dox)+const(Nut)+const(ID*Nut)+const(ID*Dox),data=tempdata,
                    cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death_SP)
#Do regression on pooled cell lines----
tempdata<-Data
tempdata$Nut=0;
tempdata$Nut[tempdata$Treatment=="nut"]=1;
tempdata$Dox=0;
tempdata$Dox[tempdata$Treatment=="dox"]=1;
tempdata$wt=0;
tempdata$wt[tempdata$TP53=="wt"]=1;
Division<-comp.risk(Surv(Age,CauseNum>0) ~ 1+wt+Dox+Nut+wt:Nut+wt:Dox,data=tempdata,
                    cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division)
Division_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(wt)+Dox+Nut+const(wt*Nut)+const(wt*Dox),data=tempdata,
                       cause=tempdata$CauseNum,causeS=1,resample.iid=1,model="additive")
summary(Division_SP)

Death<-comp.risk(Surv(Age,CauseNum>0) ~ 1+wt+Dox+Nut+wt:Nut+wt:Dox,data=tempdata,
                 cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death)

Death_SP<-comp.risk(Surv(Age,CauseNum>0) ~ 1+const(wt)+Dox+const(Nut)+const(wt*Nut)+wt:Dox,data=tempdata,
                    cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death_SP)

Death_<-comp.risk(Surv(Age,CauseNum>0) ~1+ Dox+wt:Dox,data=tempdata,
                  cause=tempdata$CauseNum,causeS=2,resample.iid=1,model="additive")
summary(Death_)
# Plot results
division_control<-predict(Division,X=c(1,0,0,0,0,0),uniform=TRUE)
division_wt<-predict(Division,X=c(1,1,0,0,0,0),uniform=TRUE)
division_wt_nut<-predict(Division,X=c(1,1,0,1,1,0),uniform=TRUE)
division_wt_dox<-predict(Division,X=c(1,1,1,0,0,1),uniform=TRUE)
division_nut<-predict(Division,X=c(1,0,0,1,0,0),uniform=TRUE)
division_dox<-predict(Division,X=c(1,0,1,0,0,0),uniform=TRUE)

plot(division_control$time,division_control$P1,
     type="s",xlab="Time to Division (Days)",ylab="Cumulative Incidence",
     ylim=c(0,1),xlim=c(0,4),lty=2,lwd=2,col="black")
lines(division_wt$time,division_wt$P1,
      type="s",lty=2,lwd=2,col="green")
lines(division_wt_nut$time,division_wt_nut$P1,
      type="s",lty=2,lwd=2,col="red")
lines(division_wt_dox$time,division_wt_dox$P1,
      type="s",lty=2,lwd=2,col="blue")
lines(division_nut$time,division_nut$P1,
      type="s",lty=2,lwd=2,col="pink")
lines(division_dox$time,division_dox$P1,
      type="s",lty=2,lwd=2,col="violetred")
legend(2.3,0.6,legend=c("wtP53","mutP53","mutP53+nut","wtP53+dox","wtP53+nut","mutP53+dox"),lty=c(2,2,2,2,2,2),
       lwd=c(2,2,2,2,2,2),col=c("green","black","pink","blue","red","violetred"),bty="n")

division_control_sp<-predict(Division_SP,X=c(1,0,0),Z=c(0,0,0),uniform=TRUE)
division_wt_sp<-predict(Division_SP,X=c(1,0,0),Z=c(1,0,0),uniform=TRUE)
division_wt_nut_sp<-predict(Division_SP,X=c(1,0,1),Z=c(1,1,0),uniform=TRUE)
division_wt_dox_sp<-predict(Division_SP,X=c(1,1,0),Z=c(1,0,1),uniform=TRUE)
division_nut_sp<-predict(Division_SP,X=c(1,0,1),Z=c(0,0,0),uniform=TRUE)
division_dox_sp<-predict(Division_SP,X=c(1,1,0),Z=c(0,0,0),uniform=TRUE)

lines(division_control_sp$time,division_control_sp$P1,
     type="s",ylim=c(0,1),xlim=c(0,4),lty=1,lwd=1,col="black")
lines(division_wt_sp$time,division_wt_sp$P1,
      type="s",lty=1,lwd=1,col="green")
lines(division_wt_nut_sp$time,division_wt_nut_sp$P1,
      type="s",lty=1,lwd=1,col="red")
lines(division_wt_dox_sp$time,division_wt_dox_sp$P1,
      type="s",lty=1,lwd=1,col="blue")
lines(division_nut_sp$time,division_nut_sp$P1,
      type="s",lty=1,lwd=1,col="pink")
lines(division_dox_sp$time,division_dox_sp$P1,
      type="s",lty=1,lwd=1,col="violetred")
## plot death
death_control<-predict(Death,X=c(1,0,0,0,0,0),uniform=TRUE)
death_wt<-predict(Death,X=c(1,1,0,0,0,0),uniform=TRUE)
death_wt_nut<-predict(Death,X=c(1,1,0,1,1,0),uniform=TRUE)
death_wt_dox<-predict(Death,X=c(1,1,1,0,0,1),uniform=TRUE)
death_nut<-predict(Death,X=c(1,0,0,1,0,0),uniform=TRUE)
death_dox<-predict(Death,X=c(1,0,1,0,0,0),uniform=TRUE)

plot(death_control$time,death_control$P1,
     type="s",xlab="Time to Death (Days)",ylab="Cumulative Incidence",
     ylim=c(0,0.5),xlim=c(0,4),lty=2,lwd=2,col="black")
lines(death_wt$time,death_wt$P1,
      type="s",lty=2,lwd=2,col="green")
lines(death_wt_nut$time,death_wt_nut$P1,
      type="s",lty=2,lwd=2,col="red")
lines(death_wt_dox$time,death_wt_dox$P1,
      type="s",lty=2,lwd=2,col="blue")
lines(death_nut$time,death_nut$P1,
      type="s",lty=2,lwd=2,col="pink")
lines(death_dox$time,death_dox$P1,
      type="s",lty=2,lwd=2,col="violetred")
legend(0,0.5,legend=c("mutP53+dox","mutP53+nut","wtP53+nut","wtP53+dox","wtP53","mutP53"),lty=c(2,2,2,2,2,2),
       lwd=c(2,2,2,2,2,2),col=c("violetred","pink","red","blue","green","black"),bty="n")

death_control_sp<-predict(Death_SP,X=c(1,0,0),Z=c(0,0,0),uniform=TRUE)
death_wt_sp<-predict(Death_SP,X=c(1,0,0),Z=c(1,0,0),uniform=TRUE)
death_wt_nut_sp<-predict(Death_SP,X=c(1,0,0),Z=c(1,1,1),uniform=TRUE)
death_wt_dox_sp<-predict(Death_SP,X=c(1,1,1),Z=c(0,0,0),uniform=TRUE)
death_nut_sp<-predict(Death_SP,X=c(1,0,0),Z=c(0,1,0),uniform=TRUE)
death_dox_sp<-predict(Death_SP,X=c(1,1,0),Z=c(0,0,0),uniform=TRUE)

lines(death_control_sp$time,death_control_sp$P1,
     type="s",lty=1,lwd=1,col="black")
lines(death_wt_sp$time,death_wt_sp$P1,
      type="s",lty=1,lwd=1,col="green")
lines(death_wt_nut_sp$time,death_wt_nut_sp$P1,
      type="s",lty=1,lwd=1,col="red")
lines(death_wt_dox_sp$time,death_wt_dox_sp$P1,
      type="s",lty=1,lwd=1,col="blue")
lines(death_nut_sp$time,death_nut_sp$P1,
      type="s",lty=1,lwd=1,col="pink")
lines(death_dox_sp$time,death_dox_sp$P1,
      type="s",lty=1,lwd=1,col="violetred")
##pie charts to show proportion of cell fates and mean cycle time comparison----
tempdata<-Data
#get lifetimes of wt p53 cells (division, death, censored lifetimes) for each treatment group
wt_cont_div<-tempdata$Age[(tempdata$CauseNum==1 & tempdata$TP53=="wt"&tempdata$Treatment=="control")] 
wt_nut_div<-tempdata$Age[(tempdata$CauseNum==1 & tempdata$TP53=="wt"& tempdata$Treatment=="nut")]
wt_dox_div<-tempdata$Age[(tempdata$CauseNum==1 & tempdata$TP53=="wt"& tempdata$Treatment=="dox")]
wt_cont_death<-tempdata$Age[(tempdata$CauseNum==2 & tempdata$TP53=="wt"&tempdata$Treatment=="control")]
wt_nut_death<-tempdata$Age[(tempdata$CauseNum==2 & tempdata$TP53=="wt"& tempdata$Treatment=="nut")]
wt_dox_death<-tempdata$Age[(tempdata$CauseNum==2 & tempdata$TP53=="wt"& tempdata$Treatment=="dox")]
wt_cont_cens<-tempdata$Age[(tempdata$CauseNum==0 & tempdata$TP53=="wt"&tempdata$Treatment=="control")]
wt_nut_cens<-tempdata$Age[(tempdata$CauseNum==0 & tempdata$TP53=="wt"& tempdata$Treatment=="nut")]
wt_dox_cens<-tempdata$Age[(tempdata$CauseNum==0 & tempdata$TP53=="wt"& tempdata$Treatment=="dox")]
wt_cont_total<-length(wt_cont_div)+length(wt_cont_death)+length(wt_cont_cens)
wt_nut_total<-length(wt_nut_div)+length(wt_nut_death)+length(wt_nut_cens)
wt_dox_total<-length(wt_dox_div)+length(wt_dox_death)+length(wt_dox_cens)

#proportion of fate outcomes 
wt_cont_div_prop<-length(wt_cont_div)/wt_cont_total
wt_cont_death_prop<-length(wt_cont_death)/wt_cont_total
wt_cont_cens_prop<-length(wt_cont_cens)/wt_cont_total
wt_nut_div_prop<-length(wt_nut_div)/wt_nut_total
wt_nut_death_prop<-length(wt_nut_death)/wt_nut_total
wt_nut_cens_prop<-length(wt_nut_cens)/wt_nut_total
wt_dox_div_prop<-length(wt_dox_div)/wt_dox_total
wt_dox_death_prop<-length(wt_dox_death)/wt_dox_total
wt_dox_cens_prop<-length(wt_dox_cens)/wt_dox_total

wt.cont_prop=c(wt_cont_div_prop,wt_cont_death_prop,wt_cont_cens_prop)
wt.nut_prop=c(wt_nut_div_prop,wt_nut_death_prop,wt_nut_cens_prop)
wt.dox_prop=c(wt_dox_div_prop,wt_dox_death_prop,wt_dox_cens_prop)

#get lifetimes of MUT p53 cells (division, death, censored lifetimes) for each treatment group
mut_cont_div<-tempdata$Age[(tempdata$CauseNum==1 & tempdata$TP53=="mut"&tempdata$Treatment=="control")]
mut_nut_div<-tempdata$Age[(tempdata$CauseNum==1 & tempdata$TP53=="mut"& tempdata$Treatment=="nut")]
mut_dox_div<-tempdata$Age[(tempdata$CauseNum==1 & tempdata$TP53=="mut"& tempdata$Treatment=="dox")]
mut_cont_death<-tempdata$Age[(tempdata$CauseNum==2 & tempdata$TP53=="mut"&tempdata$Treatment=="control")]
mut_nut_death<-tempdata$Age[(tempdata$CauseNum==2 & tempdata$TP53=="mut"& tempdata$Treatment=="nut")]
mut_dox_death<-tempdata$Age[(tempdata$CauseNum==2 & tempdata$TP53=="mut"& tempdata$Treatment=="dox")]
mut_cont_cens<-tempdata$Age[(tempdata$CauseNum==0 & tempdata$TP53=="mut"&tempdata$Treatment=="control")]
mut_nut_cens<-tempdata$Age[(tempdata$CauseNum==0 & tempdata$TP53=="mut"& tempdata$Treatment=="nut")]
mut_dox_cens<-tempdata$Age[(tempdata$CauseNum==0 & tempdata$TP53=="mut"& tempdata$Treatment=="dox")]
mut_cont_total<-length(mut_cont_div)+length(mut_cont_death)+length(mut_cont_cens)
mut_nut_total<-length(mut_nut_div)+length(mut_nut_death)+length(mut_nut_cens)
mut_dox_total<-length(mut_dox_div)+length(mut_dox_death)+length(mut_dox_cens)

#proportion of fate outcomes 
mut_cont_div_prop<-length(mut_cont_div)/mut_cont_total
mut_cont_death_prop<-length(mut_cont_death)/mut_cont_total
mut_cont_cens_prop<-length(mut_cont_cens)/mut_cont_total
mut_nut_div_prop<-length(mut_nut_div)/mut_nut_total
mut_nut_death_prop<-length(mut_nut_death)/mut_nut_total
mut_nut_cens_prop<-length(mut_nut_cens)/mut_nut_total
mut_dox_div_prop<-length(mut_dox_div)/mut_dox_total
mut_dox_death_prop<-length(mut_dox_death)/mut_dox_total
mut_dox_cens_prop<-length(mut_dox_cens)/mut_dox_total

mut.cont_prop=c(mut_cont_div_prop,mut_cont_death_prop,mut_cont_cens_prop)
mut.nut_prop=c(mut_nut_div_prop,mut_nut_death_prop,mut_nut_cens_prop)
mut.dox_prop=c(mut_dox_div_prop,mut_dox_death_prop,mut_dox_cens_prop)

#plot pie charts
par(mfcol=c(2,3))
pie(wt.cont_prop,names<-c("Division" ,"Death", "Right censored"),edges=5000,main="WT control")
pie(wt.nut_prop, names<-c("Division" ,"Death", "Right censored"),edges=5000,main="WT nut")
pie(wt.dox_prop, names<-c("Division" ,"Death", "Right censored"),edges=5000,main="WT dox")
pie(mut.cont_prop,names<-c("Division" ,"Death", "Right censored"),edges=5000,main="MUT control")
pie(mut.nut_prop, names<-c("Division" ,"Death", "Right censored"),edges=5000,main="MUT nut")
pie(mut.dox_prop, names<-c("Division" ,"Death", "Right censored"),edges=5000,main="MUT dox")

#mean cycle times
mean(wt_cont_div) # mean cycle time of control WTp53 
mean(wt_nut_div) # mean cycle time of WTp53 treated with nut
mean(wt_dox_div) # mean cycle time of WTp53 treated with dox
mean(mut_cont_div) # mean cycle time of control WTp53 
mean(mut_nut_div) # mean cycle time of WTp53 treated with nut
mean(mut_dox_div) # mean cycle time of WTp53 treated with dox

#Welch two sample t-test for difference in mean cycle time between treatment groups and genotype
cont_mut_vs_wt_ttest<-t.test(wt_cont_div,mut_cont_div)
print(cont_mut_vs_wt_ttest)
nut_mut_vs_wt_ttest<-t.test(wt_nut_div,mut_nut_div)
print(nut_mut_vs_wt_ttest)
dox_mut_vs_wt_ttest<-t.test(wt_dox_div,mut_dox_div)
print(dox_mut_vs_wt_ttest)
wt_cont_vs_nut_ttest<-t.test(wt_cont_div,wt_nut_div)
print(wt_cont_vs_nut_ttest)
wt_cont_vs_dox_ttest<-t.test(wt_cont_div,wt_dox_div)
print(wt_cont_vs_dox_ttest)