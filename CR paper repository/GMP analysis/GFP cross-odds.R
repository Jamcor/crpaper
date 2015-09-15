#Cross odds ratio for competing risks. GFP expression
# Warning. Use 

library(mets)
library(timereg)

#Perform clustering analysis on clones
#Load table called GFP regresssion
RegData=read.table("H:\\CellRxLabDoc\\Modeling\\Matlab\\Cell Branching Systems Analysis\\Analysis in R\\GFPregression.txt")

# load table called GFPPairs

data=read.table( "H:\\CellRxLabDoc\\Modeling\\Matlab\\Cell Branching Systems Analysis\\Analysis in R\\GFPPairs.txt")

# model cumulative incidence functions for GFP expression and apoptosis
# using the timereg package
# get non-parametric cumulative hazard for GFP onset

# get cells treated with G-CSF
# Cells which synchonised in the same generation have a very high COR

GCSF_ind=(RegData$Growth_Factor==1)&(RegData$Generation==1)
GCSFData=RegData[GCSF_ind,]
GFPCIF_GCSF=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(Clone),data=GCSFData,
	causeS=1,cause=GCSFData$Cause,n.sim=0,times=GCSFData$Event_Time,
	model="fg",max.clust=500)
ApoptosisCIF_GCSF=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(Clone),data=GCSFData,
	cause=GCSFData$Cause,causeS=2,n.sim=0,times=GCSFData$Event_Time,
	model="fg",max.clust=500)

# Plot estimated cumulative hazards and incidence 
X1=GFPCIF_GCSF$cum[,1];Y1=GFPCIF_GCSF$cum[,2];Y2=1-exp(-Y1)
X2=ApoptosisCIF_GCSF$cum[,1];Y3=ApoptosisCIF_GCSF$cum[,2];Y4=1-exp(-Y3)
par(mfrow=c(1,4),lty=1)
plot(X1,Y1,type="l",ylab="cumulative hazard",xlab="days",main="GCSF: GFP onset")
plot(X1,Y2,type="l",ylab="cumulative incidence",xlab="days",main="GCSF: GFP onset")
plot(X2,Y3,type="l",ylab="cumulative hazard",xlab="days",main="GCSF: Apoptosis")
plot(X2,Y4,type="l",ylab="cumulative incidence",xlab="days",main="GCSF: Apoptosis")

# This will take a long time!
COR_GFP_GFP_Clone_GCSF=cor.cif(GFPCIF_GCSF,data=GCSFData,cause1=1,cause2=1,)
summary(COR_GFP_GFP_Clone_GCSF)
par(mfrow=c(1,3),lty=1)
con=concordance(COR_GFP_GFP_Clone_GCSF,cif1=Y2);
plot(X1,con$intercept$concordance[,1],type="l")
plot(X1,con$intercept$concordance[,2],type="l")
plot(X1,con$intercept$concordance[,3],type="l")

# calculate for M-CSF
MCSF_ind=(RegData$Growth_Factor==1)&(RegData$Generation==3)
MCSF_ind=MCSF_ind
MCSFData=RegData[MCSF_ind,]
GFPCIF_MCSF=comp.risk(Event(Event_Time,Cause,cens.code=0)~+1+cluster(Clone),data=MCSFData,
	cause=1,n.sim=0,times=MCSFData$Event_Time,
	model="fg",max.clust=500)
ApoptosisCIF_MCSF=comp.risk(Event(Event_Time,Cause,cens.code=0)~+1+cluster(Clone),data=MCSFData,
	cause=2,n.sim=0,times=MCSFData$Event_Time,
	model="fg",max.clust=NULL)

# Plot estimated cumulative hazards and incidence 
X1=GFPCIF_MCSF$cum[,1];Y1=GFPCIF_MCSF$cum[,2];Y2=1-exp(-Y1)
X2=ApoptosisCIF_MCSF$cum[,1];Y3=ApoptosisCIF_MCSF$cum[,2];Y4=1-exp(-Y3)
par(mfrow=c(1,4),lty=1)
plot(X1,Y1,type="l",ylab="cumulative hazard",xlab="days",main="MCSF: GFP onset")
plot(X1,Y2,type="l",ylab="cumulative incidence",xlab="days",main="MCSF: GFP onset")
plot(X2,Y3,type="l",ylab="cumulative hazard",xlab="days",main="MCSF: Apoptosis")
plot(X2,Y4,type="l",ylab="cumulative incidence",xlab="days",main="MCSF: Apoptosis")
COR_Division_MCSF=cor.cif(GFPCIF_MCSF,data=MCSFData,cause1=1,cause2=1,sym=1)
summary(COR_Division_MCSF)
# This will take a long time!
theta.des=model.matrix(~-1+factor(MCSFData$Generation));
COR_GFP_GFP_Clone_MSCF=cor.cif(GFPCIF_MCSF,data=MCSFData,cause1=1,cause2=1,sym=1,theta.des=theta.des)
summary(COR_GFP_GFP_Clone_MSCF)








# plot coefficients (non-parametric log cumulative incidence) and cumulative incidence function
# here the link function is h(x)=1-exp(-exp(x)) and its inverse which is used to get the cause
# specific cumulative incidence is log(-log(1-g(t)))
# However the cum parameter is the cumulative hazard (called the time varying coefficient)
# P1(t) =1-exp(-A(t))
GFPCIF=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(Clone),data=Data,
	cause=Data$Cause,causeS=1,n.sim=0,times=Data$Event_Time,
	model="fg",max.clust=NULL)

X1=GFPCIF$cum[,1]
Y1=GFPCIF$cum[,2]
Y2=1-exp(-Y1)
par(mfrow=c(1,2),lty=1)
plot(X1,Y1,type="l",ylab="cumulative hazard",xlab="days",main="GFP onset")
plot(X1,Y2,type="l",ylab="cumulative incidence",xlab="days",main="GFP onset")

# do the same for apoptosis
ApoptosisCIF=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(Clone),data=Data,
	cause=Data$Cause,causeS=2,n.sim=0,times=Data$Event_Time,
	model="fg",max.clust=NULL)
X2=ApoptosisCIF$cum[,1]
Y3=ApoptosisCIF$cum[,2];
Y4=1-exp(-Y3)
par(mfrow=c(1,4),lty=1)
plot(X1,Y1,type="l",ylab="cumulative hazard",xlab="days",main="GFP onset")
plot(X1,Y2,type="l",ylab="cumulative incidence",xlab="days",main="GFP onset")
plot(X2,Y3,type="l",ylab="cumulative hazard",xlab="days",main="Apoptosis")
plot(X2,Y4,type="l",ylab="cumulative incidence",xlab="days",main="Apoptosis")
# This will take a long time!
COR_GFP_GFP_Clone=cor.cif(GFPCIF,data=Data,cause1=1,cause2=1)
summary(COR_GFP_GFP_Clone)



# SISTERS data$Relatedness=1

sister_ndx=data$Relatedness==1
SisterData=data[sister_ndx,]
GCSF_ndx=data$Growth_Factor==1
SisterGCSFData=data[sister_ndx&GCSF_ndx,]
SisterMCSFData=data[sister_ndx&(!GCSF_ndx),]

SistersGFP=comp.risk(Event(Event_Time,Cause,cens.code=0)~+1+cluster(PairID),data=SisterData,
	cause=1,n.sim=0,times=SisterData$Event_Time,
	model="fg",max.clust=20000)
SistersApoptosis=comp.risk(Event(Event_Time,Cause,cens.code=0)~+1+cluster(PairID),data=SisterData,
	cause=2,n.sim=0,times=SisterData$Event_Time,
	model="fg",max.clust=20000)
par(mfrow=c(1,2))
plot(SistersGFP,xlab="days",ylab="log cumulative hazard",main="GCSF/MCSF")
plot(SistersApoptosis,xlab="days",ylab="log cumulative hazard",main="GCSF/MCSF")

GFPgivenGFP=cor.cif(SistersGFP,data=SisterData,cause1=1,cause2=1,theta=0)
summary(GFPgivenGFP)
GFPgivenApoptosis=cor.cif(cif=SistersGFP,cif2=SistersApoptosis,data=SisterData,cause1=1,cause2=2)
summary(GFPgivenApoptosis)
ApoptosisgivenApoptosis=cor.cif(cif=SistersApoptosis,data=SisterData,cause1=2,cause2=2)
summary(ApoptosisgivenApoptosis)

# GCSF group

SistersGCSFGFP=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(PairID),data=SisterGCSFData,
	cause=SisterGCSFData$Cause,causeS=1,n.sim=0,times=SisterGCSFData$Event_Time,
	model="fg",max.clust=NULL)
SistersGCSFApoptosis=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(PairID),data=SisterGCSFData,
	cause=SisterGCSFData$Cause,causeS=2,n.sim=0,times=SisterGCSFData$Event_Time,
	model="fg",max.clust=NULL)
plot(SistersGCSFGFP,xlab="days",ylab="log cumulative hazard",sub="GCSF")
plot(SistersGCSFApoptosis,xlab="days",ylab="log cumulative hazard",sub="GCSF")

GFPgivenGFP_GCSF=cor.cif(SistersGCSFGFP,data=SisterGCSFData,cause1=1,cause2=1)
summary(GFPgivenGFP_GCSF)
GFPgivenApoptosis_GCSF=cor.cif(cif=SistersGCSFGFP,cif2=SistersGCSFApoptosis,data=SisterGCSFData,cause1=1,cause2=2)
summary(GFPgivenApoptosis_GCSF)
ApoptosisgivenApoptosis_GCSF=cor.cif(cif=SistersGCSFApoptosis,data=SisterGCSFData,cause1=2,cause2=2)
summary(ApoptosisgivenApoptosis_GCSF)

# MCSF group

SistersMCSFGFP=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(PairID),data=SisterMCSFData,
	cause=SisterMCSFData$Cause,causeS=1,n.sim=0,times=SisterMCSFData$Event_Time,
	model="fg",max.clust=NULL)
SistersMCSFApoptosis=comp.risk(Surv(Event_Time,Cause>0)~+1+cluster(PairID),data=SisterMCSFData,
	cause=SisterMCSFData$Cause,causeS=2,n.sim=0,times=SisterMCSFData$Event_Time,
	model="fg",max.clust=NULL)
plot(SistersMCSFGFP,xlab="days",ylab="log cumulative hazard",sub="MCSF")
plot(SistersMCSFApoptosis,xlab="days",ylab="log cumulative hazard",sub="MCSF")

GFPgivenGFP_MCSF=cor.cif(SistersMCSFGFP,data=SisterMCSFData,cause1=1,cause2=1)
summary(GFPgivenGFP_MCSF)
GFPgivenApoptosis_MCSF=cor.cif(cif=SistersMCSFGFP,cif2=SistersMCSFApoptosis,data=SisterMCSFData,cause1=1,cause2=2)
summary(GFPgivenApoptosis_MCSF)
ApoptosisgivenApoptosis_MCSF=cor.cif(cif=SistersMCSFApoptosis,data=SisterMCSFData,cause1=2,cause2=2)
summary(ApoptosisgivenApoptosis_MCSF)

# perform regression with respect to growth factor
theta.des=model.matrix(~-1+factor(SisterData$Growth_Factor))

GFPgivenGFP_GF=cor.cif(SistersGFP,data=SisterData,cause1=1,cause2=1,theta.des=theta.des)
summary(GFPgivenGFP_GF)
GFPgivenApoptosis_GF=cor.cif(cif=SistersGFP,cif2=SistersApoptosis,data=SisterData,
	cause1=1,cause2=2,theta.des=theta.des)
summary(GFPgivenApoptosis_GF)
ApoptosisgivenApoptosis_GF=cor.cif(cif=SistersApoptosis,data=SisterData,
	theta.des=theta.des,cause1=2,cause2=2)
summary(ApoptosisgivenApoptosis_GF)

# perform regression with respect to sisters and cousins and cousing once removed
theta.des=model.matrix(~-1+(log(1/SisterData$Likelihood_Ratio)))
GFPgivenGFP_LLR=cor.cif(SistersGFP,data=SisterData,cause1=1,cause2=1,theta.des=theta.des)
summary(GFPgivenGFP_LLR)


