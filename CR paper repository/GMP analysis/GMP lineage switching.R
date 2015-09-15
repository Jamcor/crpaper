# Performs regression modeling of the competing risk subdistribution function
# on GMP data, looking at the onset of GFP versus apoptosis. 
# Note that event times are relative to the start of the experiment 
# (time of exposure to growth factor)

# load table called "GFPregression"


fname=file.choose()
n=nchar(fname); 
data=read.table(fname)
library(cmprsk)


# Estimate cumulative incidence functions from competing risk data and test equality accross groups
group=factor(data$Growth_Factor,c(1,2),c("GCSF","MCSF"))
ftime=data$Event_Time
fstatus=data$Cause
stats=cuminc(ftime,fstatus,group)
y=timepoints(stats,seq(0,7,by=0.01))
plot(stats,xlab='days',title="Empirical cumulative incidence for onset of GFP fluorescence",
	col=1:4,lty=1,curvlab=c("GCSF GFP","MCSF GFP",
   "GCSF apoptosis","MCSF apoptosis"))

path=getwd()
fp=paste(path[1],"/GMP lineage switching stats.txt",sep="")
 
write("Effect of growth factor, K groups test",fp)
write.table(stats$Tests,fp,sep="\t",quote=TRUE,row.names=FALSE,col.names=TRUE,append=TRUE)
write("Row Names:",fp,append=TRUE)
write("Time",fp,append=TRUE)
write(row.names(y$est),fp,append=TRUE)
write("1 is GFP expression, 2 is apoptosis",fp,append=TRUE)
write("Cumulative incidence",fp,append=TRUE)
write.table(y$est,fp,sep="\t",quote=TRUE,row.names=FALSE,col.names=TRUE,append=TRUE)
write("Variance",fp,append=TRUE)
write.table(y$var,fp,sep="\t",quote=TRUE,row.names=FALSE,col.names=TRUE,append=TRUE)



# regression analysis of the effect of growth factors - GCSF versus MCSF
ftime=data$Event_Time
fstatus=data$Cause
cov1=data$Growth_Factor
regression=crr(ftime,fstatus,cov1)
v=summary(regression)
write("Regression performed on GFP expression",fp,append=TRUE)
write("Call",fp,append=TRUE);write(v$call,fp,append=TRUE)
write("Converged",fp,append=TRUE);write(v$converged,fp,append=TRUE)
write("n",fp,append=TRUE);write(v$n,fp,append=TRUE)
write("loglik",fp,append=TRUE);write(v$loglik,fp,append=TRUE)
write("coef",fp,append=TRUE)
write.table(v$coef,fp,sep="\t",quote=TRUE,append=TRUE)
write("Confidence Interval",fp,append=TRUE)
write.table(v$conf.int,fp,sep="\t",quote=TRUE,append=TRUE)
write("logtest",fp,append=TRUE)
write.table(v$logtest,fp,sep="\t",quote=TRUE,append=TRUE)

# compare proportional hazards model with nonparametric
Y=predict(regression,cov1=1:2)
plot(Y,xlab='days',ylab='Probability',main="Prediction of 
	cumulative incidence for onset of GFP fluorescence using regression model",
	col=1:2,lty=1)
legend(0,0.7, c("GCSF", "MCSF"),lty=1,col=1:2)
write("Predicted regression",fp,append=TRUE)
write("Columns: Time GCSF MCSF",fp,append=TRUE)
write.table(Y,fp,sep="\t",quote=TRUE,row.names=FALSE,col.names=FALSE,append=TRUE)

# perform regression on apoptosis

regression=crr(ftime,fstatus,failcode=2,cov1)
write("Regression performed on apoptosis",fp,append=TRUE)

v=summary(regression)
write("Call",fp,append=TRUE);write(v$call,fp,append=TRUE)
write("Converged",fp,append=TRUE);write(v$converged,fp,append=TRUE)
write("n",fp,append=TRUE);write(v$n,fp,append=TRUE)
write("loglik",fp,append=TRUE);write(v$loglik,fp,append=TRUE)
write("coef",fp,append=TRUE)
write.table(v$coef,fp,sep="\t",quote=TRUE,append=TRUE)
write("Confidence Interval",fp,append=TRUE)
write.table(v$conf.int,fp,sep="\t",quote=TRUE,append=TRUE)
write("logtest",fp,append=TRUE)
write.table(v$logtest,fp,sep="\t",quote=TRUE,append=TRUE)

# compare proportional hazards model with nonparametric
Y=predict(regression,cov1=1:2)
plot(Y,xlab='days',ylab='Probability',main="Predicted
	cumulative incidence of apoptosis",
	col=1:2,lty=1)
legend(0,0.7, c("GCSF", "MCSF"),lty=1,col=1:2)
write("Predicted regression",fp,append=TRUE)
write("Columns: Time GCSF MCSF",fp,append=TRUE)
write.table(Y,fp,sep="\t",quote=TRUE,row.names=FALSE,col.names=FALSE,append=TRUE)

# Perform competing risk regression using clustering as well
# First perform regression on GFP expression
library(crrSC)
clusteredregression=crrc(ftime,fstatus,cov1,cluster=data$Clone,failcode=1)
v=summary(clusteredregression)
write("Regression performed on mitosis including clustering on clone ID",fp,append=TRUE)
write("Call",fp,append=TRUE);write(v$call,fp,append=TRUE)
write("Converged",fp,append=TRUE);write(v$converged,fp,append=TRUE)
write("n",fp,append=TRUE);write(v$n,fp,append=TRUE)
write("loglik",fp,append=TRUE);write(v$loglik,fp,append=TRUE)
write("coef",fp,append=TRUE)
write.table(v$coef,fp,sep="\t",quote=TRUE,append=TRUE)
write("Confidence Interval",fp,append=TRUE)
write.table(v$conf.int,fp,sep="\t",quote=TRUE,append=TRUE)
write("logtest",fp,append=TRUE)
write.table(v$logtest,fp,sep="\t",quote=TRUE,append=TRUE)

clusteredregression=crrc(ftime,fstatus,cov1,cluster=data$Clone,failcode=2)
v=summary(clusteredregression)
write("Regression performed on apoptosis including clustering on clone ID",fp,append=TRUE)
write("Call",fp,append=TRUE);write(v$call,fp,append=TRUE)
write("Converged",fp,append=TRUE);write(v$converged,fp,append=TRUE)
write("n",fp,append=TRUE);write(v$n,fp,append=TRUE)
write("loglik",fp,append=TRUE);write(v$loglik,fp,append=TRUE)
write("coef",fp,append=TRUE)
write.table(v$coef,fp,sep="\t",quote=TRUE,append=TRUE)
write("Confidence Interval",fp,append=TRUE)
write.table(v$conf.int,fp,sep="\t",quote=TRUE,append=TRUE)
write("logtest",fp,append=TRUE)
write.table(v$logtest,fp,sep="\t",quote=TRUE,append=TRUE)





