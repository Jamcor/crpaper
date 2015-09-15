PlotPairedData<-function(Sisters,MothersAndDaughters,CauseID=1,TypeOfEvent="Division") {
  # plots paired data using ClusterID
  # ignores incomplete clusters
  # Returns Pearson correlation coefficient
  # Cause ID is 1, division, 2 death, 0 right censored
  ndx<-(!is.na(Sisters$SisterClusterID))
  tempdata<-Sisters[ndx,]
  x<-NULL
  y<-NULL
  for(i in 1:length(tempdata$SisterClusterID)) {
    ndx<-(tempdata$SisterClusterID[i]==tempdata$SisterClusterID)
    Event_Time<-tempdata$Event_Time[ndx]
    Cause<-tempdata$Cause[ndx]
    if ((length(Event_Time)==2)&(sum(Cause==CauseID)==2)&(tempdata$SisterClusterID[i]>0)) {
      x<-c(x,Event_Time[1])
      y<-c(y,Event_Time[2])
      tempdata$SisterClusterID[ndx]=0 # zero pair now that it has been recorded
    }
    else if (tempdata$SisterClusterID[i]>0){
      tempdata$SisterClusterID[ndx]=0 # zero pair now 
    }
  }
  SisterStats<-cor.test(x,y,method="pearson")
  
  r<-max(max(x),max(y))
  xrange<-c(min(x),max(x))
  regline<-lm(y~-(1)+x)
  yrange<-regline$coefficients[[1]]*xrange
  lx<-c(0,r)
  ly<-c(0,r)
  plot(x,y,xlab="Sibling 1 age (days)",ylab="Sibling 2 age (days)",
       type="p",pch=20,col="black",main=TypeOfEvent,,xlim=lx,ylim=ly)
  
  lines(lx,ly,type="l",lty=2,lwd=1,col="black")
  lines(xrange,yrange,type="l",lty=1,lwd=1,col="red")
  txt<-paste("r = ",format(SisterStats$estimate[[1]],digits=3),", p = ",
             format(SisterStats$p.value,digits=3))
  text(0,0.9*r,txt,pos=4)

  ndx<-(!is.na(MothersAndDaughters$ClusterID))
  tempdata<-MothersAndDaughters[ndx,]
  x<-NULL
  y<-NULL
  for(i in 1:length(tempdata$ClusterID)) {
    ndx<-(tempdata$ClusterID[i]==tempdata$ClusterID)
    Event_Time<-tempdata$Event_Time[ndx]
    Cause<-tempdata$Cause[ndx]
    Generation<-tempdata$Generation[ndx]
    if ((sum(Cause==2)==1)&(sum(Cause==1)==1)&(tempdata$ClusterID[i]>0)) { # only plot mother and daughter who divide: Cause 1 and 2
      if ((Generation[1]>0)&(Generation[2]>0)) { # only plot if generation 1 or greater
        x<-c(x,Event_Time[(Cause==2)]) # x axis is daughter
        y<-c(y,Event_Time[(Cause==1)]) # y axis in mother              
        tempdata$ClusterID[ndx]=0 # zero pair now that it has been recorded        
      }
      
    }
    else if (tempdata$ClusterID[i]>0){
      tempdata$SisterClusterID[ndx]=0 # zero pair now 
    }
  }
  MothersAndDaughtersStats<-cor.test(x,y,method="pearson")
  regline<-lm(y~-(1)+x)
  
  r<-max(max(x),max(y))
  xrange<-c(min(x),max(x))
  yrange<-regline$coefficients[[1]]*xrange
  lx<-c(0,r)
  ly<-c(0,r)
  plot(x,y,xlab="Daughter age (days)",ylab="Mother age (days)",
       type="p",pch=20,col="black",main=TypeOfEvent,,xlim=lx,ylim=ly)
  
  lines(lx,ly,type="l",lty=2,lwd=1,col="black")
  lines(xrange,yrange,type="l",lty=1,lwd=1,col="red")
  txt<-paste("r = ",format(MothersAndDaughtersStats$estimate[[1]],digits=3),", p = ",
             format(MothersAndDaughtersStats$p.value,digits=3))
  text(0,0.9*r,txt,pos=4)
  stats<-c(SisterStats,MothersAndDaughtersStats)
  return(stats)
}