getMotherAndDaughterpairs<-function(tempdata,EventID,CauseID,CauseValue) {
  v<-NULL
  ndx<-names(tempdata)==EventID # Event ID
  times<-tempdata[,ndx]
  ndx<-names(tempdata)==CauseID
  causes<-tempdata[,ndx]
  PairID<-tempdata$PairID
  is.mother<-tempdata$Mothers
  for (i in 1:length(causes)) {
    ndx<-(PairID==PairID[i])
    if (sum(causes[ndx]==CauseValue)==2) {
      causes[ndx]<-0 # don't select again
      mother<-is.mother[ndx]
      t<-times[ndx]
      # daughter is x, mother is y
      pair<-c(t[!mother],t[mother])
      v<-rbind(v,pair)  
    }
  }
  return (v)
}
