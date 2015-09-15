ClusterMothersAndDaughters<-function(P0,IsMotherAndDaughters) {
  if (IsMotherAndDaughters) {
    # create mother and daughter pairs
    k<-1
    IsMother<-NULL
    ClusterID<-NULL
    Row<-NULL
    # Find all mothers and daughters firstly with generation>0
    for (i in 1:length(P0$Progeny)) {
      if ((P0$Generation[i]>1)) {
        ClusterID<-c(ClusterID,k)
        IsMother<-c(IsMother,FALSE)
        Row<-c(Row,i)
        # find mother in same clone and group and copy label
        ndx<-(P0$Clone==P0$Clone[i])&(P0$Group==P0$Group[i])&(P0$Progeny==(P0$Progeny[i]%/%2))
        ClusterID<-c(ClusterID,k)
        IsMother<-c(IsMother,TRUE)
        Row<-c(Row,which(ndx))
        k<-k+1
      } 
    }
    # mark clusters ot Keep. No replicated rows, otherwise we are replicating data!
    Keep<-vector(length=length(Row))
    Keep[1:length(Keep)]<-NA
    for (i in 1:length(Row)) {
      if (is.na(Keep[i])) {
        # find rows for cluster Pair
        ndx<-((ClusterID==ClusterID[i])&(is.na(Keep)))
        if (sum(ndx)==2) {
          Keep[ndx]<-TRUE
          # block out all other rows with these indices
          n<-which(ndx)
          BlockNdx<-((Row[n[1]]==Row))|((Row[n[2]]==Row))
          BlockNdx[ndx]<-FALSE # don't block select positions to Keep
          Keep[BlockNdx]<-FALSE
        }
      }  
    }
    # Create a new DataFrame  for MothersAndDaughters
    P0$IsMother<-FALSE
    P0$ClusterID<-0
    P0$Keep<-FALSE
    ndx<-vector(length=length(P0$Generation))
    MothersAndDaughters<-P0[ndx,] #create empty data frame

    NewEntry<-P0[1,]
    for (i in 1:length(Row)) {
      NewEntry<-P0[Row[i],]
      NewEntry$IsMother<-IsMother[i]
      NewEntry$ClusterID<-ClusterID[i]
      NewEntry$Keep<-Keep[i]
      MothersAndDaughters<-rbind(MothersAndDaughters,NewEntry)
    }

    return(MothersAndDaughters)
  }
  
  else {
    # writes Sister Pair
    P0$SisterClusterID<-NA
    k<-1
    # label even progeny and ancestor
    for (i in 1:length(P0$Progeny)) {
      if ((P0$Progeny[i]%%2==0)) {
        P0$SisterClusterID[i]<-k   
        # find odd sister and copy label
        ndx<-(P0$Clone==P0$Clone[i])&(P0$Group==P0$Group[i])&(P0$Progeny==(P0$Progeny[i]+1))
        P0$SisterClusterID[ndx]<-k
        k<-k+1
      }    
    }   
    return(P0)
  }
  
}