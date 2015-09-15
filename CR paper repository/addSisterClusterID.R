addSisterClusterID<-function(data) {
  # adds extra field call SisterClusterID
  data$SisterClusterID<-0
  k<-1
 
  # label even progeny and ancestor
  for (i in 1:length(data$Progeny)) {
    if ((data$Progeny[i]==1)|(data$Progeny[i]%%2==0)) {
      data$SisterClusterID[i]<-k
      k<-k+1
    }    
  }
  # copy label to sister
  for (i in 1:length(data$Progeny)) {
    if (data$SisterClusterID[i]==0) {# find sister
      ndx<-(((data$Progeny[i]-1)==data$Progeny)&(data$Clone[i]==data$Clone)&(data$Group[i]==data$Group))
      if (sum(ndx)==1) {
        data$SisterClusterID[i]<-data$SisterClusterID[ndx]
      }
    }
  }
  # fill in remaining odd cells (no sisters)
  for (i in 1:length(data$Progeny)) {
    if (data$SisterClusterID[i]==0) {
      data$SisterClusterID[i]<-k
      k<-k+1
    }
  }
  return(data)
}