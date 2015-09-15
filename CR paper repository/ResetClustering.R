ResetClustering<-function(ClusterIDs) {
  mask<-is.na(ClusterIDs)
  n<-1
  for (i in 1:length(ClusterIDs)) {
    ndx<-(ClusterIDs[i]==ClusterIDs)&(!mask)
    mask[ndx]<-TRUE # dont do it again
    ClusterIDs[ndx]<-n
    if (sum(ndx)>0) n<-n+1
  }
  return(ClusterIDs)
}