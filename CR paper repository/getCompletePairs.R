getCompletePairs<-function(PairID){
  n=length(PairID)
  v<-sort(PairID)
  v_order<-order(PairID)
  b<-(PairID==0.4)
  newb<-b
  ndx<-v[1]
  for (i in 2:n){
    if (v[i]!=ndx) ndx<-v[i] 
    else {
      b[i]=TRUE
      b[i-1]=TRUE
      v[i]<-i-1
      v[i-1]<-i-1
    }
  }
  for (i in 1:n) {
    m<-v_order[i]
    newb[m]<-b[i]
    PairID[m]<-v[i]
  }
  out<-list(newb,PairID)
  return(out)
}
  