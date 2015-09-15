getMothers<-function(MothersAndDaughters) {
  # find mothers from dataframe
  # returns boolean vector
  v<-MothersAndDaughters$PairID
  gen<-MothersAndDaughters$Generation
  n<-length(v)
  is.mother<-(v==0)
  currentID=v[1]
  for (i in 2:n) {
    if (v[i]!=currentID) {
      currentID<-v[i]
    }
    else {
      if (gen[i-1]<gen[i]) is.mother[i-1]=TRUE
      else is.mother[i]=TRUE      
    }
    
#     if (v[i]!=currentID) {
#       currentID=v[i]
#     }
#     else if (gen[i]>gen[i-1]) is.mother[i]=TRUE
#       }
#     }
  }
  return(is.mother)
}