SampleClones<-function(tempdata)
  # Uses tempdata a dataset with paired Relatedness scoring to select cells pairs
  # without replication of progeny in pairs
{
  CloneIDs<-unique(tempdata$Clone)
  SampledClones<-data.frame()
  for (i in CloneIDs) {
    ndx<-(i==tempdata$Clone)
    CloneData<-tempdata[ndx,]
    print(paste("Clone ", format(i),": "))
    PairIDs<-unique(CloneData$PairID)
    PairIDs<-permute(PairIDs)
    nd<-dim(PairIDs)    
    if (length(nd)>0) {
      PairIDs<-RemoveDuplicatedProgeny(PairIDs,CloneData)
      if (length(PairIDs)>0) {

        Scored<-ScorePairIDs(PairIDs,CloneData)
        
        # find first best score
        if (sum(is.na(Scored$Entropy))==0) {
          if (sum(Scored$Entropy==0)>0) {
            ndx<-cumsum(max(Scored$N)==Scored$N)>0
            ndx<-cumsum(ndx)==1
          }                              # get ndx of first best score
          BestPairIDs<-PairIDs[ndx]    
          BestClonePairIDs<-is.element(CloneData$PairID,unlist(BestPairIDs))
          NewCloneData<-CloneData[BestClonePairIDs,]
          SampledClones<-rbind(SampledClones,NewCloneData)        
        }
      }
      
    }
  }
  return(SampledClones)
}


ScorePairIDs<-function(PairIDs,CloneData) {
  # Would like to make sure maximum number of Relatedness types 
  n<-length(PairIDs)
  Entropy<-vector(mode="numeric",length=n)
  NTypes<-vector(mode="numeric",length=n)
  N<-vector(mode="numeric",length=n)
  for (i in 1:n) {
    #calculate entropy
    if (i>length(PairIDs)) browser()
    ndx<-is.element(CloneData$PairID,PairIDs[[i]])
    Relatedness<-CloneData$Relatedness[ndx]
    types<-unique(Relatedness)
    m<-length(types)
    t<-vector(mode="numeric",m)
    for (j in 1:m) {  

      t[j]<-sum(is.element(types[j],Relatedness))

    }
    Entropy[i]<--sum((t/sum(t))*log2(t/sum(t))) # add 1 to the entropy
    NTypes[i]<-m    
    N[i]<-length(PairIDs[[i]])
    
   
  }
  Score<-list(Entropy=Entropy,NTypes=NTypes,N=N)
  return(Score)
}

RemoveDuplicatedProgeny<-function(PairIDs,CloneData)
  # Removes PairIDs from permutation which have duplicated progeny
  # Returns a list of PairIDs 
{
  n<-dim(PairIDs);nr<-n[1];nc<-n[2]
  newPairIDs<-list()
  for (i in 1:nr) {
    ndx<-rep(TRUE,nc)
    for (j in 1:(nc-1)) {
      if (ndx[j]) {        
        b1<-CloneData$PairID==PairIDs[i,j]
        p1<-CloneData$Progenitor[b1]
        for (k in (j+1):nc) {        
          b2<-CloneData$PairID==PairIDs[i,k]
          p2<-CloneData$Progenitor[b2]
          if (sum(is.element(p1,p2))>0) ndx[k]<-FALSE
        }
      }        
    }
    # only keep PairIDs with two elements
    tempPairIDs<-PairIDs[i,ndx]
    b<-rep(FALSE,length(tempPairIDs))
    for (j in length(tempPairIDs)) {
      b[j]<-(sum(CloneData$PairID==tempPairIDs[j])==2)
    }
    tempPairIDs<-tempPairIDs[b]
    if (length(tempPairIDs)>0) {
      newPairIDs[[length(newPairIDs)+1]]<-tempPairIDs
      b<-is.element(CloneData$PairID,unlist(PairIDs[i,ndx]))
      Progenitors<-CloneData$Progenitor[b]
      if (length(Progenitors)>length(unique(Progenitors))) {
        print("error")
        browser()      
      }
    }
  }
  return(newPairIDs)
}


permute<-function(sequence)
  # finds all possible permutations of sequence
  # and generates random permutations if the length of the sequence
  # is greater than 5
{
  if (length(sequence)>1) {
    n<-length(sequence)
    if (n<6) {
      temp<-matrix(ncol=1,nrow=n)
      temp[,1]<-sequence # first column
      for (i in 2:n) {
        m<-dim(temp)
        m<-m[1] # number of rows in temp
        newtemp<-NULL
        for (j in 1:m) {
          remainder<-setdiff(sequence,temp[j,])
          for (k in 1:length(remainder)) {
            d<-dim(newtemp)
            if (length(d)>0) {
              if (length(c(temp[j,],remainder[k]))!=d[2]) browser()
            }
            newtemp<-rbind(newtemp,c(temp[j,],remainder[k]))
          }        
        }
        temp<-newtemp
      }
    } # get all permutations
    else {
      m=factorial(5)
      temp=matrix(nrow=m,ncol=n)
      for (i in 1:m) {
        temp[i,]<-sample(sequence,n)
        if (sum(is.na(temp[i,]))>0) browser()
      }
    }
  } else temp=sequence  
  return(temp)
}

SamplePairs<-function(tempdata,number) {
  # just takes random sample of pairs. Assumes pairs are unique
  # used to shorten regression fitting 
  PairIDs<-sample(unique(tempdata$PairID),number)
  b<-is.element(tempdata$PairID,PairIDs)
  out<-tempdata[b,]
  return(out)
}
  
  
  
