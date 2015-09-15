getMothersGFPStatus<-function(P0) {
  # get mother's GFP status
  P0$IsMotherGFPPos<-NA
  for (i in 1:length(P0$Generation)) {
    if (P0$Generation[i]>0) {
      MotherProgenyID<-P0$Progeny[i]%/%2
      # search for mother
      ndx<-(MotherProgenyID==P0$Progeny)&(P0$Group[i]==P0$Group)&(P0$Clone[i]==P0$Clone) # match group and clone
      P0$IsMotherGFPPos[i]<-P0$GFP[ndx]
    }
  }
  return(P0)
}
