ProbandwisePlot<-function(CrossOddsRatio,marg.cif1,marg.cif2=NULL,title="Probandwise probability",pos=c(0,1),
                          eventtxt1="T1",eventtxt2="T2") {
  v<-summary(CrossOddsRatio) # stats for CrossOddsRatio
  if (is.null(marg.cif2)) {
    marg.cif2<-marg.cif1
    sym=TRUE
  }
  else {
    sym=FALSE
  }
  if (sym!=v$sym) warning(paste("sym should be ",sym))
  con<-concordance(CrossOddsRatio,cif1=marg.cif1$P1,cif2=marg.cif2$P1)

  plot(marg.cif2$time,con$intercept$probandwise[,1],type="s",
       main=title,ylim=c(0,1),lwd=2,xlim=c(0,max(marg.cif2$time)),
       xlab="Time (days)",ylab="Probability")
  lines(marg.cif2$time,con$intercept$probandwise[,2],type="s",lty=2)
  lines(marg.cif2$time,con$intercept$probandwise[,3],type="s",lty=2)
  lines(marg.cif2$time,marg.cif2$P1,type="s",col=2,lwd=2)
  
  if (sym) {
    txt1<-paste("P[(",eventtxt1,"<t) given (",eventtxt2,"<t)]")
    txt2<-paste("P[(",eventtxt1,"<t)]")
    text(pos[1],pos[2],paste("log(cor)=",format(CrossOddsRatio$theta,digits=3),"+/-",
                             format(sqrt(CrossOddsRatio$var.theta),digits=3),", symmetric,"," p=",
                             format(v$estimates[4],digits=3)),pos=4) 
  }
  else {
    txt1<-paste("P[(",eventtxt2,"<t) given (",eventtxt1,"<t)]")
    txt2<-paste("P[(",eventtxt2,"<t)]")
    text(pos[1],pos[2],paste("log(cor)=",format(CrossOddsRatio$theta,digits=3),"+/-",
                             format(sqrt(CrossOddsRatio$var.theta),digits=3),", asymmetric,"," p=",
                             format(v$estimates[4],digits=3)),pos=4)
  }
  legend(pos[1],pos[2],legend=c(txt1,"95% CI",txt2),col=c(1,1,2),
         lwd=c(2,1,2),lty=c(1,2,1),bty="n")
  
}