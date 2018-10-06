fbSmry<-function(fb){
  
  names(fb)[c(14:17)]=c("l50","l50min","l50max","a50")
  fb=fb[,c("species","order","family","linf","k","t0","a","b","a50","l50","l50min","l50max")]
  
  fb[is.na(fb["l50"]),"l50"]=(fb[is.na(fb["l50"]),"l50min"]+fb[is.na(fb["l50"]),"l50max"])/2
  fb[fb$t0>=0,"t0"]=NA
  
  fb$species=as.character(fb$species)
  fb$species=factor(as.character(fb$species))
  fb$family=factor(as.character(fb$family))
  
  fb=subset(fb,!is.na(linf))
  
  smn=ddply(fb,.(species,order),with,data.frame(
    k   =mean(k,na.rm=T),
    linf=mean(linf,na.rm=T),
    t0  =mean(t0,na.rm=T),
    a   =mean(a,na.rm=T),
    b   =mean(b,na.rm=T),
    a50 =mean(a50,na.rm=T),
    l50 =mean(l50,na.rm=T)))
  
  sn =ddply(fb,.(species,order),with,data.frame(
    k   =sum(!is.na(k)),
    linf=sum(!is.na(linf)),
    t0  =sum(!is.na(t0)),
    a   =sum(!is.na(a)),
    b   =sum(!is.na(b)),
    a50 =sum(!is.na(a50)),
    l50 =sum(!is.na(l50))))
  
  svr=ddply(fb,.(species,order),with,data.frame(
    k   =var(k,na.rm=T),
    linf=var(linf,na.rm=T),
    t0  =var(t0,na.rm=T),
    a   =var(a,na.rm=T),
    b   =var(b,na.rm=T),
    a50 =var(a50,na.rm=T),
    l50 =var(l50,na.rm=T)))
  
  pmn=FLife:::mf2FLPar(smn[,-(1:2)])
  
  dmns=dimnames(pmn)
  dmns[["quant"]]=c("mean","se","var","n")
  dmns=dmns[c(1,3,2)]
  
  fb=FLPar(array(NA,dim=unlist(laply(dmns,length)),dimnames=dmns))
  dimnames(fb)$iter=svr[,"species"]
  
  fb[,"mean"]=FLife:::mf2FLPar(smn[,-(1:2)])
  fb[,"var"]=FLife:::mf2FLPar(svr[,-(1:2)])
  fb[,"n"]  =FLife:::mf2FLPar(sn[,-(1:2)])
  fb[,"se"] =(fb[,"var"]/(fb[,"n"]-1))^0.5
  
  fbvcov=var(model.frame(fb[,"mean"])[,1:7],na.rm=TRUE)
  fbcor =FLPar(cov2cor(fbvcov))
  
  FLPars(par=fb,cor=fbcor)}
