eqSegreg<-function(object){
  
  sr=as.FLSR(object,model="segreg")
  
  lower(sr)[1:2]=c(c(min(rec(sr),na.rm=TRUE)/max(ssb(sr),na.rm=TRUE)),
                   c(min(ssb(sr),na.rm=TRUE)))
  upper(sr)[1:2]=c(c(max(rec(sr),na.rm=TRUE)/min(ssb(sr),na.rm=TRUE)),
                   c(max(ssb(sr),na.rm=TRUE)))
  sr=fmle(sr,control=list(trace=FALSE),method="L-BFGS-B")
  
  params(sr)["a"][is.na(params(sr)["a"])]=median(params(sr)["a"],na.rm=TRUE)
  params(sr)["b"][is.na(params(sr)["b"])]=median(params(sr)["b"],na.rm=TRUE)

  res=brp(FLBRP(object,sr=sr))

  return(res)}

eqBevholt<-function(object){

  sr=as.FLSR(object,model="bevholt")
  lower(sr)[1:2]=0
  sr=fmle(sr,control=list(trace=FALSE),method="L-BFGS-B")
  
  params(sr)["a"][is.na(params(sr)["a"])]=median(params(sr)["a"],na.rm=TRUE)
  params(sr)["b"][is.na(params(sr)["b"])]=median(params(sr)["b"],na.rm=TRUE)

  res=brp(FLBRP(object,sr=sr))
  
  return(res)}

eqGeomean<-function(object){
  
  sr=as.FLSR(object,model="geomean")
  sr=fmle(sr,control=list(trace=FALSE),method="L-BFGS-B")
  
  params(sr)["a"][is.na(params(sr)["a"])]=median(params(sr)["a"],na.rm=TRUE)
  
  res=brp(FLBRP(object,sr=sr))
  
  return(res)}
