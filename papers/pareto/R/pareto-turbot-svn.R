library(plyr)
library(dplyr)
library(reshape)

library(FLCore)
library(FLBRP)
library(FLasher)
library(ggplotFL)
library(FLife)
library(randtests)

library(mydas)

library(caret)
library(ggplot2)
library(mco)
library(kernlab)

source('~/Desktop/sea++/mydas/pkg/R/svr.R')

load("/home/laurence/Desktop/cx1/mydas/pareto-turbot-pm.RData")

if (!FALSE){
  library(doParallel)
  library(foreach)
  
  cl=makeCluster(3)
  registerDoParallel(cl)
  
  dirOM ="/home/laurence/Desktop/cx1/mydas"
  dirRes=dirOM
  }

nits=500

scen=expand.grid(nyr    =c(3,5,7),
                 cv     =c(0.1,0.2,0.3),
                 spp    ="turbot",
                 stringsAsFactors=FALSE)

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,
           "sprat"=52,"razor"=54, "lobster"=57)

### Random variaton for control
controls=list()
set.seed(123456)
for (j in 1:10){
  controls[[j]]=rbind(FLPar(k1   =runif(nits, 0.0, 1.5)),
                      FLPar(k2   =runif(nits, 0.0, 1.5)),
                      FLPar(gamma=rep(1, nits)))}
controls=FLPars(controls)

set.seed(1233)
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),.3,b=0.0)

pareto=ddply(pm,.(cv,nyr), with, {

  reward   =c('safety','blim','kobe','yield','yieldAav')
  predictor='k1 + k2'
  
  dat=data.frame(k1=k1,k2=k2,
                 safety=safety,kobe    =kobe,blim=blim,
                 yield =yield, yieldAav=yieldAav)                 
  
  svr=mlply(reward,function(x) fitSvr(dat, x, predictor,tuneLength=10)$finalModel)
  
  res=mco::nsga2(fn=svrFn2, svr[[1]],svr[[2]],svr[[3]],svr[[4]],svr[[5]],
             idim=2, odim=5, 
             lower.bounds=rep(0.0, 2),upper.bounds=rep(1.5, 2), 
             popsize=1000, generations=100, cprob=0.9)
  
  res=data.frame(par=res$par,value=res$value)
  names(res)=c("k1","k2",'safety','blim','kobe','yield','yieldAav')
  res})

save(pareto,file=file.path(dirRes,"pareto-turbot-pareto.RData"))


