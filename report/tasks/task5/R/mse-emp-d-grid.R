library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)

library(doParallel)
library(foreach)
registerDoParallel(4)

#library(devtools)
#install_github("lauriekell/mydas")
library(mydas)

dirMy="/home/laurence/Desktop/Dropbox/mydasOMs"
#dirMy="/ichec/home/users/laurie"
dirDat=file.path(dirMy,"data")
dirRes=file.path(dirMy,"results")

source('~/Desktop/sea++/mydas/pkg/R/mseSBTD.R')
source('~/Desktop/sea++/mydas/pkg/R/hcrSBTD.R')

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

### Stochasticity
set.seed(1234)
nits=500
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.2,b=0.0)
uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.3,b=0.0)
save(srDev,file=file.path(dirRes,"srDev.RData"))
save(uDev, file=file.path(dirRes,"uDev.RData"))

### grid for control
set.seed(123456)
control=FLPars("0"=rbind(FLPar(k1   =rep(0.0, nits)),
                         FLPar(k2   =rep(0.0, nits)),
                         FLPar(gamma=rep(1,   nits))),
               "1"=rbind(FLPar(k1   =rep(0.5, nits)),
                         FLPar(k2   =rep(0.5, nits)),
                         FLPar(gamma=rep(1,   nits))),
               "2"=rbind(FLPar(k1   =rep(0.5, nits)),
                         FLPar(k2   =rep(1.0, nits)),
                         FLPar(gamma=rep(1.0, nits))),
               "3"=rbind(FLPar(k1   =rep(1.0, nits)),
                         FLPar(k2   =rep(0.5, nits)),
                         FLPar(gamma=rep(1,   nits))),
               "4"=rbind(FLPar(k1   =rep(1.0, nits)),
                         FLPar(k2   =rep(1.0, nits)),
                         FLPar(gamma=rep(1,   nits))))

save(control,file=file.path(dirRes,"gCtrl.RData"))


scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster"),
                 iCtrl=1:5,
                 stringsAsFactors=FALSE)

fn<-function(iSpp,iCtrl){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))

  mse=mseSBTD(om,eq,
            control=control[[iCtrl]],
            srDev  =srDev,uDev=uDev,
            start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46)

  save(mse,file=file.path(dirRes,paste("gridD-",iSpp,"-",iCtrl,".RData",sep="")))
  data.frame(spp=iSpp,control=iCtrl)}

## missing
gridD<-foreach(i=c(3,10,17,24),#seq(dim(scen)[1]), 
                .combine=rbind,
                .multicombine=TRUE,
                .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                            "FLasher","FLBRP","FLife")) %dopar%{
                            
        with(scen[i,],fn(spp,control))}


## Ray with 10 year index
scen=expand.grid(spp =c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster")[3],
                 iCtrl=1:5,
                 stringsAsFactors=FALSE)

fn<-function(iSpp,iCtrl){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))
  
  mse=mseSBTD(om,eq,
              control=control[[iCtrl]],
              srDev  =srDev,uDev=uDev,
              start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46,
              nyrs=10)
  
  save(mse,file=file.path(dirRes,paste("gridD-",iSpp,"-10yr-",iCtrl,".RData",sep="")))
  data.frame(spp=iSpp,control=iCtrl)}

gridD<-foreach(i=seq(dim(scen)[1]), 
               .combine=rbind,
               .multicombine=TRUE,
               .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                           "FLasher","FLBRP","FLife")) %dopar%{
                             
                             with(scen[i,],fn(spp,iCtrl))}

## Sprat with recruit index
scen=expand.grid(spp  =c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster")[5],
                 iCtrl=1:5,
                 stringsAsFactors=FALSE)

fn<-function(iSpp,iCtrl){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))
  
  mse=mseSBTD(om,eq,
              control=control[[iCtrl]],
              srDev  =srDev,uDev=uDev,
              start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46,
              nyrs=5,
              cpueFn =function(x) apply(stock.n(x)%*%(1-mat(x))%*%stock.wt(x),c(2,6),sum))
  
  save(mse,file=file.path(dirRes,paste("gridD-",iSpp,"-uJuv-",iCtrl,".RData",sep="")))
  data.frame(spp=iSpp,control=iCtrl)}

gridD<-foreach(i=seq(dim(scen)[1]), 
               .combine=rbind,
               .multicombine=TRUE,
               .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                           "FLasher","FLBRP","FLife")) %dopar%{
                             
                             with(scen[i,],fn(spp,iCtrl))}


load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-sprat-uJuv-1.RData")
sprat=FLStocks("k1=0; k2=0"=mse)
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-sprat-uJuv-2.RData")
sprat["k1=0.5; k2=0.5"]=mse
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-sprat-uJuv-3.RData")
sprat["k1=0.5; k2=1.0"]=mse
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-sprat-uJuv-4.RData")
sprat["k1=1.0; k2=0.5"]=mse
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-sprat-uJuv-5.RData")
sprat["k1=1.0; k2=1.0"]=mse
sprat=FLStocks(llply(sprat,window,end=100))
plot(sprat)

load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-ray-10yr-1.RData")
ray=FLStocks("k1=0; k2=0"=mse)
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-ray-10yr-2.RData")
ray["k1=0.5; k2=0.5"]=mse
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-ray-10yr-3.RData")
ray["k1=0.5; k2=1.0"]=mse
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-ray-10yr-4.RData")
ray["k1=1.0; k2=0.5"]=mse
load("/home/laurence/Desktop/Dropbox/mydasOMs/results/gridD-ray-10yr-5.RData")
ray["k1=1.0; k2=1.0"]=mse
ray=FLStocks(llply(ray,window,end=104))
plot(ray)

fn<-function(om){
  as.data.frame(FLQuants(om, 
                         "ssb" = function(x) ssb(x)%/%refpts( eq)["msy","ssb"], 
                         "f" =   function(x) fbar(x)%/%refpts(eq)["msy","harvest"], 
                         "rec" = function(x) rec(x)%/%refpts( eq)["msy","rec"], 
                         "catch"=function(x) landings(x)%/%refpts(eq)["msy","yield"]))}

load(file.path(dirDat,"sprat.RData"))
sprat2=cbind(spp="sprat",ldply(sprat,fn))
sprat2=subset(sprat2,year>mseStart["sprat"])
sprat2$year=sprat2$year-mseStart["sprat"]

load(file.path(dirDat,"ray.RData"))
ray2=cbind(spp="ray",ldply(ray,fn))
ray2=subset(ray2,year>mseStart["ray"])
ray2$year=ray2$year-mseStart["ray"]

smry2=rbind(sprat2,ray2)
save(smry2,file=file.path(dirRes,"smry2.RData"))

