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
registerDoParallel(3)
b
library(devtools)
install_github("lauriekell/mydas")
library(mydas)

dirMy="/home/laurence/Desktop/Dropbox/mydasOMs"
#dirMy="/ichec/home/users/laurie"
dirDat=file.path(dirMy,"data")
dirRes=file.path(dirMy,"results")

### Stochasticity
set.seed(1234)
nits=500
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.2,b=0.0)
uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.3,b=0.0)
save(srDev,file=file.path(dirRes,"srDev.RData"))
save(uDev, file=file.path(dirRes,"uDev.RData"))

### Random variaton for control
control=list()
set.seed(123456)
for (j in 1:12){
  control[[j]]=rbind(FLPar(k1   =runif(nits, 0.0, 1.0)),
                     FLPar(k2   =runif(nits, 0.0, 1.0)),
                     FLPar(gamma=rep(1, nits)))}
control=FLPars(control)

save(control,file=file.path(dirRes,"dCtrl.RData"))

scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat"),
                 control=1:12,
                 stringsAsFactors=FALSE)

fn<-function(iSpp,iCtrl){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))
  
  mse=mseSBTD(om,eq,
              control=control[[iCtrl]],
              srDev  =srDev,uDev=uDev,
              start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46)
  
  mse=window(mse,start=mseStart[iSpp])
  
  save(mse,file=file.path(dirRes,paste("ranD-",iSpp,"-",iCtrl,".RData",sep="")))

  data.frame(spp=iSpp,control=iCtrl)}

ranD<-foreach(i=seq(dim(scen)[1]), 
               .combine=rbind,
               .multicombine=TRUE,
               .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                           "FLasher","FLBRP","FLife")) %dopar%{
                             
      with(scen[i,],fn(spp,control))}

scen=data.frame(spp    =c("brill","pollack","ray","sprat","sprat","turbot","turbot"),
           control=c(11,     11,       10,    9,     12,     9,       12),
           stringsAsFactors=FALSE)