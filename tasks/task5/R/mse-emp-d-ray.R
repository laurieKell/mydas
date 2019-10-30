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

library(devtools)
install_github("lauriekell/mydas")
library(mydas)

dirMy="/home/laurence/Desktop/Dropbox/mydasOMs"
#dirMy="/ichec/home/users/laurie"
dirDat=file.path(dirMy,"data")
dirRes=file.path(dirMy,"results")

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

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

scen=expand.grid(spp    ="ray",
                 control=1:12,
                 stringsAsFactors=FALSE)

fn<-function(iSpp,iCtrl){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))
  
  mse=mseSBTD(om,eq,
              control=control[[iCtrl]],
              srDev  =srDev,uDev=uDev,
              start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46,
              nyrs=10)
  
  mse=window(mse,start=mseStart[iSpp])
  
  save(mse,file=file.path(dirRes,paste("ranD-y10-",iSpp,"-",iCtrl,".RData",sep="")))

  data.frame(spp=iSpp,control=iCtrl)}

ranD<-foreach(i=c(5,8,11), #seq(dim(scen)[1]), 
               .combine=rbind,
               .multicombine=TRUE,
               .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                           "FLasher","FLBRP","FLife")) %dopar%{
                             
      with(scen[i,],fn(spp,control))}

