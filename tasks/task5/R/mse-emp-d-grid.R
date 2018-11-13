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

### Stochasticity
set.seed(1234)
nits=500
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.2,b=0.0)
uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.3,b=0.0)
save(srDev,file=file.path(dirRes,"srDev.RData"))
save(uDev, file=file.path(dirRes,"uDev.RData"))

### grid for control
set.seed(123456)
control=FLPars("1"=rbind(FLPar(k1   =rep(0.5, nits)),
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

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)
scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster"),
                 control=1:4,
                 stringsAsFactors=FALSE)

fn<-function(iSpp,iCtrl){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))

  mse=mseSBTD(om,eq,
            control=control[[iCtrl]],
            srDev  =srDev,uDev=uDev,
            start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46)

  save(mse,file=file.path(dirRes,paste("gridD-",iSpp,"-",iCtrl,".RData",sep="")))
  data.frame(spp=iSpp,control=iCtrl)}

gridD<-foreach(i=c(3,10,17,24),#seq(dim(scen)[1]), 
                .combine=rbind,
                .multicombine=TRUE,
                .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                            "FLasher","FLBRP","FLife")) %dopar%{
                            
        with(scen[i,],fn(spp,control))}

