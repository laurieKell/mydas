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

#source('~/Desktop/sea++/mydas/pkg/R/omOut.R')
#source('~/Desktop/sea++/mydas/pkg/R/smryStat.R')
#source('~/Desktop/sea++/mydas/pkg/R/mseSBTD.R')
#source('~/Desktop/sea++/mydas/pkg/R/hcrSBTD.R')

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

##### D ##########################################################
#### Comstant catch with k1=k2=0

scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster"),
                 stringsAsFactors=FALSE)

control=rbind(FLPar(k1   =rep(0.0, nits)),
              FLPar(k2   =rep(0.0, nits)),
              FLPar(gamma=rep(1.0, nits)))

constD=NULL
constD<-foreach(i=3, #(seq(dim(scen)[1])), 
                .combine=rbind,
                .multicombine=TRUE,
                .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                            "FLasher","FLBRP","FLife")) %dopar%{
                              
        load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
                              
        mse=mseSBTD(om,eq,
                    control=control,
                    srDev  =srDev,uDev=uDev,
                    start  =mseStart[scen[i,"spp"]]+1,end=mseStart[scen[i,"spp"]]+46,nyrs=3)
                              
        mse=mse[,ac(mseStart[scen[i,"spp"]]:(min(mseStart[scen[i,"spp"]]+46+2,dims(mse)$maxyear)))]
                              
        save(mse,file=file.path(dirRes,paste("constD-",scen[i,"spp"],".RData",sep="")))
        
        
        cbind(spp=scen[i,"spp"],omSmry(mse,eq,lh))}

save(constD,file=file.path(dirRes,"constD.RData"))
