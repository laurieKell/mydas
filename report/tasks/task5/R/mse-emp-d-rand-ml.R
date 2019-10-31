library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)

if (FALSE){
  library(devtools)
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("flr/FLasher",  force=TRUE)
  devtools::install_github("flr/FLFishery",       force=TRUE)
  devtools::install_github("lauriekell/mydas-pkg",force=TRUE)
}

library(FLasher)
library(FLBRP)
library(FLife)

library(r4ss)
library(stringr)
library(plyr)
library(dplyr)
library(mydas)

#source('~/Desktop/sea++/mydas/pkg/R/smryStat.R')
#source('~/Desktop/sea++/mydas/pkg/R/omOut.R')

##Parallel stuff
if (!FALSE){
  library(doParallel)
  library(foreach)
  
  cl=makeCluster(3)
  registerDoParallel(cl)
  
  dirMy="/home/laurence/Desktop/sea++/mydas/project/tasks/task5/results"
  dirOM="/home/laurence/Desktop/sea++/mydas/project/tasks/task4/data"
}else{
  library(doMPI)
  cl <- startMPIcluster()
  registerDoMPI(cl)
  
  dirMy="/rdsgpfs/general/user/lkell/ephemeral/mydas"
  dirOM="/rdsgpfs/general/user/lkell/home/mydas"
  }

setwd(dirMy)

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

### Stochasticity
set.seed(1233)
nits=500
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.3, b=0.0)

### Random variaton for control
control=list()
set.seed(123456)
for (j in 1:12){
    control[[j]]=rbind(FLPar(k1   =runif(nits, 0.0, 1.5)),
                       FLPar(k2   =runif(nits, 0.0, 1.5)),
                       FLPar(gamma=rep(1, nits)))}
control=FLPars(control)
  
scen=expand.grid(spp    =c("turbot",   "ray",   "pollack",   "sprat"),
                 control=1:12,
                 nyr    =c(3,5,7),
                 cv     =c(0.1,0.2,0.3),
                 stringsAsFactors=FALSE)

ranD<-foreach(i=seq(dim(scen)[1])[1], 
              .combine=rbind,
              .multicombine=TRUE,
              .export=c("dirMy","scen","srDev","uDev","control","mseStart","fn"),
              .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                          "FLasher","FLBRP","FLife","mydas")) %dopar%{
    
      fn<-function(iSpp,iCtrl,cv,nyr){
          load(file.path(dirOM,paste(iSpp,".RData",sep="")))
                              
          set.seed(1235)
          uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),cv,b=0.0)
                      
          mse=mseSBTD(om,eq,
                      control=control[[iCtrl]],
                      srDev  =srDev,uDev=uDev,
                      start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46,
                      nyrs   =nyr)
                              
          save(mse,file=file.path(dirMy,paste("ranD-",iSpp,"-",iCtrl,"-",cv,"-",nyr,".RData",sep="")))
                              
          data.frame(spp=iSpp,control=iCtrl,cv=cv,nyr=nyr)}
    
      with(scen[i,], fn(spp,control,cv,nyr))}

process=function(spp,ctrl,nyr,cv,
                 dirOM,dirRes=dirMy,
                 start=1,end=40){
  
  load(file.path(dirOM,paste(spp,".RData",sep="")))
  load(file.path(dirRes,"dCtrl.RData"))
  load(file.path(dirRes,paste("ranD-",spp,"-",ctrl,"-",cv,"-",nyr,".RData",sep="")))
  
  mse=window(mse,start=mseStart[spp]+start,end=mseStart[spp]+end)
  res=transform(omSmry(mse,eq,lh),
                iter=factor(iter,level=ac(sort(as.numeric(ac(unique(iter)))))))
  ctl=transform(model.frame(control[[ctrl]]),
                iter=factor(iter,level=ac(sort(as.numeric(ac(iter))))))
  res=transform(res,year=year-start)
  pm =ddply(res,.(iter), smryStat)
  pm=merge(ctl,pm,by="iter")
  pm=merge(pm,model.frame(lh),by="iter")
  pm}

pm=mdply(scen[1,],function(spp,control,nyr,cv) 
       process(spp,ctrl,nyr,cv,dirOM,dirMy))

save(pm,file=file.path(dirMy,"pm.RData"))
