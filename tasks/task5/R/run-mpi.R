#srun -p DevQ -N 1 -A gmlif003b -t 1:00:00 –-pty bash

library(plyr,    lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(dplyr,   lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(reshape, lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(ggplot2, lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")

library(FLCore,  lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(ggplotFL,lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(FLasher, lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(FLBRP,   lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(FLife,   lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")

library(mydas,   lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")

library(doMPI,   lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(doRNG,   lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")

dirMy="/ichec/home/users/laurie"
dirDat=file.path(dirMy,"data")
dirRes=file.path(dirMy,"results")

cl=startMPIcluster(count = 39)  # where count is number of processors -1
registerDoMPI(cl)
clusterSize(cl)                  ## this just prints to confirm your cluster is up and running


load(file.path(dirDat,"uDev.RData"))
load(file.path(dirDat,"srDev.RData"))
load(file.path(dirDat,"dCtrl.RData"))

scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat"),
                 control=1:8,
                 stringsAsFactors=FALSE)

fn<-function(iSpp,iCtrl){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))
  
  mse=mseSBTD(om,eq,
              control=control[[iCtrl]],
              srDev  =srDev,uDev=uDev,
              start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46)
  
  mse=mse[,ac(mseStart[iSpp]:min(mseStart[iSpp]+46+2),dims(mse)$maxyear)]
  
  save(mse,file=file.path(dirRes,paste("gridD-",iSpp,"-",iCtrl,".RData",sep="")))

  data.frame(spp=iSpp,control=iCtrl)}

ranD<-foreach(i=seq(dim(scen)[1]), 
               .combine=rbind,
               .multicombine=TRUE,
               .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                           "FLasher","FLBRP","FLife")) %dorng% {
                             
      with(scen[i,],fn(spp,comtrol))}

closeCluster(cl)
mpi.quit()
