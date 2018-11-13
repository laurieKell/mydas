library(plyr) #lib.loc="/ichec/home/users/laurie/R/x86_64-pc-linux-gnu-library/3.5")
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)

library(mydas)

library(doMPI)
library(doRNG)

library(rslurm)

dirMy="/ichec/home/users/laurie"
dirDat=file.path(dirMy,"data")
dirRes=file.path(dirMy,"results")

load(file.path(dirDat,"uDev.RData"))
load(file.path(dirDat,"srDev.RData"))
load(file.path(dirDat,"dCtrl.RData"))

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat"),
                 control=1:8,
                 stringsAsFactors=FALSE)

#Pass function and the parameters data frame to `slurm_apply`,
#specifiying the number of cluster nodes to use and the number of CPUs per node.
#The latter (`cpus_per_node`) determines how many processes will be forked on
#each node, as the `mc.cores` argument of `parallel::mcMap`.

fn<-function(spp,control){
  
  load(file.path(dirDat,paste(iSpp,".RData",sep="")))
  
  #mse=mseSBTD(om,eq,
  #            control=control[[iCtrl]],
  #            srDev  =srDev,uDev=uDev,
  #            start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46)
  #mse=mse[,ac(mseStart[iSpp]:min(mseStart[iSpp]+46+2),dims(mse)$maxyear)]

  mse=fwd(om,catch=catch(om)[,ac(60:100)],sr=eq,residuals=srDev)
  
  save(mse,file=file.path(dirRes,paste("gridD-",iSpp,"-",iCtrl,".RData",sep="")))
  
  mse}

#!/bin/sh
#
#SBATCH --job-name=randomD
#SBATCH --output=slurm_%a.out
#SBATCH --nodes=1
#SBATCH --time=00:20:00
#SBATCH -A gmlif003b
#SBATCH -p DevQ

#module load r

#RScript slurm_run.R

sopt=list(nodes=1, time="00:20:00") #,"A gmlif003b","p DevQ")

sjob<-slurm_apply(fn, scen, jobname='randomD',
                    nodes=1, cpus_per_node=40, submit=FALSE, slurm_options=sopt)
