library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)

library(devtools)
devtools::install_github("lauriekell/mydas", subdir="pkgs/mydas")

library(mydas)

sessionInfo()

dirMy ="/home/laurence/Desktop/sims/wklife"
dirDat=file.path(dirMy,"data")
dirRes=file.path(dirMy,"results")

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

### Stochasticity
nits=500
set.seed(1234)
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:500)),0.2,b=0.0)

### OEM
uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:500)),0.3,b=0.0)

## MSE for Derivate empirical MP
scen=expand.grid(spp=c("turbot","lobster","ray","pollack","razor","brill","sprat")[1],
                 k1=seq(1.5,2.5,0.5),k2=seq(1.5,2.5,0.5),gamma=seq(1.0,1.25,0.25),
                 stringsAsFactors=FALSE)

empD=NULL
for (i in seq(dim(scen)[1])){
  load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
  om=iter(om,seq(nits))
  eq=iter(eq,seq(nits))
  lh=iter(lh,seq(nits))
  
  res =mseSBTD(om,eq,control=with(scen[i,],c(k1=k1,k2=k2,gamma=gamma)),
               start=mseStart[scen[i,"spp"]],end=mseStart[scen[i,"spp"]]+45,srDev=srDev,uDev=uDev)
  empD=rbind(empD,cbind(scen=i,stock=scen[i,"spp"],k1r=scen[i,"k1"],k2=scen[i,"k2"],gamma=scen[i,"gamma"],omSmry(res,eq,lh)))
  
  save(empD,file=file.path(dirRes,"empD.RData"))}

## MSE for Proportion empirical MP
scen=expand.grid(spp=c("turbot","lobster","ray","pollack","razor","brill","sprat")[1],
                 k1=c(0.50,0.25),k2=c(0.50,0.25),
                 stringsAsFactors=FALSE)
empP=NULL
for (i in seq(dim(scen)[1])){
  load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
  om=iter(om,seq(nits))
  eq=iter(eq,seq(nits))
  lh=iter(lh,seq(nits))

  res =mseSBTP(om,eq,control=with(scen[i,],c(k1=k1,k2=k2)),
               start=mseStart[scen[i,"spp"]],end=mseStart[scen[i,"spp"]]+45,
               srDev=srDev,uDev=uDev,refYr=40)
  empP=rbind(empP,cbind(scen=i,spp=scen[i,"spp"],k1=scen[i,"k1"],k2=scen[i,"k2"],omSmry(res,eq,lh)))
  
  save(empP,file=file.path(dirRes,"empP.RData"))}

names(empP)[10]="catchJuv"
names(empP)[2]="spp"
pm   =ddply(empP,.(spp,iter,k1,k2),smryStat)
empP.=subset(empP,iter%in%101:500)
pm.  =subset(pm,  iter%in%101:500)

dbWriteTable(conLK, "mydas_empp",    value=empP., append=!FALSE,overwrite=!TRUE,row.names=FALSE)
dbWriteTable(conLK, "mydas_empp_pm", value=pm.,   append=!FALSE,overwrite=!TRUE,row.names=FALSE)

