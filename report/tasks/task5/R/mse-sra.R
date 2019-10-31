library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)
library(mpb)
library(mydas)

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

#source('~/Desktop/flr/mpb/R/mseMPB.R')
#source("/home/laurence/Desktop/flr/mpb/R/hcr.R")
#source('~/Desktop/flr/mpb/R/mseSRA.R')
#source('~/Desktop/flr/mpb/R/biodyn-sra.R')
#source('~/Desktop/flr/FLife/R/omOut.R')


theme_set(theme_bw())

dirMy ="/home/laurence/Desktop/sea++/mydas/tasks/task4"
dirDat=file.path(dirMy,"data")

## OM
load(file.path(dirDat,"turbot.RData"))

om=FLCore:::iter(window(om,start=20,end=90),1:500)

## MP
mp=mpb:::setMP(as(window(om,end=55),"biodyn"),
         r =median(prior["r"],na.rm=T),
         k =median(prior["v"],na.rm=T),
         b0=0.8,
         p =median(mpb:::p(prior["bmsy"]/prior["v"]),na.rm=TRUE))

nits=dims(om)$iter
set.seed(1234)
 srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.3,b=0.0)

eq=FLCore:::iter(eq,seq(nits))

#mse=mseSRA(om,eq,mp,start=55,ftar=0.5,srDev=srDev)

scen=expand.grid(stock=c("turbot","lobster","ray","pollack","razor","brill","sprat"),
                  ftar =c(1.0,0.7),
                  btrig=c(0.5,0.6))

runMPB<-function(stock,ftar,btrig) {
  load(file.path(dirDat,paste(stock,".RData",sep="")))

  om=window(om,start=20)
  nits=dim(om)[6]
  eq=FLCore:::iter(eq,seq(nits))
  
  set.seed(1234)
  srDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.3,b=0.0)
  uDev  =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.2,b=0.0)
  selDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100,age=dimnames(m(om))$age)),0.2,b=0.0)

  mp=mpb:::setMP(as(window(om,end=55),"biodyn"),
                 r =median(prior["r"],na.rm=T),
                 k =median(prior["v"],na.rm=T),
                 b0=0.8,
                 p =median(mpb:::p(prior["bmsy"]/prior["v"]),na.rm=TRUE))
  
  res=mseSRA(om,eq,mp,
              start=mseStart[stock],end=mseStart[stock]+40,
              ftar=ftar,btrig=btrig,
              srDev=srDev)
  
  res2=data.frame("stock"=stock,"ftar"=ftar,"btrig"=btrig,omSmry(res,eq,lh[c("a","b")]))
  
  res2}

srampb=NULL
for (i in seq(dim(scen)[1])){
  
  res=with(scen[i,],runMPB(stock,ftar,btrig))
  
  srampb=rbind(srampb,cbind(scen=i,res))
  
  save(srampb,file="/home/laurence/Desktop/sea++/mydas/tasks/task5/data/srampb.RData")}
