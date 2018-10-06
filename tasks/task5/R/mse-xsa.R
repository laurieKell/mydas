library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)
library(FLAssess)
library(FLXSA)
library(mpb)
library(mydas)

theme_set(theme_bw())

dirMy ="/home/laurence/Desktop/sea++/mydas/tasks/task4"
dirDat=file.path(dirMy,"data")
source('~/Desktop/sea++/mydas/pkgs/mydas/R/omOut.R')

## sets up intial MP 
xsaMP<-function(om,pg=10,ctrl=xsaControl){
  stk  =setPlusGroup(om,pg)
  idx  =FLIndex(index=stock.n(stk))
  range(idx)[c("plusgroup","startf","endf")]=c(pg,0.1,.2)
  stk+FLXSA(stk,idx,control=ctrl,diag.flag=FALSE)}

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

source('~/Desktop/sea++/mydas/pkgs/mydas/R/mseXSA.R')

runXSA<-function(stock,ftar=1,bpa=0.5,sigma=0.3) {
  load(file.path(dirDat,paste(stock,".RData",sep="")))
  load(file.path(dirDat,"xsaCtrl.RData"))
  #om  =iter(om,1:10)
  
  range(om)[c("minfbar","maxfbar")]=ceiling(mean(lh["a1"]))
  range(eq)[c("minfbar","maxfbar")]=ceiling(mean(lh["a1"]))
  
  nits=dim(om)[6]
  om  =window(om,start=20)
  eq  =iter(eq,seq(nits))
  lh  =iter(lh,seq(nits))
  
  set.seed(1234)
  srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100)),0.3,b=0.0)
  uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:100,age=dimnames(m(om))$age)),0.2,b=0.0)
  
  if (stock=="ray")
    mp=xsaMP(trim(window(om,end=60),age=3:40),ctrl=xsaCtrl[[stock]],pg=ifelse(stock=="ray",15,10))
  else
    mp=xsaMP(window(om,end=60),ctrl=xsaCtrl[[stock]],pg=ifelse(stock=="ray",15,10))
  
  res=mseXSA(om,
             eq,
             mp,control=xsaCtrl[[stock]],
             ftar=ftar,bpa=bpa,sigma=sigma,
             interval=1,start=60,end=90,
             srDev=srDev,uDev=uDev)
  
  res=data.frame("stock"=stock,"ftar"=ftar,"bpa"=bpa,omSmry(res,eq,lh[c("a","b")]))
  
  res}

scen=expand.grid(stock=c("turbot","lobster","ray","pollack","razor","brill","sprat"),
                 ftar=c(1.0),bpa=c(0.5))

xsa=NULL
for (i in seq(dim(scen)[1])){
  res=with(scen[i,],runXSA(stock,ftar,bpa))
  
  xsa=rbind(xsa,res)
  
  save(xsa,file="/home/laurence/Desktop/sea++/mydas/tasks/task5/data/xsa.RData")}


