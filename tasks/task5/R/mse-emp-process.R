#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(mydas)

dirMy="/home/laurence/Desktop/Dropbox/mydasOMs"
dirDat=file.path(dirMy, "data")
dirRes=file.path(dirDat,"results")

parD=NULL
for (i in seq(dim(scen)[1])){
  #drop_download(path=paste(file.path("mydasOMs",scen[i,"spp"]),".RData",sep=""),overwrite=T)
  load(file.path(dirRes,paste(scen[i,"spp"],".RData",sep="")))
  parD=rbind(parD,cbind(spp=scen[i,"spp"],model.frame(prior),model.frame(controlD)[,-4]))
}

empD$spp =ac(empD$spp)
empD     =subset(empD,year>=mseStart[spp]&year<=mseStart[spp]+45)
empD$year=empD$year-mseStart[empD$spp]

empd_pm=ddply(empD,.(spp,iter), smryStat)

empd_pm=merge(transform(empd_pm,iter=as.numeric(ac(iter)),spp=ac(spp)),
              transform(parD,   iter=as.numeric(ac(iter)),spp=ac(spp)),by=c("spp","iter"))
save(empD,empd_pm,controlD,file=file.path(dirRes,"empd-results.RData"))
