#!/usr/bin/env Rscript

library(plyr)
library(dplyr)
library(reshape)
library(mydas)
library(FLCore)
library(ggplotFL)
library(FLBRP)
library(FLife)

dirMy="/home/laurence/Desktop/Dropbox/mydasOMs"
dirDat=file.path(dirMy, "data")
dirRes=file.path(dirDat,"results")

source('~/Desktop/sea++/mydas/pkg/R/omOut.R')

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

fl=c("/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-1.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-2.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-3.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-4.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-5.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-6.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-7.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-8.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-9.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-10.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-11.RData",
     "/home/laurence/Desktop/Dropbox/mydasOMs/results/ranD-y10-ray-12.RData")

load("/home/laurence/Desktop/Dropbox/mydasOMs/data/dCtrl.RData")

empD=NULL
for (i in c(1:4,6:7,9:10,12)){
  load(fl[i])
  empD=rbind(empD,cbind(iCtrl=i,omSmry(mse,eq,lh),model.frame(control[[i]])))}

empD     =subset(empD,year>=mseStart["sprat"]&year<=mseStart["sprat"]+45)[,-26]
empD$year=empD$year-mseStart["sprat"]
empD=transform(empD,iter=as.numeric(ac(iter))+(iCtrl-1)*500)

empd_pm=ddply(empD,.(iter,k1,k2), smryStat)




