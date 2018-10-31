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

library(mydas)

#sessionInfo()

dirMy="/home/laurence/Desktop/Dropbox/mydasOMs"
dirDat=file.path(dirMy, "data")
dirRes=file.path(dirMy,"results")

mseStart=            c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)
scen=expand.grid(spp=c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster"),
                 j  =1:4, 
                 stringsAsFactors=FALSE)

### Stochasticity
set.seed(1234)
nits=500

srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.2,b=0.0)
uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.3,b=0.0)

##### D ##########################################################

controlD=list()
for (j in 1:4){
k1=runif(nits, 0.0, 1.0)
k2=runif(nits, 0.0, 1.0)
controlD[[j]]=rbind(FLPar(k1   =k1),
               FLPar(k2   =k2),
               FLPar(gamma=runif(nits, 1, 1)))}

empD=NULL
empD<-foreach(i=(seq(dim(scen)[1])), 
             .combine=rbind,
             .multicombine=TRUE,
             .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                         "FLasher","FLBRP","FLife")) %dopar%{

  load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
  om=FLCore:::iter(om,seq(nits))
  eq=FLCore:::iter(eq,seq(nits))
  lh=FLCore:::iter(lh,seq(nits))
  
  mse=mseSBTD(om,eq,
              control=controlD[[scen[i,"j"]]],
              srDev=srDev,uDev=uDev,
              start  =mseStart[scen[i,"spp"]]+5,end=mseStart[scen[i,"spp"]]+45,nyrs=5)
  res=cbind(spp=scen[i,"spp"],j=scen[i,"j"],omSmry(mse,eq,lh))
  
  save(res,file=file.path(dirRes,paste("empd-",i,".RData",sep="")))
  
  res}


parD=NULL
for (i in seq(dim(scen)[1])){
  #drop_download(path=paste(file.path("mydasOMs",scen[i,"spp"]),".RData",sep=""),overwrite=T)
  load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
  parD=rbind(parD,cbind(spp=scen[i,"spp"],j=scen[i,"j"],model.frame(prior),model.frame(controlD[[scen[1,"j"]]])[,-4]))
  }

empD$spp =ac(empD$spp)
empD     =subset(empD,year>=mseStart[spp]&year<=mseStart[spp]+45)
empD$year=empD$year-mseStart[empD$spp]

empd_pm=ddply(empD,.(spp,iter,j), smryStat)
empd_pm=merge(transform(empd_pm,iter=as.numeric(ac(iter)),spp=ac(spp),j=as.numeric(ac(j))),
              transform(parD,   iter=as.numeric(ac(iter)),spp=ac(spp),j=as.numeric(ac(j))),
              by=c("spp","iter","j"))

save(empD,empd_pm,controlD,file=file.path(dirRes,"empd-results.RData"))

##### P ##########################################################
set.seed(1234)
controlP=rbind(FLPar(k1=runif(nits, 0.0,1.0)),
               FLPar(k2=runif(nits, 0.0,1.0)))

empP<-foreach(i=(seq(dim(scen)[1])), 
              .combine=rbind,
              .multicombine=TRUE,
              .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                          "FLasher","FLBRP","FLife")) %dopar%{
  
  load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
  #om=FLCore:::iter(om,seq(nits))
  #eq=FLCore:::iter(eq,seq(nits))
  #lh=FLCore:::iter(lh,seq(nits))

  res =mseSBTP(om,eq,
               control =controlP,
               start   =mseStart[scen[i,"spp"]],end=mseStart[scen[i,"spp"]]+45,
               srDev   =srDev,uDev=uDev,
               refU    =39:41,
               refCatch=23:25)
  
  save(res,file=file.path(dirRes,paste("empp-",i,".RData",sep="")))
  
  cbind(spp=scen[i,"spp"],omSmry(res,eq,lh))
  }

empP=NULL
for (i in 1:7){
  load(paste("/home/laurence/Desktop/Dropbox/mydasOMs/data/",scen[i,"spp"],".RData",sep=""))
  load(paste("/home/laurence/Desktop/Dropbox/mydasOMs/results/empp-",i,".RData",sep=""))
  
  empP=rbind(empP,cbind(spp=scen[i,"spp"],omSmry(res,eq,lh)))
  }

parP=NULL
for (i in seq(dim(scen)[1])){
  #drop_download(path=paste(file.path("mydasOMs",scen[i,"spp"]),".RData",sep=""),overwrite=T)
  load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
  parP=rbind(parP,cbind(spp=scen[i,"spp"],model.frame(prior),model.frame(controlP)[,-3]))
  }

empP$spp =ac(empP$spp)
empP     =subset(empP,year>=mseStart[spp]&year<=mseStart[spp]+45)
empP$year=empP$year-mseStart[empP$spp]
empp_pm=ddply(empP,.(spp,iter), smryStat)

empp_pm=merge(transform(empp_pm,iter=as.numeric(ac(iter)),spp=ac(spp)),
              transform(parP,   iter=as.numeric(ac(iter)),spp=ac(spp)),by=c("spp","iter"))

save(empP,empp_pm,controlP,file=file.path(dirRes,"empp-results.RData"))


#https://stats.stackexchange.com/questions/160479/practical-hyperparameter-optimization-random-vs-grid-search

#There is a simple probabilistic explanation for the result: for any distribution over a sample space with a finite maximum, the maximum of 60 random observations lies within the top 5% of the true maximum, with 95% probability. That may sound complicated, but it’s not. Imagine the 5% interval around the true maximum. Now imagine that we sample points from his space and see if any of it lands within that maximum. Each random draw has a 5% chance of landing in that interval, if we draw n points independently, then the probability that all of them miss the desired interval is (1−0.05)n

#So the probability that at least one of them succeeds in hitting the interval is 1 minus that quantity. We want at least a .95 probability of success. To figure out the number of draws we need, just solve for n in the equation:
  
  1−(1−0.05)^n>0.95

We get n⩾60

. Ta-da!
  
  The moral of the story is: if the close-to-optimal region of hyperparameters occupies at least 5% of the grid surface, then random search with 60 trials will find that region with high probability.