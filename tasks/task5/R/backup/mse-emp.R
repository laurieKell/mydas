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

source('~/Desktop/sea++/mydas/pkg/R/omOut.R')
source('~/Desktop/sea++/mydas/pkg/R/smryStat.R')
source('~/Desktop/sea++/mydas/pkg/R/mseSBTD.R')
source('~/Desktop/sea++/mydas/pkg/R/hcrSBTD.R')

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

### Stochasticity
set.seed(1234)
nits=500

srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.2,b=0.0)
uDev =FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),0.3,b=0.0)

##### D ##########################################################
#### Comstant catch with k1=k2=0

scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster"),
                 stringsAsFactors=FALSE)

control=rbind(FLPar(k1   =rep(0.0, nits)),
              FLPar(k2   =rep(0.0, nits)),
              FLPar(gamma=rep(1.0, nits)))

constD=NULL
constD<-foreach(i=(seq(dim(scen)[1])), 
             .combine=rbind,
             .multicombine=TRUE,
             .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                         "FLasher","FLBRP","FLife")) %dopar%{

  load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))

  mse=mseSBTD(om,eq,
              control=control,
              srDev  =srDev,uDev=uDev,
              start  =mseStart[scen[i,"spp"]]+1,end=mseStart[scen[i,"spp"]]+46,nyrs=3)
  
  mse=mse[,ac(mseStart[scen[i,"spp"]]:(mseStart[scen[i,"spp"]]+46+2))]
  
  save(mse,file=file.path(dirRes,paste("constD-",scen[i,"spp"],".RData",sep="")))
  
  cbind(spp=scen[i,"spp"],omSmry(mse,eq,lh))}

save(constD,file=file.path(dirRes,"constD.RData"))

### Random variaton for control
control=list()
set.seed(123456)
for (j in 1:12){
  control[[j]]=rbind(FLPar(k1   =runif(nits, 0.0, 1.0)),
                     FLPar(k2   =runif(nits, 0.0, 1.0)),
                     FLPar(gamma=runif(nits, 1,   1)))}
save(controltD,file=file.path(dirRes,"dCtrl.RData"))

scen=expand.grid(spp    =c("brill",   "turbot",   "ray",   "pollack",   "sprat",   "razor",   "lobster")[4],
                 control=1:12,
                 stringsAsFactors=FALSE)

empD=NULL
empD<-foreach(i=(seq(dim(scen)[1])), 
                .combine=rbind,
                .multicombine=TRUE,
                .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                            "FLasher","FLBRP","FLife")) %dopar%{
                              
          load(file.path(dirDat,paste(scen[i,"spp"],".RData",sep="")))
                              
          mse=mseSBTD(om,eq,
                      control=control[[scen[i,"control"]]],
                      srDev  =srDev,uDev=uDev,
                      start  =mseStart[scen[i,"spp"]]+1,end=mseStart[scen[i,"spp"]]+46,nyrs=3)
                              
        mse=mse[,ac(mseStart[scen[i,"spp"]]:(mseStart[scen[i,"spp"]]+46+2))]
                              
        save(mse,file=file.path(dirRes,paste("constD-",i,".RData",sep="")))
                              
        cbind(spp=scen[i,"spp"],control=scen[i,"control"],omSmry(mse,eq,lh))}

###### Processing
empD     =res
empD$spp =ac(empD$spp)
empD     =subset(empD,year>=mseStart[spp]&year<=mseStart[spp]+45)
empD$year=empD$year-mseStart[empD$spp]

empd_pm=ddply(empD,.(spp,iter,j), smryStat)

dt=transform(empd_pm,k1=cut(k1,seq(0,1,0.03)),
                     k2=cut(k2,seq(0,1,0.03)))

dt=ddply(dt,.(k1,k2,spp),with,data.frame(
          safe =mean(safety),
          kobe =mean(kobe.n/45),
          yield=mean(yield),
          aav  =1-mean(yieldAav)))

ggplot(dt)+
  geom_tile(aes(k1,k2,fill=kobe))+
  scale_fill_gradientn(colours=c("navy","blue","cyan","lightcyan","yellow","red","red4"))+
  facet_wrap(.~spp)

empd_pm=merge(transform(empd_pm,iter=as.numeric(ac(iter)),spp=ac(spp),j=as.numeric(ac(j))),
              transform(parD,   iter=as.numeric(ac(iter)),spp=ac(spp),j=as.numeric(ac(j))),
              by=c("spp","iter","j"))

save(empD,empd_pm,controlD,file=file.path(dirRes,"empd-results-pollack.RData"))

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