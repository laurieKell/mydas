library(plyr)
library(dplyr)
library(reshape)

library(FLCore)
library(FLBRP)
library(FLasher)
library(ggplotFL)
library(FLife)
library(randtests)

library(mydas)

if (FALSE){
  library(doParallel)
  library(foreach)
  
  cl=makeCluster(3)
  registerDoParallel(cl)
  
  dirOM ="/home/laurence/Desktop/sea++/mydas/project/tasks/task4/data"
  dirRes="/home/laurence/Desktop/sea++/mydas/project/papers/pareto/data"
  }else{
  library(doMPI)
  cl <- startMPIcluster()
  registerDoMPI(cl)
  
  dirOM ="/rdsgpfs/general/user/lkell/home/mydas"
  dirRes=dirOM 
  }

nits=500

scen=expand.grid(nyr    =c(3,5,7),
                 cv     =c(0.1,0.2,0.3),
                 spp    ="turbot",
                 ctrl   =1:10,
                 stringsAsFactors=FALSE)

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,
           "sprat"=52,"razor"=54, "lobster"=57)

### Random variaton for control
controls=list()
set.seed(123456)
for (j in 1:10){
    controls[[j]]=rbind(FLPar(k1   =runif(nits, 0.0, 1.5)),
                        FLPar(k2   =runif(nits, 0.0, 1.5)),
                        FLPar(gamma=rep(1, nits)))}
controls=FLPars(controls)

set.seed(1233)
srDev=FLife:::rlnoise(nits,FLQuant(0,dimnames=list(year=1:105)),.3,b=0.0)

## Run MSE with random control
ranD<-foreach(i=seq(dim(scen[1:3,])[1]), 
              .combine=rbind,
              .multicombine=TRUE,
              .export=c("dirRes","dirOM","scen","srDev","controls","mseStart"),
              .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                          "FLasher","FLBRP","FLife","mydas")) %dopar%{
    
      fn<-function(iSpp,iCtrl,cv,nyr){
          load(file.path(dirRes,paste("base-",iSpp,".RData",sep="")))
                              
          set.seed(1235)
          uDev =FLife:::rlnoise(500,FLQuant(0,dimnames=list(year=1:105)),cv,b=0.0)
                      
          mse=mseSBTD(om,eq,
                      control=controls[[iCtrl]],
                      srDev  =srDev,uDev=uDev,
                      start  =mseStart[iSpp]+1,end=mseStart[iSpp]+46,
                      nyrs   =nyr)
                              
          save(mse,controls[[iCtrl]],file=file.path(dirRes,paste("pareto-ranD-",
                          iSpp,"-",iCtrl,"-",cv,"-",nyr,".RData",sep="")))
                              
          data.frame(spp=iSpp,control=iCtrl,cv=cv,nyr=nyr)}
    
  with(scen[i,], fn(spp,ctrl,cv,nyr))}

pm=mdply(scen,function(spp,ctrl,nyr,cv,
                       start=1,end=40){
  
  load(file.path(dirRes,"base-turbot.RData"))
  load(file.path(dirRes,paste("pareto-ranD-",
                               iSpp,"-",iCtrl,"-",cv,"-",nyr,".RData",sep="")))
  
  mse=window(mse,start=mseStart[spp]+start,end=mseStart[spp]+end)
  res=transform(omSmry(mse,eq,lh),
                iter=factor(iter,level=ac(sort(as.numeric(ac(unique(iter)))))))
  res=transform(res,year=year-start)
  pm =ddply(res,.(iter), smryStat)
  pm=merge(pm,transform(model.frame(control[[ctrl]]),
                        iter=factor(iter,level=ac(sort(as.numeric(ac(iter)))))),
           by="iter")
  pm})

pm$kobe    =pm$kobe.n/40
pm$yieldAav=(1-pm$yieldAav)

save(pm,file=file.path(dirRes,"pareto-turbot-pm.RData"))

