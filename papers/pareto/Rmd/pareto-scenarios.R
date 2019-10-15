dirOM ="/home/laurence/Desktop/sea++/mydas/project/data/OM"
dirRes="/home/laurence/Desktop/sea++/mydas/project/papers/pareto/data"

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)

library(FLife)
library(mydas)
library(randtests)

library(caret)
library(ggplot2)
library(mco)
library(kernlab)

source('~/Desktop/sea++/mydas/pkg/R/svr.R')

scen=expand.grid(nyr    =c(3,5,7)[2],
                 cv     =c(0.1,0.2,0.3)[2],
                 spp    ="turbot",
                 ctrl   =1:9,
                 stringsAsFactors=FALSE)

### Random variaton for control
nits=500
scen=expand.grid(nyr    =c(3,5,7),
                 cv     =c(0.1,0.2,0.3),
                 spp    ="turbot",
                 ctrl   =1:9,
                 stringsAsFactors=FALSE)

### Random variaton for control
controls=list()
set.seed(123456)
for (j in 1:9){
  controls[[j]]=rbind(FLPar(k1   =runif(nits, 0.0, 1.5)),
                      FLPar(k2   =runif(nits, 0.0, 1.5)),
                      FLPar(gamma=rep(1, nits)))}
controls=FLPars(controls)

## Performance measures
scen=expand.grid(nyr    =c(3,5,7),
                 cv     =c(0.1,0.2,0.3),
                 stringsAsFactors=FALSE)

pms<-foreach(i=seq(dim(scen)[1]), 
             .combine=rbind,
             .multicombine=TRUE,
             .export=c("dirRes","dirOM","ks"),
             .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                         "FLasher","FLBRP","FLife","randtests")) %dopar%{
                           
     source('~/Desktop/sea++/mydas/pkg/R/smryStat.R')
     source('~/Desktop/sea++/mydas/pkg/R/omOut.R')
                           
     fn<-function(cv,nyr){
             spp="turbot"
             load(file.path(dirRes,paste("base-",spp,".RData",sep="")))
                             
             load(file.path(dirRes,paste("ranD-final-",
                                     spp,"-",cv,"-",nyr,".RData",sep="")))
             mse=window(mse,start=80,end=100)
             res=transform(omSmry(mse,eq,lh),
                     iter=factor(iter,level=ac(sort(as.numeric(ac(unique(iter)))))))
             res=transform(res,year=year-79)
             pm =ddply(res,.(iter), smryStat)
                             
             pm}
                           
 with(scen[i,], cbind(cv=cv,nyr=nyr,fn(cv,nyr)))}

save(pms,file=file.path(dirRes,"pms.RData"))

svrs<-foreach(i=seq(dim(scen)[1]), 
             .combine=rbind,
             .multicombine=TRUE,
             .export=c("dirRes","dirOM","scen"),
             .packages=c("plyr","dplyr","reshape","ggplot2","FLCore","ggplotFL",
                         "FLasher","FLBRP","FLife","randtests",
                         "caret","ggplot2","mco","kernlab")) %dopar%{
                           
#for (i in seq(dim(scen)[1])){    
      source('~/Desktop/sea++/mydas/pkg/R/omOut.R')
                           
      predictors='k1 + k2'
      rewards=c('yield','yieldAav','safety','blim')
                           
                           
      svr=mlply(rewards,function(x) fitSvr(pm., x, predictors,
                       tuneLength=10)$finalModel)
      names(svr)=rewards
                           
      pareto=mco::nsga2(fn=svrFn2, svr[[1]],svr[[2]], svr[[3]],svr[[4]],idim=2, odim=4, 
                      lower.bounds=rep(0.0, 2),upper.bounds=rep(1.5, 2), 
                      popsize=1000, generations=100, cprob=0.9)
                           
      dat=data.frame(par=pareto$par,value=pareto$value)
      names(dat)=c("k1","k2",rewards)
                           
      pts=rbind(cbind("What"="good",subset(pm,plim>=0.95)[,rewards]),
                cbind("What"="bad", subset(pm,plim< 0.95)[,rewards]),
                cbind("What"="fits", dat[,rewards]))
                
      pts}                   

save(svrs,file=file.path(dirRes,"svrs.RData"))
