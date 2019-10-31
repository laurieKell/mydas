dirOM ="/home/laurence/Desktop/sea++/mydas/project/data/OM"
dirRes="/home/laurence/Desktop/sea++/mydas/project/papers/pareto/data"

library(ggplot2)
library(plyr)
library(dplyr)
library(reshape)

library(FLife)
library(mydas)
library(randtests)

source('~/Desktop/sea++/mydas/pkg/R/smryStat.R')
source('~/Desktop/sea++/mydas/pkg/R/omOut.R')

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
nits    =500
controls=list()
set.seed(123456)
for (j in 1:9){
  controls[[j]]=rbind(FLPar(k1   =runif(nits, 0.0, 1.5)),
                      FLPar(k2   =runif(nits, 0.0, 1.5)),
                      FLPar(gamma=rep(1, nits)))}
controls=FLPars(controls)

## Performance measures
pm=mdply(scen,function(spp,ctrl,nyr,cv,
                       start=76,end=100){
  
  load(file.path(dirRes,"base-turbot.RData"))
  load(file.path(dirRes,paste("base-ranD-",
                              spp,"-",ctrl,"-",cv,"-",nyr,".RData",sep="")))
  
  mse=window(mse,start=start,end=end)
  res=transform(omSmry(mse,eq,lh),
                iter=factor(iter,level=ac(sort(as.numeric(ac(unique(iter)))))))
  res=transform(res,year=year-start)
  pm =ddply(res,.(iter), smryStat)
  pm=merge(pm,transform(model.frame(controls[[ctrl]]),
                        iter=factor(iter,level=ac(sort(as.numeric(ac(iter)))))),
           by="iter")
  pm})

pm$kobe    =pm$kobe.n/25
#pm$yieldAav=(1-pm$yieldAav)

save(pm,file=file.path(dirRes,"base-turbot-pm.RData"))

###### Pareto Frontiers
predictors='k1 + k2'
rewards=c('yield','yieldAav','kobe','safety','blim','plim')

svr=mlply(rewards,function(x) fitSvr(pm, x, predictors,
                                     tuneLength=10)$finalModel)
names(svr)=rewards

pareto=mco::nsga2(fn=svrFn2, svr[[1]],svr[[2]], svr[[3]],svr[[4]],
                  idim=2, odim=4, 
                  lower.bounds=rep(0.0, 2),upper.bounds=rep(1.5, 2), 
                  popsize=1000, generations=100, cprob=0.9)

save(pareto,svr,file=file.path(dirRes,"base-turbot-pareto.RData"))
