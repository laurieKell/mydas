library(FLCore)
library(ggplotFL)
library(FLBRP)
library(FLasher)
library(FLife)
library(mydas)
library(popbio)
library(spatstat)

library(plyr)
library(dplyr)
library(reshape)
library(GGally)

library(LBSPR)

library(doParallel)
library(foreach)

registerDoParallel(3)

dirDat="/home/laurence/Desktop/sea++/mydas/project/papers/ROC/data"

load("/home/laurence/Desktop/sea++/mydas/project/papers/ROC/data/lhs.RData")

## Scenarios
design=data.frame(s      =c( 0.9, 0.7, 0.9, 0.9, 0.9, 0.9),
                  sel2   =c( 1.0, 1.0,5000, 1.0, 1.0, 1.0),
                  sel3   =c(5000,5000,5000,  50,5000,5000),
                  nsample=c( 500, 500, 500, 500, 250, 500),
                  m      =c(rep("gislason",5),"constant"))
design=mdply(expand.grid(Stock=names(lhs),CV=c("0.3","0.5","AR"),stringsAsFactors=FALSE), 
             function(Stock,CV) design)

f=FLQuant(c(rep(0.1,60),seq(0.1,2.5,length.out=40)[-40],
            seq(2.5,1,length.out=11),rep(1,20)))

## Stochasticity
nits=100
set.seed(234)
srDev=FLQuants(NULL)
srDev[["0.3"]]=rlnoise(nits,f%=%0,0.3,0.0);set.seed(234)
srDev[["0.5"]]=rlnoise(nits,f%=%0,0.5,0.0);set.seed(234)
srDev[["AR"]] =rlnoise(nits,f%=%0,0.3,0.7)

## Runs
res<-foreach(i=dimnames(design)[[1]][(13:18)], 
              .combine=rbind,
              .multicombine=TRUE,
              .export=c("dirDat","design","srDev","lhs","f"),
              .packages=c("FLCore","ggplotFL","FLBRP","FLasher",
                          "FLife","mydas","popbio","spatstat",
                          "plyr","dplyr","reshape","GGally","LBSPR")) %dopar% {
                            
#for (i in dimnames(design)[[1]][-(13:18)]){
                            
  source('~/Desktop/sea++/mydas/pkg/R/lbspr.R')
  source('~/Desktop/sea++/mydas/pkg/R/oemLn.R')
                            
  ## modify parameters                                                  
  par        =lhs[[design[i,"Stock"]]]
  par["s"]   =design[i,"s"]
  par["sel2"]=design[i,"sel2"]
  par["sel3"]=design[i,"sel3"]
  
  ## simulated stock
  if (design[i,"m"]=="gislason")
     eq=lhEql(par,m=function(x,params) {
           x=wt2len(stock.wt(x),params)
           gislason(x,params)})
  else
     eq=lhEql(par,m=function(x,params) {
                           x=stock.wt(x)
                           x[]=params["linf"]
                           gislason(x,params)})
  
  fbar(eq)=f%*%refpts(eq)["msy","harvest"]
  om      =as(eq,"FLStock")
  om      =propagate(om,nits)
  om      =fwd(om,fbar=fbar(eq)[,2:130],
                  sr  =eq,residuals=srDev[[design[i,"CV"]]])
  
  ## length frequencies
  ak   =invAlk(par,cv=0.1)  
  lopt =c(2/3*par["linf"])
  prior=popdyn(par)
  
  ## Indicators
  ## catch
  lfd=lenSample(catch.n(om)[,50:125],ak,nsample=design[i,"nsample"])
  
  ### model based
  tmp=capture.output({
    lb=lbspr(lfd,prior)
    save(lb,file=file.path(dirDat,"sims","catch",paste("lb",i,"RData",sep=".")))
  },type="message")

  ### empirical
  ind =ddply(subset(as.data.frame(lfd,drop=TRUE),data>0), .(year,iter), with, 
             lenInd(len,data,lopt))
  save(ind,file=file.path(dirDat,"sims","catch",paste("ind",i,"RData",sep=".")))
  
  ## Stock
  lfd=lenSample(stock.n(om)[,50:125],ak,nsample=design[i,"nsample"])
  
  ### model based
  tmp=capture.output({
    lb=lbspr(lfd,prior)
    save(lb,file=file.path(dirDat,"sims","stock",paste("lb",i,"RData",sep=".")))
  },type="message")
  
  ### empirical
  ind =ddply(subset(as.data.frame(lfd,drop=TRUE),data>0), .(year,iter), with, 
             lenInd(len,data,lopt))
  save(ind,file=file.path(dirDat,"sims","stock",paste("ind",i,"RData",sep=".")))

  design[i,]} 
