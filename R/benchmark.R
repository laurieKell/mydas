library(FLCore)
library(FLAssess)
library(FLXSA)

library(plyr)

library(doParallel)
library(foreach)

load("/home/laurence/Desktop/tmp/om.RData")

ctrl=FLXSA.control(tol    =1e-09,maxit   =50, 
                   min.nse=0.3,  fse     =0.1, 
                   rage   =0,       qage    =10, 
                   shk.n  =TRUE,   shk.f   =TRUE, 
                   shk.yrs=4,    shk.ages=5, 
                   window =100,   tsrange =10, 
                   tspower=0,
                   vpa    =FALSE)

xsa<-function(om,pg=10){
  stk=setPlusGroup(om,pg)
  idx=FLIndex(index=stock.n(stk))
  range(idx)[c("plusgroup","startf","endf")]=c(pg,0.1,.2)
  
  xsa=stk+FLXSA(stk,idx,control=ctrl,diag.flag=FALSE)
  
  xsa}

cat(NULL,file="bmk.txt")
for (iClust in c(1,4,8,16,24,48,72)){

  registerDoParallel(iClust)

  start_time <- Sys.time()
  rmse<-foreach(i=rev(seq(72)), 
              .combine=rbind,
              .multicombine=TRUE,
              .packages=c("FLCore","FLAssess","FLXSA","plyr")) %dopar%{
                
          stk=window(om,end=100-i)
          res=xsa(stk)
                  
          mean(((computeStock(res)-stock(stk))/stock(stk))^2)
          }
    cat(iClust,Sys.time()-start_time,"\n",append=TRUE,file="bmk.txt")  
    }
