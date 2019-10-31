library(fishmethods)
library(FLCore)
library(ggplotFL)

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/brill.RData")
om1=iter(om,1)

catch =as.data.frame(om1@catch)
a50   =ceiling(iter(lh,1)["a50"][[1]])
k     =refpts(eq)["virgin","biomass",1]
b1k   =1
btk   =c(ssb(om)[,75,,,,1]/k)
fmsysm=log(c(prior["fm",1]))
bmsyk =c(prior["bmsy",1]/prior["v",1])
M     =c(log(apply(m(eq),6,mean)[[1]]))

res=dbsra(
      year   =catch$year, 
      catch  =catch$data, 
      catchCV=0.3, 
      catargs=list(dist="none",low=0,up=Inf,unit="MT"),
      agemat=a50, 
      k     =list(low=500,up=2500,tol=0.01,permax=100), 
      b1k   =list(dist="none", low=0.8,  up=1.0,  mean=b1k,   sd=0.01),
      btk   =list(dist="beta", low=0.1,  up=0.5,  mean=btk,   sd=0.11, refyr=75),
      fmsym =list(dist="lnorm",low=0.5,  up=1.5,  mean=fmsym, sd=0.12),
      bmsyk =list(dist="beta", low=0.10, up=0.5,  mean=bmsyk, sd=0.05),
      M     =list(dist="lnorm",low=0.001,up=1,    mean=M,     sd=0.04),
      nsims =10000)

res$Initial
res$Estimates

prior[c("msy","bmsy","fmsy","v"),1]
ssb(om)[,75,,,,1]



