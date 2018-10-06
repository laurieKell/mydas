library(plyr)
library(dplyr)
library(ggplot2)
library(reshape)

library(FLCore)
library(FLBRP)
library(FLasher)
library(ggplotFL)

# data(ple4)
# harvest( ple4)=harvest( ple4)%=%apply(harvest( ple4),1,mean)[drop=T]
# stock.wt(ple4)=stock.wt(ple4)%=%apply(stock.wt(ple4),1,mean)[drop=T]
# catch.wt(ple4)=catch.wt(ple4)%=%apply(catch.wt(ple4),1,mean)[drop=T]
# catch.n(ple4) =catch.n(ple4)%=%apply(catch.n(ple4),1,mean)[drop=T]
# landings.n(ple4)=catch.n(ple4)
# discards.n(ple4)[]=0
# 
# om=ple4

#harvest.spwn(om)=0

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/turbot.RData")

om=iter(om,2)
range(om)[c("minfbar","maxfbar")]=10

sr=as.FLSR(om,model="bevholt")

sr=fmle(sr,
         control=list(trace=TRUE),
         method="L-BFGS-B")

eq=FLBRP(om,sr=sr)
fbar(eq)[]=computeRefpts(eq)["msy","harvest"]

fbar=FLQuant(rep(c(computeRefpts(eq)["msy","harvest"]),each=dim(m(om))[2]-1),
             dimnames=dimnames(fbar(om)[,-1]))
om  =fwd(om,fbar=fbar,sr=eq)

plot(window(ssb(om),start=45)%/%computeRefpts(eq)["msy","ssb"])+
  geom_hline(aes(yintercept=1),col="red")

computeRefpts(eq)
fbar(eq)

