#source('~/Desktop/flr/FLife/R/omOut.R') 

library(reshape)
library(FLCore)
library(ggplotFL)
library(FLife)
library(MLZ)
library(mydas)

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/turbot.RData")

source('~/Desktop/sea++/mydas/pkgs/mydas/R/mlz.R')
source('~/Desktop/sea++/mydas/pkgs/mydas/R/omSmry.R')

ts=omSmry(om,eq,lh)

object=subset(ts,year%in%50:60)
res=mdply(sort(as.numeric(unique(ts$iter))), function(i,object,params){
    res=mlzFn(subset(object,iter==i),params[,i])
    cbind(var=c("z1","z2","year","sigma"),res)
    },object=object,params=prior)


ggplot(subset(res,!(var%in%c("sigma","year"))))+
  geom_boxplot(aes(x=var,as.numeric(as.character(Estimate))),outlier.colour = NA)+
  xlab("")+ylab("Z")+
  scale_y_continuous(limits=c(0,1))

ggplot(subset(res,(var%in%c("year"))))+
  geom_histogram(aes(as.numeric(as.character(Estimate))))+
  xlab("Year")+ylab("Count")