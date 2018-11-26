library(quantreg)
library(ggplot2)

gis=read.csv("/home/laurence/Desktop/sea++/mydas/project/papers/FLife/inputs/gislason-et-al-2010.csv",sep="\t")

m.50=rq(log(m)~log(l)+log(k)+log(linf),data=gis,tau=0.5)
m.25=rq(log(m)~log(l)+log(k)+log(linf),data=gis,tau=0.25)
m.75=rq(log(m)~log(l)+log(k)+log(linf),data=gis,tau=0.75)

ggplot(rbind(data.frame(length=15:50,m=exp(predict(m.25,newdata=data.frame(k=0.3,linf=50,l=15:50))),quantile="25th"),
      data.frame(length=15:50,m=exp(predict(m.50,newdata=data.frame(k=0.3,linf=50,l=15:50))),quantile="50th"),
      data.frame(length=15:50,m=exp(predict(m.75,newdata=data.frame(k=0.3,linf=50,l=15:50))),quantile="75th")))+
  geom_line(aes(length,m,col=quantile))

