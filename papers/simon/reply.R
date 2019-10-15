library(FLife)
ls=expand.grid(linf=seq(25,100,length.out=10),s=seq(.55,1,length.out=10))
par=lhPar(FLPar(linf=ls$linf))
par["s"   ][]=ls[,"s"]
eql=lhEql(par)

ls$k   =c(par["k"])
ls$msy =c(refpts(eql)["msy","yield"])
ls$fmsy=c(refpts(eql)["msy","harvest"])
ls$bmsy=c(refpts(eql)["msy","ssb"])
names(ls)[2]="h"

ggplot(melt(ls,id=c("linf","Steepness","k")))+
  geom_line(aes(k,value,col=Steepness))+
  facet_grid(variable~.,scale="free")

ggplot(melt(ls,id=c("linf","Steepness","k")))+
  geom_line(aes(as.numeric(Steepness),value,col=ac(k)))+
  facet_grid(variable~.,scale="free")+xlab("Steepness")

fbar(eql)=fbar(eql)[,1]
fbar(eql)[]=refpts(eql)["msy","harvest"]*2

ls$cmsy =c(catch(eql)%/%refpts(eql)["msy","yield"])
ls$bbmsy=c(ssb(eql)%/%refpts(eql)["msy","ssb"])

ggplot(melt(ls,id=c("linf","s","k")))+
  geom_line(aes(k,value,col=ac(s)))+
  facet_wrap(~variable,scale="free",ncol=2)
