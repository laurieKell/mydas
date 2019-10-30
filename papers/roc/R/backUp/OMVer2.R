library(FLCore)
library(FLBRP)
library(FLasher)
library(FLife)
library(mydas)
library(popbio)

library(plyr)
library(GGally)
library(spatstat)
library(mvtnorm)
library(reshape)

source('~/Desktop/sea++/mydas/pkg/R/oemLn.R')

simCovar<-function(x,cv,cor,nits=100){
  
  mn=aaply(x,1,mean,na.rm=TRUE)
  y =exp(rmvnorm(nits,       log(mn[dimnames(cor)[[1]]]),
                 cor2cov(cor,log(mn[dimnames(cor)[[1]]])*cv*cv)))
  
  res=FLPar(t(array(unlist(c(y)),c(nits,dim(cor)[1]),
                    dimnames=list(iter=seq(100),params=dimnames(cor)[[1]]))))

  mn=propagate(FLPar(mn),nits)
  mn[dimnames(cor)[[1]]]=res
  mn["l50"]=res["linf"]*res["l50linf"]
  lhPar(mn)}

roc<-function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), 
             FPR=cumsum(!labels)/sum(!labels),
             labels,
             reference=sort(scores))}

## load life history data from fishbase
#load("/home/laurence/Desktop/papers/generic/fishnets-master/data/fishbase-web/fishbase-web.RData")

load(url("https://github.com//fishnets//fishnets//blob//master//data//fishbase-web//fishbase-web.RData?raw=True"))

lh=subset(fb,species%in%c("Pollachius pollachius","Psetta maxima","Scophthalmus rhombus","Raja clavata",
                          "Sprattus sprattus","Sprattus sprattus sprattus")) 

names(lh)[c(14:17)]=c("l50","l50min","l50max","a50")
lh=lh[,c("species","linf","k","t0","a","b","a50","l50","l50min","l50max")]

lh[is.na(lh["l50"]),"l50"]=(lh[is.na(lh["l50"]),"l50min"]+lh[is.na(lh["l50"]),"l50max"])/2
lh=lh[,-(9:10)]

lh[lh$t0>=0,"t0"]=NA
lh=rbind(lh[1,],lh)
lh[1,-c(1,7:8)]=c(84.6,0.19,-0.94,0.01017,2.98)
lh[1,"l50"]=0.72*lh[1,"linf"]^0.93

fctr=ddply(lh,.(species),with,mean(linf))
fctr=as.character(fctr[order(fctr$V1),"species"])

lh$species=factor(as.character(lh$species),levels=fctr)

lh=transform(lh,l50linf=l50/linf)
cor=cor(model.frame(lh)[,c("linf","k","l50","t0","a","b","l50linf")],use="pairwise.complete.obs")

mp=ddply(lh,.(species), with, data.frame(linf=mean(linf,na.rm=T),
                                         k   =mean(k,   na.rm=T),
                                         t0  =mean(t0,  na.rm=T),
                                         l50 =mean(l50, na.rm=T),
                                         a   =mean(a,   na.rm=T),
                                         b   =mean(b,   na.rm=T)))
mp=as(mp[,-1],"FLPar")
#Linf*(3/(3+M/K)) 
mp=rbind(mp,FLPar("lopt" =unlist(2/3*mp["linf"])))
mp=rbind(mp,FLPar("lmega"=unlist(1.1*mp["lopt"])))
attributes(mp)$species=fctr

my_smooth <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_point(...,size=.5)+
    geom_smooth(...,method="lm",se=FALSE)}

my_density <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_histogram(...,lwd=1)}

theme_set(theme_bw(base_size=20))

ggpairs(transform(lh[,c(1:3,8,6)],linf=log(linf),k=log(k),l50=log(l50)),
          mapping = ggplot2::aes(color=species),
          lower = list(continuous = wrap(my_smooth)),
          diag=list(continuous=wrap(my_density,alpha=0.2)),
          title = "")+
  theme(legend.position ="none",
        panel.grid.major =element_blank(),
        axis.ticks       =element_blank(),
        axis.text.x      =element_blank(),
        axis.text.y      =element_blank(),
        panel.border     =element_rect(linetype = 1, colour="black", fill=NA))


lhs=dlply(lh,.(species), with, FLPar(linf=linf,k=k,t0=t0,l50=l50,a=a,b=b,l50linf=l50linf))
lhs=llply(lhs, simCovar,cv=0.1,cor=cor[c("linf","k","l50linf"),
                                       c("linf","k","l50linf")],nits=100)

save(lhs,file="/home/laurence/tmp/lhs.RData")

f=FLQuant(c(rep(0.1,60),seq(0.1,2.5,length.out=40)[-40],
            seq(2.5,1,length.out=11),rep(1,20)))

nits=100
set.seed(234)
srDev=rlnoise(nits,f%=%0,0.3,0)

oms=llply(lhs, function(x,f,srDev){
  eq      =lhEql(x)
  fbar(eq)=f%*%refpts(eq)["msy","harvest"]
  
  ##check fmsy is estimated
  if (any(is.na(refpts(eq)["msy","harvest"])))
    fbar(eq)[,,,,,is.na(refpts(eq)["msy","harvest"])]=
      f%*%refpts(eq)["f0.1","harvest",is.na(refpts(eq)["msy","harvest"])]
  
  om      =as(eq,"FLStock")
  om      =fwd(om,fbar=fbar(eq)[,-1],
               sr=eq,residuals=srDev)
  om},f=f,srDev=srDev)

## LFDs
alks =llply(lhs,alk,cv=0.1)

inds=mdply(seq(length(oms)), function(iSpp){
  
  print(iSpp)
  
  om=window(oms[[iSpp]],start=2)
  
  if (any(any(catch.n(om)<0)))
    catch.n(om)[catch.n(om)<0]=0.0
  
  #lfd=maply(data.frame(iter=seq(dim(om)[6])), function(iter) 
  #  lenSample(iter(catch.n(om),iter),alks[[iSpp]][,,iter],nsample=500))
  lfd=lenSample(catch.n(om),alks[[iSpp]],nsample=500)
  lfd=transform(melt(lfd),data=value)[,-4]
  lfd=subset(as.data.frame(lfd),data>0)
  lfd=cbind(lfd,lopt=c(mp["lopt",iSpp]))
  
  ddply(lfd, .(year,iter), with, {
                       lenInd(len,data,unique(lopt))})
  })

ind=subset(inds[,-6],year%in%70:100)
ind=transform(ind,overfished=!(as.numeric(year)%in%90:105))
ind=melt(ind[,-2],id=c("X1","iter","overfished"))
ind=ddply(ind,.(X1,variable),with, roc(overfished,value))

ggplot(ind)+ 
  geom_line(aes(FPR,TPR,group=variable,col=as.character(variable)))+
  geom_line(aes(x,y),data.frame(x=seq(0,1,0.01),y=seq(0,1,0.01)))+
  xlab("False Negative Rate")+ylab("True Positive Rate")+
  #scale_color_manual("Indicator",values=rainbow(5))+
  theme_bw(12)+
  theme(legend.position="bottom")+
  facet_wrap(~fctr[X1])

## OMs
design=data.frame(s   =c( 0.9, 0.7, 0.9, 0.9, 0.9, 0.9),
                  sel2=c( 1.0, 1.0,5000, 1.0, 1.0, 1.0),
                  sel3=c(5000,5000,5000,  50,5000,5000),
                  m1  =c(0.55,0.55,0.55,0.55,0.05,0.55),
                  m2  =c(-1.6,-1.6,-1.6,-1.6,-1.6,-1.0))
design=mdply(expand.grid(Stock=seq(5),CV=c("0.3","0.5")), function(Stock,CV) design)

## Stochasticity
set.seed(234)
srDev0.3=rlnoise(nits,f%=%0,0.3,0)
set.seed(234)
srDev0.5=rlnoise(nits,f%=%0,0.5,0)
srDev=FLQuants("0.5"=srDev0.5,"0.3"=srDev0.3)
rm(srDev0.5,srDev0.3)

## Runs
om2=FLStocks(list(NULL))

for (i in seq(60)){
  par        =simPar[[design[i,"Stock"]]]
  par["s"]   =design[i,"s"]
  par["sel2"]=design[i,"sel2"]
  par["sel3"]=design[i,"sel3"]
  par["m1"]  =design[i,"m1"]
  par["m2"]  =design[i,"m2"]
  
  eq      =lhEql(lhPar(par))
  fbar(eq)=f%*%refpts(eq)["msy","harvest"]
  om      =as(eq,"FLStock")
  om2[[i]]=fwd(om,fbar=fbar(eq)[,2:130],
                  sr  =eq,residuals=srDev[[design[i,"CV"]]])
  } 

## summary stuff


rs =ldply(res2,rbind)
key=cbind(i=seq(60)[!laply(res2,is.null)],ldply(res2,dim)[1])
tmp=mdply(key,function(i,V1) melt(rep(i,V1)))[,3]
rs =cbind(rs,i=tmp)




