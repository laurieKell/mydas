library(FLCore)
library(FLBRP)
library(FLasher)
library(FLife)
library(mydas)
library(popbio)

library(plyr)
library(GGally)
library(spatstat)

source('~/Desktop/sea++/mydas/pkg/R/oemLn.R')

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


## create lh parameter objects
# create covariance
cor=cor(model.frame(lh)[,c("linf","k","l50","b")],use="pairwise.complete.obs")

simCovar<-function(x,niters=500,se.missing=0.3,cv=0.1,cor=NULL){

  mn=aaply(x,1,mean, na.rm=TRUE)
  sd=aaply(x,1,var,  na.rm=TRUE)^0.5
  n =aaply(x,1,function(x) sum(!is.na(x)))
  se=sd/n^0.5
  
  if (any(is.na(se))) se[is.na(se)]=se.missing
  
  if (is.null(cor)){
    cor=array(0,rep(length(mn),2))
    diag(cor)=1}
  
  y=rmvnorm(niters,mn,cor2cov(cor,mn*cv*cv))
  
  res=FLPar(array(unlist(c(t(y))),c(dim(x)[1],niters)))
  
  dimnames(res)$params=names(mn)
  
  res}

par =as(lh[,-c(1,7)],"FLPar")
rng=0.2

nits=100
simPar=dlply(data.frame(species=lh$species,iter=seq(dim(lh)[1])),.(species), with,{
  set.seed(1234)  
  
  x=par[,iter]
  x=FLPar(apply(x,1,mean,na.rm=T))
  x=propagate(x,nits)
  x["linf"]=runif(nits, x["linf"]*(1-rng),x["linf"]*(1+rng))
  x["k"]   =runif(nits, x["k"]*(1-rng),   x["k"]*(1+rng))
  x["l50"] =runif(nits, x["l50"]*(1-rng), x["l50"]*(1+rng))
  lhPar(x)})

simPar[[3]]["a50"][is.na(simPar[[3]]["a50"])]=4
simPar[[3]]["sel1"][is.na(simPar[[3]]["sel1"])]=4

lopts=llply(simPar,function(x) mean(lopt(x)))
alks =llply(simPar,alk,cv=0.1)

save(simPar,alks,lopts,file="/home/laurence/Desktop/sea++/mydas/project/tasks/task6/wklifeIX/simPar.RData")

f=FLQuant(c(rep(0.1,60),seq(0.1,2,length.out=40)[-40],
            seq(2,1,length.out=11),rep(1,20)))

## Run OMs
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
# bug
#pds =llply(simPar,popdyn)


res2=list(NULL)
vals=seq(60)[!laply(oms,is.null)]
21, 26, 50:51, 56

for (i in vals[vals>56]){ 
     #seq(60)[!laply(oms,is.null)]){

  iStk=as.numeric(design[i,"Stock"])
  
  ## Simulate length frequencies by year
  lfd=lenSample(catch.n(oms[[i]]),alks[[iStk]],nsample=250)
  lfd=subset(as.data.frame(lfd),data>0)

  lop=mean(lopts[[iStk]])
  inds=ddply(lfd, .(year,iter), with, lenInd(len,data,lopt=lop))
  
  rf=FLQuants(  
    lbar =as.FLQuant(transmute(inds,year=year,iter=iter,data=lbar))/mean(simPar[[iStk]]["l50"]),
    lfm  =as.FLQuant(transmute(inds,year=year,iter=iter,data=lfm)),
    pmega=as.FLQuant(transmute(inds,year=year,iter=iter,data=pmega)),
    l95  =as.FLQuant(transmute(inds,year=year,iter=iter,data=l95/mean(simPar[[iStk]]["linf"]))),
    l25  =as.FLQuant(transmute(inds,year=year,iter=iter,data=l25/mean(simPar[[iStk]]["linf"]))))

  res2[[i]]=model.frame(rf,drop=T)
  
  save(res2,file="/home/laurence/Desktop/sea++/mydas/project/tasks/task6/wklifeIX/res2.RData")
  }

rs =ldply(res2,rbind)
key=cbind(i=seq(60)[!laply(res2,is.null)],ldply(res2,dim)[1])
tmp=mdply(key,function(i,V1) melt(rep(i,V1)))[,3]
rs =cbind(rs,i=tmp)


roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), 
             FPR=cumsum(!labels)/sum(!labels),
             labels,
             reference=sort(scores))}

ind=subset(rs,year%in%70:105)
ind=transform(ind,overfished=as.numeric(!(year%in%90:105)))
ind=melt(ind[,-1],id=c("i","iter","overfished"))
ind=ddply(ind,.(i,variable),with, roc(overfished,value))

ind=cbind(ind,design[ind$i,])
ggplot(ind)+ 
  geom_line(aes(FPR,TPR,group=i), size=0.2)+ #,col=as.character(CV)),size=.2)+
  geom_line(aes(x,y),data.frame(x=seq(0,1,0.01),y=seq(0,1,0.01)))+
  xlab("False Negative Rate")+ylab("True Positive Rate")+
  #scale_color_manual("Indicator",values=rainbow(5))+
  theme(legend.position="bottom")+
  facet_wrap(~variable)


