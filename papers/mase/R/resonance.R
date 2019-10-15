library(FLCore)
library(FLBRP)
library(FLasher)
library(FLife)
library(mydas)
library(ggplotFL)
library(scales)
library(plyr)
library(dplyr)
library(grid)
library(reshape)
library(popbio)
library(magrittr)
library(broom)
library(GGally)

source('~/Desktop/flr/FLife/R/omOut.R')


smry<-function(par,f,srDev,burn=20,m=function(length,params)
  exp(0.55)*(length^-1.61)%*%(params["linf"]^1.44)%*%params["k"],
  fbar=srDev%=%1){
  ##set up equilibrium object
  
  if ("numeric"%in%is(f))
    f=FLPar(f=f,c(1,length(f)))
    
  ## need to add interactions for f and par
  if (dim(par)[2]>1&dim(f)[2]>1){
    npar=dim(par)[2]
    par=as(mdply(seq(dim(f)[2]), with, cbind(model.frame(par)))[,-c(1,dim(par)[1]+2)],"FLPar")
    f  =rep(c(f),each=npar)
    f  =FLPar(array(f,c(1,length(f))))
    }
  
  eql=lhEql(par,m=m)
  
  ## convert to FLStock with constant F
  eq=eql

  fbar(eq)=fbar
  mlt=FLPar(f=array(c(f)*c(eq@refpts["msy","harvest"]),c(1,length(c(f)*c(eq@refpts["msy","harvest"])))))
  fbar(eq)=fbar(eq)%*%mlt
  
  stk=fwd(as(eq,"FLStock"),fbar=fbar(eq)[,-1],
          sr=eq,residuals=srDev)
  
  ## Summary stats
  srr=model.frame(FLQuants(eql,"ssb"=ssb,"rec"=rec),drop=TRUE)
  srp=model.frame(FLQuants("rp"=setPlusGroup(stock.n(eq)[,1]*stock.wt(eq)[,1]*eq@mat[,1],12)),drop=T)
  ts =model.frame(FLQuants(stk,"ssb"=ssb,"biomass"=stock,"rec"=rec,"catch"=catch,
                               "dev"=function(x) propagate(srDev,dim(x)[6])),drop=T)
  
  ts=subset(ts,year>burn)
  ts$year=ts$year-burn

  ind=omSmry(stk,eql,par)

  key=cbind(model.frame(par),f=c(f))
  
  list(srr=srr,srp=srp,ts=ts,ind=ind,key=key)}

par=propagate(lhPar(FLPar(linf=100,s=.9)),16)
dat=expand.grid(bg=c(3,4),sr=c(5000,5),s=c(0.75,0.9),k=c(0.1653,0.1653*2))

par["bg"]=dat$bg
par["sr"]=dat$sr
par["s" ]=dat$s
par["k" ]=dat$k

f  =FLPar(f=array(c(0.1,1,3),c(1,3)))

set.seed(234)
srDev=rlnoise(1,FLQuant(0,dimnames=list(year=1:1021)),0.3,0)
gis.1=smry(par,f,srDev)
m2.1 =smry(par,f,srDev,m=0.2)

set.seed(234)
srDev=rlnoise(1,FLQuant(0,dimnames=list(year=1:1021)),0.3,0.6)
gis.2=smry(par,f,srDev)
m2.2 =smry(par,f,srDev,m=0.2)

set.seed(234)
srDev=rlnoise(1,FLQuant(0,dimnames=list(year=1:1021)),0.5)
gis.3=smry(par,f,srDev)
m2.3 =smry(par,f,srDev,m=0.2)

ts=rbind(
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.1$key,gis.1$ts,by="iter")),
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="0.2",     merge( m2.1$key, m2.1$ts,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="Gislasen",merge(gis.2$key,gis.2$ts,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="0.2",     merge( m2.2$key, m2.2$ts,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.3$key,gis.3$ts,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="0.2",     merge( m2.3$key, m2.3$ts,by="iter")))

srr=rbind(
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.1$key,gis.1$srr,by="iter")),
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="0.2",     merge( m2.1$key, m2.1$srr,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="Gislasen",merge(gis.2$key,gis.2$srr,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="0.2",     merge( m2.2$key, m2.2$srr,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.3$key,gis.3$srr,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="0.2",     merge( m2.3$key, m2.3$srr,by="iter")))

srp=rbind(
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.1$key,gis.1$srp,by="iter")),
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="0.2",     merge( m2.1$key, m2.1$srp,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="Gislasen",merge(gis.2$key,gis.2$srp,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="0.2",     merge( m2.2$key, m2.2$srp,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.3$key,gis.3$srp,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="0.2",     merge( m2.3$key, m2.3$srp,by="iter")))

ind=rbind(
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.1$key,gis.1$ind,by="iter")),
  cbind(CV=0.3,AR=0.0,deviates="0.3",M="0.2",     merge( m2.1$key, m2.1$ind,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="Gislasen",merge(gis.2$key,gis.2$ind,by="iter")),
  cbind(CV=0.3,AR=0.6,deviates="0.3",M="0.2",     merge( m2.2$key, m2.2$ind,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="Gislasen",merge(gis.3$key,gis.3$ind,by="iter")),
  cbind(CV=0.5,AR=0.0,deviates="0.3",M="0.2",     merge( m2.3$key, m2.3$ind,by="iter")))

gis.prior=priors(par)
m2.prior =priors(par,eq=lhEql(lhPar(par),m=0.2))

prior=rbind(
  cbind("M"="0.2",      model.frame(gis.prior)),
  cbind("M"="Gislasin", model.frame( m2.prior)))

smry=list(ts=ts,srr=srr,srp=srp,ind=ind,key=gis.1$key,params=par,priors=prior)

save(smry,file="/home/laurence/Desktop/sea++/mydas/project/papers/mase/data/smry.RData",compress="xz")


#### Production
library(GGally)

expt=c("lopt","lc","spr0","mk","fmsy","msy","bmsy","r","rc","fm")

points.matrix=scale(smry$priors[,expt])

my_smooth <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_point(...,size=.5)+
    geom_smooth(...,method="lm",se=FALSE)}

my_density <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_density(...,lwd=1)}

ggpairs(as.data.frame(points.matrix)[,5:10],
        lower = list(continuous = wrap(my_smooth)),
        diag=list(continuous=wrap(my_density,alpha=0.2)),
        title = "")+
  theme(legend.position ="none",
        panel.grid.major =element_blank(),
        axis.ticks       =element_blank(),
        axis.text.x      =element_blank(),
        axis.text.y      =element_blank(),
        panel.border     =element_rect(linetype = 1, colour="black", fill=NA))

ggpairs(as.data.frame(points.matrix)[, c(1,3,4,5,6,7)],
        lower = list(continuous = wrap(my_smooth)),
        diag=list(continuous=wrap(my_density,alpha=0.2)),
        title = "")+
  theme(legend.position ="none",
        panel.grid.major =element_blank(),
        axis.ticks       =element_blank(),
        axis.text.x      =element_blank(),
        axis.text.y      =element_blank(),
        panel.border     =element_rect(linetype = 1, colour="black", fill=NA))

pellat <- function(biomass, params){
  a=sweep(biomass,6,params["r"]/params["p"],"*")
  b=sweep(biomass,6,params["k"],"/")
  c=sweep(b,      6,params["p"],"^")
  a*(1-c)}


kclusts <- data.frame(k=2:20) %>% 
  group_by(k) %>% 
  do(kclust=kmeans(points.matrix, .$k, nstart=25))

clusters <- kclusts %>% group_by(k) %>% do(tidy(.$kclust[[1]]))

assignments <- kclusts %>% group_by(k) %>% do(augment(.$kclust[[1]], points.matrix))
clusterings <- kclusts %>% group_by(k) %>% do(glance(.$kclust[[1]]))
ggplot(clusterings, aes(k, tot.withinss)) + 
  geom_line()+
  xlab("Number of Clusters")+ylab("Within SS")+
  theme_bw()

names(assignments)=tolower(names(assignments))
names(clusters)   =tolower(names(clusters))

ggplot(subset(assignments, k%in% 5:16), aes(x1, x2))+ 
  geom_point(aes(color=.cluster),size=0.2)+ 
  geom_point(data=subset(clusters, k%in% 5:16), size=2, shape=21,size=2,col="grey5")+
  stat_ellipse(aes(x1, x2,group=.cluster,col=.cluster))+    
  facet_wrap(~ k)+ 
  theme_bw()+
  theme(legend.position="none")

dt2=ddply(assignments, .(k), with,
          data.frame(ncluster=k,cluster=rep(.cluster,each=11),dat))

ggplot(subset(dt2,ncluster%in%seq(5,20,5)))+
  geom_point(aes(lag,acf))+
  geom_smooth(aes(lag,acf))+
  facet_grid(ncluster~cluster)

## Fournier ####################################################
load("/home/laurence/Desktop/sea++/mydas/project/papers/mase/data/smry.RData")

dat=ddply(smry$ts, .(CV,AR,M,linf,k,f,s,sr,bg), with, 
                as.data.frame(spec.pgram(ssb/mean(ssb),taper=0,log="no",plot=FALSE)[c("freq","spec")]))
dat=transform(dat,wavelen=1/freq)

ggplot(dat)+
  geom_point(aes(wavelen,spec,col=as.character(k)))+
  facet_grid(f~CV+AR,scale="free")

######################
ggplot(smry$srr)+
  geom_line(aes(ssb,rec))+
  geom_line(aes(x,rec,col=iter),
            data=transform(smry$ts,rec=rec/dev,
                                  x  =rescale(year,to=c(-200,0.8*max(ssb)))))+
  geom_path(aes(ssb,y,col=iter),
            data=transform(smry$ts,
                           y=(diags:::minMax(year)*0.95*mean(rec))))+
  scale_x_continuous(limits=c(-2e2,1e3))

ggplot(ddply(smry$srp,.(iter), transform, rp=rp/max(rp)))+
  geom_line(aes(age,rp,col=iter))

ggplot(ddply(smry$ts,.(iter), transform, ssb=ssb-mean(ssb)))+
  geom_line(aes(year,ssb,col=iter))

ggplot(ddply(smry$ts,.(iter), transform, rec=rec-mean(rec)))+
  geom_line(aes(year,rec,col=iter))

key=data.frame(m=rep(rep(c("Gislasen","0.2"),each=3),4),F=rep(c(0.1,1,3),8),model.frame(par[c("bg","sr")]))[,-5]
dat$iter=as.numeric(dat$iter)
dat=ddply(smry$ts, .(iter), with, acfFn(ssb,10))
dat=cbind(dat,key[dat$iter,])
       
ggplot(dat)+
    geom_line(aes(lag,acf,col=m,linetype=as.character(bg)))+
    xlab("Year")+ylab("")+
    facet_grid(F~sr)
  
ggplot(ddply(smry$ts, .(iter), with, acfFn(ssb,10)))+
  geom_errorbar(aes(lag,ymax=acf,ymin=0,col=iter),width=0)+
  xlab("Year")+ylab("")+
  facet_grid(iter~.)

thm=theme(legend.position = "none", 
      axis.title.y = element_text(colour='NA'), 
      axis.text.y  = element_text(colour="NA", angle=90), 
      axis.ticks.y = element_line(colour="NA"),
      axis.ticks =   element_line(colour="NA"),
      
      axis.title.x = element_blank(), 
      axis.text.x  = element_blank(), 
      axis.ticks.x = element_blank(), 
      
      plot.margin = unit(c(0, 0, 0, 1), "lines"),
      panel.background = element_rect(fill   ="NA", colour="NA"), 
      panel.border     = element_rect(fill   ="NA", colour="NA"), 
      panel.grid.major = element_line(colour ="NA"), 
      panel.grid.minor = element_line(colour ="NA")              )


grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(12, 12)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(sr+thm,     vp=define_region( 1:6,  2:10))

print(relEgg+thm, vp=define_region( 8:10, 5:9))

print(rP+thm,     vp=define_region( 7:9,  1:4))
print(sP+thm,     vp=define_region( 7:9,  9:12))

print(acfR+thm,     vp=define_region( 10:12,  1:4))
print(acfS+thm,     vp=define_region( 10:12,  9:12))

(max(rec(stk)%/%srDev)-min(rec(stk)%/%srDev))/min(rec(stk)%/%srDev)
(max(rec(stk))-min(rec(stk)))/min(rec(stk))

var(rec(stk))^0.5/mean(rec(stk))
var(ssb(stk))^0.5/mean(ssb(stk))


