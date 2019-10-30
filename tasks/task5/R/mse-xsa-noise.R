library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)
library(FLAssess)
library(FLXSA)
library(mpb)
library(mydas)

#source('~/Desktop/flr/mpb/R/hcr.R')
#source('~/Desktop/flr/mpb/R/mseXSA.R')
#source('~/Desktop/flr/FLife/R/omOut.R')
#source('~/Desktop/flr/mse/R/msy.R')

theme_set(theme_bw())

dirMy ="/home/laurence/Desktop/sea++/mydas/tasks/task4"
dirDat=file.path(dirMy,"data")

xsa<-function(om,pg=10,ctrl=xsaControl){
  stk=setPlusGroup(om,pg)
  idx=FLIndex(index=stock.n(stk))
  range(idx)[c("plusgroup","startf","endf")]=c(pg,0.1,.2)
  stk+FLXSA(stk,idx,control=ctrl,diag.flag=FALSE)}

#load(file.path(dirDat,"brill.RData"))
#load(file.path(dirDat,"ray.RData"))
load(file.path(dirDat,"sprat.RData"))
range(om)[c("minfbar","maxfbar")]=ceiling(mean(lh["a1"]))
range(eq)[c("minfbar","maxfbar")]=ceiling(mean(lh["a1"]))

##OM
om=window(om,start=25)
om=iter(om,1:10)
eq=iter(eq,1:10)

##MP
xsaControl=FLXSA.control(tol    =1e-09, maxit   =150, 
                         min.nse=0.3,   fse     =1.0, 
                         rage   =1,     qage    =6, 
                         shk.n  =TRUE,  shk.f   =TRUE, 
                         shk.yrs=1,     shk.ages=4, 
                         window =10,    tsrange =10, 
                         tspower= 0,
                         vpa    =FALSE)
mp=xsa(window(om,end=60),ctrl=xsaControl,pg=10)
#mp=xsa(window(trim(om,age=3:20),end=75),ctrl=xsaControl,pg=20)
mp=xsa(window(om,end=60),ctrl=xsaControl,pg=5)

plot(FLStocks(list("xsa"=mp,"om"=om)))

nits=dim(mp)[6]
set.seed(4321)
srDev=FLife:::rlnoise(nits,rec(    om)[,,,,,1]%=%0,0.3,b=0.0)
uDev =FLife:::rlnoise(nits,stock.n(om)[,,,,,1]%=%0,0.2,b=0.0)

mseRay=mseXSA(om,
            eq,
            mp,control=xsaControl,
            ftar=1.0,
            interval=1,start=60,end=90,
            srDev=srDev,uDev=uDev)

##OM
omYr=om
set.seed(1234)
m(omYr)=m(om)%*%rlnoise(nits,iter(m(om),1)%=%0,sd=0.0,b=0.0,what="year")
set.seed(4321)
srDev=FLife:::rlnoise(nits,rec(     om)[,,,,,1]%=%0,0.3,b=0.0)

omYr=fwd(omYr,fbar=fbar(om)[,-1],sr=eq,residuals=srDev)


omYc=om
set.seed(1234)
m(omYc)=m(om)%*%rlnoise(nits,iter(m(om)[1,],1)%=%0,sd=0.05,b=0.6,what="year")

omYc=fwd(omYc,fbar=fbar(omYc)[,-1],sr=eq,residuals=srDev)

plot(FLStocks("Year"=omYr[,,,,,c((stock(omYr)[,ac(55)]-stock(om)[,ac(55)])/stock(om)[,ac(55)])>-0.5],
              "AR"  =omYc[,,,,,c((stock(omYc)[,ac(55)]-stock(om)[,ac(55)])/stock(om)[,ac(55)])>-0.5],
              "Recruitment"=om))

mse2=mseXSA(window(omYr,start=25),
            eq,
            mp,control=xsaControl,
            ftar=1.0,
            interval=1,start=50,end=80,
            srDev=srDev,uDev=uDev)

mse3=mseXSA(window(omYc,start=25),
            eq,
            mp,control=xsaControl,
            ftar=1.0,
            interval=1,start=50,end=80,
            srDev=srDev,uDev=uDev)

plot(FLStocks(llply(FLStocks("1"=mse1, "2"=mse2, "3"=mse3),window,start=20,end=80)))+
  facet_grid(qname~stock,scale="free")

save(mse1,mse2,mse3,file="/home/laurence/Desktop/mse.RData")
