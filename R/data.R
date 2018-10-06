library(FLCore)
library(plyr)

load(url("https://github.com//fishnets//fishnets//blob//master//data//fishbase-web//fishbase-web.RData?raw=True"))

fbSmry<-function(fb){
  
  names(fb)[c(14:17)]=c("l50","l50min","l50max","a50")
  fb=fb[,c("species","order","family","linf","k","t0","a","b","a50","l50","l50min","l50max")]
  
  fb[is.na(fb["l50"]),"l50"]=(fb[is.na(fb["l50"]),"l50min"]+fb[is.na(fb["l50"]),"l50max"])/2
  fb[fb$t0>=0,"t0"]=NA
  
  fb$species=as.character(fb$species)
  fb$species=factor(as.character(fb$species))
  fb$family=factor(as.character(fb$family))
  
  fb=subset(fb,!is.na(linf))
  
  smn=ddply(fb,.(species,order),with,data.frame(
    k   =mean(k,na.rm=T),
    linf=mean(linf,na.rm=T),
    t0  =mean(t0,na.rm=T),
    a   =mean(a,na.rm=T),
    b   =mean(b,na.rm=T),
    a50 =mean(a50,na.rm=T),
    l50 =mean(l50,na.rm=T)))
  
  sn =ddply(fb,.(species,order),with,data.frame(
    k   =sum(!is.na(k)),
    linf=sum(!is.na(linf)),
    t0  =sum(!is.na(t0)),
    a   =sum(!is.na(a)),
    b   =sum(!is.na(b)),
    a50 =sum(!is.na(a50)),
    l50 =sum(!is.na(l50))))
  
  svr=ddply(fb,.(species,order),with,data.frame(
    k   =var(k,na.rm=T),
    linf=var(linf,na.rm=T),
    t0  =var(t0,na.rm=T),
    a   =var(a,na.rm=T),
    b   =var(b,na.rm=T),
    a50 =var(a50,na.rm=T),
    l50 =var(l50,na.rm=T)))
  
  pmn=FLife:::mf2FLPar(smn[,-(1:2)])
  
  dmns=dimnames(pmn)
  dmns[["quant"]]=c("mean","se","var","n")
  dmns=dmns[c(1,3,2)]
  
  fb=FLPar(array(NA,dim=unlist(laply(dmns,length)),dimnames=dmns))
  dimnames(fb)$iter=svr[,"species"]
  
  fb[,"mean"]=FLife:::mf2FLPar(smn[,-(1:2)])
  fb[,"var"]=FLife:::mf2FLPar(svr[,-(1:2)])
  fb[,"n"]  =FLife:::mf2FLPar(sn[,-(1:2)])
  fb[,"se"] =(fb[,"var"]/(fb[,"n"]-1))^0.5
  
  fbvcov=var(model.frame(fb[,"mean"])[,1:7],na.rm=TRUE)
  fbcor =FLPar(cov2cor(fbvcov))
  
  FLPars(par=fb,cor=fbcor)}

saury   ="Cololabis saira"
mackerel=c("Scomber scombrus","Scomber colias","Scomber japonicus","Scomber australasicus","Scomber indicus")
herring =c("Clupea harengus harengus","Strangomera bentincki","Clupea harengus pallasii") 
anchovy =c("Engraulis encrasicolus","Engraulis anchoita","Engraulis mordax","Engraulis japonicus","Engraulis ringens","Engraulis capensis")
sardine =c("Sardina pilchardus","Sardinops sagax","Sardinops melanostictus","Sardinops caeruleus","Sardinops ocellatus","Sardinella lemuru","Sardinella brasiliensis","Sardinella zunasi","Sardinella longiceps","Sardinella gibbosa","Sardinella aurita","Sardinella maderensis","Dussumieria acuta")
tuna    =c("Thunnus alalunga","Thunnus maccoyii","Thunnus orientalis","Thunnus thynnus","Thunnus obesus","Thunnus albacares","Katsuwonus pelamis")

npfc=fbSmry(subset(fb,species%in%c(mackerel,saury,tuna,sardine,herring,anchovy)))

save(npfc,compress="xz",
     file="/home/laurence/Desktop/sea++/mydas/R/shiny/data/npfc.RData")

mydas=subset(fb,species%in%c("Scophthalmus rhombus","Raja clavata","Psetta maxima",
                             "Pollachius pollachius","Pollachius virens",
                             "Sprattus sprattus sprattus"))

mydas=fbSmry(mydas)

save(mydas,compress="xz",
     file="/home/laurence/Desktop/sea++/mydas/R/shiny/data/mydas.RData")

grps=ac(unique(read.csv("/home/laurence/Desktop/sea++/mydas/R/shiny/data/Fishdata.csv")[,5]))
grps=subset(fb,species%in%grps)
grps=fbSmry(grps)

save(grps,compress="xz",
     file="/home/laurence/Desktop/sea++/mydas/R/shiny/data/grps.RData")

load("/home/laurence/Desktop/sea++/mydas/papers/FLife/data/scombrid.RData")
scombrid=scombrid[,c("name","taxon","k","linf","t0","a","b","a50","amat","l50")]

names(scombrid)[1]="species"
scombrid=subset(scombrid,!is.na(linf))

smn=ddply(scombrid,.(species,taxon),with,data.frame(
  k   =mean(k,na.rm=T),
  linf=mean(linf,na.rm=T),
  t0  =mean(t0,na.rm=T),
  a   =mean(a,na.rm=T),
  b   =mean(b,na.rm=T),
  a50 =mean(a50,na.rm=T),
  amat=mean(amat,na.rm=T),
  l50 =mean(l50,na.rm=T)))

sn =ddply(scombrid,.(species,taxon),with,data.frame(
  k   =sum(!is.na(k)),
  linf=sum(!is.na(linf)),
  t0  =sum(!is.na(t0)),
  a   =sum(!is.na(a)),
  b   =sum(!is.na(b)),
  a50 =sum(!is.na(a50)),
  amat=sum(!is.na(amat)),
  l50 =sum(!is.na(l50))))

svr=ddply(scombrid,.(species,taxon),with,data.frame(
  k   =var(k,na.rm=T),
  linf=var(linf,na.rm=T),
  t0  =var(t0,na.rm=T),
  a   =var(a,na.rm=T),
  b   =var(b,na.rm=T),
  a50 =var(a50,na.rm=T),
  amat=var(amat,na.rm=T),
  l50 =var(l50,na.rm=T)))

pmn=FLife:::mf2FLPar(smn[,-(1:2)])

dmns=dimnames(pmn)
dmns[["quant"]]=c("mean","se","var","n")
dmns=dmns[c(1,3,2)]

scom=FLPar(array(NA,dim=unlist(laply(dmns,length)),dimnames=dmns))
dimnames(scom)$iter=svr[,"species"]

scom[,"mean"]=FLife:::mf2FLPar(smn[,-(1:2)])
scom[,"var"]=FLife:::mf2FLPar(svr[,-(1:2)])
scom[,"n"]  =FLife:::mf2FLPar(sn[,-(1:2)])
scom[,"se"] =(scom[,"var"]/(scom[,"n"]-1))^0.5

vr=var(model.frame(scom[,"mean"])[,1:8],na.rm=TRUE)
cr=FLPar(cov2cor(vr))

scom=FLPars(par=scom,cor=cr)
save(scom,compress="xz",
     file="/home/laurence/Desktop/sea++/mydas/R/shiny/data/scom.RData")

corrplot(grps[[2]][drop=T])

library(ggplot2)
library(GGally)

p <- ggpairs(log(model.frame(grps[[1]][,"mean"]))[,c(1:2,6:7)]))+ theme_bw()
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07","blue")) +
      scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07","blue"))  
  }
}


