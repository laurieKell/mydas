cod=transform(cod,len=(data/a)^(1/b))
codPar=FLPar(c(optim(par=c(t0=0,k=.2,linf=10), fn,dat=cod)$par,a=0.018,b=3))
#bluefin
bftPar=FLPar(linf=318.85*.96,k=0.093,t0=-0.97-.16667,a=1.96e-5,b=3.0092,ar=4.0)
pars=FLPar(c(cbind(herPar,codPar,bftPar[dimnames(herPar)$params])),
dimnames=list(params=dimnames(herPar)$params,stock=c("Herring","Cod","Bluefin"),iter=1),units="NA")
pars=cbind(gislasim(FLPar(pars[,"Herring",drop=T])),
gislasim(FLPar(pars[,"Cod",drop=T])),
gislasim(FLPar(pars[,"Bluefin",drop=T])))
pars=FLPar(c(pars),
dimnames=list(params=dimnames(pars)$params,stock=c("Herring","Cod","Bluefin"),iter=1),units="NA")
parLh=pars
save(parLh,file=paste(dirMy,"data/parLh.RData",sep="/"))
load("/home/laurie/Ubuntu One/papers/submitted/3stocks/data/cod4.RData")
load("/home/laurie/Ubuntu One/papers/submitted/3stocks/data/her4.RData")
load("/home/laurie/Desktop/scrs2014/scrs2014-20/data/OM.RData")
names(cod4)[2]="WG"
names(her4)[2]="WG"
her=her4
cod=cod4
rm(her4,cod4)
bft=OM
rm(OM)
bft[[1]]=iter(bft[[1]],1)[,ac(1971:1990)]
bft[[2]]=iter(bft[[2]],1)[,ac(1971:1990)]
bft[[1]]=bft[[2]]
names(bft)=names(her)
mat(bft[["biol"]])=mat(bft[["biol"]])*fec(lFn(stock.wt(bft[["biol"]]),pars[,"Bluefin"]))/stock.wt(bft[["biol"]])
save(cod,her,bft,file=paste(dirMy,"data/stks.RData",sep="/"))
atage=rbind(cbind(stock="Bluefin",Assumptions="WG",  as.data.frame(bft[["WG"]][[  "mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Bluefin",Assumptions="Life History",as.data.frame(bft[["biol"]][["mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Cod",    Assumptions="WG",  as.data.frame(cod[["WG"]][[  "mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Cod",    Assumptions="Life History",as.data.frame(cod[["biol"]][["mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Herring",Assumptions="WG",  as.data.frame(her[["WG"]][[  "mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Herring",Assumptions="Life History",as.data.frame(her[["biol"]][["mat","m","stock.wt","ep"]],drop=T)))
atage=transform(atage,qname=factor(qname,labels=c("Mass","Maturity","Egg Production","Natural Mortality"),
#levels=c("Mass","Maturity","Egg Production","Natural Mortality")
levels=c("stock.wt","mat","ep","m")
))
atage[with(atage,stock=="Bluefin" & qname=="Maturity" & Assumptions=="Life History"),"data"]=
atage[with(atage,stock=="Bluefin" & qname=="Maturity" & Assumptions=="WG"),"data"]
save(atage,file=paste(dirMy,"data/atage.RData",sep="/"))
ggplot(ddply(atage,.(stock,Assumptions,qname),transform,data=data/max(data,na.rm=T)))+
geom_point( aes(age,data,col=Assumptions),size=.75)+
geom_smooth(aes(age,data,group=Assumptions,col=Assumptions))+
facet_grid(qname~stock,scale="free")+
xlab("Age")+ylab("")+theme_bw()
logisticFn <- function(params,data) { #x, a50, ato95){
pow <-function(a,b) a^b
func <- function(x,a50,ato95,asym){
if ((a50-x)/ato95 > 5)
return(0)
if ((a50-x)/ato95 < -5)
return(asym)
return(asym/(1.0+pow(19.0,(a50-x)/ato95)))}
sapply(data,func,params["a50"],params["ato95"],params["asym"])}
fn <- function(x,dat) {
res=sum((dat$data-logisticFn(FLPar(a50=x[1],ato95=x[2],asym=1),dat$age))^2)
res}
for (spp in c("Herring","Cod","Bluefin")){
dat=subset(atage, qname=="Maturity" & stock==spp & Assumptions=="WG")
pars[c("a50","ato95"), spp,1][]=optim(par=c(a50  =pars["a50",  spp],
ato95=pars["ato95",spp]), fn,  dat=dat)$par
}
source('~/Desktop/flr/git/FLBRP/R/lh.R', echo=TRUE)
library(FLBRP)
dirMy="/home/laurie/Ubuntu One/papers/Journal/elasticity"
dir3S ="/home/laurie/Ubuntu One/papers/submitted/3stocks"
#len cm
#Farley JH, Davis TL (1998) Reproductive dynamics of southern bluefin tuna, Thunnus maccoyii. Fishery Bulletin 96: 223–236.
fec=function(len) (4.78242*10e7)*len^7.530
lFn=function(wt,par) (wt/par["a"])^(1/par["b"])
#Herring
load(paste(dir3S,"data/her4new.RData", sep="/"))
herLW=read.csv(paste(dir3S,"/inputs/herLW.csv",sep=""),sep="\t")
herL =melt(read.csv(paste(dir3S,"/inputs/herL.csv", sep="")),id="year")
names(herL)[2:3]=c("age","len")
her=transform(merge(herL,herLW,by="year"),wt=a*len^b,age=as.numeric(age)-1)
fn <- function(x,dat)
sum((dat$len-vonB(FLPar(t0=x[1],k=x[2],linf=x[3]),as.numeric(ac(dat$age))))^2)
herPar=FLPar(c(optim(par=c(t0=0,k=.2,linf=50), fn,dat=her)$par,a=mean(herLW$a)*.001,b=mean(herLW$b)))
#cod
load(paste(dir3S,"data/cod4.RData", sep="/"))
cod=cbind(as.data.frame(stock.wt(cod4[["ices"]]),drop=T),a=0.018,b=3)
cod=transform(cod,len=(data/a)^(1/b))
codPar=FLPar(c(optim(par=c(t0=0,k=.2,linf=10), fn,dat=cod)$par,a=0.018,b=3))
#bluefin
bftPar=FLPar(linf=318.85*.96,k=0.093,t0=-0.97-.16667,a=1.96e-5,b=3.0092,ar=4.0)
pars=FLPar(c(cbind(herPar,codPar,bftPar[dimnames(herPar)$params])),
dimnames=list(params=dimnames(herPar)$params,stock=c("Herring","Cod","Bluefin"),iter=1),units="NA")
pars=cbind(gislasim(FLPar(pars[,"Herring",drop=T])),
gislasim(FLPar(pars[,"Cod",drop=T])),
gislasim(FLPar(pars[,"Bluefin",drop=T])))
pars=FLPar(c(pars),
dimnames=list(params=dimnames(pars)$params,stock=c("Herring","Cod","Bluefin"),iter=1),units="NA")
parLh=pars
save(parLh,file=paste(dirMy,"data/parLh.RData",sep="/"))
pars
load("/home/laurie/Ubuntu One/papers/submitted/3stocks/data/cod4.RData")
load("/home/laurie/Ubuntu One/papers/submitted/3stocks/data/her4.RData")
load("/home/laurie/Desktop/scrs2014/scrs2014-20/data/OM.RData")
names(cod4)[2]="WG"
names(her4)[2]="WG"
her=her4
cod=cod4
rm(her4,cod4)
bft=OM
rm(OM)
bft[[1]]=iter(bft[[1]],1)[,ac(1971:1990)]
bft[[2]]=iter(bft[[2]],1)[,ac(1971:1990)]
bft[[1]]=bft[[2]]
names(bft)=names(her)
mat(bft[["biol"]])=mat(bft[["biol"]])*fec(lFn(stock.wt(bft[["biol"]]),pars[,"Bluefin"]))/stock.wt(bft[["biol"]])
save(cod,her,bft,file=paste(dirMy,"data/stks.RData",sep="/"))
atage=rbind(cbind(stock="Bluefin",Assumptions="WG",  as.data.frame(bft[["WG"]][[  "mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Bluefin",Assumptions="Life History",as.data.frame(bft[["biol"]][["mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Cod",    Assumptions="WG",  as.data.frame(cod[["WG"]][[  "mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Cod",    Assumptions="Life History",as.data.frame(cod[["biol"]][["mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Herring",Assumptions="WG",  as.data.frame(her[["WG"]][[  "mat","m","stock.wt","ep"]],drop=T)),
cbind(stock="Herring",Assumptions="Life History",as.data.frame(her[["biol"]][["mat","m","stock.wt","ep"]],drop=T)))
atage=transform(atage,qname=factor(qname,labels=c("Mass","Maturity","Egg Production","Natural Mortality"),
#levels=c("Mass","Maturity","Egg Production","Natural Mortality")
levels=c("stock.wt","mat","ep","m")
))
atage[with(atage,stock=="Bluefin" & qname=="Maturity" & Assumptions=="Life History"),"data"]=
atage[with(atage,stock=="Bluefin" & qname=="Maturity" & Assumptions=="WG"),"data"]
save(atage,file=paste(dirMy,"data/atage.RData",sep="/"))
ggplot(ddply(atage,.(stock,Assumptions,qname),transform,data=data/max(data,na.rm=T)))+
geom_point( aes(age,data,col=Assumptions),size=.75)+
geom_smooth(aes(age,data,group=Assumptions,col=Assumptions))+
facet_grid(qname~stock,scale="free")+
xlab("Age")+ylab("")+theme_bw()
logisticFn <- function(params,data) { #x, a50, ato95){
pow <-function(a,b) a^b
func <- function(x,a50,ato95,asym){
if ((a50-x)/ato95 > 5)
return(0)
if ((a50-x)/ato95 < -5)
return(asym)
return(asym/(1.0+pow(19.0,(a50-x)/ato95)))}
sapply(data,func,params["a50"],params["ato95"],params["asym"])}
fn <- function(x,dat) {
res=sum((dat$data-logisticFn(FLPar(a50=x[1],ato95=x[2],asym=1),dat$age))^2)
res}
for (spp in c("Herring","Cod","Bluefin")){
dat=subset(atage, qname=="Maturity" & stock==spp & Assumptions=="WG")
pars[c("a50","ato95"), spp,1][]=optim(par=c(a50  =pars["a50",  spp],
ato95=pars["ato95",spp]), fn,  dat=dat)$par
}
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
fn <- function(x,len,m)
sum((m-fnM(FLPar(M1=x[1],M2=x[2]),len))^2)
for (spp in c("Herring","Cod","Bluefin")){
len=wt2len(pars[,spp],subset(atage, qname=="Mass" & stock==spp & Assumptions=="WG" & age>0 )$data)
m  =subset(atage, qname=="Natural Mortality" & stock==spp & Assumptions=="wg"  & age>0 )$data
pars[c("M1","M2"), spp,1][]=optim(par=c(M1=pars["M1",spp],
M2=pars["M2",spp]), fn,  len=len, m=m)$par
}
parWg=pars
save(parWg,file=paste(dirMy,"data/parWg.RData",sep="/"))
dmns=dimnames(parLh)
dmns[["Assumptions"]]=c("WG","Life History")
dmns=dmns[c(1:2,4,3)]
par=FLPar(aperm(array(c(parWg,parLh),laply(dmns,length),dimnames=dmns),c(1,3,2,4)))
save(par,file=paste(dirMy,"data/par.RData",sep="/"))
par
hat=function(stock,par){
ages=0:10
ages=FLQuant(ages,dimnames=list(age=ages))
#growth
as.data.frame(FLQuants("Mass"             =len2wt(  par[,stock],vonB(par[,stock],ages)),
"Maturity"         =logistic(par[,stock],ages),
"Egg Production"   =logistic(par[,stock],ages)*len2wt(  par[,stock],vonB(par[,stock],ages)),
"Natural Mortality"=fnM(     par[,stock],vonB(par[,stock],ages))),drop=T)
}
wg=mdply(data.frame(stock=c("Herring","Cod","Bluefin"),stringsAsFactors=F),hat,par=parWg)
lh=mdply(data.frame(stock=c("Herring","Cod","Bluefin"),stringsAsFactors=F),hat,par=parLh)
ggplot(subset(atage,qname=="Mass"))+
geom_point( aes(age,data,col=Assumptions),size=1)+
geom_smooth(aes(age,data,group=Assumptions,col=Assumptions),span=.3)+
facet_wrap(~stock,scale="free",ncol=1)+
geom_line(aes(age,data),col="green", data=subset(lh,qname=="Mass"))+
geom_line(aes(age,data),col="yellow",data=subset(wg,qname=="Mass"))+
xlab("Age")+ylab("")+theme_bw()
ggplot(subset(atage,qname=="Natural Mortality"))+
geom_point( aes(age,data,col=Assumptions),size=1)+
geom_smooth(aes(age,data,group=Assumptions,col=Assumptions),span=.3)+
facet_wrap(~stock,scale="free",ncol=1)+
geom_line(aes(age,data),col="green", data=subset(lh,qname=="Natural Mortality"))+
geom_line(aes(age,data),col="yellow",data=subset(wg,qname=="Natural Mortality"))+
xlab("Age")+ylab("")+theme_bw()
ggplot(subset(atage,qname=="Maturity" & Assumptions=="WG"))+
geom_point( aes(age,data,col=Assumptions),size=1)+
geom_smooth(aes(age,data,group=Assumptions,col=Assumptions),span=.3)+
facet_wrap(~stock,scale="free",ncol=1)+
geom_line(aes(age,data),col="green", data=subset(lh,qname=="Maturity"))+
geom_line(aes(age,data),col="yellow",data=subset(wg,qname=="Maturity"))+
xlab("Age")+ylab("")+theme_bw()
ggplot(subset(atage,qname=="Egg Production" & Assumptions=="WG"))+
geom_point( aes(age,data,col=Assumptions),size=1)+
geom_smooth(aes(age,data,group=Assumptions,col=Assumptions),span=.3)+
facet_wrap(~stock,scale="free",ncol=1)+
geom_line(aes(age,data),col="green", data=subset(lh,qname=="Egg Production"))+
geom_line(aes(age,data),col="yellow",data=subset(wg,qname=="Egg Production"))+
xlab("Age")+ylab("")+theme_bw()
library(numDeriv)
func=function(x,dmns,pr,what) {
pr[dmns] =exp(x)
pr["t0"] =pr["t0"]
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(pr,fnM=fnM,sr="bevholt")
fbar(res)=seq(0,0.75,length.out=51)*refpts(res)["crash","harvest"]
#refpts(res)["crash.75","harvest"]=refpts(res)["crash","harvest"]*0.75
refpts(res)["crash","harvest"]=refpts(res)["crash","harvest"]*.75
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
log(res)}
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
dmns=dimnames(pr)$params
what="msy"
x=c(pr)
func(x,dmns,pr,what)
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
dmns=dimnames(pr)$params
what="msy"
x=c(pr)
pr[dmns] =exp(x)
pr["t0"] =pr["t0"]
x
pr[dmns] =exp(x)
pr["t0"] =pr["t0"]
pr
pr=pars[,1,1]
pr
c(pr)
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
pr
x=log(c(pr))
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
dmns=dimnames(pr)$params
what="msy"
x=log(c(pr))
func(x,dmns,pr,what)
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
dmns=dimnames(pr)$params
what="msy"
x=log(c(pr))
x
pr=pars[,1,1]
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
dmns=dimnames(pr)$params
what="msy"
x=log(c(pr))
pr[dmns] =exp(x)
pr["t0"] =pr["t0"]
p
pr
pr["t0"] =pr["t0"]
pr["t0"] =-pr["t0"]
pr
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(pr,fnM=fnM,sr="bevholt")
plot(res)
fbar(res)=seq(0,0.75,length.out=51)*refpts(res)["crash","harvest"]
#refpts(res)["crash.75","harvest"]=refpts(res)["crash","harvest"]*0.75
refpts(res)["crash","harvest"]=refpts(res)["crash","harvest"]*.75
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
log(res)
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
what="msy"
jbn=jacobian(func,log(c(pr)),dmns=dimnames(pr)$params,par=pr,what=what)
pr=pars[,1,1]
pr[1,1]=-pr[1,1]
what="msy"
log(c(pr))
dimnames(pr)$params
pr
what
func=function(x,dmns,pr,what) {
pr[dmns] =exp(x)
pr["t0"] =-pr["t0"]
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(pr,fnM=fnM,sr="bevholt")
fbar(res)=seq(0,0.75,length.out=51)*refpts(res)["crash","harvest"]
#refpts(res)["crash.75","harvest"]=refpts(res)["crash","harvest"]*0.75
refpts(res)["crash","harvest"]=refpts(res)["crash","harvest"]*.75
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
log(res)}
jbn=jacobian(func,log(c(pr)),dmns=dimnames(pr)$params,pr=pr,what=what)
jbn
par=pars[,1,1]
x  =log(abs(par))
nc=function(x,what,par,dmns=dimnames(par)$params,sr="bevholt") {
par.=par
par[dmns] = exp(x)
if (any(par.<0)) par[par.<0] = -par[par.<0]
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(par,fnM=fnM,sr=sr)
fbar(res)=seq(0,0.75,length.out=51)*refpts(res)["crash","harvest"]
#refpts(res)["crash.75","harvest"]=refpts(res)["crash","harvest"]*0.75
refpts(res)["crash","harvest"]=refpts(res)["crash","harvest"]*.75
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
log(res)}
func=function(x,what,par,dmns=dimnames(par)$params,sr="bevholt") {
par.=par
par[dmns] = exp(x)
if (any(par.<0)) par[par.<0] = -par[par.<0]
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(par,fnM=fnM,sr=sr)
fbar(res)=seq(0,0.75,length.out=51)*refpts(res)["crash","harvest"]
#refpts(res)["crash.75","harvest"]=refpts(res)["crash","harvest"]*0.75
refpts(res)["crash","harvest"]=refpts(res)["crash","harvest"]*.75
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
log(res)}
par=pars[,1,1]
what="msy"
jbn=jacobian(func,log(abs(par)),par=par,dmns=dimnames(par)$params,what=what,sr="bevholt")
jbn
doIt=function(x,y){
func=function(x,what,par,dmns=dimnames(par)$params,sr="bevholt",
fbarFn=function(x,nyr=51) seq(0,0.75,length.out=nyr)*refpts(x)["crash","harvest"]*.75) {
par.=par
par[dmns] = exp(x)
if (any(par.<0)) par[par.<0] = -par[par.<0]
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(par,fnM=fnM,sr=sr)
fbar(res)=fbarFn(res)
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
log(res)}
}
func=function(x,what,par,dmns=dimnames(par)$params,sr="bevholt",
fbarFn=function(x,nyr=51) seq(0,0.75,length.out=nyr)*refpts(x)["crash","harvest"]*.75) {
par.=par
par[dmns] = exp(x)
if (any(par.<0)) par[par.<0] = -par[par.<0]
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(par,fnM=fnM,sr=sr)
fbar(res)=fbarFn(res)
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
log(res)}
jbn=jacobian(func,log(abs(par)),par=par,dmns=dimnames(par)$params,what=what,sr="bevholt")
jbn
what=c("msy","f0,1")
jbn=jacobian(func,log(abs(par)),par=par,dmns=dimnames(par)$params,what=what,sr="bevholt")
stk
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(ssb(    stk)/refpts(res)[what,"ssb"]),
c(stock(  stk)/refpts(res)[what,"biomass"]),
c(catch(  stk)/refpts(res)[what,"harvest"]),
c(catch(  stk)/refpts(res)[what,"yield"]),
c(rec(    stk)/refpts(res)[what,"rec"]))
res=FLBRP(stk)
refpts(res)[what,"ssb"])
c(ssb(    stk)/refpts(res)[what,"ssb"])
refpts(res)[what,"ssb"]
refpts(res)[what[1],"ssb"]
what
what=c("msy","f0.1")
what
refpts(res)[what,"ssb"]
ssb(    stk)/refpts(res)[what,"ssb"]),
ssb(    stk)/refpts(res)[what,"ssb"]
ssb(    stk)%/%refpts(res)[what,"ssb"]
maply(what, function(x) ssb(stk)/refpts(res)[x,"ssb"])
mdply(what, function(x) ssb(stk)/refpts(res)[x,"ssb"])
maply(what, function(x) ssb(stk)/refpts(res)[x,"ssb"])
c(maply(what, function(x) ssb(stk)/refpts(res)[x,"ssb"]))
func=function(x,what,par,dmns=dimnames(par)$params,sr="bevholt",
fbarFn=function(x,nyr=51) seq(0,0.75,length.out=nyr)*refpts(x)["crash","harvest"]*.75) {
par.=par
par[dmns] = exp(x)
if (any(par.<0)) par[par.<0] = -par[par.<0]
fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))
res=lh(par,fnM=fnM,sr=sr)
fbar(res)=fbarFn(res)
refpts(res)["crash",2:8]=NA
dimnames(refpts(res))$refpt[5]=".crash"
res=brp(res)
dimnames(refpts(res))$refpt[5]="crash"
stk=fwd(brp(res),maxF=3)
stk=fwd(stk,f=fbar(res)[,-1],sr=res)
res=c(c(ssb(    stk)),
c(stock(  stk)),
c(catch(  stk)),
c(catch(  stk)),
c(rec(    stk)),
c(maply(what, function(x) ssb(stk)/refpts(res)[x,"ssb"])),
c(maply(what, function(x) ssb(stk)/refpts(res)[x,"biomass"])),
c(maply(what, function(x) ssb(stk)/refpts(res)[x,"harvest"])),
c(maply(what, function(x) ssb(stk)/refpts(res)[x,"yield"])),
c(maply(what, function(x) ssb(stk)/refpts(res)[x,"rec"])))
log(res)}
jbn=jacobian(func,log(abs(par)),par=par,dmns=dimnames(par)$params,what=what,sr="bevholt")
jbn
length(jbn)
func
savehistory("~/Desktop/U1/papers/Journal/elasticity/R/exStk_.R")
