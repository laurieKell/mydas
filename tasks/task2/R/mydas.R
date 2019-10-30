#R version 3.5.1, 64 bit

#library(devtools)
#install.packages(c("FLCore","FLFishery","FLasher","FLBRP","mpb","FLife"), repos="http://flr-project.org/R")
#install_github(c("flr/FLCore", "flr/FLFishery", "flr/kobe", "flr/diags"))
#install_github("lauriekell/mpb")
#install_github("lauriekell/mydas-pkg")

library(plyr)
library(FLife)
library(FLBRP)

load(url("https://github.com//fishnets//fishnets//blob//master//data//fishbase-web//fishbase-web.RData?raw=True"))

sp <- sort(unique(as.character(fb$species)))
sp[substring(sp,1,3)=='Lop']

#lh = subset(fb,species == "Psetta maxima")
lh = subset(fb,species == "Lophius budegassa")

names(lh)[c(14:17)] = c("l50","l50min","l50max","a50")
lh = lh[,c("species","linf","k","t0","a","b","a50","l50","l50min","l50max")]
lh = apply(lh[,-1],2, mean, na.rm = T)

lh = FLPar(lh)
par = lhPar(lh) 
# not all parameters are explained in the help file for lhPar: asym, bg, m1, m2, a1

eq = lhEql(par) #sr=bevholt by default, help file does not list other options, need to reference FLCore SRModels? 
plot(eq)

eq = lhEql(par, sr='segreg') #generates lots of missing values
plot(eq)

eq = lhEql(par, sr='ricker') #generates some missing values
plot(eq)

library(FLasher)
gTime = round(FLife:::genTime(FLPar(par)))

eq@fbar = refpts(eq)["msy","harvest"]%*%FLQuant(c(rep(.1,19),
                                              seq(.1,2,length.out = 30),
                                              seq(2,1.0,length.out = gTime)[-1],
                                              rep(1.0,61)))[,1:105]

plot(eq@fbar)
om = as(eq,"FLStock")

om = fwd(om,fbar = fbar(om)[,-1], sr = eq)
plot(om)


library(FLXSA)
library(mydas)
xsa = function(om,pg = 10,ctrl = xsaControl){
  stk = setPlusGroup(om,pg)
  idx = FLIndex(index=stock.n(stk))
  range(idx)[c("plusgroup","startf","endf")] = c(pg,0.1,.2)
  stk + FLXSA(stk, idx, control = ctrl, diag.flag = FALSE)
}


# this only works if you have data for a50 in lh
range(om)[c("minfbar","maxfbar")] = ceiling(mean(lh["a50"]))
range(eq)[c("minfbar","maxfbar")] = ceiling(mean(lh["a50"]))

range(om)[c("minfbar","maxfbar")] = ceiling(par["a50"])
range(eq)[c("minfbar","maxfbar")] = ceiling(par["a50"])



om = window(om, start = 25)
om = iter(om,1:1)
eq = iter(eq,1:1)


xsaControl = FLXSA.control(tol    =1e-09, maxit   =150, 
                         min.nse=0.3,   fse     =1.0, 
                         rage   =1,     qage    =3, 
                         shk.n  =TRUE,  shk.f   =TRUE, 
                         shk.yrs=1,     shk.ages=4, 
                         window =10,    tsrange =10, 
                         tspower= 0,
                         vpa    =FALSE)


mp = xsa(window(om,start = 25, end = 75), ctrl = xsaControl,pg = 10)

plot(FLStocks(list("xsa"= mp,"om" = om)))

nits = dim(mp)[6]
set.seed(4321)
srDev = FLife:::rlnoise(nits,rec(    om)[,,,,,1] %=% 0, 0.3, b = 0.0)
uDev = FLife:::rlnoise(nits,stock.n(om)[,,,,,1] %=% 0, 0.2, b = 0.0)

# function not in namespace, no help file yet
mseXSA=mydas:::mseXSA

mse1=mseXSA(om,
            eq,
            mp,control = xsaControl,
            ftar = 1.0,
            interval = 1,start = 54,end = 80,
            srDev = srDev, uDev = uDev)

plot(mse1)

###########
### Stochasticity
library(popbio)
nits = 100
set.seed(1234)
srDev = FLife:::rlnoise(nits,FLQuant(0, dimnames = list(year = 1:100)), 0.2, b = 0.0)
plot(iter(srDev,1))

#srDev = FLife:::rlnoise(1,FLQuant(0, dimnames = list(year = 1:100)), 1.0, b = 0.0)
#plot(iter(srDev,1))

### OEM
uDev = FLife:::rlnoise(nits, FLQuant(0, dimnames = list(year = 1:100)), 0.3, b  = 0.0)

## MSE for Derivative empirical MP
scen = expand.grid(stock = c("turbot")[1],
                 k1 = seq(0.2,0.8,0.2),k2 = seq(0.2,0.8,0.2), gamma = seq(0.75,1.25,0.25),
                 stringsAsFactors = FALSE)
library(reshape)
empD=NULL
for (i in seq(dim(scen)[1])){
  
  om=iter(om,seq(nits))
  eq=iter(eq,seq(nits))
  lh=iter(lh,seq(nits))
  
res = mydas:::mseSBTD(om,eq, control = with(scen[i,],
     FLPar(c(k1=k1, k2=k2, gamma=gamma))),
     start = 60, end = 100, srDev = srDev, uDev = uDev)

empD = rbind(empD, 
 cbind(scen=i,stock=scen[i,"stock"],
 k1 = scen[i,"k1"], k2 = scen[i,"k2"], 
 gamma = scen[i,"gamma"], 
 mydas:::omSmry(res, eq, lh)))
}
