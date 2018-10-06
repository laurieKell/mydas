library(FLBRP)

plusGrp  =40
m        =FLQuant(c(0.49,0.24,0.24,0.24,0.24,0.20,0.175,0.125,rep(0.1,32)),dimnames=list(age=1:plusGrp))
mat      =FLQuant(c(0,0,0,0.5,rep(1,36)),dimnames=list(age=1:plusGrp))

sex      =0.5
par      =gislasim(FLPar(c(linf=314.9,k=0.089,t0=-1.13,a=1.96e5,b=3.00,s=.99,v=10e6,spr0=470)))

range    =c(min=1,max=plusGrp,minfbar=1,maxfbar=plusGrp,plusgroup=plusGrp)
bft=lh(par,m=m,mat=mat,range=range)
fbar(bft)=fbar(bft)*.1
refpts(bft)=refpts(bft)[-3,]