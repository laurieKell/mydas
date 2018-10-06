library(FLBRP)
library(numDeriv)

dirMy="/home/laurie/Ubuntu One/papers/Journal/elasticity"

par =gislasim(FLPar(linf=100,t0=.1,M1=2.1104327,M2=1.7023068,s=0.9,v=1000,fec=1))

#M
fnM=function(par,len,T=290){
  a=FLPar(c(a=-par["M1"],b=-par["M2"],c=1.5067827,d=0.9664798,e=763.5074169))
  exp(a[1]+a[2]*log(len) + a[3]*log(par["linf"]) + a[4]*log(par["k"]) + a[5]/T)}

fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))

#Fbar
fnFbar=function(x,nyr=51) seq(0,0.75,length.out=nyr)*2*refpts(x)["msy","harvest"]


x=c(par)
func=function(x,rp,qnt,yr,nyr,
              par,dmns=dimnames(par)$params,
              fnFbar=function(x,nyr) seq(0,0.75,length.out=nyr)*2*refpts(x)["msy","harvest"],
              fnM   =function(par,len) exp(par["M1"]-par["M2"]*log(len)))
              {
  
    par[dmns]=x
    
    brp=lh(par,fnM=fnM)
    fbar(brp)=fnFbar(brp,nyr)[1:yr]
 
    stk=fwd(brp(brp),maxF=3)
    stk=fwd(stk,f=fbar(brp)[,-1],sr=brp)
    
    cat(1)
    
    stk[[qnt]][[1]][,yr]/refpts(brp)[rp,qnt] 
    }

func(c(par),"msy","ssb",21,51,par)

par =gislasim(FLPar(linf=100,t0=.1,M1=2.1104327,M2=1.7023068,s=0.9,v=1000,fec=1))
hss=hessian(func,c(par),rp="msy",qnt="ssb",nyr=51,yr=20,par=par)
  