load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/turbot.RData")
om1=iter(om,1:100)
eq1=iter(eq,1:100)

set.seed(1233)
srDev=FLife:::rlnoise(100,FLQuant(0,dimnames=list(year=1:105)),.3,b=0.0)

omf=fwd(om1,fbar =fbar( om1)[,-1],sr=eq1,residuals=srDev)
omc=fwd(om1,catch=catch(om1)[,-1],sr=eq1,residuals=srDev)

omm=om1
m(omm)=m(omm)%*%rlnoise(100,iter(m(omm)[1,],1)%=%0,sd=0.05,b=0.6)
omm=fwd(omm,catch=catch( om1)[,-1],sr=eq1,residuals=srDev)


omm2=om1
m(omm2)=m(omm2)%*%rlnoise(100,iter(m(omm2),1)%=%0,sd=0.05,b=0.0)
omm2=fwd(omm2,catch=catch( om1)[,-1],sr=eq1,residuals=srDev)

plot(FLStocks("Catch"=omc,"F"=omf,"1"=om1,"m"=omm,"m2"=omm2))

xsa<-function(om,pg=10,ctrl=xsaControl){
  stk=setPlusGroup(om,pg)
  idx=FLIndex(index=stock.n(stk))
  range(idx)[c("plusgroup","startf","endf")]=c(pg,0.1,.2)
  stk+FLXSA(stk,idx,control=ctrl,diag.flag=FALSE)}

xsaControl=FLXSA.control(tol    =1e-09, maxit   =150, 
                         min.nse=0.3,   fse     =1.0, 
                         rage   =1,     qage    =6, 
                         shk.n  =TRUE,  shk.f   =TRUE, 
                         shk.yrs=1,     shk.ages=4, 
                         window =10,    tsrange =10, 
                         tspower= 0,
                         vpa    =FALSE)
omm.=omm
m(omm.)=m(om1)

mp=xsa(window(omm.,end=75),ctrl=xsaControl,pg=10)

plot(FLStocks(list("xsa"=mp,"om"=omm)))

