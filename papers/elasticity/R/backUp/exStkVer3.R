library(FLAdvice)
library(numDeriv)

dirMy="/home/lkell/Desktop/Stuff/Publications/gbyp-sam/papers/journal/elastic"

doIt=function(what,par){
  
  func=function(x,dmns,par,what) {
    
    par[dmns] =exp(x)
    par["t0"] =-par["t0"]
        
    fnM=function(par,len,T=290,a=FLPar(c(a=-par["M1"],b=-par["M2"],c=1.5067827,d=0.9664798,e=763.5074169),iter=dims(par)$iter))
                     exp(a[1]+a[2]*log(len) + a[3]*log(par["linf"]) + a[4]*log(par["k"]) + a[5]/T)
   
    if ("vb" %in% dimnames(par)$params)
       res=lh(par,fnM=fnM,sr=list(model="bevholt",s=par["s"],   v=par["vb"]))
    else  
       res=lh(par,fnM=fnM,sr=list(model="bevholt",alpha    =par["alpha"],beta   =par["beta"]))
             
    fbar(res)=FLQuant(0.07,dimnames=list(year=1:41))
    prj=fwd(prj,catch=FLQuant(seq(40,80,length.out=40),dimnames=list(year=2:41)),sr=res)
    fbar(res)=FLQuant(seq(0.07,0.375, length.out=41))
    prj=fwd(brp(res))
    
    res=c(c(ssb(    prj)),
          c(stock(  prj)),
          c(catch(  prj)),
          c(catch(  prj)),
          c(rec(    prj)),
          c(ssb(    prj)/refpts(res)[what,"ssb"]),
          c(stock(  prj)/refpts(res)[what,"biomass"]),
          c(catch(  prj)/refpts(res)[what,"harvest"]),
          c(catch(  prj)/refpts(res)[what,"yield"]),
          c(rec(    prj)/refpts(res)[what,"rec"]))
      
    log(res)}
  
    jbn=jacobian(func,log(c(par)),dmns=dimnames(par)$params,par=par,what=what)
      
    res=data.frame(Year    =rep(seq(41),10),
                   Quantity=rep(rep(c("SSB","Biomass","Harvest","Yield","Recruits"),each=41),2),
                   Type    =rep(c("Absolute","Relative"),each=205),jbn)
    res=melt(res,id.var=c("Year","Quantity","Type"))
    res$Parameter      =factor(dimnames(par)$params[res$variable])
        
    p.=data.frame(Parameter=c("linf",  "t0",    "M1","M2","s",  "vb", "a",     "b",     "bg",      "k",       "ato95",      "sl",         "sr",      "a50",     "asym",       "a1",         "fec"),
                  Process  =c("Growth","Growth","M", "M", "SRR","SRR","Growth","Growth","Maturity","Growth",  "Maturity","Selectivity","Selectivity","Maturity","Maturity","Selectivity","Maturity"))
      
    res=merge(res,p.)
    
    return(res[,c("Year","Quantity","Type","Parameter","Process","value")])}
    
br =lh(gislasim(FLPar(linf=100)))
fig1=ggplot(br[[c("m","stock.wt","mat","catch.sel")]]) + theme_ms(8) +
    geom_line(aes(age,data))+facet_wrap(~qname,scale="free")+scale_x_continuous(limits=c(0,20))+
    xlab("Age")+ylab("")
fig1$data$qname=factor(fig1$data$qname,levels=c("stock.wt","m","mat","catch.sel"),labels=c("Mass","M","Proportion Mature","Selectivity"))

fbar(br)=seq(0,1,length.out=101)*refpts(br)["crash","harvest"]
refpts(br)=refpts(br)[c(4,1,5)]
dimnames(refpts(br))$refpt=c("MSY","F0.1","Fcrash")
fig2=plot(br) + theme_ms(8) + xlab("") + ylab("")

fbar(br)=rep(0.7,41)*refpts(br)["MSY","harvest"]
stk=fwd(brp(br))
stk=fwd(stk,catch=FLQuant(seq(1,1.25,length.out=31)*catch(stk)[,1,drop=T],dimnames=list(year=11:41)),sr=br)
fig3=plot(stk) + theme_ms(8)

fig4=kobe(data.frame(ssb=c(ssb(stk)/refpts(br)["MSY","ssb"]),f=c(fbar(stk)/refpts(br)["MSY","harvest"]),year=1:41, Decade=10*(1:41 %/% 10)))+
      geom_path( aes(ssb,f))                            +
      geom_point(aes(ssb,f,colour=Decade),size=2.25)    +
      scale_colour_gradient2(low="blue",high="red", midpoint=40)+theme_flr(size=20) +
      scale_x_continuous(limits=c(0,4))+scale_y_continuous(limits=c(0,4))           +
      xlab(expression(Stock/B[MSY])) + xlab("Year") + theme_ms(8)                               
 
ggsave(fig1,       filename=paste(dirMy,"/tex/fig1.png",sep=""), height=5,width=8)
ggsave(fig2,       filename=paste(dirMy,"/tex/fig2.png",sep=""), height=5,width=8)
ggsave(fig3,       filename=paste(dirMy,"/tex/fig3.png",sep=""), height=5,width=8)
ggsave(fig4,       filename=paste(dirMy,"/tex/fig4.png",sep=""), height=3,width=4)
         
par =gislasim(FLPar(linf=100,t0=.1,M1=2.1104327,M2=1.7023068,s=0.9,vb=1000,fec=1))
  
fec=mdply(c("f0.1","msy","crash"),doIt,par)
   
fec$BRP     =factor(c("MSY","F0.1","Fcrash")[fec$X1], levels=c("MSY","F0.1","Fcrash"))
fec$Process =factor(fec$Process,  levels=c("Growth","Maturity","M","SRR","Selectivity"))
fec$Quantity=factor(fec$Quantity, levels=c("Biomass","SSB","Harvest","Yield","Recruits"))
hvt=fec   
save(hvt,file=paste(dirMy,"/data/hvt.RData",sep=""))

load(paste(dirMy,"/data/fec.RData",sep=""))


cf=transform(subset(fec, Type=="Relative" & Parameter %in% c("M2","s") & BRP=="F0.1" & Quantity=="SSB"),Parameter=ac(Parameter))

cf=rbind(cf[,c("Year","Parameter","value","Type")],
      transform(as.data.frame((stk[["ssb"]][[1]]-refpts(br)["MSY","ssb"])/refpts(br)["MSY","ssb"]), value=data, Year=year, Parameter="SSB", Type="Absolute")[,c("Year","Parameter","value","Type")])
cf$What="Elasticity"
cf[cf$Type=="Absolute","What"]="Time Series"

fig5=ggplot(cf)+geom_line(aes(Year,value,group=Parameter,colour=Parameter))+facet_wrap(~What,ncol=1,scale="free")+
  theme_ms(8)+
  geom_hline(yintercept=0) +
    geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red") +
  ylab("")
ggsave(fig5,       filename=paste(dirMy,"/tex/fig5.png",sep=""), height=4,width=8)
 
fig6=ggplot(transform(subset(fec,Year %in% 5:40 & Type=="Relative" & Quantity=="SSB"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_hline(aes(yintercept=0), colour="black")+
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)    

fig7=ggplot(transform(subset(fec,Year %in% 5:40 & Type=="Relative" & Quantity=="Harvest"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_hline(aes(yintercept=0), colour="black")+
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)    

fig8=ggplot(transform(subset(fec,Year %in% 5:40 & Type=="Relative" & Quantity=="Biomass"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_hline(aes(yintercept=0), colour="black")+
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)    

grob <- ggplotGrob(fig6)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig6,       filename=paste(dirMy,"/tex/fig6.png",sep=""), height=4,width=4)
grob <- ggplotGrob(fig7)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig7,       filename=paste(dirMy,"/tex/fig7.png",sep=""), height=4,width=4)
