library(FLBRP)
library(numDeriv)

dirMy="/home/laurie/Ubuntu One/papers/Journal/elasticity"

load(paste(dirMy,"data/par.RData",sep="/"))

what=c("msy","f0.1")
nyrs=51
#par=FLPar(par[,1,1,drop=T])

par =gislasim(FLPar(linf=100,t0=.1,M1=2.1104327,M2=1.7023068,s=0.9,v=1000,fec=1))

#M
fnM=function(par,len,T=290){
  a=FLPar(c(a=-par["M1"],b=-par["M2"],c=1.5067827,d=0.9664798,e=763.5074169))
  exp(a[1]+a[2]*log(len) + a[3]*log(par["linf"]) + a[4]*log(par["k"]) + a[5]/T)}

fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))

#Fbar
fnFbar=function(x,nyr=51) seq(0,0.75,length.out=nyr)*2*refpts(x)["msy","harvest"]

elasticity=function(what,par,nyr=51,refYrs=21:25,fnFbar=fnFbar,fnM=fnM){
  
  func=function(x,par,dmns=dimnames(par)$params,
                what=c("msy","fmax","f0.1","spr.30","crash"),
                nyr=51,refYrs=21:25,
                fnFbar=fnFbar,fnM=fnM) {
  
    ## reset any negative values
    par.     =par
    par[dmns]=exp(x)h
    if (any(par.<0)) par[par.<0]=-par[par.<0]
        
    res=lh(par,fnM=fnM)
    fbar(res)=fnFbar(res) 
 
    stk=fwd(brp(res),maxF=3)
    stk=fwd(stk,f=fbar(res)[,-1],sr=res)
 
    res=c(c(ssb( stk)),
          c(stock(stk)),
          c(fbar( stk)),
          c(catch(stk)),
          c(rec(  stk)),
    
          c(ssb(  stk)/mean(ssb(  stk)[,refYrs])),
          c(stock(stk)/mean(stock(stk)[,refYrs])),
          c(fbar( stk)/mean(fbar( stk)[,refYrs])),
          c(catch(stk)/mean(catch(stk)[,refYrs])),
          c(rec(  stk)/mean(rec(  stk)[,refYrs])),
          
          unlist(mdply(what, function(x) c(ssb(  stk)/refpts(res)[x,"ssb"]    ))[,-1]),
          unlist(mdply(what, function(x) c(stock(stk)/refpts(res)[x,"biomass"]))[,-1]),
          unlist(mdply(what, function(x) c(fbar( stk)/refpts(res)[x,"harvest"]))[,-1]),
          unlist(mdply(what, function(x) c(catch(stk)/refpts(res)[x,"yield"]  ))[,-1]),
          unlist(mdply(what, function(x) c(rec(  stk)/refpts(res)[x,"rec"]    ))[,-1])
          )
      
    log(res)}
  
#     res=func(log(abs(par)),par,dmns=dimnames(par)$params,
#              what,
#              nyr=51,refYrs=21:25,
#              fnFbar=function(x,nyr=51) seq(0,0.75,length.out=nyr)*2*refpts(x)["msy","harvest"],
#              fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len)))
    
    jbn=data.frame(
        jacobian(func,log(abs(par)),par=par,what=what,
                 nyr=51,refYrs=21:25,
                 fnFbar=function(x,nyr=51) seq(0,0.75,length.out=nyr)*2*refpts(x)["msy","harvest"],
                 fnM=function(par,len) exp(par["M1"]-par["M2"]*log(len))))
    names(jbn)=dimnames(par)$params

    Process       =c("Growth","Growth","M", "M", "SRR","SRR","Growth","Growth","Maturity","Growth","Maturity","Selectivity","Selectivity","Maturity","Maturity","Selectivity","Maturity")
    names(Process)=c("linf",  "t0",    "M1","M2","s",  "v",  "a",     "b",     "bg",      "k",     "ato95",   "sl",         "sr",         "a50",     "asym",    "a1",         "fec")

    factors=expand.grid(Refpt    =c("Absolute","Relative",what),
                        Quantity =c("SSB","Biomass","Harvest","Yield","Recruits"),
                        Year     =1:nyrs,
                        Parameter=dimnames(par)$params)
  
    res=data.frame(factors,jbn)
    res=melt(res,id.var=c("Year","Quantity","Refpt","Parameter"))
    res=cbind(res,Process=factor(Process[res$Parameter]))
    
    return(res[,c("Year","Quantity","Refpt","Parameter","Process","value")])}
  
par =gislasim(FLPar(linf=100,t0=.1,M1=2.1104327,M2=1.7023068,s=0.9,v=1000,fec=1))
j1  =elasticity(what,par,nyr=51,refYrs=21:25,fnFbar=fnFbar,fnM=fnM)
  

res         =mdply(c("f0.1","spr.30","fmax","msy","crash"),elasticity,par)
res$BRP     =factor(c("F0.1","spr%30","Fmax","MSY","Fcrash")[res$X1], levels=c("MSY","F0.1","Fmax","spr%30","Fcrash"))
res$Process =factor(res$Process,  levels=c("Growth","Maturity","M","SRR","Selectivity"))
res$Quantity=factor(res$Quantity, levels=c("Biomass","SSB","Harvest","Yield","Recruits"))

save(res,file=paste(dirMy,"/data/exStk.RData",sep=""))

br =lh(gislasim(FLPar(linf=100)))
fig1=ggplot(br[[c("m","stock.wt","mat","catch.sel")]]) + theme_ms(8) +
    geom_line(aes(age,data))+facet_wrap(~qname,scale="free")+scale_x_continuous(limits=c(0,20))+
    xlab("Age")+ylab("")
fig1$data$qname=factor(fig1$data$qname,levels=c("stock.wt","m","mat","catch.sel"),labels=c("Mass","M","Proportion Mature","Selectivity"))

fbar(br)=seq(0,1,length.out=101)*refpts(br)["crash","harvest"]
dimnames(refpts(br))$refpt=c("F0.1","Fmax","spr30%","MSY","Fcrash")
fig2=plot(br) + theme_ms(8) + xlab("") + ylab("")

fbar(br)=seq(0,0.75,length.out=51)*refpts(br)["Fcrash","harvest"]
stk=fwd(brp(br),maxF=3)
fig3=plot(stk) + theme_ms(8)

fig4=kobe(data.frame(ssb=c(ssb(stk)/refpts(br)["MSY","ssb"]),f=c(fbar(stk)/refpts(br)["MSY","harvest"]),year=1:51, Decade=10*(1:51 %/% 10)))+
      geom_path( aes(ssb,f))                            +
      geom_point(aes(ssb,f,colour=Decade),size=2.25)    +
      scale_colour_gradient2(low="blue",high="red", midpoint=40)+theme_flr(size=20) +
      scale_x_continuous(limits=c(0,4))+scale_y_continuous(limits=c(0,4))           +
      xlab(expression(Stock/B[MSY])) + xlab("Year") + theme_ms(8)                               


fig5=ggplot(transform(subset(res,Year %in% 1:40 & Process=="Selectivity" & Type=="Relative"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Quantity~BRP,scale="free")  + 
  ylab("") +theme_ms(8)    

fig6=ggplot(transform(subset(res,Year %in% 1:40 & Process=="Growth" & Type=="Relative"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Quantity~BRP,scale="free") + 
  ylab("") +theme_ms(8)    

fig7=ggplot(transform(subset(res,Year %in% 1:40 & Process=="Maturity" & Type=="Relative"), Parameter=ac(Parameter)))+ 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Quantity~BRP,scale="free") + 
  ylab("") +theme_ms(8)    

fig8=ggplot(transform(subset(res,Year %in% 1:40 & Process=="M" & Type=="Relative"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Quantity~BRP,scale="free")  + 
  ylab("") +theme_ms(8)    

fig9=ggplot(transform(subset(res,Year %in% 1:40 & Process=="SRR" & Type=="Relative"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Quantity~BRP,scale="free")  + 
  ylab("") +theme_ms(8)    

fig10=ggplot(transform(subset(res,Year %in% 1:40 & Type=="Absolute"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Quantity~Process,scale="free") + 
  ylab("") +theme_ms(8)    

fig11=ggplot(transform(subset(res,Year %in% 1:40 & Type=="Absolute"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~Quantity,scale="free") + 
  ylab("") +theme_ms(8)    

fig12=ggplot(transform(subset(res,Year %in% 1:40 & Type=="Relative" & BRP=="MSY"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Quantity~Process,scale="free") + 
  ylab("") +theme_ms(8)    

fig13=ggplot(transform(subset(res,Year %in% 1:40 & Type=="Relative" & BRP=="MSY"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red")  +
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~Quantity,scale="free") + 
  ylab("") +theme_ms(8)    

redo=function(fig,flnm,height,width){
  
   grob <- ggplotGrob(fig)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
   }

ggsave(fig1,       filename=paste(dirMy,/tex/fig1.png",sep=""), height=5,width=8)
ggsave(fig2,       filename=paste(dirMy,/tex/fig2.png",sep=""), height=5,width=8)
ggsave(fig3,       filename=paste(dirMy,/tex/fig3.png",sep=""), height=5,width=8)
ggsave(fig4,       filename=paste(dirMy,/tex/fig4.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig5)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig5,       filename=paste(dirMy,/tex/fig5.png",sep=""),sep=""), height=6,width=10)
grob <- ggplotGrob(fig6)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig6,       filename=paste(dirMy,/tex/fig6.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig7)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob) 
ggsave(fig7,       filename=paste(dirMy,/tex/fig7.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig8)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig8,       filename=paste(dirMy,/tex/fig8.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig9)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig9,       filename=paste(dirMy,/tex/fig9.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig10)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig10,      filename=paste(dirMy,/tex/fig10.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig11)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig11,      filename=paste(dirMy,/tex/fig11.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig12)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig12,      filename=paste(dirMy,/tex/fig12.png",sep=""), height=6,width=10)
grob <- ggplotGrob(fig13)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Max]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[4], label=expression(F["SPR%30"]))
   grob <- geditGrob(grob, strip_elem[5], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig13,      filename=paste(dirMy,/tex/fig13.png",sep=""), height=6,width=10)

cf=transform(subset(res, Type=="Relative" & Parameter %in% c("M1","s") & BRP=="MSY" & Quantity=="SSB"),Parameter=ac(Parameter))

cf=rbind(cf[,c("Year","Parameter","value","Type")],
      transform(as.data.frame((stk[["ssb"]][[1]]-refpts(br)["MSY","ssb"])/refpts(br)["MSY","ssb"]), value=data, Year=year, Parameter="SSB", Type="Absolute")[,c("Year","Parameter","value","Type")])
cf$What="Elasticity"
cf[cf$Type=="Absolute","What"]="Time Series"

fig14=ggplot(cf)+geom_line(aes(Year,value,group=Parameter,colour=Parameter))+facet_wrap(~What,ncol=1,scale="free")+
  theme_ms(8)+
  geom_hline(yintercept=0) +
    geom_vline(aes(xintercept=7), colour="green")+
  geom_vline(aes(xintercept=16),colour="red") +
  ylab("")

ggsave(fig14,      filename=paste(dirMy,/tex/fig14.png",sep=""), height=6,width=10)



  
