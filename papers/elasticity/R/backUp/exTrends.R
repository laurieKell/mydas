library(FLAdvice)
library(numDeriv)

dirMy="/home/lkell/Desktop/Stuff/Publications/gbyp-sam/papers/journal/elastic"

    
br =lh(gislasim(FLPar(linf=100)))
figDyn1=ggplot(br[[c("m","stock.wt","mat","catch.sel")]]) + theme_ms(8) +
    geom_line(aes(age,data))+facet_wrap(~qname,scale="free")+scale_x_continuous(limits=c(0,20))+
    xlab("Age")+ylab("")
figDyn1$data$qname=factor(figDyn1$data$qname,levels=c("stock.wt","m","mat","catch.sel"),labels=c("Mass","M","Proportion Mature","Selectivity"))

fbar(br)=seq(0,1,length.out=101)*refpts(br)["crash","harvest"]
refpts(br)=refpts(br)[c(4,1,5)]
dimnames(refpts(br))$refpt=c("MSY","F0.1","Fcrash")
figDyn2=plot(br) + theme_ms(8) + xlab("") + ylab("")

fbar(br)=FLQuant(c(seq(.0,.75,length.out=21),seq(.75,.0,length.out=21)[-1]),dimnames=list(year=1:41))
fbar(br)=c(fbar(br)[,1:37],rep(refpts(br)["MSY","harvest"]*.60,19))
stk=fwd(brp(br))
figDyn3=plot(stk) + theme_ms(8)

figDyn4=kobe(data.frame(ssb=c(ssb(stk)/refpts(br)["MSY","ssb"]),f=c(fbar(stk)/refpts(br)["MSY","harvest"]),year=1:56, Decade=10*(1:56 %/% 10)))+
      geom_path( aes(ssb,f))                            +
      geom_point(aes(ssb,f,colour=Decade),size=2.25)    +
      scale_colour_gradient2(low="blue",high="red", midpoint=40)+theme_flr(size=20) +
      scale_x_continuous(limits=c(0,5))+scale_y_continuous(limits=c(0,5))           +
      xlab(expression(Stock/B[MSY])) + xlab("Year") + theme_ms(8)                               
 
ggsave(figDyn1,       filename=paste(dirMy,"/tex/figDyn1.png",sep=""), height=5,width=8)
ggsave(figDyn2,       filename=paste(dirMy,"/tex/figDyn2.png",sep=""), height=5,width=8)
ggsave(figDyn3,       filename=paste(dirMy,"/tex/figDyn3.png",sep=""), height=5,width=8)
ggsave(figDyn4,       filename=paste(dirMy,"/tex/figDyn4.png",sep=""), height=3,width=4)
         
par =gislasim(FLPar(linf=100,t0=.1,M1=2.1104327,M2=1.7023068,s=0.75,vb=1000,fec=1))

resDyn=mdply(c("msy","crash"),doIt,par,dynamic=TRUE,fbar=fbar(br))   
resDyn=transform(resDyn, BRP     =factor(c("MSY","Fcrash")[X1], levels=c("MSY","Fcrash")),
                         Process =factor(Process,               levels=c("Growth","Maturity","M","SRR","Selectivity")),
                         Quantity=factor(Quantity,              levels=c("Biomass","SSB","Harvest","Yield","Recruits")))
save(resDyn,file=paste(dirMy,"/data/resDyn.RData",sep=""))

load(paste(dirMy,"/data/resDyn.RData",sep=""))

fec=resDyn
figDyn6=ggplot(transform(subset(fec,Year %in% 1:56 & Type=="Relative" & Quantity=="SSB"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=6.5), colour="green")  +
  geom_vline(aes(xintercept=13),  colour="red")  +
  geom_vline(aes(xintercept=35),  colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)    

figDyn7=ggplot(transform(subset(fec,Year %in% 1:56 & Type=="Relative" & Quantity=="Harvest"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=6.5), colour="green")  +
  geom_vline(aes(xintercept=13),  colour="red")  +
  geom_vline(aes(xintercept=35),  colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)    

figDyn8=ggplot(transform(subset(fec,Year %in% 1:56 & Type=="Relative" & Quantity=="Biomass"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=6.5), colour="green")  +
  geom_vline(aes(xintercept=13),  colour="red")  +
  geom_vline(aes(xintercept=35),  colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)
                 
figDyn9=ggplot(transform(subset(fec,Year %in% 1:56 & Type=="Relative" & Quantity=="Yield"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=6.5), colour="green")  +
  geom_vline(aes(xintercept=13),  colour="red")  +
  geom_vline(aes(xintercept=35),  colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(Year,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)

grob <- ggplotGrob(figDyn6)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(figDyn6,       filename=paste(dirMy,"/tex/figDyn6.png",sep=""), height=4,width=4)
grob <- ggplotGrob(figDyn7)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(figDyn7,       filename=paste(dirMy,"/tex/figDyn7.png",sep=""), height=4,width=4)
grob <- ggplotGrob(figDyn8)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(figDyn8,       filename=paste(dirMy,"/tex/figDyn8.png",sep=""), height=4,width=4)
grob <- ggplotGrob(figDyn9)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(figDyn9,       filename=paste(dirMy,"/tex/figDyn9.png",sep=""), height=4,width=4)


sn=sweep(stock.n(stk),2,apply(stock.n(stk),2,sum),"/")
sn=sweep(sn,1,apply(sn,1,mean),"/")
ggplot(sn[1:12])+geom_point(aes(year,age,size=data))


