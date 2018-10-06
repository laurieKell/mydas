library(FLAdvice)
library(numDeriv)

#dirMy="/home/lkell/Desktop/Stuff/Publications/gbyp-sam/papers/journal/elastic"
dirMy="c:/Projects/gbyp-sam/trunk/papers/journal/elastic"
    
par=gislasim(FLPar(linf=100))

br =lh(par)
fig1=ggplot(br[[c("m","stock.wt","mat","catch.sel")]]) + theme_ms(8) +
    geom_line(aes(age,data))+facet_wrap(~qname,scale="free")+scale_x_continuous(limits=c(0,20))+
    xlab("Age")+ylab("")
fig1$data$qname=factor(fig1$data$qname,levels=c("stock.wt","m","mat","catch.sel"),labels=c("Mass","M","Proportion Mature","Selectivity"))

refpts(br)=refpts(br)[c(4,1,5)]
dimnames(refpts(br))$refpt=c("MSY","F0.1","Fcrash")
fig2=plot(br) + theme_ms(8) + xlab("") + ylab("")
 
fbar(br)=seq(0,1,length.out=20)*refpts(br)["Fcrash","harvest"]
br=brp(br)
fig3=kobe(data.frame(ssb=c(ssb(br)/refpts(br)["MSY","ssb"]),f=c(fbar(br)/refpts(br)["MSY","harvest"])))+
      geom_path( aes(ssb,f))                            +
      geom_point(aes(ssb,f),size=2.25)    +
      scale_colour_gradient2(low="blue",high="red", midpoint=40)+
      scale_x_continuous(limits=c(0,4))+scale_y_continuous(limits=c(0,4))           +
      xlab(expression(Stock/B[MSY])) + theme_ms(8)                               
 
ggsave(fig1,       filename=paste(dirMy,"/tex/fig1.png",sep=""), height=5,width=8)
ggsave(fig2,       filename=paste(dirMy,"/tex/fig2.png",sep=""), height=5,width=8)
ggsave(fig3,       filename=paste(dirMy,"/tex/fig3.png",sep=""), height=5,width=8)
         
#par =gislasim(FLPar(linf=100))#,t0=-.1,fec=1))
# To use doIt function, need to change sign of the negative parameters
# This is so logs of the magnitude can be calculated.
par["t0"]=-1*par["t0"]
par["M2"]=-1*par["M2"]

resEql=mdply(c("f0.1","msy","crash"),doIt,par,fbar=fbar(br))
save(resEql, file="resEqlraw.Rdata")
   
resEql=transform(resEql, BRP     =factor(c("MSY","F0.1","Fcrash")[X1], levels=c("MSY","F0.1","Fcrash")),
                         Process =factor(Process,                      levels=c("Growth","Maturity","M","SRR","Selectivity")),
                         Quantity=factor(Quantity,                     levels=c("Biomass","SSB","Harvest","Yield","Recruits")))
fb=as.data.frame(fbar(br),drop=T)
names(fb)=c("Year","F")
resEql=merge(resEql,fb)                 
save(resEql,file=paste(dirMy,"/data/resEql.RData",sep=""))

fec=resEql
fig4=ggplot(transform(subset(fec,Year %in% 1:15 & Type=="Relative" & Quantity=="SSB"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=3.3332e-01), colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(F,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)    

fig5=ggplot(transform(subset(fec,Year %in% 1:15 & Type=="Relative" & Quantity=="Harvest"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=3.3332e-01), colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(F,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)    

fig6=ggplot(transform(subset(fec,Year %in% 1:15 & Type=="Relative" & Quantity=="Biomass"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=3.3332e-01), colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(F,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)
                 
fig7=ggplot(transform(subset(fec,Year %in% 1:15 & Type=="Relative" & Quantity=="Yield"), Parameter=ac(Parameter))) + 
  geom_vline(aes(xintercept=3.3332e-01), colour="red")  +
  geom_hline(aes(yintercept=0),  colour="black")+
  geom_line(aes(F,value,colour=Parameter,group=paste(Parameter,Type)))+facet_grid(Process~BRP,scale="free")  + 
  ylab("") +theme_ms(5)

grob <- ggplotGrob(fig4)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig4,       filename=paste(dirMy,"/tex/fig4.png",sep=""), height=6,width=6)
grob <- ggplotGrob(fig5)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig5,       filename=paste(dirMy,"/tex/fig5.png",sep=""), height=6,width=6)
grob <- ggplotGrob(fig6)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig6,       filename=paste(dirMy,"/tex/fig6.png",sep=""), height=6,width=6)
grob <- ggplotGrob(fig7)
   strip_elem <- grid.ls(getGrob(grob, "strip.text.x", grep=TRUE, global=TRUE))$name
   grob <- geditGrob(grob, strip_elem[1], label=expression(F[MSY]))
   grob <- geditGrob(grob, strip_elem[2], label=expression(F[0.1]))
   grob <- geditGrob(grob, strip_elem[3], label=expression(F[Crash])) 
   grid.draw(grob)
ggsave(fig7,       filename=paste(dirMy,"/tex/fig7.png",sep=""), height=6,width=6)
