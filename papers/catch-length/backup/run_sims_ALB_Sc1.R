setwd("/home/laurence/Desktop/sea++/mydas/project/papers/maite")
source("sim.pop.R")
library("dplyr")

Nyears<-20
incr<-seq(0,by=0.25,length=9)
  stb<-rep(incr[9],5)
  decr<-seq(stb[4],by=-0.25,length=6)
Fdynamics<-c(incr,stb,decr) # Scenario 1, Dep=0.2


for(i in 1:100){
  set.seed(i)
  Outs<-sim_pop(Fdynamics=Fdynamics)
  saveRDS(object = Outs,file = paste("data/Outs_sim",i,".RDS",sep="")) #your folder name
}

### Plots
pdf("data/ALB_Sc1_Dep0.2.pdf",height = 8, width = 6)
par(mfrow=c(4,2),mar=c(0,4,0,1),oma=c(4,1,3,1))
Sim1<-readRDS(paste("data/Outs_sim",1,".RDS",sep=""))
plot(Sim1$Outs$Year,Sim1$Outs$Cw_t,type="n",lwd=3,xaxt="n",ylab="Catch (t)",ylim=c(0,7000))
tmp.mean<-data.frame(Sim1$Outs$Year,Sim1$Outs$Cw_t)
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$Cw_t,lwd=1,col="grey")
  tmp.mean[,i]<-Outs$Outs$Cw_t
}
lines(rowMeans(tmp.mean[,-1]),lwd=3)


plot(Sim1$Outs$Year,Fdynamics,type="n",lwd=3,ylab="Fishing mortality",xaxt="n",ylim=c(0,3))
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$F_t,lwd=1,col="grey")
}
lines(Sim1$Outs$Year,Fdynamics,lwd=3)

tmp.mean<-data.frame(Sim1$Outs$Year,Sim1$Outs$TB_t)
plot(Sim1$Outs$Year,Sim1$Outs$TB_t,type="n",lwd=3,ylab="Total Biomass (t)",xaxt="n",ylim=c(1000,35000))
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$TB_t,lwd=1,col="grey")
  tmp.mean[,i]<-Outs$Outs$TB_t
}
lines(rowMeans(tmp.mean[,-1]),lwd=3)

tmp.mean<-data.frame(Sim1$Outs$Year,Sim1$Outs$Dep)
plot(Sim1$Outs$Year,Sim1$Outs$Dep,type="n",lwd=3,ylab="Depletion",xaxt="n",ylim=c(0,1.1))
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$Dep,lwd=1,col="grey")
  tmp.mean[,i]<-Outs$Outs$Dep
}
lines(rowMeans(tmp.mean[,-1]),lwd=3)
abline(h=0.2,lty=2,lwd=2)

tmp.mean<-data.frame(Sim1$Outs$Year,Sim1$Outs$Mean.length)
plot(Sim1$Outs$Year,Sim1$Outs$Mean.length,type="n",lwd=3,ylab="Mean length (cm)",xaxt="n",ylim=c(70,120))
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$Mean.length,lwd=1,col="grey")
  tmp.mean[,i]<-Outs$Outs$Mean.length
}
lines(rowMeans(tmp.mean[,-1]),lwd=3)

tmp.mean<-data.frame(Sim1$Outs$Year,Sim1$Outs$SPR_t)
plot(Sim1$Outs$Year,Sim1$Outs$SPR_t,type="n",lwd=3,ylab="SPR",xaxt="n",ylim=c(0,1))
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$SPR_t,lwd=1,col="grey")
  tmp.mean[,i]<-Outs$Outs$SPR_t
}
lines(rowMeans(tmp.mean[,-1]),lwd=3)
abline(h=0.4,lty=2)

tmp.mean<-data.frame(Sim1$Outs$Year,Sim1$Outs$RecDev)
plot(Sim1$Outs$Year,Sim1$Outs$RecDev,type="n",lwd=3,ylab="Recruitment deviations",ylim=c(-2,2))
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$RecDev,lwd=1,col="grey")
  tmp.mean[,i]<-Outs$Outs$RecDev
}
lines(rowMeans(tmp.mean[,-1]),lwd=3)

abline(h=0,lty=2)

tmp.mean<-data.frame(Sim1$Outs$Year,Sim1$Outs$BBmsy)
plot(Sim1$Outs$Year,Sim1$Outs$BBmsy,type="n",lwd=3,ylab="B/Bmsy",ylim=c(0,3.3))
for(i in 2:100){
  Outs<-readRDS(paste("data/Outs_sim",i,".RDS",sep=""))
  lines(Outs$Outs$Year,Outs$Outs$BBmsy,lwd=1,col="grey")
  tmp.mean[,i]<-Outs$Outs$BBmsy
}
lines(rowMeans(tmp.mean[,-1]),lwd=3)
abline(h=1,lty=2)
mtext(text = "Year",side = 1,line = 2.5,outer = T)
mtext(text = "Albacore_Harvest Scenario 1_Dep=0.2",side = 3,line = 1,outer = T)
dev.off()
