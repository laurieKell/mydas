library(plyr)
library(reshape)
library(ggplot2)
library(dplyr)

library(FLCore)
library(ggplotFL)
library(FLRP)
library(FLife)

library(LBSPR)

lbsprWrap<-function(len,params,species="",units="cm"){

  pars        =new("LB_pars")
  pars@Linf   =c(params["linf"]) 
  pars@L50    =vonB(c(params["a50"]),params) 
  pars@L95    =pars@L50+vonB(c(params["ato95"]),params)
  pars@MK     =c(params["mk"])
  pars@Species=species
  pars@L_units=units
  
  #labs=dimnames(len)[[1]]
  #brks=cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
  #           upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  #mid=aaply(brks,1,mean)
  
  LBlen       =new("LB_lengths")
  LBlen@LMids =as.numeric(dimnames(len)[[1]])
  LBlen@LData =len
  LBlen@Years =as.numeric(dimnames(len)[[2]])
  LBlen@NYears=dim(len)[2] 
  
  LBSPRfit(pars,LBlen)}

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/brill.RData")

om1 =window(iter(om,1),start=21,end=100)

n   =catch.n(om1)
ln  =vonB(ages(catch.n(om1)),iter(lh,1))
sd  =ln*0.2

bin =0:ceiling(max(ln)*1.10)+0.5

# Length frequency distribution
### Create lengths-at-age
lfq=ddply(model.frame(FLQuants(ln=ln,sd=sd,n=n)),.(age,year,unit,season,area,iter), 
           with, data.frame(length=bin,data=dnorm(bin,ln,sd)*n))

### sum up over ages 
lfq=ddply(lfq,.(length,year,unit,season,area,iter), 
           with, data.frame(freq=sum(data)))

### plot
ggplot(lfq)+
  geom_histogram(aes(length,weight=freq),binwidth=1)+
  facet_wrap(~year,scale="free")+
  xlab("Length (cm)")+ylab("Frequency")

# run LBSPR
len=t(daply(lfq,.(year,length), with, sum(freq)))

brl=lbsprWrap(len,prior[,1])

brl10=lbsprWrap(len[,ac(21,91,10)],prior[,1])

### Compare results
ggplot(data.frame("Year"=21:100,
               "MP"=brl@Ests[,"FM"],
               "OM"=c(apply(harvest(om[,ac(21:100)])%/%m(om[,ac(21:100)]),c(2),mean))))+
  geom_line(aes(Year,OM))+
  geom_line(aes(Year,MP))

