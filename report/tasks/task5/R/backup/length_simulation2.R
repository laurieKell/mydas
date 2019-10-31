library(plyr)
library(reshape)
library(ggplot2)
library(dplyr)

library(FLCore)
library(ggplotFL)
library(FLRP)
library(FLife)

library(LBSPR)

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/brill.RData")

om1 =window(iter(om,1:2),start=61,end=70)
ln  =vonB(ages(stock.n(om1)),iter(lh,1:2))
sd  =ln*0.2
n   =stock.n(iter(om1,1))
bin =0:ceiling(max(ln)*1.10)+0.5

res2=ddply(model.frame(FLQuants(ln=ln,sd=sd,n=n)),.(age,year,unit,season,area,iter), 
           with, data.frame(length=bin,data=dnorm(bin,ln,sd)*n))

res3=ddply(res2,.(length,year,unit,season,area,iter), 
           with, data.frame(freq=sum(data)))

ggplot(res3)+
  geom_histogram(aes(length,weight=freq),binwidth=1)+
  facet_grid(year~iter,scale="free")+
  xlab("Length (cm)")+ylab("Frequency")


le = 1:c(prior["linf"]) 
#integrate rather than using a midpoint
leA = seq(1.5,c(prior["linf"])+ 0.5 ,1)  # upper point
leB = seq(0.5,c(prior["linf"])- 0.5 ,1)  # lower point

# the proportion at age assuming no fishing and M = 0.3
result=list()
mylist=list()

for(j in 1:100){


N = 1
#put noise on natural mortality
M = rlnorm(j,log(c(prior["k"])/c(prior["mk"])),0.2) #default 0.3
for(i in 1:length(a)-1){
  N = c(N, exp(-M[j])*N[i])
}

Ntot = sum(N)
p = (N/Ntot)*c(om1@stock.n[,j,,,,]*10000)


# the length distribution * proportion at age

ldist = NULL
for(i in a) {
  d2 = (pnorm((leA-mu[i])/sdev[i]) - pnorm((leB-mu[i])/sdev[i]))*p[i]
  ldist = rbind(ldist, d2)
}
#transpose matrices
tmatrix = t(ldist)
propatatlen = setNames(melt(tmatrix), c('length', 'x1', 'prop'))
outlen  = ddply(propatatlen, .(length), summarise, prop=sum(prop))
mylist = data.frame(year=paste(j), length=outlen[,1], prop=outlen[,2])
result = data.frame(rbind(result,mylist))
rownames(result) = 1:nrow(result)
}

dist=cast(result, length~year,sum)
len=data.matrix(dist, rownames.force = NA)

matplot(dist$length, dist[,2], type="l")
# these can be plotted individually
matplot(le, t(ldist), type="l")


  pars        =new("LB_pars")
  pars@Linf   =c(prior["linf"]) 
  pars@L50    =vonB(c(prior["a50"]),prior) 
  pars@L95    =pars@L50+vonB(c(prior["ato95"]),prior)
  pars@MK     =c(prior["mk"])
  pars@Species="brill"
  pars@L_units="cm"
  pars@Walpha = c(prior["a"])
  pars@Wbeta  = c(prior["b"])

  mid=as.numeric(dimnames(len)[[1]])+0.5
  
  LBlen       =new("LB_lengths")
  LBlen@LMids =mid
  LBlen@LData =len[,2:101]
  LBlen@Years =as.numeric(dimnames(len)[[2]][2:101])
  LBlen@NYears=length(LBlen@Years)
  
  brl = LBSPRfit(LB_pars=pars, LB_lengths=LBlen)
  
  plotSize(brl)

  plotMat(brl)

  plotEsts(brl)
  
  pars@SL50=brl@SL50[1]
  pars@SL95=brl@SL95[1]
  pars@SPR =0.3
  pars@FM =numeric()
  sim =LBSPRsim(pars, Control=list(modtype=c("GTG")))
  plotSim(sim)
 