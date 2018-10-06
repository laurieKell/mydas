library(LBSPR)

lbsprWrap<-function(len,params,species="",units="cm"){

  pars        =new("LB_pars")
  pars@Linf   =c(params["linf"]) 
  pars@L50    =vonB(c(params["a50"]),params) 
  pars@L95    =pars@L50+vonB(c(params["ato95"]),params)
  pars@MK     =c(params["mk"])
  pars@Species=species
  pars@L_units=units
  
  labs=dimnames(len)[[1]]
  brks=cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
             upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
  mid=aaply(brks,1,mean)

  LBlen       =new("LB_lengths")
  LBlen@LMids =mid
  LBlen@LData =len
  LBlen@Years =as.numeric(dimnames(len)[[2]])
  LBlen@NYears=dim(len)[2] 
  
  LBSPRfit(pars,LBlen)}

# len=t(daply(data.frame(year=51:55,len=20),.(year), with, table(cut(rlnorm(200,log(len),.2),seq(1,40,1)))))
# 
# brl=lbsprWrap(len,prior[,1])
# 
# plotSize(brl)
# 
# plotMat(brl)
# 
# plotEsts(brl)
# 
# brl@Ests