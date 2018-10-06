#######################################################
## LBSPR
#######################################################

library(LBSPR)
         
minVal <- 0
maxVal <- 600
mn <- (maxVal - minVal)/2   

   LB_lengths       =new("LB_lengths")
   LB_lengths@LMids =seq(2.5,150,5)
   LB_lengths@LData =t(daply(data.frame(year=51:55,len=12),.(year), with, table(cut(rlnorm(200,log(len),.2),seq(1,27,1)))))
   LB_lengths@Years =51:55
   LB_lengths@NYears=5           

plotSize(LB_lengths)

## Example from package
datdir <- DataDir()
MyPars <- new("LB_pars")
MyPars@Species <- "MySpecies"
MyPars@Linf <- 100 
MyPars@L50 <- 66 
MyPars@L95 <- 70
MyPars@MK <- 1.5 
MyPars@L_units <- "mm"
Len1 <- new("LB_lengths", LB_pars=MyPars, file=paste0(datdir, "/LFreq_MultiYr.csv"), 
            dataType="freq")

plotSize(Len1)


   LB_pars <- new("LB_pars")
   LB_pars@Species <- "Brill"
   LB_pars@MK <- mean(m(eq))/FLPar(aaply(lh,1,mean))["k"][[1]]
   LB_pars@Linf <- mean(lh$linf)
   LB_pars@L50 <- mean(lh$l50)
   LB_pars@L95 <- mean(lh$linf)
   #LB_pars@SL50 <- 50
   #LB_pars@SL95 <- 60
   #LB_pars@SPR <- 0.4
   LB_pars@Walpha <- mean(lh$a)
   LB_pars@Wbeta <-  mean(lh$b)
   LB_pars@BinWidth <- 5 
   #LB_pars@R0 <- 1
   #LB_pars@Steepness <- ifelse(lh$h==1, 0.99, lh$h)

   lbspr_res <- LBSPRfit(LB_pars=LB_pars, LB_lengths=LB_lengths, Control=list(modtype=c("GTG")))


