#### Simulated population using an age structude model #################

sim_pop<-function(Nyears=20,Fdynamics=Fdynamics,
                  SigmaR=0.5,
                  rho=0, # no autocorrelation
                  SigmaF=0.2,
                  SigmaI=0.2,
                  SigmaC=0.1,
                  R0=1000,
                  AgeMax = 15,
                  start_ages=0,
                  # including age 0 and MaxAge = 15
                  binwidth=2,
                  linf=122,
                  vbk=0.209,
                  t0=(-1.338),
                  lwa=1.34e-05,
                  lwb=3.1066,
                  CVlen=0.1,
                  S50=85,
                  S95=90,
                  M50=90,
                  M95=100,
                  M=0.3,
                  h=0.7,
                  qcoef=0.00001,
                  Nl_sample=1000){

    ages <- seq(start_ages, to=AgeMax, by=1) 
    
    
    ## recruitment deviations
    RecDev <- c(0,rnorm(Nyears-1, mean = 0, sd = SigmaR))
    #plot(RecDev,type="l")
    ## autocorrelated recruitment deviations
    RecDev_AR <- rep(NA, length(RecDev))
    RecDev_AR[1] <- RecDev[1]
    for (t in 2:length(RecDev)) {
      RecDev_AR[t] <-
        RecDev_AR[t - 1] * rho + sqrt(1 - rho ^ 2) * RecDev[t]
    }
    #plot(RecDev_AR,type="l")
    
    ## fishing mortality deviations
    FishDev <-  c(1,rnorm(Nyears-1, mean = -(SigmaF ^ 2) / 2, sd = SigmaF))
    
    ## abundance index observation error
    IndexDev <-  rnorm(Nyears, mean = -(SigmaI ^ 2) / 2, sd = SigmaI)
    
    ## catch observation error
    CatchDev <- rnorm(Nyears, mean = -(SigmaC ^ 2) / 2, sd = SigmaC)
    
    # Rdynamics=="BevertonHolt"
    R_t <- rep(NA, Nyears)
    R_t[1] <- R0 * exp(RecDev_AR[1]) # RecDev the first year is 0 for now so the first year is equal to 1
    
    # Fdynamics
    Fdynamics <- Fdynamics
    F_t <- Fdynamics * exp(FishDev)
    
    #####################################################################################################  
    ## length bins
    mids <- seq((binwidth/2), linf*1.3, by=binwidth) 
    highs <- mids + (binwidth/2)
    lows <- mids - (binwidth)/2
    ## growth at age
    L_a <- linf*(1-exp(-vbk*(ages - t0)))
    W_a <- lwa*L_a^lwb  
    W_l <- lwa*mids^lwb
    
    ##probability of being in a length bin given age
    lbprobs <- function(mnl,sdl) return(pnorm(highs,mnl,sdl)-pnorm(lows,mnl,sdl))
    vlprobs <- Vectorize(lbprobs,vectorize.args=c("mnl","sdl"))
    plba <- t(vlprobs(L_a,L_a*CVlen))
    plba <- plba/rowSums(plba)
    dim(plba) #age-length prob matrix
    
    ## maturity 
    ML50 <- M50
    ML95 <- M95
    Mat_l <- (1 /(1 + exp(-log(19)*(mids-ML50)/(ML95-ML50)))) # Maturity at length
    Mat_a <- apply(t(plba)*Mat_l, 2, sum) # Maturity at age
    
    # selectivity
    SL50 <- S50
    SL95 <- S95
    
    S_l<- (1 /(1 + exp(-log(19)*(mids-SL50)/(SL95-SL50))))
    S_a <- apply(t(plba)*S_l, 2, sum) 
    S_a[1] <- 1e-6 # it should be zero for age 0
    
    ## fishing mortality associated with each age (selectivity)
    F_at <- matrix(NA, nrow=length(ages), ncol=Nyears)
    for(t in 1:Nyears){
      for(a in 1:length(ages)){
        F_at[a,t] <- F_t[t] * S_a[a]
      }
    }
    
    ##########################
    ## Data objects
    ##########################
    ## year 1 abundance at age
    N_at <- N_at0 <- matrix(0, nrow = length(L_a), ncol = Nyears) # 
    
    for (a in 1:length(L_a)) {
      if (a == 1) {
        N_at[a,1] <- R_t[1]
        N_at0[a,1] <- R_t[1]
      }
      if (a > 1 & a < length(L_a)) {
        N_at[a,1] <- N_at[a - 1, 1] * exp(-M - F_at[a-1,1]) # with fishing
        N_at0[a,1] <- N_at0[a - 1, 1] * exp(-M)
      }            
      if (a == length(L_a)) {
        N_at[a, 1] <- (N_at[a - 1, 1] * exp(-M - F_at[a-1,1])) / (1 - exp(-M - F_at[a-1,1]))
        N_at0[a, 1] <- (N_at0[a - 1, 1] * exp(-M)) / (1 - exp(-M))        
      }
    }
    
    ## year 1 biomass quantities
    TB_t <- SB_t <- rep(NA, Nyears)
    TB_t[1] <- sum(N_at[, 1] * W_a) # kilos 
    SB_t[1] <- sum(N_at[, 1] * W_a * Mat_a) # kilos 
    
    ## year 1 catch
    Cn_at <- Cw_at <- array(NA, c(length(L_a), Nyears))
    Cn_at[,1] <- N_at[,1] * (1 - exp(-M - F_at[,1])) * (F_at[,1] / (M + F_at[a,1])) # in numbers
    Cw_at[,1] <- Cn_at[,1] * W_a  # in weight kilos
    
    # to Calculate SB0
    calc_equil_abund <- function(ages, M, F, S_a, R0){  
      
      N_a <- rep(NA, length(ages))
      N_a[1] <- R0
      Fvec <- S_a*F
      
      for(i in 2:length(ages)){
        if(i<length(ages)) N_a[i] <- N_a[i-1]*exp(-(Fvec[i-1]+M))
        if(i==length(ages)) N_a[i] <- N_a[i-1]*(exp(-(Fvec[i-1]+M)))/(1-exp(-(Fvec[i-1]+M)))
      }
      return(N_a) # numbers at age at equilibruim
    }
    
    Na0<-calc_equil_abund (ages=ages, M=M, F=0, S_a=S_a, R0=R0)
    ## unfished spawning biomass
    B0<-sum(Na0 * W_a)  
    SB0 <- sum(Na0 * W_a * Mat_a)  
    
    # MSY calculations
    calc_msy <- function(F, ages, M, R0, W_a, S_a){
      Nage <- calc_equil_abund(ages=ages, M=M, F=F, R0=R0, S_a=S_a)
      F_a<-F*S_a
      C_a<- Nage * F_a/(F_a+M) * (1-exp(-(F_a+M))) #Baranov equation 
      YPR <- sum(C_a*W_a) # Yield per recruit
      return(YPR)
    }
    
    #calc_msy (ages=ages, M=M, F=7, S_a=S_a, R0=R0, W_a=W_a)
    
    vec<-seq(0,9,0.01)
    Fs<-rep(NA,length(vec))
    for(f in 1:length(vec)){
      Fs[f]<-calc_msy(F=vec[f],ages=ages, M=M, R0=R0, W_a=W_a, S_a=S_a)
    }
    #plot(vec,Fs,ylab="Yield")
    Fmsy <- optimize(calc_msy, ages=ages, M=M, R0=R0, W_a=W_a, S_a=S_a, lower=0, upper=10, maximum=TRUE)$maximum
    
    msy <- calc_msy(F=Fmsy, ages=ages, M=M, R0=R0, W_a=W_a, S_a=S_a)
    TBmsy = sum(calc_equil_abund(ages=ages, M=M, F=Fmsy, S_a=S_a, R0=R0)*W_a)
    
    
    for (t in 2:Nyears) {
      ## fishing effort based on spawning biomass
      R_t[t] <- ((4 * h * R0 * SB_t[t-1]) / (SB0 * (1-h) + SB_t[t-1] * (5*h-1))) * exp(RecDev_AR[t])
      
      ## age-structured dynamics
      for (a in 1:length(L_a)) {
        if (a == 1) {
          N_at[a, t] <- R_t[t]
          N_at0[a, t] <- R_t[t]
        }
        
        if (a > 1 & a < length(L_a)) {
          N_at[a, t] <- N_at[a - 1, t - 1] * exp(-M - F_at[a-1,t-1])
          N_at0[a, t] <- N_at0[a - 1, t - 1] * exp(-M)
        }              
        
        if (a == length(L_a)) {
          N_at[a, t] <- (N_at[a - 1, t - 1] * exp(-M - F_at[a-1,t-1])) + (N_at[a, t - 1] * exp(-M - F_at[a,t-1]))
          N_at0[a, t] <- (N_at0[a - 1, t - 1] * exp(-M)) + (N_at0[a, t - 1] * exp(-M))
        }
      }
      ## spawning biomass
      SB_t[t] <- sum((N_at[, t] * W_a * Mat_a))
      TB_t[t] <- sum(N_at[, t] * W_a)
      
      ## catch
      Cn_at[,t] <- N_at[,t] * (1 - exp(-M - F_at[,t])) * (F_at[,t] / (M + F_at[,t]))
      Cw_at[,t] <- Cn_at[,t] * W_a
    }
    
    ################################################### funtion needed for calc_equil
    calc_SPR <- function(ages, Mat_a, W_a, M, F, S_a){
      ## calculate spawning biomass per recruit in fished and unfished conditions
      Na0 <- calc_equil_abund(ages=ages, M=M, F=0, S_a=S_a, R0=1)
      Naf <- calc_equil_abund(ages=ages, M=M, F=F, S_a=S_a, R0=1)
      
      Nofish <- sum(Na0*Mat_a*W_a)
      Fish <- sum(Naf*Mat_a*W_a)
      
      ## SPR
      ratio <- Fish/Nofish
      return(ratio)
    }
    
    SPR_t <-
      sapply(1:length(F_t), function(x)
        calc_SPR(
          ages = ages,
          Mat_a = Mat_a,
          W_a = W_a,
          M = M,
          S_a = S_a,
          F = F_t[x]
        ))
    SPR <- SPR_t[length(SPR_t)]
    
    #plot(SPR_t)
    
    Cn_t <- colSums(Cn_at)
    Cw_t <- colSums(Cw_at)
    
    
    ## relative spawning biomass (depletion)
    D_t <- SB_t / SB0 
    #plot(D_t)
    
    ## abundance index
    I_t <- qcoef * TB_t * exp(IndexDev)
    Effort_t <- Cw_t/(qcoef*TB_t)
    
    ######### Age-converted-to-length structure for compositional data
    ################################################
    ## Probability being in harvested at an age
    ################################################
    page <- matrix(ncol = dim(plba)[1], nrow = Nyears)
    for (y in 1:Nyears){
      page[y,] <- N_at[,y] * S_a
    }
    page <- page / rowSums(page)
    
    ################################################
    ## Probability of sampling a given length bin
    ################################################
    plb <- matrix(ncol = length(highs), nrow = Nyears)
    for (y in 1:Nyears)
      plb[y, ] <- page[y, ] %*% plba
    plb <- plb / rowSums(plb)
    
    #######################
    ## Length frequencies
    #######################
    comp_sample<-rep(Nl_sample,Nyears)
    LF <- array(0, dim = dim(plb))
    rownames(LF) <- 1:Nyears
    for (y in 1:Nyears) {
      LF[y, ] <- rmultinom(n = 1,size = comp_sample[y],
                           prob = plb[y, ])
      
    }
    rowSums(LF) 
    colnames(LF) <- mids
    LF_df <- reshape2::melt(LF)
    names(LF_df) <- c("Year", "Length", "Numbers")
    
    LF_df <- LF_df %>%
      dplyr::group_by(Year, Length,Numbers)
    LF_df$Year <- factor(LF_df$Year)
    LF_df$Length <- as.numeric(LF_df$Length)
    bins <- as.numeric(colnames(LF))
    bw <- bins[1]
    years <- unique(as.numeric(LF_df$Year))[order(unique(as.numeric(LF_df$Year)))] #unique(LF_df$Year)[order(as.numeric(unique(LF_df$Year)))]
    nyears <- length(years)
    # p <- ggplot(data=LF_df, aes(x=Length, y=Numbers)) +
    #   geom_bar(stat="identity",fill="steelblue")+
    #   facet_wrap(Year~., ncol=5, dir="v") +
    #   xlab("Length bin (cm)") + ylab("Numbers")
    # p
    
    ### Plots
    # par(mfrow=c(4,2),mar=c(0,4,0,1),oma=c(4,1,1,1))
    # plot(years,TB_t/1000,type="l",lwd=3,xaxt="n",ylab="Total Biomass (t)")
    # plot(years,SPR_t,type="l",lwd=3,ylim=c(0,1),xaxt="n",ylab="SPR")
    # abline(h=0.4,lty=2)
    # plot(years,I_t,type="l",lwd=3,xaxt="n",ylab="Index of abundance")
    # plot(years,Cw_t/1000,type="l",lwd=3,xaxt="n",ylab="Catch (t)")
    # plot(years,F_t,type="l",lwd=3,ylab="Fishing mortality",xaxt="n")
    # 
    # head(LF_df)
    
    ### Mean length 
    data.by.fish<-data.frame(Year=NA,Length=NA)
    Sizes<-mids #mid 
    
    for(j in 1:length(Sizes)){ 
      for(i in 1:nrow(LF)){ 
        if(!is.na(LF[i,j])==TRUE & LF[i,j]>0){
          Year<-rep(years[i],LF[i,j])
          Length<-rep(Sizes[j],LF[i,j])
          Length<-Length+(runif(length(Length), min=-1, max=1))
          tmp<-data.frame(Year,Length)
          data.by.fish<-rbind(data.by.fish,tmp)
        }
      }
    } 
    data.by.fish<-data.by.fish[-1,]
    data.by.fish<-data.by.fish[order(data.by.fish$Year),]
   # head(data.by.fish)
    
    summary(data.by.fish)
    data.by.fish$Year<-as.factor(data.by.fish$Year)
    Mean.length<-aggregate(x = data.by.fish$Length,by=list(data.by.fish$Year),FUN=mean)
    Mean.length<-Mean.length$x
    # boxplot(data.by.fish$Length~data.by.fish$Year,ylab="Length (cm)",xaxt="n")
    # abline(h=90,lty=2)
    # mtext(text = "Year",side = 1,line = 2.5,outer = T)
    # 
    # Reference points
    FFmsy<-F_t/Fmsy
    BBmsy<-TB_t/TBmsy
    # plot(years,FFmsy,type="l",lwd=3,ylab="FFmsy",ylim=c(0,5))
    # abline(h=1,lty=2)
    # plot(years,BBmsy,type="l",lwd=3,ylab="BBmsy",ylim=c(0,3.2))
    # abline(h=1,lty=2)
    # 
    Dep<-TB_t/B0
    
    return(list(Outs=data.frame(Year=years,Cw_t,CatchDev,I_t,IndexDev,F_t,FishDev,TB_t,SB_t,Dep,SPR_t,RecDev,RecDev_AR,BBmsy,FFmsy,Mean.length),
           der.quants=c(Fmsy,TBmsy,B0,msy),LF=data.frame(LF),Bio_a=data.frame(Age=ages,Length=L_a,Weight=W_a,Sel=S_a,Mat=Mat_a),Bio_l=data.frame(Length=bins,Weight=W_l,Sel=S_l,Mat=Mat_l)))
}


