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

# MSY calculations
calc_msy <- function(F, ages, M, R0, W_a, S_a){
  Nage <- calc_equil_abund(ages=ages, M=M, F=F, R0=R0, S_a=S_a)
  F_a<-F*S_a
  C_a<- Nage * F_a/(F_a+M) * (1-exp(-(F_a+M))) #Baranov equation 
  YPR <- sum(C_a*W_a) # Yield per recruit
  
  return(YPR)}

Nyears=20
Fdynamics=c(0.00, 0.25, 0.50, 0.75, 1.00, 
            1.25, 1.50, 1.75, 2.00, 2.00, 
            2.00, 2.00, 2.00, 2.00, 2.00,
            1.75, 1.50, 1.25, 1.00, 0.75)
SigmaR=0.5
rho=0
SigmaF=0.2
SigmaI=0.2
SigmaC=0.1
R0=1000
AgeMax = 15
start_ages=0
binwidth=2
linf=122
vbk=0.209
t0=(-1.338)
lwa=1.34e-05
lwb=3.1066
CVlen=0.1
S50=85
S95=90
M50=90
M95=100
M=0.3
h=0.7
qcoef=0.00001
Nl_sample=1000
  

FishDev <-  c(1,rnorm(Nyears-1, mean = -(SigmaF ^ 2) / 2, sd = SigmaF))


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

yld=mdply(data.frame(F=seq(0,10,0.1)), function(F) 
          calc_msy(F, ages, M, R0, W_a, S_a))
ggplot(yld)+
  geom_line(aes(F,V1))+
  geom_point(aes(F,V1),data=yld[max(yld$V1)==yld$V1,])


  