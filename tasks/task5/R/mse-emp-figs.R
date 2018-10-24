library(ggpmisc)

formula <- y ~ poly(x, 2, raw = TRUE)

### Derivative #################################################################
load("/home/laurence/Desktop/Dropbox/mydasOMs/empd-results.RData")

empd_pm=transform(empd_pm,kobe.p=kobe.n/45,yieldAav=pmin(0.5,yieldAav))

dat=melt(empd_pm[,c("safety","kobe.year","kobe.p","yield",
    "yieldAav","spp","k1","k2","gamma")], 
         id=c("spp","k1","k2","gamma"))
colnames(dat)[5:6] = c("objective","objectval")
dat=melt(dat,id=c("spp","objective","objectval"))
colnames(dat)[4:5] = c("param","paramval")

dat=subset(dat,objective!="kobe.year"&param!="gamma")
ggplot(aes(objectval,paramval, col=spp),data=dat)+
  geom_point(size=.05)+
  geom_smooth(se=FALSE,span=1)+
  #stat_poly_eq(formula = formula, 
  #             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)+
  facet_grid(param~objective,scale="free")+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Performance Measure")+ylab("Hyperparameter")


### Derivative
load("/home/laurence/Desktop/Dropbox/mydasOMs/empd-results.RData")

empd_pm=transform(empd_pm,kobe.p=kobe.n/45,yieldAav=pmin(0.5,yieldAav))

dat=melt(empd_pm[,c("safety","kobe.year","kobe.p","yield",
                    "yieldAav","spp","k1","k2","gamma")], 
         id=c("spp","k1","k2","gamma"))
colnames(dat)[5:6] = c("objective","objectval")
dat=melt(dat,id=c("spp","objective","objectval"))
colnames(dat)[4:5] = c("param","paramval")

dat=subset(dat,objective!="kobe.year"&param!="gamma")
ggplot(aes(objectval,paramval, col=spp),data=dat)+
  geom_point(size=.05)+
  geom_smooth(se=FALSE,span=1)+
  #stat_poly_eq(formula = formula, 
  #             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)+
  facet_grid(param~objective,scale="free")+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Performance Measure")+ylab("Hyperparameter")


### Proportional ###############################################################
load("/home/laurence/Desktop/Dropbox/mydasOMs/empp-results.RData")

empp_pm=transform(empp_pm,kobe.p=kobe.n/45,yieldAav=pmin(0.5,yieldAav))

dat=melt(empp_pm[,c("safety","kobe.year","kobe.p","yield",
                    "yieldAav","spp","k1","k2")], 
         id=c("spp","k1","k2"))
colnames(dat)[4:5] = c("objective","objectval")
dat=melt(dat,id=c("spp","objective","objectval"))
colnames(dat)[4:5] = c("param","paramval")

dat=subset(dat,objective!="kobe.year")
ggplot(aes(objectval,paramval, col=spp),data=dat)+
  geom_point(size=.05)+
  geom_smooth(se=FALSE,span=1)+
  #stat_poly_eq(formula = formula, 
  #             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)+
  facet_grid(param~objective,scale="free")+
  theme_bw()+
  theme(legend.position="bottom")+
  xlab("Performance Measure")+ylab("Hyperparameter")

For empd, the trends are the same for the stocks, but values of safety vary

for k1

effective parameter values are from 1 to 4
k1= 4 gives high safety and probability of being in the green kobe quadrant (kope.p)
K1=2 gives high yield and low AAV

for k2

effective values are around 2
this has a big effect on safety and kobe.p but not yield or AAV

For empp

for k1

This has little effect, maninly as the stock never declines

for k2

low values of k1 give high safety and kobe.p
yield is less affected than AAV

Thinking of the next steps, how about? 

Compare with XSA, by adding vertical lines for median XSA results for safety, ...
Contour safety, kobe.p, yield and AAV by k1 and k2,  to see where the maximum lies, then compare this to xsa
create a ulility function by combining the 4 performance measures and then contour by k1 and k2