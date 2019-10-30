library(RPostgreSQL)
library(DBI)

drv =dbDriver("PostgreSQL")
con1=dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                    dbname='FLRout',
                    port = 5432,
                    user = 'MydasAdmin',
                    password = 'datapoor1!')

dbListTables(con1)
# "albnsbt1"  "albnsbt2" "albnwg"  "albnpella"

# only get 1st 10 to demonstrate
alb1=dbGetQuery(con1,"SELECT * from albnsbt1")
alb2=dbGetQuery(con1,"SELECT * from albnsbt2")
alb3=dbGetQuery(con1,"SELECT * from albnpella")
alb4=dbGetQuery(con1,"SELECT * from albnrobust")

# eval of right hand side of formula 
#eval(mod[[3]], c(as(par, 'list'), list(ssb=ssb(psr))))
# predict(predictModel)
# predict(predictModel(model=mod, params=par), ssb=ssb(psr))

dat=transmute(subset(alb1,scen==1),
          scen=paste(k1,k2,refYr),    
          iter=iter,
          year=year,
          bmsy=ssb/msy_ssb,
          fmsy=fbar/msy_harvest,
          msy =catch/msy_yield,
          rec =rec/virgin_rec)


library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)
library(FLAssess)
library(FLXSA)
library(mpb)
library(DBI)
library(RPostgreSQL)

drv  = dbDriver("PostgreSQL")
con1 = dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                 dbname='FLRout',
                 port = 5432,
                 user = 'MydasAdmin',
                 password = 'datapoor1!')

mpb=dbGetQuery(con1, paste("select* from mpb"))

ss_labels <- c("0.7"="Ftar=0.7", "1"="Ftar=1")
ggplot(subset(mpb, spp=="turbot" & year < 94),aes(as.factor(year), ssb/msy_ssb,fill=as.factor(btrig)))+
  geom_boxplot(outlier.size=0.1, position=position_dodge(1),width=0.8, lwd=0.05, notch=TRUE)+
  stat_summary(fun.y=mean, geom="line", aes(group=1))+
  geom_hline(aes(yintercept=1), size=0.75, colour= "red", linetype="dashed")+ 
  facet_wrap(~ftar,ncol=1, labeller = labeller(ftar = ss_labels),scale="free_y") + theme_bw() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position="bottom") + 
  scale_colour_manual(values=c("white", "#56B4E9"), labels=c("0.5","0.6"))  + 
  scale_fill_manual(name=expression(B[lim]),values=c("white", "#56B4E9"), labels=c("0.5","0.6")) + 
  scale_x_discrete(breaks = c(50,60,70,80,90))+
  scale_y_continuous(breaks = c(0:5))+
  xlab("year")+ylab(expression(B/B[MSY]))

ggplot(subset(mpb, spp=="turbot" & year < 94),aes(as.factor(year), catch/msy_yield,fill=as.factor(btrig)))+
  geom_boxplot(outlier.size=0.1, position=position_dodge(1),width=0.8, lwd=0.05, notch=TRUE)+
  stat_summary(fun.y=mean, geom="line", aes(group=1))+
  geom_hline(aes(yintercept=1), size=0.75, colour= "red", linetype="dashed")+ 
  facet_wrap(~ftar,ncol=1, labeller = labeller(ftar = ss_labels),scale="free_y") + theme_bw() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position="bottom") + 
  scale_colour_manual(values=c("white", "#56B4E9"), labels=c("0.5","0.6"))  + 
  scale_fill_manual(name=expression(B[lim]),values=c("white", "#56B4E9"), labels=c("0.5","0.6")) + 
  scale_x_discrete(breaks = c(50,60,70,80,90))+
  scale_y_continuous(breaks = c(0:5))+
  xlab("year")+ylab(expression(catch/catch[MSY]))

ggplot(subset(mpb, spp=="turbot" & year < 94),aes(as.factor(year), fbar/msy_harvest,fill=as.factor(btrig)))+
  geom_boxplot(outlier.size=0.1, position=position_dodge(1),width=0.8, lwd=0.05, notch=TRUE)+
  stat_summary(fun.y=mean, geom="line", aes(group=1))+
  geom_hline(aes(yintercept=1), size=0.75, colour= "red", linetype="dashed")+ 
  facet_wrap(~ftar,ncol=1, labeller = labeller(ftar = ss_labels),scale="free_y") + theme_bw() +
  theme(panel.grid.major = element_blank(),
        text = element_text(size=14),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position="bottom") + 
  scale_colour_manual(values=c("white", "#56B4E9"), labels=c("0.5","0.6"))  + 
  scale_fill_manual(name=expression(B[lim]),values=c("white", "#56B4E9"), labels=c("0.5","0.6")) + 
  scale_x_discrete(breaks = c(50,60,70,80,90))+
  scale_y_continuous(breaks = c(0:5))+
  xlab("year")+ylab(expression(f/f[MSY]))




