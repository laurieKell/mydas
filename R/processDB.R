library(ggplot2)
library(plyr)
library(reshape)

library(RPostgreSQL)
library(DBI)

library(FLCore)
library(kobe)
library(diags)

mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)

drv  =dbDriver("PostgreSQL")
conAl=dbConnect(drv, 
               host    = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
               dbname  ='FLRout',
               port    =5432,
               user    ='MydasAdmin',
               password='datapoor1!')
dbListTables(conAl)

conLK=dbConnect(drv, 
               host    ='wklife.csrzweaa3tbm.eu-west-2.rds.amazonaws.com',
               dbname  ='wklife',
               port    =5432,
               user    ='mydas',
               password='Yes_Meski')
dbListTables(conLK)

om  =dbGetQuery(conAl,"SELECT * from om")
xsa =dbGetQuery(conAl,"SELECT * from xsanew")
mpb =dbGetQuery(conAl,"SELECT * from mpb")
sbt1=dbGetQuery(conAl,"SELECT * from sbt1")
sbt2=dbGetQuery(conAl,"SELECT * from sbt2")

sbt1=sbt1[,-1]
sbt2=sbt2[,-1]
mpb =mpb[ ,-c(1:2)]
xsa =xsa[ ,-1]

sbt1=subset(sbt1,year>=51&gamma!=0.75)
sbt2=subset(sbt2,year>=51)
mpb =subset(mpb, year>=51)
xsa =subset(xsa, year>=51)

sbt1=subset(sbt1,year<mseStart[ac(spp)]+36)
sbt2=subset(sbt2,year<mseStart[ac(spp)]+36)
mpb =subset(mpb, year<mseStart[ac(spp)]+36)
xsa =subset(xsa, year<mseStart[ac(spp)]+36)

ggplot(ddply(sbt1,.(spp,year), with, data.frame(ssb=median(ssb/msy_ssb))))+
  geom_line(aes(year,ssb,group=spp))

ggplot(ddply(sbt2,.(spp,year), with, data.frame(ssb=median(ssb/msy_ssb))))+
  geom_line(aes(year,ssb,group=spp))

ggplot(ddply(mpb,.(spp,year), with, data.frame(ssb=median(ssb/msy_ssb))))+
  geom_line(aes(year,ssb,group=spp))

ggplot(ddply(xsa,.(spp,year), with, data.frame(ssb=median(ssb/msy_ssb))))+
  geom_line(aes(year,ssb,group=spp))

sbt1PM=ddply(sbt1,.(spp,iter,k1,k2,gamma), smryStat)
sbt2PM=ddply(sbt2,.(spp,iter,k1,k2),       smryStat)
mpbPM =ddply(mpb, .(spp,iter,ftar,btrig),  smryStat)
xsaPM =ddply(xsa, .(spp,iter),             smryStat)

ggplot(melt(subset(sbt1PM,spp=="pollack"),id=c("spp","iter","k1","k2","gamma")))+
  geom_boxplot(aes(ac(k2),as.numeric(value)),outlier.size=0.2)+
  facet_grid(variable~ac(paste(k1,gamma,sep=";")),scale="free_y",space="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Gain Term When Stock Increasing")+ylab("")+
  theme_bw()
b1=melt(subset(sbt1PM,spp=="pollack"&k2==2&k1==1.5&gamma==1))

ggplot(melt(subset(sbt2PM,spp=="pollack"),id=c("spp","iter","k1","k2")))+
  geom_boxplot(aes(ac(k2),as.numeric(value)),outlier.size=0.2)+
  facet_grid(variable~ac(paste(k1)),scale="free_y",space="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Gain Term When Stock Increasing")+ylab("")+
  theme_bw()
b2=melt(subset(sbt2PM,spp=="pollack"&k2==0.5&k1==0.25))

ggplot(melt(subset(mpbPM,spp=="pollack"),id=c("spp","iter","ftar","btrig")))+
  geom_boxplot(aes(ac(ftar),as.numeric(value)),outlier.size=0.2)+
  facet_grid(variable~ac(paste(btrig)),scale="free_y",space="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("F times FMSY")+ylab("")+
  theme_bw()
b3=melt(subset(mpbPM,spp=="pollack"&ftar==0.7))

ggplot(melt(subset(xsaPM,spp=="pollack"),
            id=c("spp","iter")))+
  geom_boxplot(aes(1,as.numeric(value)),outlier.size=0.2)+
  facet_grid(variable~.,scale="free_y",space="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("")+theme_bw()

compare=rbind.fill(cbind(MP="Age Based HCR",  melt(subset(xsaPM,spp=="pollack"),id=c("spp","iter"))),
                   cbind(MP="Biomass HCR",    b3),
                   cbind(MP="SBT Relative",   b2),
                   cbind(MP="SBT Trend",      b1))
compare=subset(compare,!(variable%in%c("k1","k2","gamma","ftar","btrig")))

ggplot(subset(compare,MP%in%c("Age Based HCR","Biomass HCR")))+
  geom_boxplot(aes(y=as.numeric(value)),outlier.size=0.2)+
  facet_grid(variable~MP,scale="free_y",space="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("MP")+ylab("")+theme_bw()

ggplot(compare)+
  geom_boxplot(aes(y=as.numeric(value)),outlier.size=0.2)+
  facet_grid(variable~MP,scale="free_y",space="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("MP")+ylab("")+theme_bw()

dbWriteTable(conLK, "mydas_empd_pm", value=sbt1PM, append=FALSE,overwrite=TRUE,row.names=FALSE)
dbWriteTable(conLK, "mydas_empp_pm", value=sbt2PM, append=FALSE,overwrite=TRUE,row.names=FALSE)
dbWriteTable(conLK, "mydas_mpb_pm",  value=mpbPM,  append=FALSE,overwrite=TRUE,row.names=FALSE)
dbWriteTable(conLK, "mydas_xsa_pm",  value=xsaPM,  append=FALSE,overwrite=TRUE,row.names=FALSE)

dbWriteTable(conLK, "mydas_empd",   value=sbt1,   append=FALSE,overwrite=TRUE,row.names=FALSE)
dbWriteTable(conLK, "mydas_empp",   value=sbt2,   append=FALSE,overwrite=TRUE,row.names=FALSE)
dbWriteTable(conLK, "mydas_mpb",    value=mpb,    append=FALSE,overwrite=TRUE,row.names=FALSE)
dbWriteTable(conLK, "mydas_xsa",    value=xsa,    append=FALSE,overwrite=TRUE,row.names=FALSE)


