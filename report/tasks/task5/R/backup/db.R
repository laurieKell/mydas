library(DBI)
library(RPostgreSQL)
library(dplyr)
library(plyr)
library(reshape)

drv  = dbDriver("PostgreSQL")
con1 = dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                 dbname='FLRout',
                 port = 5432,
                 user = 'MydasAdmin',
                 password = 'datapoor1!')
##write
res$om="Turbot"
dbWriteTable(con1, "data_om", value = res, append=FALSE,overwrite=FALSE,row.names=FALSE)

##read
omtest = dbGetQuery(con1, paste0("SELECT *from data_om"))

scenarios=expand.grid(stock=c("brill","turbot","ray","pollack","sprat","lobster","razor"),
                      mp   =c("xsa","mpb","sra","lbspr","2/3","pid","sbt1","sbt2","irate"))
xsaOptions=expand.grid(mp="xsa",ftar=c(0.7,1,1,2),interval=1:3)
scenarios=merge(scenarios,xsaOptions,all=TRUE)

mse =mseXSA(om,eq,
            mp,control=xsaControl,
            ftar=1.0,
            interval=1,start=0,end=100,
            srDev=srDev,uDev=uDev)


runs=foreach(i=seq(dim(scenarios)[1]), 
               .combine=c,
               .multicombine=TRUE,
               .packages=c("FLCore","FLAssess","FLXSA","FLasher","FLRP","plyr")) %dopar%{
                 
                res=mse(scenario[i,"stock"],scenario[i,"mp"],..)
                
                writeDB(cbind(scenario=i,omSmry(res)))
                }
