library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)

library(FLCore)
library(ggplotFL)
library(FLasher)
library(FLBRP)
library(FLife)
library(mpb)

library(DBI)
library(RPostgreSQL)

library(doParallel)
library(foreach)

dirMy ="/home/laurence/Desktop/sea++/mydas/tasks/task4"
dirDat=file.path(dirMy,"data")

drv =dbDriver("PostgreSQL")

registerDoParallel(4)

## parallel
res<-foreach(i=(seq(dim(scen)[1])[1:4]), 
              .combine=c,
              .multicombine=TRUE,
              .packages=c("FLCore","FLasher","FLBRP","FLife","plyr","reshape","DBI","RPostgreSQL")) %dopar%{

    con1=dbConnect(drv, host    ='postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                        dbname  ='FLRout',
                        port    =5432,
                        user    ='MydasAdmin',
                        password='datapoor1!')
          
    ##run
    
    dbWriteTable(con1,"empsbt2rev",value=rs2,append=TRUE,overwrite=FALSE,row.names=FALSE)
    }
            
dbGetQuery(con1, "SELECT * FROM emp")
















