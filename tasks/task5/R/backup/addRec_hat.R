library(FLBRP)
library(RPostgreSQL)
library(DBI)

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/turbot.RData")
srr=cbind("spp"="turbot",model.frame(params(eq)))

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/brill.RData")
srr=rbind(srr,cbind("spp"="brill",model.frame(params(eq))))

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/ray.RData")
srr=rbind(srr,cbind("spp"="ray",model.frame(params(eq))))

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/pollack.RData")
srr=rbind(srr,cbind("spp"="pollack",model.frame(params(eq))))

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/sprat.RData")
srr=rbind(srr,cbind("spp"="sprat",model.frame(params(eq))))

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/razor.RData")
srr=rbind(srr,cbind("spp"="razor",model.frame(params(eq))))

load("/home/laurence/Desktop/sea++/mydas/tasks/task4/data/lobster.RData")
srr=rbind(srr,cbind("spp"="lobster",model.frame(params(eq))))

save(srr,file="/home/laurence/Desktop/sea++/mydas/tasks/task5/data/srr.RData",compress="xz")

drv =dbDriver("PostgreSQL")
con1=dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
               dbname='FLRout',
               port = 5432,
               user = 'MydasAdmin',
               password = 'datapoor1!')

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/sbt1.RData")
names(sbt1)[c(2,9)]=c("spp","stock")
sbt1=merge(sbt1,srr,by=c("spp","iter"))
sbt1=transform(sbt1,rec_hat=a*ssb/(b+ssb))
sbt1=sbt1[,c(1:23,28,24:25)]
save(sbt1,file="/home/laurence/Desktop/sea++/mydas/tasks/task5/data/sbt1.RData",compress="xz")
dbWriteTable(con1, "sbt1", sbt1, append=!TRUE, overwrite=TRUE)

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/sbt2.RData")
names(sbt2)[c(2,8)]=c("spp","stock")
sbt2=merge(sbt2,srr,by=c("spp","iter"))
sbt2=transform(sbt2,rec_hat=a*ssb/(b+ssb))
sbt2=sbt2[,c(1:22,27,23:24)]
save(sbt2,file="/home/laurence/Desktop/sea++/mydas/tasks/task5/data/sbt2.RData",compress="xz")
dbWriteTable(con1, "sbt2", sbt2, append=!TRUE, overwrite=TRUE)
dbDisconnect(con1)

xsa=dbGetQuery(con1, "SELECT * from xsa")
xsa=xsa[,c(21,1:20)]

nms=c("sprat","pollack","brill","turbot","lobster","razor")
names(nms)=c("spr","pol","bll","tur","lbe","raz")
xsa$spp=as.character(nms[xsa$spp])

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/srr.RData")
srr$spp=as.character(srr$spp)

xsa=merge(xsa,srr)
xsa=transform(xsa,rec_hat=a*ssb/(b+ssb))
dbWriteTable(con1, "xsa", xsa, append=!TRUE, overwrite=TRUE)

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/mpb.RData")
names(mpb)[c(2,8)]=c("spp","stock")
dbWriteTable(con1, "mpb", mpb, append=!TRUE, overwrite=TRUE)

load("/home/laurence/Desktop/sea++/mydas/tasks/task5/data/srampb.RData")
names(srampb)[c(2,8)]=c("spp","stock")
dbWriteTable(con1, "srampb", srampb, append=!TRUE, overwrite=TRUE)


