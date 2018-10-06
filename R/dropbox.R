library(FLCore)
library(ggplotFL)

library(rdrop2)
token<-drop_auth()

saveRDS(token, "/home/laurence/Dropbox/mydas/token.RDS")

drop_get(path='mydas/sprat.RData', local_file = "sprat.RData")

