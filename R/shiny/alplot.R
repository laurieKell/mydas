
library(ggplot2)
library(GGally)
lhdat=read.csv("/Users/alextidd/Documents/functionalgroups/lifedata.csv")
lhdat$logwinf=log(lhdat$Winfinity)
lhdat$logK=log(lhdat$K)
lhdat$logLoo=log(lhdat$Loo)
p <- ggpairs(lhdat[c(1,12:14)], aes(color = Family, alpha=0.6))+ theme_bw()
for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + 
      scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07","blue")) +
      scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07","blue"))  
  }
}
p