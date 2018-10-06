load("data/shiny.RData")
names(priors)[1]="OM"
iters=unique(mseIts$iter)

library(plyr)
library(dplyr)
library(ggplot2)