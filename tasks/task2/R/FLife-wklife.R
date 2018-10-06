## ----knitr_init, echo=FALSE, results="hide"------------------------------
library(knitr)
## Global options
opts_chunk$set(cache     =!TRUE,
               echo      =TRUE,
               eval      =TRUE,
               prompt    =FALSE,
               comment   =NA,
               message   =FALSE,
               warning   =FALSE,
               tidy      =FALSE,
               fig.height=6,
               fig.width =8)

iFig=0

## ---- pkgs, message=FALSE------------------------------------------------
library(ggplot2)
library(FLife)
library(plyr)
library(reshape)

## ---- theme, echo=FALSE--------------------------------------------------
theme_set(theme_bw())
options(digits=3)

## ---- data---------------------------------------------------------------
data(wklife)

wklife

## ---- fig.height=8, echo=FALSE-------------------------------------------
library(GGally)

my_smooth <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
  geom_point(...,size=.5)+
  geom_smooth(...,method="lm",se=FALSE)}

my_density <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
  geom_density(...,lwd=1)}

ggpairs(transform(wklife[,-c(1:5)],linf=log(linf),k=log(k),l50=log(l50)),
  lower = list(continuous = wrap(my_smooth)),
  diag=list(continuous=wrap(my_density,alpha=0.2)),
  title = "")+
  theme(legend.position ="none",
  panel.grid.major =element_blank(),
  axis.ticks       =element_blank(),
  axis.text.x      =element_blank(),
  axis.text.y      =element_blank(),
  panel.border     =element_rect(linetype = 1, colour="black", fill=NA))

## ----FLPar---------------------------------------------------------------
wkpar=as(wklife[,6:13],"FLPar")
attributes(wkpar)[names(wklife)[1:5]]=wklife[,1:5]

## ----m-gislason----------------------------------------------------------
par=lhPar(wkpar)

## ----eqls----------------------------------------------------------------
library(FLBRP)

eql=lhEql(par)

## ----vectors-------------------------------------------------------------
ggplot(FLQuants(eql,"m","catch.sel","mat","catch.wt"))+
  geom_line(aes(age,data,col=attributes(wkpar)$name[iter]))+
  facet_wrap(~qname,scale="free")+
  scale_x_continuous(limits=c(0,15))+ 
  guides(colour=guide_legend(title="Species",title.position="top"))

## ----eql-----------------------------------------------------------------
plot(iter(eql,7))

## ----eql-lmsl------------------------------------------------------------
lmsl=as(iter(eql,7),"FLStock")

plot(lmsl)

