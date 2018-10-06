####################################
## read in life history data      ##
## Three sets                     ##
## i)   ERA                       ##
## ii)  scombridbrid              ##
## iii) scombridbrid fecundity    ##
## iv)  Teleost                   ##
## v)   Herring                   ##
####################################

library(FLCore)
library(plyr)

dirInp="/home/laurence/Desktop/sea++/mydas/papers/FLife/inputs"
dirDat="/home/laurence/Desktop/sea++/mydas/papers/FLife/data"

## i) ERA ##
pel=read.table(file=file.path(dirInp,"era.csv"),row.names=NULL,sep=",")
nms=tolower(names(pel))
pel=pel[,-60]
names(pel)=nms[-1]
pel=pel[,-c(1,9,11:17,22:23,25:30)]
pel$ocean=ifelse(pel$ocean=="IND","Indian","Atlantic")

save(pel,file=file.path(dirDat,"pel.RData"),compress="xz")

## iii) scombrid
scombrid=read.csv(file=file.path(dirInp,"scombrid.csv"))

names(scombrid)=tolower(names(scombrid))
fecEqn     =scombrid[!is.na(scombrid[,"fecundity_equation"]),c("sc_name","fecundity_equation")]
scombrid       =scombrid[,c(1:9,12:14,16:19,23:25,27,29:33,42,46,47,48,49)]

nm1=c("sc_name",         "c_name",          "code",            "taxo_group",     
      "tribe",           "genera",          "stock_name",      "ocean",          
      "ocean_region",    "area",            "lat",             "lon",      
      "startyear_study", "finalyear_study", "total_size",      "sex",             
      "linf",            "k",               "to",              "tmax",            
      "maturity_size",   "lm",              "lm50",            "tm",              
      "tm50",            "spawning_season", "wl_a",            "wl_b",            
      "wl_units",        "source")                     

nm2=c("name",            "common",          "code",            "taxon",     
      "tribe",           "genus",           "stock",           "ocean",          
      "region",          "area",            "lat",             "lon",      
      "minyear",         "maxyear",         "n",               "sex",             
      "linf",            "k",               "t0",              "amax",            
      "nmat",            "lmat",            "l50",             "amat",             
      "a50",             "spwn",            "a",               "b",               
      "units",           "cite")                     

names(scombrid)=nm2

save(scombrid,file=file.path(dirDat,"scombrid.RData"),compress="xz")


## iii)  scombridbrid fecundity       
# Length in cm
# Weight in Kg
# Batch fecundity in millions

scombrid$spwn=gsub("([]])", " ", scombrid$spwn)
scombrid$spwn=gsub("([\\[])", " ", scombrid$spwn)
scombrid$spwn=gsub(',', "-", scombrid$spwn)
scombrid$spwn=gsub('-{2}', "-", scombrid$spwn)
scombrid$spwn=gsub('-{2}', "-", scombrid$spwn)

fec<-function(length,params)
            params["a"]+params["b"]*(params["c"]*L)^params["d"]

#                                               BF 0.545+0.0372*OFBW (ovary free body weight in kg)
# BF -948860+246.964*length(cm);                BF= -18800 +0.14448 *W (Kg)
# BF -509775+144.792*length(cm);                BF=8663+0.09485*weight(Kg)
# BF -265749+98.667*length(cm);                 BF=60212+0.077.51*weight(Kg)
# BF 6184.5*length(cm)+21191;                   BF=59.657*GonadFreeweight(g)+198509; BF=5499.7*OW(OvaryWeight,g)-98360
# BF 22304*length(cm)\x8a_\xea661391;           BF = 86.3*SW(somatic weight,g)+225940
# BF (in thousands)=0.0013*length(cm)^3.0029;   BF (in thousands)=35.973*BW(g)^1.3417; BF (in thousands)=1.7369*GW(g)^0.9622
# BF(in thousands)=0.0000007*length(cm)^4.7119; BF(in thousands)=15.667*weight(g)^1.5141; BF(in thousands)=0.4885*GonadWeight(g)^1.2525;
# BF 0.00000003*length(cm*0.1)^4.88
# BF 5991.3 length(cm)-146389;                  BF=82.012*weight(g) + 15786
# BF 0.00054*length(cm*0.1)^3.05
# BF 0.0011* length(cm*0.1)^2.896;              BF=31087*WW(whole weight, kg)^1.384
# BF  -159197.73+61.017*length(cm);             BF= -8210.81+160.33* SW (somatic weight, g)
#                                               BF 56.7weight(g)-0.000571
#                                               BF (in millions) = 0.00334 * MG (fresh gonad mass in g) \x8a_\xea 0.323
# BF 0.2934*length(cm)^3.2673;                  BF=62173*weight(kg) + 225310
#                                               BF 163048*weight(kg)-5062591
# BF 1.1015^-8*length(cm*0.1)^4.679
# BF 0.0003747*length(cm)^3.180758;             BF=150400 + 62941* weight(kg)
# BF (in millions)=4.46x 10^-9*length(cm)^4.16; BF(in millions)=1.18x10^-2*weight(kg)^1.43
# BF (in millions)=4.78242X10^-17*L(cm)^7.530
# BF 0.0058*length(cm)^3.994 (off Java);        BF=0.0018* length(cm)^4.175 (off Hawaii)
# BF 8.815x 10^4* length^4.419;                 BF= 6.153x 10^3* weight^1.543;
# BF (in millions)= 3.2393*10^5 length(cm) \x8a_\xea 5.2057*10^7

val=c(NA,     NA,     NA,
      -948860,246.964,  1,
      -509775,144.792,  1,
      -265749,98.667,   1,
      21191, 6184.5,   1,
      22304, 661391,   1,
      0,     1.3,      3.0029,
      0,     0.0007,   4.7119,
      0,     3.315e-8, 4.88,
      -146389,5991.3,   1,
      0,     4.812755e-7, 3.05,
      0,     1.397632e-6, 2.896,              
      -159197.73,61.017, 1,             
      NA,NA,NA,
      NA,NA,NA,
      0,  0.2934,3.2673,                 
      NA, NA, NA,
      0, 9.663267e-06,4.679,
      0, 0.0003747,3.180758,             
      0,4.46e-3,4.16,
      0, 4.78242e-17,7.530,
      0, 0.0058,3.994,        
      0,8.815e4,4.419,                
      3.2393e11, 5.2057e7, 1)

id=c("20","192","193","194","257" ,"278" ,"291","292","561","656" ,"868","955","1034","1105",
     "1253","1336","1337","1351","1353","1355","1420","1498","1511","1521")

fec=FLPar(array(val,dim=c(3,length(id)),
      dimnames=list(params=c("a","b","c"),id=id)))

attributes(fec)$name=scombrid[dimnames(fec)$id,"name"]

save(fec,file=file.path(dirDat,"fec.RData"),compress="xz")

## iv)  Teleost 

load("/home/laurence/Desktop/papers/elasticity/data/fishbase-web.RData")

names(fb)[14]="l50"

lhPar=names(fb)[c(9:11,14,20:21)]

#get reference set where all data present 
good=fb[apply(fb[,lhPar],1,function(x) !any(is.na(x))),
        names(fb)[c(1:11,14,20,21,22:33)]]

#select families 
good=subset(good,family%in% 
              c("Clupeidae","Engraulidae","Gadidae","Scombridae",
                "Merlucciidae","Sciaenidae","Sparidae","Lutjanidae",
                "Mullidae","Carangidae","Mugilidae","Sebastidae",
                "Serranidae","Scophthalmidae","Pleuronectidae","Soleidae",
                "Cyprinidae","Salmonidae"))

#average lh params by species
good=ddply(good,.(species),with, data.frame("linf"=mean(linf),
                                            "k"   =mean(k),
                                            "t0"  =mean(t0),
                                            "l50" =mean(l50),
                                            "a"   =mean(a),
                                            "b"   =mean(b)))
#taxa stuff
taxa=subset(fb,species%in%good$species)
taxa=taxa[,names(fb)[c(1,2,5,6,13,25,32)]]
taxa=taxa[!duplicated(taxa$species)&taxa$species%in%good$species,]
taxa=transform(taxa,species=ac(species))
good=transform(good,species=ac(species))
full=merge(good,taxa)

full$t0=ifelse(full$t0==0,-0.1,full$t0)

teleost=full

save(teleost,file=file.path(dirDat,"teleost.RData"),compress="xz")

## v)  Herring
load("/home/laurence/Desktop/papers/elasticity/data/fishbase-web.RData")
names(fb)[14]="l50"
nms=names(fb)[c(9:11,14,20:21)]

herring=subset(fb,species==sort(ac(unique(fb$species)))[381])[,nms]

save(herring,file=file.path(dirDat,"herring.RData"),compress="xz")

wklife<-read.csv("/home/laurence/Desktop/flr/FLife/mystuff/datasets/inputs/wklife.csv")
save(wklife,file=file.path(dirDat,"wklife.RData"),compress="xz")

rays=data.frame(t(array(c(
  "Amblyraja radiata",1.8195,-0.6383,
  "Amblyraja radiata",1.8506,-1.0000,
  "Amblyraja radiata",2.0792,-0.8861,
  "Amblyraja radiata",2.1038,-0.9586,
  "Dipturus batis",2.4048,-1.2441,
  "Dipturus innominatus",2.1790,-1.0223,
  "Dipturus pullopunctatus",1.8871,-1.3010,
  "Dipturus pullopunctatus",2.1239,-1.0969,
  "Leucoraja erinacea",1.7168,-0.5376,
  "Leucoraja erinacea",1.7218,-0.4559,
  "Leucoraja naevus",1.8762,-0.7959,
  "Leucoraja naevus",1.9619,-0.9666,
  "Leucoraja ocellata",2.0569,-0.8416,
  "Leucoraja ocellata",2.0864,-1.1308,
  "Leucoraja ocellata",2.1367,-1.2291,
  "Leucoraja wallacei",1.6253,-0.5850,
  "Okamejei hollandi",1.6580,0.1106,
  "Okamejei hollandi",1.6990,0.1271,
  "Okamejei hollandi",1.7559,-0.2832,
  "Okamejei hollandi",1.7987,-0.3893,
  "Raja brachyura",2.0719,-0.7212,
  "Raja brachyura",2.1430,-0.9208,
  "Raja clavata",1.9325,-0.6778,
  "Raja clavata",2.0043,-0.7447,
  "Raja clavata",2.0043,-0.6676,
  "Raja clavata",2.0212,-0.6676,
  "Raja clavata",2.0253,-0.8697,
  "Raja clavata",2.0569,-0.6882,
  "Raja clavata",2.0671,-0.9747,
  "Raja clavata",2.0719,-0.8539,
  "Raja clavata",2.0719,-0.7959,
  "Raja clavata",2.1021,-1.0088,
  "Raja clavata",2.1072,-1.0458,
  "Raja clavata",2.1206,-1.0969,
  "Raja clavata",2.1430,-1.0458,
  "Raja clavata",2.1461,-1.0315,
  "Raja clavata",2.1553,-1.0458,
  "Raja eglanteria",1.8949,-0.6990,
  "Raja microocellata",2.1367,-1.1308,
  "Raja miraletus",1.9440,-0.7144,
  "Raja miraletus",1.9633,-0.7696,
  "Raja montagui",1.8370,-0.7212,
  "Raja montagui",1.8597,-0.5171,
  "Raja montagui",1.8621,-0.7447,
  "Raja montagui",1.8943,-0.5918,
  "Raja montagui",1.8987,-0.6778,
  "Raja montagui",1.9903,-0.8297,
  "Raja rhina",1.9854,-0.6021,
  "Raja rhina",2.0294,-0.7959,
  "Raja undulata",2.0374,-0.9586,
  "Raja undulata",2.0492,-1.0000,
  "Zearaja nasuta",1.9605,-0.7959),c(3,53)))) 

rays=data.frame(spp=as.character(rays[,1]),linf=10^as.numeric(as.character(rays[,2])),
                k   =10^as.numeric(as.character(rays[,3])))

mylo=data.frame(t(array(c(
  "Myliobatis californica",2.0000,-0.6402,
  "Myliobatis californica",2.2014,-1.0022),c(3,2)))) 

mylo=data.frame(spp=as.character(mylo[,1]),linf=10^as.numeric(as.character(mylo[,2])),
                k   =10^as.numeric(as.character(mylo[,3])))
rays=rbind(rays,mylo)

save(rays,file="/home/laurence/Desktop/papers/elasticity/data/rays.RData")


load("/home/laurence/Desktop/papers/generic/fishnets-master/data/fishbase-web/fishbase-web.RData")

names(fb)[c(14,17)]=c("l50","a50")
fb=fb[,c("species","linf","k","t0","a","b","a50","l50")]

ray    =subset(fb,species=="Raja clavata")
sprat  =subset(fb,species=="Sprattus sprattus sprattus")
bigeye =subset(fb,species=="Thunnus obesus")
turbot =subset(fb,species=="Scophthalmus rhombus")
brill  =subset(fb,species=="Psetta maxima")
pollack=subset(fb,species%in%c("Pollachius pollachius","Pollachius virens"))  
plaice=subset(fb,species=="Pleuronectes platessa")

saury   ="Cololabis saira"
mackerel=c("Scomber scombrus","Scomber colias","Scomber japonicus","Scomber australasicus","Scomber indicus")
herring =c("Clupea harengus harengus","Strangomera bentincki","Clupea harengus pallasii") 
anchovy =c("Engraulis encrasicolus","Engraulis anchoita","Engraulis mordax","Engraulis japonicus","Engraulis ringens","Engraulis capensis")
sardine =c("Sardina pilchardus","Sardinops sagax","Sardinops melanostictus","Sardinops caeruleus","Sardinops ocellatus","Sardinella lemuru","Sardinella brasiliensis","Sardinella zunasi","Sardinella longiceps","Sardinella gibbosa","Sardinella aurita","Sardinella maderensis","Dussumieria acuta")
tuna    =c("Thunnus alalunga","Thunnus maccoyii","Thunnus orientalis","Thunnus thynnus","Thunnus obesus","Thunnus albacares","Katsuwonus pelamis")
