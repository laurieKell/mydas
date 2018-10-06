library(DBI)
library(RPostgreSQL)
library(dplyr)
library(plyr)
library(reshape)
library(maptools)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
options(scipen = 999)


drv  = dbDriver("PostgreSQL")

con  = dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                 dbname='mydasDev',
                 port = 5432,
                 user = 'MydasApplication',
                 password = 'gmit2017!')

stecf            = dbGetQuery(con,"SELECT * FROM data_stecflandings")
land2016         = subset(stecf, year %in% c(2016))
#remove area 4bc etc
land2016$flag    = ifelse(land2016$species %in% "LIN" & land2016$latitude <57.5 & land2016$area %in% "4", 1, 0)
land2016         = subset(land2016, flag==0)
land2016         = subset(land2016, !(area %in% "6B RFMO"))
area             = dbGetQuery(con,"SELECT * FROM div_area")

#ices division isnt broken down by division in North sea
land2016$division = ifelse(land2016$area %in% "4", "4A", land2016$division)
#convert in order to join and sum with landings
area$division     = ifelse(area$division %in% c("12A","12B","12C"), "12",
                    ifelse(area$division %in% c("14A","14B"), "14", area$division))
#areas not in euro zone so remove
area              = subset(area, !(area_27 %in% c("8.d.1","7.k.1","7.c.1","6.b.1","7.j.1")))

totarea           = ddply(area, .(division), summarise, totareakm=sum(area_km2))

allareas          = inner_join(land2016, totarea)

allland           = ddply(allareas, .(country, gear, mesh, stock, speciesgp, length), summarise, totland=sum(landings))
allland           = subset(allland, totland >0)
alllandarea       = ddply(allareas, .(country, gear, mesh, stock, speciesgp, length,  ices_rectangle), summarise, uniarea=unique(area_km2))
allareatot        = ddply(alllandarea, .(country, gear, mesh, stock, speciesgp, length), summarise, totfleetarea=sum(uniarea))
combi             = inner_join(allland, allareatot)

divarea           = ddply(allareas, .(stock, division), summarise, area=unique(totareakm))
allstockarea      = ddply(divarea, .(stock), summarise, stockarea=sum(area))

overlap           = inner_join(combi, allstockarea)
overlap$olap_percent   = (overlap$totfleetarea/overlap$stockarea)*100

overlap$price     =  ifelse(overlap$speciesgp %in% c("BLL"), 5.62,
                     ifelse(overlap$speciesgp %in% c("TUR"), 8.54,   
                     ifelse(overlap$speciesgp %in% c("LIN"), 1.22,  
                     ifelse(overlap$speciesgp %in% c("JOD"), 4.38,
                     ifelse(overlap$speciesgp %in% c("SKA"), 1.22,
                     ifelse(overlap$speciesgp %in% c("SPR"), 0.20,
                     ifelse(overlap$speciesgp %in% c("POL"), 1.76,
                     ifelse(overlap$speciesgp %in% c("POK"), 1.10,
                     ifelse(overlap$speciesgp %in% c("GUG"), 1.45,
                     0)))))))))
#management categories - NEED to check, 6 = completely datapoor, 3 has LPUE time series etc for TAC
#overlap$category  =  ifelse(overlap$speciesgp %in% c("BLL"), 3, 
#                     ifelse(overlap$speciesgp %in% c("POL"), 4,                         
#                     ifelse(overlap$speciesgp %in% c("GUG"), 6,                             
#                     ifelse(overlap$speciesgp %in% c("SPR"), 3, NA))))
#horizontal overlap categories 3 highest 1 lowest
overlap$score_olap = ifelse(overlap$olap_percent >30, 3, 
                     ifelse(overlap$olap_percent >10 & overlap$olap_percent <30, 2,
                     ifelse(overlap$olap_percent <10, 1,  0)))
#price scoring ategories 3 highest 1 lowest
overlap$score_price = ifelse(overlap$price >1.25, 3, 
                     ifelse(overlap$price >0.8 & overlap$price <1.25, 2,
                     ifelse(overlap$price <0.8, 1,  0)))

#catchability groupings 3 high, 2 medium 1 low
overlap$score_catch = ifelse(overlap$gear %in% c("BEAM") & overlap$speciesgp %in% c("BLL","TUR","GUG","SKA"), 3,
                     ifelse(overlap$gear %in% c("OTTER") & overlap$speciesgp %in% c("BLL","TUR", "SKA", "JOD","LIN","POK","POL"), 3,   
                     ifelse(overlap$gear %in% c("OTTER") & overlap$speciesgp %in% c("GUG"),2,
                     ifelse(overlap$gear %in% c("GILL") & overlap$speciesgp %in% c("POK","POL"), 3, 
                     ifelse(overlap$gear %in% c("GILL") & overlap$speciesgp %in% c("LIN", "TUR"), 2,        
                     ifelse(overlap$gear %in% c("LONGLINE") & overlap$speciesgp %in% c("LIN"), 3,  
                     ifelse(overlap$gear %in% c("LONGLINE") & overlap$speciesgp %in% c("POL"), 2,  
                     ifelse(overlap$gear %in% c("PEL_TRAWL") & overlap$speciesgp %in% c("SPR"), 3,   
                     ifelse(overlap$gear %in% c("GILL") & overlap$speciesgp %in% c("POK","POL"), 3, 1)))))))))  


#Determination of susceptibility scores, adopted from Hobday et al. (2011)
#Evidence of post-capture release and survival =1, discarded but survivorship unknown =2, majority dead or retained =3
overlap$score_postc = ifelse(overlap$speciesgp %in% c("SKA", "POK"), 1,
                     ifelse(overlap$speciesgp %in% c("BLL","TUR","POL", "LIN", "JOD"), 2,  3))
                     
#calculate susceptibility

overlap$S          = (((overlap$score_postc*overlap$score_olap*overlap$score_price*overlap$score_catch)-1)/40)+1
#create metier

overplot           = ddply(overlap, .(gear, stock), summarise, S=mean(S))

ggplot(overplot,aes(x=gear,y=S,fill=factor(stock)))+
   geom_bar(stat="identity",position="dodge")+theme_bw()+scale_fill_brewer(palette="Spectral")+  coord_cartesian(ylim = c(1, 3)) +
   theme( text = element_text(size=16), strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16),
   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Susceptibility") + guides(fill=guide_legend("stock")) 


###################### productivity analysis #########################
#write.csv(suscept, file="/Volumes/Alex/Mydas/annual_catch/tbl_era_scores.csv")

# 3 low productivity 1 high productivity
#< 5 years 5-15 years > 15 years
overlap$tm_score    = ifelse(overlap$speciesgp %in% c("GUG", "POL", "POK", "TUR", "BLL", "SPR", "JOD"), 1,  
                      ifelse(overlap$speciesgp %in% c("SKA","LIN"), 2, 3))      
#> 20,000 eggs per year 100 – 20,000 eggs per year < 100 eggs per year                  
overlap$fec_score   = ifelse(overlap$speciesgp %in% c("SKA"), 3, 
                      ifelse(overlap$speciesgp %in% c("SPR"), 2, 1))
#Broadcast spawner Demersal egg layer Live bearer
overlap$repro_score = ifelse(overlap$speciesgp %in% c("SKA"), 2, 1)
#< 2.75 2.75 – 3.25 > 3.25
overlap$troph_score = ifelse(overlap$speciesgp %in% c("SPR"), 2, 3)
#< 40 cm 40-200 cm > 200 cm
overlap$lmat_score  = ifelse(overlap$speciesgp %in% c("SPR", "BLL", "JOD", "GUG"), 1, 2)
#< 100 cm 100-300 cm > 300 cm
overlap$linf_score  = ifelse(overlap$speciesgp %in% c("SPR", "BLL", "JOD", "GUG", "TUR", "POL"), 1, 2)
#calculate productivity http://www.montereybayaquarium.org/-/m/C3EE8C68DA2A47B18A64BE6DBA72F76F.pdf
overlap$P           = ((overlap$tm_score+overlap$fec_score+overlap$repro_score+overlap$troph_score+overlap$lmat_score+overlap$linf_score)/6)

overlap$V           = sqrt(overlap$P^2 + overlap$S^2)    
#remove species barely caught by the gears, in this case less than 1 ton

overlap$speciesgp   = tolower(overlap$speciesgp)

#create contour matrix

grid = seq(1,3,length.out=100)
a=expand.grid(grid,grid)
b = matrix(sqrt((a[,1]-1)^2 + (a[,2]-1)^2),100)
riskcontour = data.frame(cbind(a,as.vector(b)))
names(riskcontour) = c("P","S","R")
head(riskcontour)

#colorRampPalette(brewer.pal(6, "Paired"))(100))
mycolour = ggplot(riskcontour) + geom_raster(aes(x=P,y=S,fill=R)) + scale_fill_gradientn("Risk",colours=colorRampPalette(brewer.pal(9,"Blues"))(100)) + xlab("Productivity") + ylab("Susceptibility") + 
  theme(
    text=element_text(size=14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA,colour = "black",size=2)
  )

summary=ddply(overlap, .(stock, gear), summarise, meanp = mean(P), 
                                         means = mean(S), 
                                         sds   = sd(S)
)


summary[is.na(summary)] = 0

summary$ymin = summary$means - summary$sds
summary$ymax = summary$means + summary$sds
summary      = subset(summary, !(gear %in% "NONE"))
# set desired dodge width
pd = position_dodge(width = 0.2)

mycolour + stat_contour(aes(x = P, y = S, z = R),linetype="dotted",colour="darkslategrey")+ 
         geom_point(data=summary, aes(meanp,means, colour=stock), position = pd)  + geom_errorbar(data = summary, aes(meanp,ymin = ymin, ymax =ymax, colour=stock), size=1, position = pd)+ 
  scale_color_manual(values=c("white","#000000", "#E69F00", "#56B4E9", "#009E73",  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  facet_wrap(~gear, nrow=2)


