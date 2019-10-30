library(plyr)
library(dplyr)
library(rgdal)
library(ggplot2)
library(DBI)
library(RPostgreSQL)
library(reshape)

options(scipen=999)
drv  = dbDriver("PostgreSQL")

con  = dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                 dbname = 'mydasDev',
                 port = 5432,
                 user = 'MydasApplication',
                 password = 'gmit2017!')
#survdat   = dbGetQuery(con, paste0("select * from data_surveybio_rich"))
survdat   = dbGetQuery(con, paste0("select * from data_surveybio"))
rects     = dbGetQuery(con, paste0("select * from data_icesrects"))
surdense  = survdat %>% group_by(speciesgp, year, ices_rectangle, ices_division) %>% 
                        summarise(cabund= sum(densbiom_kg_sqkm/1000)) #in tonnes per km2
                       
#per km square           

#number of hauls in square

surv_eff = dbGetQuery(con, paste0("SELECT *
                                      FROM data_surveystns
                                  "))
#surv_eff = dbGetQuery(con, paste0("SELECT *
#                                      FROM data_surveystns_rich
#                                  "))

#Calculate total area/total number of stations the survey covered by division and merge with catch tonnes perkm2
nstations = surv_eff %>% group_by(year, ices_rectangle,ices_division) %>%
            summarise(totn=length(unique(haulid))) %>% inner_join(rects)%>%
            group_by(year, ices_division)%>%summarise(totarea=sum(area_km2), totstn=sum(totn)) %>% inner_join(surdense)

#Biomass = totalsurveyarea*1/total number of stations within the survey area * sum of catch
nstations$relabund = ((nstations$totarea))*(1/nstations$totstn)*nstations$cabund
nstations2 = nstations %>%group_by(year, speciesgp) %>% summarise(totabund= sum(relabund), variance= var(relabund))

stecfdat  = dbGetQuery(con, paste0("select * from data_stecf_aer_cpuedays"))
#assuming cpue=qN

catchq    = stecfdat %>% inner_join(nstations2)  %>% mutate(q=(totctch/totdays)/totabund) %>%
            group_by(year, gear_type, speciesgp) %>% summarise(mn_q=mean(q),  y10=quantile(q, 0.10), y90=quantile(q,0.90), y25= quantile(q, 0.25), y75=quantile(q, 0.75)) 


ggplot(subset(catchq, gear_type %in% c("TBB","OTB","LLS","GND","GNS") & speciesgp %in% c("TUR","SKA","LIN", "POL","JOD")), aes(year,mn_q))+
  geom_ribbon(aes(x=year, ymin = y10, ymax = y90),fill = "#009E73", alpha = 0.35)+ 
  #geom_ribbon(aes(x = year, ymin = y25, ymax = y75), fill ="#009E73", alpha = 0.45)+
  geom_line()+
  facet_grid(speciesgp~gear_type, scale="free_y")+theme_bw()

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
###################Take this dataframe for q analysis #######################################
catchq2    = stecfdat %>% inner_join(nstations2)  %>% mutate(q=(totctch/totdays)/totabund) 

tgc <- summarySE(catchq2, measurevar="q", groupvars=c("year","gear_type","speciesgp"))
pd <- position_dodge(0.3) # move them .05 to the left and right
ggplot(subset(tgc, gear_type %in% c("TBB","OTB","LLS","GND","GNS") & speciesgp %in% c("TUR","SKA","LIN", "POL","JOD")), aes(x=year, y=q, colour=gear_type)) + 
  geom_errorbar(aes(ymin=q-se, ymax=q+se), width=.1, position=pd) +
  geom_line(size=1,position=pd) +
  geom_point(position=pd, size=3, aes(shape=gear_type), fill="white")+
  facet_wrap(~speciesgp, scale="free_y")+theme_bw()+   scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","black", "green"))

#################### histograms #####################
ggplot(subset(tgc, gear_type %in% c("TBB","OTB","LLS","GND","GNS") & speciesgp %in% c("TUR","SKA","LIN", "POL","JOD")))+geom_histogram(aes(log(q), fill=gear_type), bin=100)+facet_wrap(~speciesgp, scale="free_y")


test= nstations2 %>% 
  group_by(speciesgp) %>% 
  mutate(value=(totabund)/(totabund[1])) %>%
  arrange(speciesgp, year)

######################## survey abundance #################################
ggplot(subset(test, year %in% c(2008:2016) & speciesgp %in% c("TUR","SKA","LIN", "POL","JOD")), aes(year,totabund))+geom_line()+geom_point()+
  facet_wrap(~speciesgp, scale="free_y")+expand_limits(y=0)+theme_bw()+ylab("Abundance in tonnes")
  
  





tester = cast(subset(catchq,gear_type %in% c("TBB","OTB","LLS")), year+~speciesgp, mean, value = 'mn_q')
#tester = subset(tester, gear_type %in% c("TBB", "OTB"))

dist.mat <- dist(tester, method= "euclidean", diag=F, upper=F)

cluster1 <- hclust(dist.mat, method= "ward")
plot(cluster1,hang = -0.1, labels = FALSE )
category       = cutree(cluster1, k=3)

rect.hclust(cluster1, k=3, border="red")
defined = cbind(tester, category)

