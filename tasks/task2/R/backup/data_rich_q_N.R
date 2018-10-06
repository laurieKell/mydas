library(plyr)
library(dplyr)
library(rgdal)
library(ggplot2)
library(DBI)
library(RPostgreSQL)
library(reshape)

#stock assessments
load("/home/laurence/Desktop/sea++/mydas/data/stkTs.RData")

stkTs = rename(stkTs, c(stock="abundance"))
stkTs$stock = ifelse(stkTs$.id %in% c("codns"), "nseacod",
              ifelse(stkTs$.id %in% c("hadns"), "nseahad", 
              ifelse(stkTs$.id %in% c("whgns"), "nseawhg",      
              ifelse(stkTs$.id %in% c("plens"), "nseaple", 
              ifelse(stkTs$.id %in% c("solns"), "nseasol",  
              ifelse(stkTs$.id %in% c("sains"), "nseapok", 
              ifelse(stkTs$.id %in% c("codcs"), "cseacod",
              ifelse(stkTs$.id %in% c("hadcs"), "cseahad",  
              ifelse(stkTs$.id %in% c("whgcs"), "cseawhg","none")))))))))       

nms=c("codns","hadns","whgns","plens","solns","sains","codcs","hadcs","whgcs")
names(nms)=c("nseacod","nseahad","nseawhg","nseaple","nseasol","cseapok","cseacod","cseahad","cseawhg")

options(scipen=999)
drv  = dbDriver("PostgreSQL")

con  = dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                 dbname = 'mydasDev',
                 port = 5432,
                 user = 'MydasApplication',
                 password = 'gmit2017!')
#survdat   = dbGetQuery(con, paste0("select * from data_surveybio_rich"))
survdat   = dbGetQuery(con, paste0("select * from data_surveybio_rich"))
rects     = dbGetQuery(con, paste0("select * from data_icesrects"))
surdense  = survdat %>% group_by(speciessciname, year, ices_rectangle, ices_division) %>% 
  summarise(cabund= sum(densbiom_kg_sqkm/1000)) #in tonnes per km2

surv_eff = dbGetQuery(con, paste0("SELECT *
                                  FROM data_surveystns_rich
                                  "))

#Calculate total area/total number of stations the survey covered by division and merge with catch tonnes perkm2
nstations = surv_eff %>% group_by(year, ices_rectangle,ices_division) %>%
  summarise(totn=length(unique(haulid))) %>% inner_join(rects)%>%
  group_by(year, ices_division)%>%summarise(totarea=sum(area_km2), totstn=sum(totn)) %>% inner_join(surdense)

nstations$stock=ifelse(nstations$speciessciname %in% c("Gadus morhua") & nstations$ices_division %in% c("4A","4B","4C", "7D", "3A", "4"), "nseacod",
                ifelse(nstations$speciessciname %in% c("Melanogrammus aeglefinus") & nstations$ices_division %in% c("4A","4B","4C","3A","4"), "nseahad",        
                ifelse(nstations$speciessciname %in% c("Merlangius merlangus") & nstations$ices_division %in% c("4A","4B","4C","4","7D"), "nseawhg",        
                ifelse(nstations$speciessciname %in% c("Pollachius virens") & nstations$ices_division %in% c("4A","4B","4C","4","6A","6B","3A"), "nseapok",        
                ifelse(nstations$speciessciname %in% c("Pleuronectes platessa") & nstations$ices_division %in% c("4A","4B","4C","4", "3A"), "nseaple",        
                ifelse(nstations$speciessciname %in% c('Solea solea') & nstations$ices_division %in% c("4A","4B","4C","4"), "nseasol",
                ifelse(nstations$speciessciname %in% c("Gadus morhua") & nstations$ices_division %in% c( "7E", "7F","7H", "7G","7J", "7K"), "cseacod",
                ifelse(nstations$speciessciname %in% c("Melanogrammus aeglefinus") & nstations$ices_division %in% c("7B","7C", "7D","7E", "7F","7H", "7G","7J","7K"), "cseahad",
                ifelse(nstations$speciessciname %in% c("Merlangius merlangus") & nstations$ices_division %in% c("7B","7C","7H", "7E", "7F","7H", "7G","7J","7K"), "cseawhg", 0)))))))))


nstations$scale = 1
#scale the data
nstations$scale2 = ifelse(nstations$stock %in% c("nseawhg"), 1, 
                   ifelse(nstations$stock %in% c("nseacod"), 0.45,     
                   ifelse(nstations$stock %in% c("nseahad"), 0.5,   
                   ifelse(nstations$stock %in% c("nseasol"), 0.25, 
                   ifelse(nstations$stock %in% c("nseapok"), 1.63,        
                   ifelse(nstations$stock %in% c("nseaple"), 1.95, nstations$scale))))))
#change scale here
nstations$relabund = ((nstations$totarea/nstations$scale2))*(1/nstations$totstn)*nstations$cabund

nstations2 = nstations %>%group_by(year, speciessciname, stock) %>% summarise(totabund= sum(relabund), variance= var(relabund))

stecfdat  = dbGetQuery(con, paste0("select * from data_stecf_aer_cpuedays_rich"))
#assuming cpue=qN


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
catchq3 = catchq2 %>% 
  group_by(year, country_code, stock, gear_type, vessel_length) %>%
  summarise (c = sum(totctch)) %>% group_by(year, stock, country_code) %>%
  mutate(prop = c /sum(c)) %>% filter(prop > 0.10) %>%inner_join(catchq2)

tgc <- summarySE(catchq3, measurevar="q", groupvars=c("year","gear_type","stock"))

#################### histograms #####################
#histogram q from resulting estimated abundance
ggplot(tgc)+
  geom_histogram(aes(log(q), fill=gear_type), bins=20)+facet_wrap(~stock, scale="free_y")+theme_bw()+ 
  theme(legend.position="none")

test= nstations2 %>% 
  group_by(stock) %>% 
  mutate(value=(totabund)/(totabund[1])) %>%
  arrange(stock, year)

######################## survey abundance #################################
ggplot(test, aes(year,totabund))+
  geom_line()+geom_point()+
  facet_wrap(~stock, scale="free_y")+
  expand_limits(y=0)+theme_bw()+ylab("Abundance in tonnes")

#commercial
stocks = inner_join(stkTs, nstations2) 
stockhad = subset(stocks, stock %in% "nseahad")
stockhad$abundance = stockhad$abundance*10000
stocks = subset(stocks, !(stock %in% "nseahad"))
stocks = rbind(stocks, stockhad)

stocks= rename(stocks, c(totabund = "survB", abundance="assB"))  %>% filter(!stock  %in% c("cseahad", "cseawhg", "cseacod") )

biom   = melt(stocks, id.vars=c("year", "stock"), measure.vars=c("survB","assB"))
ggplot(biom, aes(year, value, group=variable, colour=variable))+ 
  geom_point()+ 
  geom_line()+
  facet_wrap(~stock, scale = "free_y")

#estimate q from assessment abundance by gear

catchq3_1    = stecfdat %>% inner_join(stocks)  %>% mutate(q=(totctch/totdays)/assB)
catchq4 = catchq3_1 %>% 
  group_by(year, country_code, stock, gear_type, vessel_length) %>%
  summarise (c = sum(totctch)) %>% group_by(year, stock, country_code) %>%
  mutate(prop = c /sum(c)) %>% filter(prop > 0.10)  %>%inner_join(catchq2)

tgc2 <- summarySE(catchq4, measurevar="q", groupvars=c("year","gear_type","stock"))
ggplot(tgc2)+geom_histogram(aes(log(q), fill=gear_type), bins=20)+facet_wrap(~stock, scale="free_y")+theme_bw()+ theme(legend.position="none")
