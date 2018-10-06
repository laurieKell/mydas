library(dplyr)
library(reshape2)
library(ggplot2)
library(FLife)
library(DBI)
library(RPostgreSQL)
library(FLBRP)


drv  = dbDriver("PostgreSQL")

con  = dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
                 dbname='mydasDev',
                 port = 5432,
                 user = 'MydasApplication',
                 password = 'gmit2017!')

#brill
econ = dbGetQuery(con, paste0("SELECT *
                              FROM data_stecf_aer_econ
                              WHERE    year = 2015
                              "))
#BLL and TUR
econ$maingear = ifelse(econ$country_name %in% c ('United Kingdom', 'Ireland', 'Belgium') & econ$fishing_tech == 'TBB',"TBB",
                       ifelse(econ$country_name %in% c ('United Kingdom', 'Ireland', 'France') & econ$fishing_tech == 'DTS',"DTS",       
                              ifelse(econ$country_name %in% c ('United Kingdom', 'France') & econ$fishing_tech == 'DFN',"DFN", "OTH")))       

varcosts = econ %>%
  filter(variable_name %in% c('Energy costs','Wages and salaries of crew','Other variable costs'))%>%
  group_by(year,country_code, maingear, vessel_length)%>%summarise(varcosts=sum(value))

fixcosts =econ %>% 
  filter(variable_name %in% c('Repair & maintenance costs','Other non-variable costs','Annual depreciation costs','Rights costs'))%>%
  group_by(year,country_code, maingear, vessel_length)%>%summarise(fixcosts=sum(value))

days = econ %>% 
  filter(variable_name %in% 'Days at sea') %>% 
  group_by(year,country_code, maingear, vessel_length)%>%
  summarise(days=sum(value))

costs=inner_join(fixcosts, varcosts)
#f=qe assume q = 0.01 and effort = 1
f=0.005
costperday = inner_join(days, costs) %>% 
  mutate(varcost_per_day=varcosts/days/f, fixcost_per_day=fixcosts/days/f )%>%
  group_by(year) %>% filter(maingear %in% c("TBB","DTS","DFN"))%>%
  summarise(mnvarcost=mean(varcost_per_day), mnfixcost=mean(fixcost_per_day))

ct_ef   = dbGetQuery(con, paste0("SELECT * FROM data_stecf_aer_cpuedays where year = 2015"))

ct_ef2  = ct_ef %>% mutate(val=totval/totctch) %>% filter(speciesgp=="TUR") %>% 
  summarise( y10=quantile(val, 0.10), y20= quantile(val, 0.20), y30=quantile(val, 0.3), y40=quantile(val, 0.4), y50=quantile(val, 0.5), y60=quantile(val, 0.60), y70=quantile(val, 0.70), y80=quantile(val, 0.80),  y90=quantile(val,0.90), y100=quantile(val, 0.90))

#priceatage
price = melt(ct_ef2)

load(url("https://github.com//fishnets//fishnets//blob//master//data//fishbase-web//fishbase-web.RData?raw=True"))

lh=subset(fb,species=="Psetta maxima")

names(lh)[c(14:17)]=c("l50","l50min","l50max","a50")
lh=lh[,c("species","linf","k","t0","a","b","a50","l50","l50min","l50max")]

lh=apply(lh[,-1],2,mean,na.rm=T)

lh=FLPar(lh)
par=lhPar(lh)


eq=lhEql(par)
price(eq) = c(price$value,rep(price[10,2],31))*1000 #per tonne

price(eq)@units = "euro/tonnes"

vcost(eq) = costperday$mnvarcost
vcost(eq)@units = "euro/unit F"

fcost(eq) = costperday$mnfixcost
fcost(eq)@units = "euro/unit F"

brptur = brp(eq)
refpts(brptur)

plot(brptur)

dimnames(refpts(brptur))$refpt[5]="mey"
brptur=brp(brptur)
plot(brptur)

chmod 770 participants
