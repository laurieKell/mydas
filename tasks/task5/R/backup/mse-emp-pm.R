
empD$spp =ac(empD$spp)
empD     =subset(empD,year>=mseStart[spp]&year<=mseStart[spp]+45)
empD$year=empD$year-mseStart[empD$spp]

empd_pm=ddply(empD,.(spp,iter,j), smryStat)

dt=transform(empd_pm,k1=cut(k1,seq(0,1,0.03)),
             k2=cut(k2,seq(0,1,0.03)))

dt=ddply(dt,.(k1,k2,spp),with,data.frame(
  safe =mean(safety),
  kobe =mean(kobe.n/45),
  yield=mean(yield),
  aav  =1-mean(yieldAav)))

ggplot(dt)+
  geom_tile(aes(k1,k2,fill=kobe))+
  scale_fill_gradientn(colours=c("navy","blue","cyan","lightcyan","yellow","red","red4"))+
  facet_wrap(.~spp)
