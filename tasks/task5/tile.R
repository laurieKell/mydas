mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
scale_fill_gradientn(colours = mycol)+
  
dt1=merge(empd_pm,model.frame(control),by="iter")

dt2=ddply(dt1,.(k1,k2,spp),transform,
          safe =mean(safety)
          kobe =mean(kobe.n/45),
          yield=mean(yield),
          aav  =1-mean(yieldAav),
          k1.=cut(k1,seq(0,1,0.05)),
          k2.=cut(k2,seq(0,1,0.05)))


geom_contour(color = "white", alpha = 0.5) +
  
ggplot(dt2)+
  geom_tile(aes(k1,k2,fill=z))+
  scale_fill_gradientn(colours = mycol)+
  facet_wrap(.~spp)
