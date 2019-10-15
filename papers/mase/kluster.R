library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(tidyr)

set.seed(27)

centers <- tibble(
  cluster = factor(1:3), 
  num_points = c(100, 150, 50),  # number points in each cluster
  x1 = c(5, 0, -3),              # x1 coordinate of cluster center
  x2 = c(-1, 1, -2)              # x2 coordinate of cluster center
)

labelled_points <- centers %>%
  mutate(
    x1 = map2(num_points, x1, rnorm),
    x2 = map2(num_points, x2, rnorm)
  ) %>% 
  select(-num_points) %>% 
  unnest(x1, x2)

ggplot(labelled_points, aes(x1, x2, color = cluster)) +
  geom_point()

points <- labelled_points %>% 
  select(-cluster)

kclust <- kmeans(points, centers = 3)
kclust
summary(kclust)

kclusts <- tibble(k = 1:9) %>%
  mutate(
    kclust = map(k, ~kmeans(points, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, points)
  )

kclusts

clusters <- kclusts %>%
  unnest(tidied)

assignments <- kclusts %>% 
  unnest(augmented)

clusterings <- kclusts %>%
  unnest(glanced, .drop = TRUE)


p1 <- ggplot(assignments, aes(x1, x2)) +
  geom_point(aes(color = .cluster)) + 
  facet_wrap(~ k)
p1

p2 <- p1 + geom_point(data = clusters, size = 10, shape = "x")
p2

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line()

#########################
load("/home/laurence/Desktop/sea++/mydas/project/papers/mase/data/refcase.RData")

ssb=ddply(refcase$ts, .(f,k,M,s,bg,sel3,CV,AR,deviates), with,
          as.data.frame(spec.ar(ssb/mean(ssb),plot=FALSE)[c("freq","spec")]))
ssb=ddply(ssb,       .(f,k,M,s,bg,sel3,CV,AR,deviates), transform,
          spc=spec/max(spec),
          wavelen=1/freq)

ggplot(subset(ssb,freq<=0.2))+
  geom_line(aes(freq,spc,col=as.character(f),group=paste(AR,f,sel3),
                linetype=as.character(sel3)))+
  xlab("r")+ylab("")+
  facet_grid(M*s*CV~bg*k)


ssb.matrix=scale(cast(ssb,f+k+M+s+bg+sel3+CV+AR+deviates~wavelen,value="spec")[,-seq(9)])

points=aaply(ssb.matrix,1,scale)

kclusts <- tibble(k=4:20) %>%
  mutate(
    kclust = map(k, ~kmeans(points, .x, nstart=100, iter.max=100)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, points)
  )

clusters <- kclusts %>%
  unnest(tidied)
assignments <- kclusts %>% 
  unnest(augmented)
clusterings <- kclusts %>%
  unnest(glanced, .drop = TRUE)

ggplot(clusterings, aes(k, tot.withinss)) +
  geom_line()


p1 <- ggplot(assignments, aes(X1, X2)) +
  geom_point(aes(color = .cluster)) + 
  facet_wrap(~ k)
p1

p2 <- p1 + geom_point(data = clusters, size = 10, shape = "x")
p2


kclusts <- data.frame(k=4:25) %>% 
  group_by(k) %>% 
  do(kclust=kmeans(points.matrix, .$k, nstart=20,  iter.max=20))

clusters    <- kclusts %>% group_by(k) %>% do(tidy(.$kclust[[1]]))
assignments <- kclusts %>% group_by(k) %>% do(augment(.$kclust[[1]], points.matrix))
clusterings <- kclusts %>% group_by(k) %>% do(glance(.$kclust[[1]]))

ggplot(clusterings, aes(k, tot.withinss)) + 
  geom_line()+
  xlab("Number of Clusters")+ylab("Within SS")+
  theme_bw()

####################
load("/home/laurence/Desktop/sea++/mydas/project/papers/mase/data/refcase.RData")

ssb=ddply(refcase$ts, .(f,k,M,s,bg,sel3,CV,AR,deviates), with,
          as.data.frame(spec.ar(ssb/mean(ssb),plot=FALSE)[c("freq","spec")]))
ssb=ddply(ssb,       .(f,k,M,s,bg,sel3,CV,AR,deviates), transform,
          spc=spec/max(spec),
          wavelen=1/freq)

ggplot(subset(ssb,wavelen<=100))+
  geom_line(aes(wavelen,spc,col=as.character(f),group=paste(AR,f,sel3),
                linetype=as.character(sel3)))+
  xlab("Wave Length")+ylab("")+
  facet_grid(bg*CV*M~k*s)

ssb.df=cast(subset(ssb,wavelen<=100),
                    f+k+M+s+bg+sel3+CV+AR+deviates~wavelen,value="spc")
ssb.matrix=t(apply(ssb.df[,-seq(9)],1,scale))
key=ssb.df[,1:9]

kclusts <- data.frame(ncluster=4:20) %>% 
  group_by(ncluster) %>% 
  do(kclust=kmeans(ssb.matrix, .$ncluster, nstart=20,  iter.max=20))

clusters   <-kclusts %>% group_by(ncluster) %>% do(tidy(   .$kclust[[1]]))
assignments<-kclusts %>% group_by(ncluster) %>% do(augment(.$kclust[[1]], key))
clusterings<-kclusts %>% group_by(ncluster) %>% do(glance( .$kclust[[1]]))

ggplot(clusterings, aes(ncluster, tot.withinss)) + 
  geom_line()+
  xlab("Number of Clusters")+ylab("Within SS")+
  theme_bw()




kclusts <- tibble(k=4:20) %>%
  mutate(
    kclust = map(k, ~kmeans(ssb.matrix, .x, nstart=20, iter.max=20)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, ssb.matrix))

dat=cbind(ssb,"cluster"=rep(unlist(c(subset(assignments,ncluster==12)[,".cluster"])),each=500))

ggplot(subset(dt2,ncluster%in%seq(4,12,2)), aes(wavelen, spec))+  
  geom_point()+
  facet_grid(cluster~ncluster,scale="free")+ 
  theme_bw()+
  theme(legend.position="none")+
  xlab("Wave Length")+ylab("")
