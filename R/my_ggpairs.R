my_density <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_density(...,lwd=1)}

my_smooth <- function(data,mapping,...){
  ggplot(data=data,mapping=mapping)+
    geom_smooth(...,method="lm",se=FALSE)+
    geom_point(...)}
