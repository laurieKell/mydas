#https://magesblog.com/post/2015-01-06-kalman-filter-example-visualised-with-r/
#https://gist.github.com/mages/52cf15bf8c9563c3518f

library(FLCore)
library(kobe)

library(mvtnorm)

xhat=c(1.5, 0.7)
cor =matrix(c(1.0, -0.9, 
             -0.9,  1.0), ncol=2)
cov=cor2cov(cor,xhat*0.1)

prior=transform(expand.grid(x=xhat[1]*seq(0,2,length.out=101),
                          y=xhat[2]*seq(0,2,length.out=101)),
              z=dmvnorm(cbind(x,y),mean=xhat,sigma=cov))

bias=c(-0.5,0.5)
sensor=transform(expand.grid(x=(xhat[1]+bias[1])*seq(0,4,length.out=101),
                            y=(xhat[2]+bias[2])*seq(0,4,length.out=101)),
              z=dmvnorm(cbind(x,y),mean=xhat+bias,sigma=cor2cov(cor,(xhat)*0.1)))

xhatf=xhat+cov%*% t(I) %*% solve(I %*% cov %*% t(I) + R) %*% (bias)
posterior=transform(expand.grid(x=(xhatf[1])*seq(0,4,length.out=101),
                           y=(xhatf[2])*seq(0,4,length.out=101)),
                  z=dmvnorm(cbind(x,y),mean=xhatf,sigma=cor2cov(cor,(xhatf)*0.1)))

A <- matrix(c(1.2, 0,
              0, -0.2), ncol=2)
Q <- 0.3 * Sigma
K <- A %% cov%% t(I) %% solve(I%% cov %% t(I) + R)
xhatnew <- A %% xhat + K %% (y - I %% xhat)
Sigmanew <- A %% cov %% t(A) - K %% I %% cov %*% t(A) + Q
z4 <- outer(prior$x,prior$y,fn, mean=c(xhatnew), varcov=Sigmanew)

dat=rbind(cbind("What"="Prior",    prior),
          cbind("What"="Sensor",   sensor),
          cbind("What"="Posterior",posterior)) 
kobePhase(xlim=c(0,3),ylim=c(0,2))+
  geom_contour(aes(x=x, y=y,z=z, colour=stat(level)),data=dat)+
  facet_wrap(~What,ncol=2)
