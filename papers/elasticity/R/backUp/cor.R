lf=FLPar(rlnorm(100,log(FLPar(linf=300)),0.2),dimnames=list(params="linf",iter=1:100))
pars=gislasim(lf)
plotmatrix(data.frame(t(pars[c("linf","k","M1","a50"),drop=T])))
