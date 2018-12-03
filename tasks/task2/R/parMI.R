parMI=list(
"brill"  =FLPar(c(linf= 41.775,k=0.44, t0=-0.93,a=0.02225,b=2.92,a50=NA,  l50=25.0), units="NA"),
"turbot" =FLPar(c(linf= 59.1,  k=0.28, t0=-0.4, a=0.01111,b=3.15,a50=4.0, l50=43.25),units="NA"),
"pollack"=FLPar(c(linf= 84.6,  k=0.19, t0=-0.94,a=0.0107, b=2.98,a50=4.0, l50=NA),   units="NA"),
"ray"    =FLPar(c(linf=118.25, k=0.142,t0=-1.07,a=0.00279,b=3.23,a50=5.88,l50=68.84),units="NA"))
      
    
 load("/home/laurence/Desktop/Dropbox/mydasOMs/data/brill.RData")
 parFB=FLPar(brill=lh)
 load("/home/laurence/Desktop/Dropbox/mydasOMs/data/turbot.RData")
 parFB=FLPar(turbot=lh)
 load("/home/laurence/Desktop/Dropbox/mydasOMs/data/pollack.RData")
 parFB=FLPar(pollack=lh)
 load("/home/laurence/Desktop/Dropbox/mydasOMs/data/ray.RData")
 parFB=FLPar(ray=lh)
 