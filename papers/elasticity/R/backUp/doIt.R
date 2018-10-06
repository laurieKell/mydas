
doIt=function(what,par,dynamic=FALSE,fbar=FLQuant(c(seq(0,.75,length.out=21),seq(.75,.0,length.out=21)[-1]),dimnames=list(year=1:41))){
  
  func=function(x,dmns,par,what) {
   
    par[dmns] =exp(x)
    par["t0"] =-par["t0"]
    par["M2"] =-par["M2"]
     
    res=lh(par,model="bevholt")
             
    fbar(res)=fbar
    res      =brp(res)
    rp       =refpts(res)
    if (dynamic) res=fwd(res)

    smy=c(c(rp[what,"ssb"]),
          c(rp[what,"biomass"]),
          c(rp[what,"harvest"]),
          c(rp[what,"yield"]),
          c(rp[what,"rec"]),
          c(ssb(    res)),
          c(stock(  res)),
          c(fbar(   res)),
          c(catch(  res)),
          c(rec(    res)),
          c(ssb(    res)/rp[what,"ssb"]),
          c(stock(  res)/rp[what,"biomass"]),
          c(fbar(   res)/rp[what,"harvest"]),
          c(catch(  res)/rp[what,"yield"]),
          c(rec(    res)/rp[what,"rec"]))

    return(log(smy))
    }

# hand cranked version
# func returns already logged values
#browser()
#smallmult <- 1 + 1e-4
#test <- func(log(c(par)),dimnames(par)$params,par,"msy")
#par2 <- par
#par2[1,] <- par2[1,] * smallmult
#test2 <- func(log(c(par2)),dimnames(par2)$params,par2,"msy")
#elas <- (test2[106:125] - test[106:125]) / (log(100*smallmult) - log(100))
## matches output of jacobian
#
#browser()
    jbn=jacobian(func,log(c(par)),dmns=dimnames(par)$params,par=par,what=what)
                                         
    res=data.frame(Year    =c(rep(1,5),rep(seq(dims(fbar)$year),10)),
                   Quantity=c(c("SSB","Biomass","Harvest","Yield","Recruits"),rep(rep(c("SSB","Biomass","Harvest","Yield","Recruits"),each=dims(fbar)$year),2)),
                   Type    =c(rep("Reference Point",5),rep(c("Absolute","Relative"),each=dims(fbar)$year*5)))
  
    res=cbind(res,jbn)
    res=melt(res,id.var=c("Year","Quantity","Type"))
    res$Parameter      =factor(dimnames(par)$params[res$variable])
        
    p.=data.frame(Parameter=c("linf",  "t0",    "M1","M2","s",  "vb", "a",     "b",     "bg",      "k",       "ato95",      "sl",         "sr",      "a50",     "asym",       "a1",         "fec"),
                  Process  =c("Growth","Growth","M", "M", "SRR","SRR","Growth","Growth","Maturity","Growth",  "Maturity","Selectivity","Selectivity","Maturity","Maturity","Selectivity","Maturity"))
      
    res=merge(res,p.)
    
    return(res[,c("Year","Quantity","Type","Parameter","Process","value")])}
