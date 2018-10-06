xsaTest<-function(om,pg=10,ctrl=FLXSA.control()){
  
  stk=setPlusGroup(om,pg)
  idx=FLIndex(index=stock.n(stk))
  range(idx)[c("plusgroup","startf","endf")]=c(pg,0.01,.1)
  
  stk+FLXSA(stk,idx,control=ctrl,diag.flag=FALSE)}