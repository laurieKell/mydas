library(plyr)
library(dplyr)

plotMSE=function(var,mseIle,mseIts,prjIle,OM,p,ftar,btrig,blim){
  
  mseIle=mseIle[mseIle$variable==var   &
                mseIle$OM      ==OM    &
                mseIle$p       ==p     &
                mseIle$ftar    ==ftar  &
                mseIle$btrig   ==btrig &
                mseIle$blim    ==blim,]
         
  mseIts=mseIts[mseIts$variable==var   &
                mseIts$OM      ==OM    &
                mseIts$p       ==p     &
                mseIts$ftar    ==ftar  &
                mseIts$btrig   ==btrig &
                mseIts$blim    ==blim,]
  
  prjIle=prjIle[prjIle$variable==var   &
                prjIle$OM      ==OM    &
                prjIle$ftar    ==ftar  &
                prjIle$btrig   ==btrig &
                prjIle$blim    ==blim,]
  print(dim(mseIle))
  print(dim(mseIts))
  print(dim(prjIle))
  
  dat=rbind.fill(cbind(what="1",mseIle),
                 cbind(what="2",prjIle))
  
  ggplot(dat)+
    geom_ribbon(aes(x=year,ymin=ci.lower,ymax=ci.upper,fill=what), alpha=0.2)+ 
    geom_line(  aes(x=year,median,col=what))                                          +
    geom_line(  aes(x=year,value, col=iter),data=mseIts,size=0.3)                     +
    theme(legend.position="none")+xlab("Year")+ylab(var)+
    scale_y_continuous(limits=c(0,max(mseIle[,"ci.upper"])*1.5))
}

#with(scenarios,plotMSE(var,mseIle,mseIts,prjIle,OM[1],p[1],ftar[1],btrig[1],blim[1]))
  
plotMP=function(var,mp,priors,OM,p,ftar,btrig,blim,iters){
  
  df=mp[mp$p.x  ==p     &
        mp$OM   ==OM    &
        mp$ftar ==ftar  &
        mp$btrig==btrig &
        mp$blim ==blim,]
  
  names(df)[names(df)==paste(var,".x",sep="")]="var"
  dat=ddply(df,.(year), with, quantile(var, probs=c(.25,.5,.75)))
  names(dat)[2:4]=c("p25","p50","p75")
  
  names(priors)[tolower(names(priors))==paste(var,sep="")]="var"
  
  ggplot(dat)+
    geom_ribbon(aes(x=year,ymin=p25,ymax=p75), alpha=0.2,fill="red") +
    geom_line(  aes(x=year,p50), alpha=0.4,col="red") +
    geom_line(  aes(x=year,var,group=iter,col=iter),data=subset(df,iter %in% iters)) +
    theme(legend.position="none")+xlab("Year")+ylab(var)+
    scale_y_continuous(limits=c(0,max(dat[,4])*2))+
    geom_hline(aes(yintercept=var),data=subset(priors,OM==scenarios$OM[1]))}


plotKobe=function(mse,prj,OM,p,ftar,btrig,blim,yr){
  
  df1=mse[mse$p    ==p     &
          mse$OM   ==OM    &
          mse$ftar ==ftar  &
          mse$btrig==btrig &
          mse$blim ==blim,]
  dat1=ddply(df1,.(year), with, 
             data.frame(harvest=quantile(harvestRel,prob=c(0.5)),
                        stock  =quantile(ssbRel,    prob=c(0.5))))
  
  df2=prj[prj$OM   ==OM    &
          prj$ftar ==ftar  &
          prj$btrig==btrig &
          prj$blim ==blim,]
  
  dat2=ddply(df2,.(year), with, 
             data.frame(harvest=quantile(harvestRel,prob=c(0.5)),
                        stock  =quantile(ssbRel,    prob=c(0.5))))
  
  kobePhase(subset(dat1,year<=yr))+
    geom_path( aes(stock,harvest)) +
    geom_path( aes(stock,harvest),data=subset(dat2,year<=yr),col="blue") +
    geom_point(aes(ssbRel,harvestRel),data=subset(df1,year==yr))+
    geom_point(aes(ssbRel,harvestRel),data=subset(df2,year==yr),col="blue")}

plotKobeMar=function(mse,prj,OM,p,ftar,btrig,blim,yr){
  
  df1=mse[mse$p    ==p     &
          mse$OM   ==OM    &
          mse$ftar ==ftar  &
          mse$btrig==btrig &
          mse$blim ==blim,]
  dat1=ddply(df1,.(year), with, 
             data.frame(harvest=quantile(harvestRel,prob=c(0.5)),
                        stock  =quantile(ssbRel,    prob=c(0.5))))
  
  df2=prj[prj$OM   ==OM    &
          prj$ftar ==ftar  &
          prj$btrig==btrig &
          prj$blim ==blim,]
  dat2=ddply(df2,.(year), with, 
             data.frame(harvest=quantile(harvestRel,prob=c(0.5)),
                        stock  =quantile(ssbRel,    prob=c(0.5))))
  
  pts =rbind.fill(cbind(run="MSE",       subset(df1,year==yr)[,c("ssbRel","harvestRel")]),
                  cbind(run="Projection",subset(df2,year==yr)[,c("ssbRel","harvestRel")]))
  names(pts)[2:3]=c("stock","harvest")
  trks=rbind.fill(data.frame(run="MSE",       subset(dat1,year<=yr)),
                  data.frame(run="Projection",subset(dat2,year<=yr)))
  
  kobePhaseMar(pts,trks)}

#with(scenarios,plotMSE(var,mse,prj,OM[1],p[1],ftar[1],btrig[1],blim[1]))


# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {

  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  formulaText <- reactive({
    paste("Operating Model ~", input$OM)
  })

  # Return the formula text for printing as a caption
  output$caption <- renderText({
    formulaText()
  })

  output$msePlot <- renderPlot({
      p.=ifelse(input$p,"p = 1","p known")
      p=plotMSE(input$ts,mseIle,mseIts,prjIle,input$OM,p.,input$ftar,input$btrig,input$blim)
  
      print(p)})
  
   output$mpPlot <- renderPlot({
     p.=ifelse(input$p,"p = 1","p known")
     
     p=plotMP(input$var,mp,priors,input$OM,p.,input$ftar,input$btrig,input$blim,iters=iters)
     print(p)
  })
  
#   output$kobePlot <- renderPlot({
#     p=with(scenarios,plotKobe(mse,prj,OM,p,ftar,btrig,blim,yr=year))
#     print(p)
#   })
  
  })
