library(shiny)
library(plyr)
library(reshape)
library(RPostgreSQL)
library(DBI)
library(ggplot2)

#mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)
#mseStart['pollack']
#mseStart["pollack"][[1]]
#mpb =dbGetQuery(con,paste0("select * from sbt1 where spp = 'pollack' and year > 56 "))

drv  = dbDriver("PostgreSQL")

con=dbConnect(drv, host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com',
              dbname='FLRout',
              port = 5432,
              user = 'MydasAdmin',
              password = 'datapoor1!')

ui <- bootstrapPage(
  sidebarLayout(
    sidebarPanel(width=3,
                 
                 div(style="display: inline-block;vertical-align:bottom; width: 180px;", selectInput(inputId = "mp",
                                                                                                     label = "Choose a management procedure",
                                                                                                     choices = c("mpbsum","sbt1sum", "sbt2sum"),
                                                                                                     multiple = FALSE, selected='mpb')),
                 br(),
                 div(style="display: inline-block;vertical-align:bottom; width: 180px;",selectInput(inputId = "spp",
                                                                                                    label = "Choose a mydas species",
                                                                                                    choices = c("brill","turbot","ray","pollack","sprat","razor","lobster"),
                                                                                                    multiple = FALSE, selected="pollack")),
                 
                 sliderInput("discount", "Yield NPV rate:",
                             min = 0, max = 1,
                             value = 0.05, step = 0.05)
                 
    ),
    
    mainPanel( 
      wellPanel(
        tags$body(  
          
          h4( "Management Procedures"),
          h5(" Fishing targets"))),
      plotOutput("mp", width = "95%", height = "800px")
      
    )
    
  )
)


server <- shinyServer(function(input, output, session) {
  datMP = reactive({ 
    withProgress(message = 'Getting data',
                 detail = 'Please be patient... ', value = 0, {
                   for (i in 1:15) {
                     incProgress(1/15)
                     Sys.sleep(0.25)
                   }
                   mseStart=c("brill"=54,"turbot"=54,"ray"=60,"pollack"=56,"sprat"=52,"razor"=54,"lobster"=57)
                   yr = mseStart[input$spp][[1]]
                   
                   if(input$mp=="mpbsum"){
                     mp = dbGetQuery(con,paste0("select * from ",input$mp," 
                                                where spp = '",input$spp,"' and 
                                                year > ",yr,""))
                   }
                   if(input$mp=="sbt1sum"){
                     mp = dbGetQuery(con,paste0("select * from ",input$mp," 
                                                where spp = '",input$spp,"' and 
                                                year > ",yr," and 
                                                gamma >= 1"))
                   }
                   if(input$mp=="sbt2sum"){
                     mp = dbGetQuery(con,paste0("select * from ",input$mp," 
                                                where spp = '",input$spp,"' and 
                                                year > ",yr,""))
                   }
                   
                   
                   
                   dRate=function (x, r, wtAv = FALSE) 
                   {
                     if (wtAv) 
                       return(sum(x/(1 + r)^(0:(length(x) - 1)))/sum(1/(1 +  r)^(0:(length(x) - 1))))
                     return(sum(x/(1 + r)^(0:(length(x) - 1))))
                   }
                   
                   avFn<-function(object){
                     
                     o1 =object[-length(object)]
                     o2 =object[-1]
                     
                     return(abs((o2-o1)/o1))}
                   
                   green<-function(flag,year){
                     flag=flag>=1
                     year=as.numeric(year)
                     
                     if (sum(flag)==0){
                       y=max(year)+1
                       p=0
                       n=0
                     }else{
                       y   =year[flag][1]
                       flag=flag[year>y]
                       n   =length(flag[flag])
                       p   =mean(flag)}
                     
                     return(data.frame(year=y,n=n,p=p))}
                   
                   udFn<-function(object){
                     
                     o1 =object[-1]
                     o2 =object[-length(object)]
                     
                     return(sum(o1<o2))/(length(object)-1)}
                   
                   smryStat<-function(mp,dr=input$discount){
                     with(mp, data.frame(safety  =min(rec_hat/virgin_rec,na.rm=T),
                                         kobe    =green((ssb/msy_ssb>1)&(fbar/msy_harvest<1),year)[1:2],
                                         yield   =dRate(catch/msy_yield,dr)/length(catch),
                                         yieldAav=mean(avFn(pmax(catch,msy_yield*0.1)))))}
                   
                   
                   })
    
    
    
    if(input$mp=="mpbsum"){
      sm=ddply(mp,.(iter,ftar,btrig), smryStat)
    }else{
      if(input$mp=="sbt1sum"){ 
        sm=ddply(mp,.(iter,k1, k2, gamma), smryStat) 
      }else{
        if(input$mp=="sbt2sum"){  
          sm = ddply(mp,.(iter,k1, k2), smryStat)
        }
      }
    }      
    return(sm)   
                   })
  
  
  
  output$mp <-renderPlot({
    req(input$spp)
    req(input$mp)
    smplot = datMP()
    if (length(smplot) == 0){
      return(NULL)}
    if (input$mp== "mpbsum") {
      smplot$ftar = as.character(smplot$ftar)
      ggplot(melt(smplot,id=c("iter","ftar","btrig")))+
        geom_boxplot(aes(ftar, as.numeric(value), fill=ftar), notch=FALSE)+
        facet_grid(variable~btrig, scales="free")+theme_bw()+
        theme(panel.grid.major = element_blank(),
              text = element_text(size=18),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              panel.border = element_rect(colour = "black"),
              legend.position="none",
              legend.key.size =  unit(0.5, "in"),
              plot.title = element_text(hjust = 0.5)) +
        scale_fill_manual(values=c("#E69F00", "#56B4E9", "grey"))  + ylab(" ") +xlab("ftar")+ ggtitle("btrig")
    }else{
      if (input$mp== "sbt1sum") { 
        ggplot(melt(smplot,id=c("iter","k1","k2","gamma")))+
          geom_boxplot(aes(factor(gamma), as.numeric(value), fill=factor(gamma)), notch=FALSE)+
          facet_grid(variable~k1+k2, scales="free")+theme_bw()+
          theme(panel.grid.major = element_blank(),
                text = element_text(size=18),
                panel.grid.minor = element_blank(),
                strip.background = element_blank(),
                panel.border = element_rect(colour = "black"),
                legend.position="none",
                legend.key.size =  unit(0.5, "in"),
                plot.title = element_text(hjust = 0.5)) +
          scale_fill_manual(values=c("#E69F00", "#56B4E9", "grey"))  + ylab(" ") +xlab("gamma")+ ggtitle("k1/k2")
        
      }else{
        if (input$mp== "sbt2sum") { 
          ggplot(melt(smplot,id=c("iter","k1","k2")))+
            geom_boxplot(aes(factor(k1), as.numeric(value), fill=factor(k1)), notch=FALSE)+
            facet_grid(variable~k2, scales="free")+theme_bw()+
            theme(panel.grid.major = element_blank(),
                  text = element_text(size=18),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black"),
                  legend.position="none",
                  legend.key.size =  unit(0.5, "in"),
                  plot.title = element_text(hjust = 0.5)) +
            scale_fill_manual(values=c("#E69F00", "#56B4E9", "grey"))  + ylab(" ") +xlab("k1")+ ggtitle("k2")
          
        }
      }} 
    
  })
  
  
  })

shinyApp(ui, server)