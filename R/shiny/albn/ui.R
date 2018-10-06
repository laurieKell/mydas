library(shiny)
library(plyr)
library(dplyr)
library(ggplot2)

# Define UI
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Management Strategy Evaluator"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("OM", "Operating Model",
                list("Base Case"     = 108, 
                     "M"             = 17, 
                     "Steepness"     = 22, 
                     "Lorenzen"      = 23, 
                     "Mat Selection" = 29, 
                     "Dome Selection"= 53, 
                     "Juve Selection"= 85),
                selected="Data Poor"),
    
    selectInput("ts", "Time Series",
                list("SSB"           = "ssbRel",
                     "Catch"         = "catchRel",
                     "Harvest Rate"  = "harvestRel",
                     "Recruitment"   = "recRel")),
    
    selectInput("var", "Reference Point",
                list("MSY"       = "msy",
                     "FMSY"      = "fmsy",
                     "BMSY"      = "bmsy")),
 
    selectInput("ftar", "Target F",
                list("FMSY"             = 1,
                     "FMSY times 0.75"  = .75,
                     "FMSY times 0.50"  = .5)),
    
    selectInput("blim", "Biomass Limit",
                list("Blim times 0.4"      = .4,
                     "Blim times 0.3"      = .3)),
    
    selectInput("btrig", "Biomass Trigger",
                list("BTrigger times 0.8"  = .8,
                     "BTrigger times 0.6"  = .6)),
    

    checkboxInput("p", "Logistic Production Function", FALSE), 
    
    helpText("Note: ...")
  ),

  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    h3(textOutput("caption")),
    
    plotOutput("msePlot"),
    plotOutput("mpPlot"))
))
