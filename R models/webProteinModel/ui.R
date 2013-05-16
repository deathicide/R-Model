library(shiny)

#make a button that increments when clicked
incrementButton <- function(inputId, value = 0) {
  tagList(
    singleton(tags$head(tags$script(src = "js/increment.js"))),
    tags$button(id = inputId,
                class = "increment btn",
                type = "button",
                as.character(value))
  )
}

# Define UI for web protein modelling
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Protein"),
  
  #make a panel with a dropdown box containing choices for models
  sidebarPanel(
    tags$style(type="text/css",'#dv_select(max-width: 185px)'),
    p("Please select a file containing protein model information"),
    p("The file needs to be in the format of number,x,y,z"),
    fileInput("atomsFile","Open File"),
    #ask the user to either upload a file containing tmdNumber,startAtom,endAtom
    #or allow the user to enter the values in manually
    radioButtons("tmdAtomMethod","Select method to determine TMD atom numbers:",
                 c("File" = "auto","Manual" = "manual")),
    #panel with a file upload button
    conditionalPanel(
      condition = "input.tmdAtomMethod == 'auto'",p("Please select a file in the format: tmdNumber,start,end"),
        fileInput("tmdLocationNum","Open File")),
    #panel with all of the manual input boxes
    conditionalPanel(
      condition = "input.tmdAtomMethod == 'manual'",
      p("TMD1 atom numbers"),
      tags$style(type="text/css",'#tmd1s {width: 40px}'),
      tags$style(type="text/css",'#tmd1e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd1s","Start")),
          div(class = "span2", textInput("tmd1e","End"))
      ),
      p("TMD2 atom numbers"),
      tags$style(type="text/css",'#tmd2s {width: 40px}'),
      tags$style(type="text/css",'#tmd2e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd2s","Start")),
          div(class = "span2", textInput("tmd2e","End"))
      ),
      p("TMD3 atom numbers"),
      tags$style(type="text/css",'#tmd3s {width: 40px}'),
      tags$style(type="text/css",'#tmd3e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd3s","Start")),
          div(class = "span2", textInput("tmd3e","End"))
      ),
      p("TMD4 atom numbers"),
      tags$style(type="text/css",'#tmd4s {width: 40px}'),
      tags$style(type="text/css",'#tmd4e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd4s","Start")),
          div(class = "span2", textInput("tmd4e","End"))
      ),
      p("TMD5 atom numbers"),
      tags$style(type="text/css",'#tmd5s {width: 40px}'),
      tags$style(type="text/css",'#tmd5e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd5s","Start")),
          div(class = "span2", textInput("tmd5e","End"))
      ),
      p("TMD6 atom numbers"),
      tags$style(type="text/css",'#tmd6s {width: 40px}'),
      tags$style(type="text/css",'#tmd6e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd6s","Start")),
          div(class = "span2", textInput("tmd6e","End"))
      ),
      p("TMD7 atom numbers"),
      tags$style(type="text/css",'#tmd7s {width: 40px}'),
      tags$style(type="text/css",'#tmd7e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd7s","Start")),
          div(class = "span2", textInput("tmd7e","End"))
      ),
      p("TMD8 atom numbers"),
      tags$style(type="text/css",'#tmd8s {width: 40px}'),
      tags$style(type="text/css",'#tmd8e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd8s","Start")),
          div(class = "span2", textInput("tmd8e","End"))
      ),
      p("TMD9 atom numbers"),
      tags$style(type="text/css",'#tmd9s {width: 40px}'),
      tags$style(type="text/css",'#tmd9e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd9s","Start")),
          div(class = "span2", textInput("tmd9e","End"))
      ),
      p("TMD10 atom numbers"),
      tags$style(type="text/css",'#tmd10s {width: 40px}'),
      tags$style(type="text/css",'#tmd10e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd10s","Start")),
          div(class = "span2", textInput("tmd10e","End"))
      ),
      p("TMD11 atom numbers"),
      tags$style(type="text/css",'#tmd11s {width: 40px}'),
      tags$style(type="text/css",'#tmd11e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd11s","Start")),
          div(class = "span2", textInput("tmd11e","End"))
      ),
      p("TMD12 atom numbers"),
      tags$style(type="text/css",'#tmd12s {width: 40px}'),
      tags$style(type="text/css",'#tmd12e {width: 40px}'),
      div(class = "row", div(class = "span1",""),
          div(class = "span3", textInput("tmd12s","Start")),
          div(class = "span2", textInput("tmd12e","End"))
      )),
    incrementButton("btn1")
  ),
  #write the main panel
  mainPanel(
    tabsetPanel(
      tabPanel("Input",
        h3(textOutput("caption"))
        #plotOutput(outputId = "mainPlot")
               ),
      tabPanel("Inverted Model")
      )
    )
))