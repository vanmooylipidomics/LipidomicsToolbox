#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

### cruSTATION Standard Checker ###

### Load Packages ###

library(shiny)

library(xcms)

### Define the Table Creating Function ###

#Funtion to create the table
makeStandardTable = function(obj, mz, ppm, rtlow, rthigh) {
  
  #turn our minutes into seconds
  seclow <- (rtlow*60) 
  sechigh <- (rthigh*60)
  
  #create our m/z range based on the ppm
  mzrange <- mz*(0.000001*ppm)
  mzlow <- (mz-mzrange)
  mzhigh <- (mz+mzrange)
  
  #create extract a list of lists of peaks that fit our parameters 
  peaks <- chromPeaks(object = obj, 
                      mz = c(mzlow, mzhigh), 
                      rt = c(seclow, sechigh), 
                      bySample = TRUE)
  
  #
  peaksmz <- sapply(peaks, function(x) x[1])
  peaksrt <- sapply(peaks, function(x) (x[4]/60))
  peaksintensity <- sapply(peaks, function(x) x[7])
  data.frame(name = mzXMLfiles, mz = peaksmz, rt = peaksrt, intensity = peaksintensity)
  
}

#Function to create the graph
makeStandardGraph = function(obj, mz, ppm, rtlow, rthigh) {
  
  seclow <- (rtlow*60) 
  sechigh <- (rthigh*60)
  
  mzrange <- mz*(0.000001*ppm)
  mzlow <- (mz-mzrange)
  mzhigh <- (mz+mzrange)
  
  chroms <- chromatogram(object = obj, mz = c(mzlow,mzhigh),rt = c(seclow, sechigh))
  
  plot(chroms)
  
}

#Create a data.frame  for each standard
rownames <- c("mz","ppm","rtlow","rthigh")
DNPPE <- c(875.550487, 2.5, 14, 17)
DGTSd9 <- c(721.66507, 2.5, 15, 19)



#Make a list for dropdown
dropdown <- data.frame(DNPPE, DGTSd9, row.names = rownames )


### Set Up the UI ###

ui <- shinyUI(fluidPage(
  
  titlePanel("Standard Checker"),
  
  sidebarLayout(
    
    sidebarPanel(
      #Create a dropdown menu of standards
      
      selectInput(inputId = "list",
                  label = "Existing Standard",
                  choices = colnames(dropdown)),
      
      actionButton(inputId = "runexisting",
                   "Run Using Selected"),
      
      # Input: MZ
      numericInput(inputId = "mz", 
                   label = "m/z", 
                   value = NULL, 
                   min = 300, 
                   max = 2000, 
                   step = 0.00001, 
                   width = NULL),
      
      # Input: PPM
      numericInput(inputId = "ppm",
                   label = "Mass Tolerance (ppm)",
                   value = NULL,
                   min = 1,
                   max = 100,
                   step = 1,
                   width = NULL),
      
      # Input: RT Low
      numericInput(inputId = "rtmin",
                   label = "Retention Time - Min",
                   value = NULL,
                   min = 0,
                   max = 30,
                   step = 0.01,
                   width = NULL),
      
      # Input: RT High
      numericInput(inputId = "rtmax",
                   label = "Retention Time - Max",
                   value = NULL,
                   min = 0,
                   max = 30,
                   step = 0.01,
                   width = NULL),
      
      
      actionButton(inputId = "runnew",
                   "Run Using Inputs")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Table", tableOutput("table")), 
        tabPanel("Summary", verbatimTextOutput("summary")), 
        tabPanel("Plot", plotOutput("plot"))
      )
    )
  )
))

### server - The code behind the UI
server <- function(input, output) {

### Code for Exisiting imputs ###

  observeEvent(input$runexisting, {
    frame <- makeStandardTable(obj = centWave,
                               mz = dropdown[1,as.vector(input$list)],
                               ppm = dropdown[2,as.vector(input$list)],
                               rtlow = dropdown[3,as.vector(input$list)],
                               rthigh = dropdown[4,as.vector(input$list)])
    output$table <- renderTable(frame)
  })
  
  observeEvent(input$runexsiting, {
    output$plot <- renderPlot(makeStandardGraph(obj = centWave,
                                                mz = dropdown[1, as.vector(input$list)],
                                                ppm = dropdown[2,as.vector(input$list)],
                                                rtlow = dropdown[3,as.vector(input$list)],
                                                rthigh = dropdown[4,as.vector(input$list)])
    )
  }
  )
  
### Code for novel Inputs ###

  # Create a table based on the inputs when you press the buttom
  observeEvent(input$runnew, {
    frame <- makeStandardTable(obj = centWave,
                    mz = input$mz,
                    ppm = input$ppm,
                    rtlow = input$rtmin,
                    rthigh = input$rtmax)
    output$table <- renderTable(frame)
  })
  
  # Create a plot based on the inputs when you press the buttom
  observeEvent(input$runnew, {
    output$plot <- renderPlot(makeStandardGraph(obj = centWave,
                               mz = input$mz,
                               ppm = input$ppm,
                               rtlow = input$rtmin,
                               rthigh = input$rtmax)
                              )
                          }
              )
}


# Run the application 
shinyApp(ui = ui, server = server)

