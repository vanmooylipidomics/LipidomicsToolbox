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
  
  #make a data frame of our sample names
  
  samplenames <- gsub(chosenFileSubset,"",mzXMLfiles)
  
  samplenamesframe <- data.frame(samplenames,samplenumber =
                                   seq(from=1, to=length(mzXMLfiles)))
  
  #create + extract a lists of peaks that fit our parameters 
  peaks <- chromPeaks(object = centWave, 
                      mz = c(mzlow, mzhigh),
                      rt = c(seclow, sechigh))
  
  #turn our matrix into a dataframe
  peaksframe <- as.data.frame(peaks)
  
  #pull out the columns we want
  peaksnumber <- peaksframe[["sample"]]
  peaksmz <- peaksframe[["mz"]]
  peaksrt <- (peaksframe[["rt"]]/60)
  peaksintensity <-format(peaksframe[["into"]], scientific = TRUE)
  
  #make them into another dataframe
  samplevalues <- data.frame(name = peaksnumber,
                             mz = peaksmz,
                             rt = peaksrt,
                             intensity = peaksintensity)
  
  #Add the sample names back in
  merged <- merge(samplevalues, samplenamesframe, by.x="name", by.y= "samplenumber")
  
  #Reorder our coulmns so sample name comes seconds
  reordered <- merged[c(1,5,2,3,4)]
  
  #Make everything a character so we can add a page break in <- MAYBE ADD SWITCH FOR THIS
  
  ascharacters <- as.data.frame(lapply(reordered, as.character), stringsAsFactors = FALSE)
  
  Done <- head(do.call(rbind, by(ascharacters, reordered$name, rbind, "")), -1 )
  
  # Done <- reordered
  
  colnames(Done)[1] <- "Sample Number"
  colnames(Done)[2] <- "Sample Name"
  colnames(Done)[3] <- "m/z"
  colnames(Done)[4] <- "Retention Time"
  colnames(Done)[5] <- "Intensity"
  
  print(Done)
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

#Create a value for each standard
rownames <- c("mz","ppm","rtlow","rthigh")
DNPPE <- c(875.550487, 2.5, 14, 17)
DGTSd9 <- c(721.66507, 2.5, 15, 19)
OleicAcidd9 <- c(290.30509, 2.5, 8.26, 12.26)
ArachidonicAcidd11 <- c(314.30200, 2.5, 7.13,11.13)
d7MG181 <- c(381.37042, 2.5,7.42,11.42)
LysoPEd7181 <- c(487.3524, 2.5,6.14,10.14)
waxesterd5 <- c(505.56839, 2.5,17.24,21.24)
d7LysoPC181 <- c(529.39935, 2.5,6.05,10.05)
d7DG15181 <- c(605.58444, 2.5,15.57,19.57)
d7PE15181 <- c(711.56642, 2.5,14.71,18.71)
GluCerad5 <- c(733.63488, 2.5,15.11, 18.99)
PCd715181 <- c(753.61337, 2.5, 14.59,18.59)
d5PG16181 <- c(771.59064, 2.5, 13.12, 17.12)
d5TG <- c(857.83285, 2.5, 21.26, 25.26)

#Make a list for dropdown
dropdown <- data.frame("DNPPE"= DNPPE,
                       "DGTSd9"= DGTSd9, 
                       "OleicAcidd9" = OleicAcidd9,
                       "ArachidonicAcidd11" = ArachidonicAcidd11,
                       "18_1_d7_MG" = d7MG181,
                       "18_1_d7_Lyso_PE" = LysoPEd7181,
                       "wax_ester_d5" = waxesterd5,
                       "18_1_d7_Lyso_PC" = d7LysoPC181,
                       "15_0_18_1_d7_DG" = d7DG15181,
                       "15_0_18_1_d7_PE" = d7PE15181,
                       "C18_Glucosyl_Ceramide_d5" = GluCerad5,
                       "15_0_18_1_d7_PC" = PCd715181,
                       "16_0_18_1_D5_PG" = d5PG16181,
                       "16_0_18_0_16_0_D5_TG" = d5TG,
                       row.names = rownames)


### Set Up the UI ###

ui <- shinyUI(fluidPage(
  
  titlePanel("Standard Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      #Create a dropdown menu of standards
      
      selectInput(inputId = "list",
                  label = "Existing Standard",
                  choices = colnames(dropdown)),
      
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
      
      #GOOOOOOO!
      actionButton(inputId = "runtest",
                   "Run"),
      verbatimTextOutput("parameters")
      
    ),
    
    #Add click button here.
    
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
  
  ### Combinded Code ###
  
  observeEvent(eventExpr = input$runtest, {
    mz <- if(is.na(input$mz) == TRUE){
      dropdown["mz",as.vector(input$list)]
    } else {
      input$mz
    }
    ppm <- if(is.na(input$ppm) == TRUE){
      dropdown["ppm",as.vector(input$list)]
    } else {
      input$ppm
    }
    rtlow <- if(is.na(input$rtmin) == TRUE){
      dropdown["rtlow",as.vector(input$list)]
    } else {
      input$rtmin
    }
    rthigh <- if(is.na(input$rtmax) == TRUE){
      dropdown["rthigh",as.vector(input$list)]
    } else {
      input$rtmax
    }
    
    output$table <- renderTable(makeStandardTable(obj = centWave,
                                                  mz = mz,
                                                  ppm = ppm,
                                                  rtlow = rtlow,
                                                  rthigh = rthigh
    ),
    striped = TRUE,
    align = 'l',
    width = 400,
    digits = 5
    )
    
    output$plot <- renderPlot(makeStandardGraph(obj = centWave,
                                                mz = mz,
                                                ppm = ppm,
                                                rtlow = rtlow,
                                                rthigh = rthigh
    )
    )
    output$parameters <- renderText(c('Current Settings','\nm/z =',mz,
                                      '\nppm =',ppm,
                                      '\nrtlow =',rtlow,
                                      '\nrthigh =',rthigh))
  }
  )
  
}


# Run the application 
shinyApp(ui = ui, server = server)

