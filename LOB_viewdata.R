### View Data Function

LOB_viewdata <- function(LOBpeaklist){
  
  #Make sure we have our librarys loaded 
  library(shiny)
  library(ggplot2)
  library(RColorBrewer)
  
  #Rename our peak list so we can modify it and keep the complete one
  run <- LOBpeaklist
  
  # Set up our large color pallete
  palette <-c("#E41A1C","#377EB8","#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
              "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
              "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A",
              "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8",
              "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C",
              "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
              "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
              "#A65628", "#F781BF", "#999999")
  
  # Define the app
  app=shinyApp(

    ui = fluidPage(
      
      title = "Lipid Data Viewer",
      
      plotOutput('plot'),
      
      hr(),
      tabsetPanel(
        tabPanel(title = "Settings",
                 column(3,
                        h4("Lipid Data Viewer"),
                        sliderInput('rt', 'X-Axis Limits (Retention Time)', 
                                    min=floor(min(run$peakgroup_rt)), 
                                    max=ceiling(max(run$peakgroup_rt)),
                                    value = c(floor(min(run$peakgroup_rt)),ceiling(max(run$peakgroup_rt))), 
                                    step=10, round=1),
                        br(),
                        checkboxInput('text', 'Display Names'),
                        checkboxInput('oxy', 'Toggle Oxidized Compounds ')
                 ),
                 column(4,
                        selectInput('class', 'Select Lipid Class', c("All",as.character(unique(run$species))),multiple = TRUE,selected = "All"),
                        selectInput('color', 'Point Color', c('None','Carbon','Double Bonds','lpSolve Fitted','Lipid Class')),
                        selectInput('color', 'Color', c('None', names(run)))
                 ),
                 column(4,
                        selectInput('facet_row', 'Facet Row', c(None='.', names(run))),
                        selectInput('facet_col', 'Facet Column', c(None='.', names(run)))
                 )
        ),
        
        tabPanel(title = "Stats")
      )
    ),
    
    
    
    # Define server logic to draw our plot
    server = function(input, output) {
      
      # Will update as varibles change
      output$plot <- renderPlot({
        
        data <- run #so we dont change our intial data 
        
        # To plot all data
        if("All"%in%input$class!=TRUE){
          data <- data[data$species==input$class,]
        }
        
        # To elimate oxy compounds if desired 
        if(input$oxy){
          data <- data[data$degree_oxidation=="0",]
        }
        
        # Construct inital plot with limits and points
        g <- ggplot(data = data,mapping = aes(x = peakgroup_rt, y = LOBdbase_mz)) +
          geom_point() +
          xlim(c(input$rt[1],input$rt[2])) +
          xlab("Retention Time (sec)") +
          ylab("m/z")
        
        # Add colors for carbon number
        if(input$color=="Carbon"){
          g <- g + geom_point(aes(color=as.character(FA_total_no_C))) +
            scale_color_manual(values = palette)
          
        }
        
        # Add color for DB number 
        if(input$color=="Double Bonds"){
          g <- g + geom_point(aes(color=as.character(FA_total_no_DB))) +
            scale_color_manual(values = palette)
        }
        

        
        # Add colors for lpSolve solutions
        if(input$color=="lpSolve Fitted"){
          g <- g + geom_point(aes(color=as.character(lpSolve))) +
            scale_color_manual(values = c("#E8DA1E","#FF3030","#B6EEA6"))
        }
        
        # Add colors for classes
        if(input$color=="Lipid Class"){
          g <- g + geom_point(aes(color=as.character(species))) +
            scale_color_manual(values = palette)
        }
        
        # Add compound names 
        if (input$text){
          g <- g + geom_text(aes(label=compound_name),nudge_y = 5,size=2)
        }

        print(g)
        
      })
    }
  )
  runApp(app)
}


LOB_viewdata(original_data)
  