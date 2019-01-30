#prepOrbidata_v2.R
# Tyrone, we should add Jaimes preamble back at the top here when we are done. - Henry
#requires: XCMS 3.0.0
#          BioParallel
#          snowParallel
#          MSnbase

################ Initial setup and variable definition #############

# run this section before specifying user settings, below

# clear workspace
#
# WS = c(ls())
# rm(WS, list = WS)

# load required packages

library(tools)

library(xcms)

library(CAMERA)

library(rsm)

library(BiocParallel)

library(snow)

library(shiny)

library(DT)

library(ggplot2)

## Use socket based parallel processing on Windows systems
if (.Platform$OS.type == "unix") {
    register(bpstart(MulticoreParam()))
} else {
    register(bpstart(SnowParam()))
} 

#timer start
ptm <- proc.time()

#USE THIS FLAG TO LAUNCH GUI
use_gui = TRUE

# ******************************************************************
################ Basic user begin editing here #############
# ******************************************************************

#### Shiny Code for File Selection UI. If you want to set your directories in the scirpt scroll down to line 135 #####

doshiny_files <- function() {
  app=shinyApp(
    ui <- fluidPage(
      # Application title
      titlePanel(
        h1("Welcome to LOBSTAHS")
      ),
      # Sidebar with Input Buttons
      sidebarLayout(
        sidebarPanel(
          h5("Select the working directory that contains the folders for both positive and negative mzXML files"),
          actionButton(inputId = 'directory',
                       label = "Select Working Directory"
          ),
          
          h5("Select the subset of files you wish to run this time"),
          actionButton(inputId = 'mzXMLdirs',
                       label = "mzXML Sub Directory"),
          
          h5("You can also  edit your file paths in the textboxs to the left. When everything looks good press done."),
          actionButton("ending","Done")
        ),
        
        # Display the main panel with selections and info
        mainPanel(
          tabsetPanel(
            tabPanel("Selections",
                     h4("Working Directory"),
                     textInput(inputId = "textdirectory", label = NULL),
                     
                     h4("File Subset"),
                     textInput(inputId = "textmzXMLdirs", label = NULL),
                     
                     h4("NO SPACES IN FILE NAMES"),
                     h5("Make sure that your file names do not include spaces as it interferes with creation of the the file lists. Don't worry about NA at the beginning of your file path.")
            )
            
            
          )
        )
      )
    ),
    
    server <- function(input, output, session) {
      
      #dir_out <- reactiveVal(value = NULL)
      
      observeEvent(eventExpr = input$directory, {
        dir_out <- choose.dir()
        dir_out<-gsub("\\\\","/",dir_out)
        updateTextInput(session, "textdirectory", value = print(dir_out))
      })
      
      observeEvent(eventExpr = input$mzXMLdirs, {
        sub_dir_out <- choose.dir()
        sub_dir_out<-gsub("\\\\","/",sub_dir_out)
        updateTextInput(session, "textdirectory", value = print(sub_dir_out))
      )}
      
      
      #end app and send outputs the Global Enviroment
      observeEvent(input$ending, {
        fpl <- list(input$textdirectory, input$textmzXMLdirs)
        
        stopApp(returnValue = fpl)
      }
      
      )
      
    }
  )
  runApp(app)
  }

################ User: define locations of data files and database(s) #############

working_dir = "/Users/TSQ/Desktop/" # specify working directory
setwd(working_dir) # set working directory to working_dir

# specify directories subordinate to the working directory in which the .mzXML files for xcms can be found; per xcms documentation, use subdirectories within these to divide files according to treatment/primary environmental variable (e.g., station number along a cruise transect) and file names to indicate timepoint/secondary environmental variable (e.g., depth)

mzXMLdirs = c("/Users/TSQ/Desktop/Pos/")

# specify which of the directories above you wish to analyze this time through

chosenFileSubset = "/Users/TSQ/Desktop/Pos/"

# specify the ID numbers (i.e., Orbi_xxxx.mzXML) of any files you don't want to push through xcms (e.g., blanks); note that the blanks for the Pt H2O2 dataset (Orbi_0481.mzXML and Orbi_0482.mzXML) have already been removed

#IMPORTANT: ***NO SPACES IN FILE NAMES*** Make sure that your file names do not include spaces as it interferes with creation of the the file lists.

#excluded.mzXMLfiles = c("0475","0474") # specifying removal of Orbi_0475.mzXML and Orbi_0474.mzXML since chromatography was screwy, to the point that weird things started to happen when I used retcor() on them
#excluded.mzXMLfiles = c("") #excluded blank this run
# if planning to use IPO, specify the ID numbers (i.e., Orbi_xxxx.mzXML) of the files you'd like to use for optimization; otherwise, IPO will try to use the entire dataset and you'll probably wind up with errors

# if you aren't planning on running IPO to optimize centWave and/or group/retcor parameters this session, but you have some parameter values from an earlier IPO run saved in a .csv file, you can specify the file paths below. you will still be given the option later to choose which parameters to use.

if (use_gui==TRUE){
  fpl <- doshiny_files()
  working_dir <- gsub(pattern = "NA", replacement = "", fpl[[1]])
  mzXMLdirs <- gsub(pattern = "NA", replacement = "", fpl[[2]])
  chosenFileSubset <- gsub(pattern = "NA", replacement = "", fpl[[2]])
} 


#####shiny code
doshiny_cent <- function() {
  app=shinyApp(
    ui = fluidPage(
      tabsetPanel(
        tabPanel("Parameters",
      column(
        width = 6,
        numericInput('ppm', 'ppm',"2.5",step = 0.1),
        numericInput('min_peakwidth', 'min_peakwidth',"20"),
        numericInput('max_peakwidth', 'max_peakwidth',"150"),
        numericInput('mzdiff', 'mzdiff',"0.005",step = .005),
        numericInput('snthresh', 'snthresh',"50"),
        numericInput('prefilter_k', 'prefilter k',"3"),
        numericInput('prefilter_I', 'prefilter I',"100"),
        numericInput('noise', 'noise',"500"),
        selectInput('mzCenterFunc', 'mzCenterFunction2',c("Intensity Weighted Mean" = "wMean",
                                                          "Intensity Mean" = "mean",
                                                          "Peak Apex" = "apex",
                                                          "Weighted Mean Apex" = "wMeanApex3",
                                                          "Mean Apex" = "meanApex3")),
        #selectInput('fitgauss', 'fit gauss',c("true" = TRUE, "false" = FALSE)),
        #selectInput('verbose_columns', 'verbose_Columns',c("true" = TRUE, "false" = FALSE)),
        selectInput('intergrate', 'intergate',c("Mexican Hat function" = 1, "raw function" = 2)),
        checkboxInput('fitgauss', 'Fit Gauss', value = TRUE, width = 5),
        checkboxInput('verbosecol', 'Verbose Columns', value = TRUE, width = NULL)
      ),
      column(
        width = 6,
        tags$p("Centwave Params :", tags$span(id = "valueA", "")),
        tags$script(
          "$(document).on('shiny:inputchanged', function(event) {
          if (event.name === 'a') {
            $('#valueA').text(event.value);
          }
        });
        "
        ),
        tableOutput('show_inputs'),
        textOutput('list_inputs'),
        actionButton("ending","Done")
      )),
      tabPanel("Info",
        HTML("
             <html>
              <body>
             
                <h2>Parameter Info</h2>
             
                  <p>Here I have put together a guide of what I understand these parameters to really mean and how you should think about setting them. Quoted text is from the XCMS manual. If you have thoughts on this or think something should be added/deleted let Henry know.</p>
             
                <h3><p>min_peakwidth and max_peakwidth</p></h3>
                   <p>Defaults to <q>20</q> and <q>150</q></p>
             
                  <p><q><i>numeric(2) with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds.</i></q></p>
             
                  <p>This essentially the widest and thinnest your peaks could possibly be. It is set in seconds. The defaults were lowered from a range of 10-60 secs (Patti et al.) to a range of 10-45 sec based on our own Hummel Chromotography. If you are running a new type on chromatography it might be a good idea to see how wide the widest peaks are.</p>
             
                <h3><p>ppm</p></h3>
                    <p>Defaults to 2.5</p>
                  <p><q><i>numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.</i></q></p>
             
                  <p>This is simply the mass tolerence in ppm. For our system this is always 2.5 ppm.</p>
             
                <h3><p>mzdiff</p></h3>
                    <p>Defaults to 0.005</p>
             
                  <p><q><i>numeric(1) representing the minimum difference in m/z dimension for peaks with overlapping retention times; can be negatove to allow overlap.</i></q></p>
             
                  <p>Essentially saying how different do two peaks at the same retention time need to be in their mass before they are identifed as two different peaks. Shouldnt change too much with good mass accuracy. DNPPE varies across samples by ~0.00015 but for smaller compounds that could be higher.</p>
             
                <h3><p>snthresh</p></h3>
                    <p>Defaults to 50</p>
             
                  <p><q><i>numeric(1) defining the signal to noise ratio cutoff.</i></q></p>
             
                  <p>Peak must be <q>X</q> times the noise to be a peak.</p>
             
                <h3><p>prefilter_k and prefilter_I</p></h3>
                    <p>Defaults to k=3 and I=100</p>
             
                  <p><q><i>numeric(2): c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I.</i></q></p>

                  <p>Before the program detects peaks, it detects <q>ROIs</q> which stands for regions of interest. It will then look in these ROIs to find peaks. It will skip a ROI if it doesnt have <q>k</q> peaks with a certian intensity <q>I</q></p>
             
                <h3><p>noise</p></h3>
                    <p>Defaults to 500</p>
                 <p><q><i>numeric(1) allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity noise are omitted from ROI detection).</i></q></p>

                <p>This value define the hieght a peak must be in order to be considered about the noise. We have set it at 500 becuase we can sometimes see useful peaks that are e^4 but anything e^3 is liekly noise. This could be set heigher for really rich sample.</p>
             
                <h3><p>fitgauss</p></h3>
                  <p>Defaults to True</p>
                  <p><q><i>logical(1) whether or not a Gaussian should be fitted to each peak.</i></q></p>

                <p>Whether or not you want to smooth your peak shapes with a filter.</p>
             
                <h3><p>verbose.columns</p></h3>
                  <p>Defaults to True</p>
                  <p><q><i>logical(1) whether additional peak meta data columns should be returned.</i></q></p>

                <p>Whether to include Metadata. Keep true.</p>
             
             </body>
             </html>")
      )
      )
    ),
    server = function(input, output, session) {
      
      AllInputs <- reactive({
        x <- reactiveValuesToList(input)
        data.frame(
          names = names(x),
          values = unlist(x, use.names = FALSE)
        )
      })
      
      output$show_inputs <- renderTable({
        AllInputs()
      })
      
      output$list_inputs <- renderText({ 
        paste("CentWaveParam(snthresh = ", input$fitgauss, ", noise = ", input$noise, ", ppm= ",input$ppm," mzdiff = centW.mzdiff, 
                     prefilter = centW.prefilter, peakwidth = c(centW.min_peakwidth,centW.max_peakwidth), fitgauss = centW.fitgauss,
                     mzCenterFun = centW.mzCenterFun, verboseColumns = centW.verbose.columns, integrate = centW.integrate
                     )")
        paste("is fitgauss logical? ", is(input$fitgauss, "logical"))
      })
      observeEvent(input$ending, {
        cwp <- CentWaveParam(snthresh = input$snthresh, noise = input$noise, ppm= input$ppm, mzdiff = centW.mzdiff,
                             prefilter = c(input$prefilter_k, input$prefilter_I), peakwidth = c(input$min_peakwidth,input$max_peakwidth), fitgauss = input$fitgauss,
                             mzCenterFun = input$mzCenterFunc, verboseColumns = input$verbosecol, integrate = input$intergrate)
        stopApp(cwp)
      })
    }
  )
  runApp(app)
}

###shiny grouping (MUSTRUN PEAK PICKING FIRST)
doshiny_group <- function() {
  app=shinyApp(
    ui = fluidPage(
      column(
        width = 6,
        numericInput('bw', 'Bandwidth',"5",step = 1),
        numericInput('max', 'Max peakgroups as feature',"50"),
        numericInput('minfrac', 'Minimum fraction to be feature',"0.25", step = 0.01),
        numericInput('minsamp', 'Minimum samples in group to be feature',"1"),
        numericInput('mzwid', 'Width of overlapping slices',"0.015")
      ),
      column(
        width = 6,
        tags$p("Centwave Params :", tags$span(id = "valueA", "")),
        tags$script(
          "$(document).on('shiny:inputchanged', function(event) {
          if (event.name === 'a') {
          $('#valueA').text(event.value);
          }
});
          "
      ),
      tableOutput('show_inputs'),
      textOutput('list_inputs'),
      actionButton("ending","Done")
      )
      
    ),
    server = function(input, output, session) {
      
      AllInputs <- reactive({
        x <- reactiveValuesToList(input)
        data.frame(
          names = names(x),
          values = unlist(x, use.names = FALSE)
        )
      })
      
      output$show_inputs <- renderTable({
        AllInputs()
      })
      observeEvent(input$ending, {
        pdp <- PeakDensityParam(sampleGroups = centWave$sampleNames,
                                bw = input$bw, minFraction = input$minfrac, minSamples = input$minsamp, 
                                binSize = input$mzwid, maxFeatures = input$max)
        stopApp(pdp)
      })
      session$onSessionEnded(function() {
        stopApp()
      })
    }
    )
  runApp(app)
}



######### PREPROCESSING PARAMETERS ##############

## PEAK PICKING
#centWave parameters
centW.min_peakwidth = 10
centW.max_peakwidth = 45 # lowered from Patti et al. recommended HPLC setting of 60 based on visual inspection of a single sample with plotPeaks
centW.ppm = 2.5   #defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.
centW.mzdiff = 0.005  # the minimum difference in m/z dimension for peaks with overlapping retention times; can be negatove to allow overlap.
centW.snthresh = 10   #defining the signal to noise ratio cutoff.
centW.prefilter = c(3,7500) # 3.5k recommended by Patti et al. appears to be too low
centW.noise = 500  #allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity < noise are omitted from ROI detection).
#centWave constant parameters
centW.fitgauss = TRUE
#centW.sleep = 1 depreciated
centW.mzCenterFun = c("wMean")
centW.verbose.columns = TRUE
centW.integrate = 1
#centW.profparam = list(step=0.001) # not available in xcms3
#centW.nSlaves = 4 # depreciated
plotcentwave =TRUE

##GROUPING
#density
# settings for group.density below are based on the recommended HPLC/Orbitrap settings from Table 1 of Patti et al., 2012, 
#"Meta-analysis of untargeted metabolomic data from multiple profiling experiment," Nature Protocols 7: 508-516
density.bw = 5 # 15?
density.max = 50
density.minfrac = 0.25
density.minsamp = 1 ### minsamp 2 does not work with groupchrom
density.mzwid = 0.015 # 0.001?

##RETCOR
# retcor.loess settings below are the function defaults
lowess.minfrac =0.9 #numeric(1) between 0 and 1 defining the minimum required fraction of samples in which peaks for the peak group were identified. 
#Peak groups passing this criteria will aligned across samples and retention times of individual spectra will be adjusted based on this alignment. 
#For minFraction = 1 the peak group has to contain peaks in all samples of the experiment.
loess.extra = 1     #numeric(1) defining the maximal number of additional peaks for all samples to be assigned to a peak group (i.e. feature) for retention time correction. For a data set with 6 samples, extraPeaks = 1 uses all peak groups with a total peak count <= 6 + 1. 
#The total peak count is the total number of peaks being assigned to a peak group and considers also multiple peaks within a sample being assigned to the group.
loess.smoothing = "loess"
loess.span = c(0.2) #defining the degree of smoothing
loess.family = "gaussian" # want to leave outliers in for the time being


## ONCE EDITS ARE COMPLETED, SAVE AND HIT "SOURCE" IN RSTUDIO TO RUN THE REST OF THE SCRIPT

################# Define functions; run me first #############

##################run centparam tweaker

# readinteger: for a given prompt, allows capture of user input as an integer; rejects non-integer input

readinteger = function(prompttext) {
  
  n = readline(prompt=prompttext)
  
  if (!grepl("^[0-9]+$", n)) {
    
    return(readinteger(prompttext))
    
  }
  
  as.integer(n)
  
}

# readyesno: for a given prompt, allows capture of user input as y or n; rejects other input

readyesno = function(prompttext) {
  
  n = readline(prompt=prompttext)
  
  if (!grepl("y|n", n)) {
    
    return(readyesno(prompttext))
    
  }
  
  as.character(n)
  
}


# verifyFileIonMode: return the ion mode of data in a particular mzXML file, by examining "polarity" attribute of each scan in the file

verifyFileIonMode = function(mzXMLfile) {
  
  rawfile = xcmsRaw(mzXMLfile) # create an xcmsraw object out of the first file
  
  # determine ion mode by examining identifier attached to scan events
  
  if (table(rawfile@polarity)["negative"]==0 & (table(rawfile@polarity)["positive"]==length(rawfile@scanindex))) { # can deduce that the file contains positive mode data
    
    filepolarity = 1 # positive
    
  } else if (table(rawfile@polarity)["positive"]==0 & (table(rawfile@polarity)["negative"]==length(rawfile@scanindex))) { # probably negative mode data
    
    filepolarity = -1 # negative
    
  } else if (table(rawfile@polarity)["positive"]>=1 & table(rawfile@polarity)["negative"]>=1) { # scans of both mode present in the file; the original .raw files weren't split by mode during initial .mzXML conversion, or something else is wrong
    
    stop("At least one file in the current dataset contains scans of more than one ion mode. Please ensure data for different ion modes have been extracted into separate files. Stopping...") # stop script if this is the case
    
  } else if (table(rawfile@polarity)["positive"]==0 & table(rawfile@polarity)["negative"]==0) {
    
    stop("Can't determine ion mode of data in the first file. Check manner in which files were converted. Stopping...") # stop script if this is the case
    
  }
  
  filepolarity
  
}

# getSubsetIonMode: return the ion mode of a subset of files, using sapply of verifyFileIonMode

getSubsetIonMode = function(mzXMLfilelist) {
  
  ionmodecount = sum(sapply(mzXMLfilelist, verifyFileIonMode)) # get sum of ion mode indicators for the files in the subset
  
  if (ionmodecount==length(mzXMLfilelist)) { # can conclude that all files contain positive mode data
    
    subset.polarity = "positive"
    
  } else if (ionmodecount==-length(mzXMLfilelist)) { # can conclude that all files contain negative mode data
    
    subset.polarity = "negative"
    
  }
  
  subset.polarity
  
}

# selectXMLSubDir: allows user to choose which subset of files to process

selectXMLSubDir = function(mzXMLdirList) {
  
  print(paste0("mzXML files exist in the following directories:"))
  
  for (i in 1:length(mzXMLdirList)) {
    
    # get number of mzXML files in this directory
    numGoodFiles = length(list.files(mzXMLdirList[i], recursive = TRUE, full.names = TRUE, pattern = "*(.mzXML|.mzxml)"))
    
    if (numGoodFiles>0) { # there are .mzXML data files in this directory
      
      print(paste0(i, ". ", numGoodFiles," .mzXML files in directory '",mzXMLdirList[i],"'"))
      
    }
    
  }
  
  processDecision = readinteger("Specify which subset you'd like to process, using integer input: ")
  
  mzXMLdirList[processDecision]
  
}

# getFNmatches: returns index(es) of file names in a given file list containing the ID numbers in a match list

getFNmatches = function(filelist,IDnumlist) {
  
  unique(grep(paste(IDnumlist,collapse="|"),filelist, value=FALSE))
  
}

# genTimeStamp: generates a timestamp string based on the current system time

genTimeStamp = function () {
  
  output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
  output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
  output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)
  
}

# check to make sure user has specified at least something in mzXMLdirs

if (!exists("mzXMLdirs")) {
  
  stop("User has not specified any directories containing mzXML files. Specify a value for mzXMLdirs.")
  
}

################# Load in mzXML files, get xcms settings from IPO or user input #############

cat(sprintf("PrepOrbit.R running on maximum %d cores for defined dataset",
            ifelse(length(detectCores()), 
                   cores <- detectCores()[1],
                   cores <- 1)))

# load selected subset for processing

mzXMLfiles.raw = list.files(chosenFileSubset, recursive = TRUE, full.names = TRUE)

# verify the ion mode of the data in these files

subset.polarity = getSubsetIonMode(mzXMLfiles.raw)

# provide some feedback to user

print(paste0("Loaded ",length(mzXMLfiles.raw)," mzXML files. These files contain ",subset.polarity," ion mode data. Raw dataset consists of:"))

print(mzXMLfiles.raw)

# check whether user has elected to exclude any files, and exclude them if they happen to be in this subset

if (exists("excluded.mzXMLfiles") & length("excluded.mzXMLfiles")>0) {
  
  excludedfiles = getFNmatches(IDnumlist = excluded.mzXMLfiles, filelist = mzXMLfiles.raw) # index files to be excluded
  
  print(paste0("The following files will be excluded from processing based on user's input:"))
  print(mzXMLfiles.raw[excludedfiles])
  
  mzXMLfiles = mzXMLfiles.raw[-excludedfiles] # exclude the files from mzXMLfiles
  
} else {
  
  mzXMLfiles = mzXMLfiles.raw
  
}


###############################################
############NEW STUFF STARTS HERE##############
###############################################


##### PEAK PICKING

#read in only msLevel1
rawSpec   <- MSnbase::readMSData(mzXMLfiles, centroided=TRUE, mode="onDisk", msLevel = 1)
rawSpecSave <-  MSnbase::readMSData(mzXMLfiles, centroided=TRUE, mode="onDisk")

if (use_gui==TRUE){
  cwp<-doshiny_cent()
  
} else {
  print(paste0("Using DEFAULT values of centWave parameters for peak picking..."))
  #format centwave parameters
  cwp <- CentWaveParam(snthresh = centW.snthresh, noise = centW.noise, ppm= centW.ppm, mzdiff = centW.mzdiff, 
                       prefilter = centW.prefilter, peakwidth = c(centW.min_peakwidth,centW.max_peakwidth), fitgauss = centW.fitgauss,
                       mzCenterFun = centW.mzCenterFun, verboseColumns = centW.verbose.columns, integrate = centW.integrate
  )
  
}
#find peaks
centWave <- findChromPeaks(rawSpec, param = cwp)


print(paste0("Peak picking completed"))

print(centWave)

if (plotcentwave == TRUE) {
  plotChromPeakImage(centWave) 
}

#print timer
print(proc.time() - ptm)


#### GROUPING AND RETCOR

#Next we group identified chromatographic peaks across samples. We use the peak density method [@Smith:2006ic] specifying 
#that a chromatographic peak has to be present in at least 2/3 of the samples within each group to be combined to a mz-rt feature.

print(paste0("Using peak density for sample groups"))

if (use_gui==TRUE){
  pdp<-doshiny_group()
  print(paste0("Using UI formatted values for grouping..."))
} else {
  print(paste0("Using DEFAULT values for grouping..."))
  #format peak density grouping parameters
  pdp <- PeakDensityParam(sampleGroups = SAMPTEST,
                          bw = density.bw, minFraction = density.minfrac, minSamples = density.minsamp, 
                          binSize = density.mzwid, maxFeatures = density.max)
}

x_density <- groupChromPeaks(centWave, param = pdp)

## Obiwarp
#print(paste0("Doing the obiwarp alignment using the default settings...."))
#rt_adjusted <- adjustRtime(x_density, param = ObiwarpParam())

## Obiwarp
#print(paste0("Doing the obiwarp alignment using the default settings...."))
#rt_adjusted <- adjustRtime(x_density, param = ObiwarpParam())

## Loess
print(paste0("Doing the loess method RT alignment using the defined settings...."))
pgp<-PeakGroupsParam(minFraction = lowess.minfrac, extraPeaks = loess.extra, smooth = loess.smoothing,
                span = loess.span, family = loess.family)

rt_adjusted <-adjustRtime(x_density, param = pgp)

## plotting difference between RT adjustment and original results
## Calculate the difference between the adjusted and the raw retention times.
xod <- rt_adjusted
diffRt <- rtime(xod) - rtime(xod, adjusted = FALSE)

## By default, rtime and most other accessor methods return a numeric vector. To
## get the values grouped by sample we have to split this vector by file/sample
diffRt <- split(diffRt, fromFile(xod))

boxplot(diffRt, main = "alignment results", ylab = "adjusted - raw rt") 

#Second round of "grouping" after RT correction
print(paste0("Performing second peak grouping after application of retcor..."))
x_2density <- groupChromPeaks(rt_adjusted, param = pdp)

#fill peaks
print(paste0("Filling peaks... NOTE: SERIAL PROCESSING ONLY"))
x_filled <- fillChromPeaks(x_2density, BPPARAM = SnowParam(workers = 7))

#This seems like it works in parallel now on small data sets.

## convert to xset and correct for missing values
xset <- x_filled
xset <-as(xset, "xcmsSet")

## important! you might want to set/adjust the 'sampclass' of the returned xcmSet object before proceeding with the analysis.
## XCMSnExp saves this as "sampleNames" so we will copy from that
xset$class <- centWave$sampleNames

#run fill peaks on xset, use all cores -2 and with memory allocation from multiple of 4
#BPPARAM_fillpeaks <- SnowParam(min(detectCores()-2,4), progressbar = TRUE)
#xset_fill <- do.call(fillPeaks,list(object = xset,BPPARAM = BPPARAM_fillpeaks))



#####################################################################################
##### Isotope peak identification, creation of xsAnnotate object using CAMERA #######
#####################################################################################


print(paste0("Applying CAMERA to identify isotopic peaks, create xsAnnotate object, and create CAMERA pseudospectra using correlation of xcms peak groups between and within samples. These pseudospectra are the groups within which the adduct hierarchy and retention time screening criteria will be applied using LOBSTAHS"))

# first, a necessary workaround to avoid a import error; see https://support.bioconductor.org/p/69414/
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

# create annotated xset using wrapper annotate(), allowing us to perform all CAMERA tasks at once

xset_a = annotate(xset,
                  
                  quick=FALSE, # set to FALSE because we want to run groupCorr; will also cause CAMERA to run adduct annotation. while LOBSTAHS will do its own adduct identification later, it doesn't hurt to do this now if it lets CAMERA create better pseudospectra
                  sample=NA, # use all samples
                  nSlaves=4, # use 4 sockets
                  
                  # group FWHM settings
                  # using defaults for now
                  
                  sigma=6,
                  perfwhm=0.6,
                  
                  # groupCorr settings
                  # using defaults for now
                  
                  cor_eic_th=0.75,
                  graphMethod="hcs",
                  pval=0.05,
                  calcCiS=TRUE,
                  calcIso=TRUE,
                  calcCaS=FALSE, # weird results with this set to TRUE
                  
                  # findIsotopes settings
                  
                  maxcharge=4,
                  maxiso=4,
                  minfrac=0.5, # 0.25?
                  
                  # adduct annotation settings
                  
                  psg_list=NULL,
                  rules=NULL,
                  polarity=subset.polarity,
                  multiplier=3,
                  max_peaks=100,
                  
                  # common to multiple tasks
                  
                  intval="into",
                  ppm=2.5,
                  mzabs=0.0015
                  
)

cleanParallel(xset_a) # kill sockets

# at this point, should have an xsAnnotate object called "xset_a" in hand, which will serve as the primary input to the main screening and annotation function "doLOBscreen" in LOBSTAHS
print(paste0("xsAnnotate object 'xset_a' has been created. User can now use LOBSTAHS to perform screening..."))

print(xset_a)

## Issue#31 fix
xset_a@groupInfo[is.na(xset_a@groupInfo)] <- 0

#print timer
print(proc.time() - ptm)


