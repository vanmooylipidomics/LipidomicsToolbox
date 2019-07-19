#prepOrbidata_v3.R
# Add Jaimes preamble back at the top here when we are done. - Henry
# Requires: XCMS 3.0.0
#           BioParallel
#           snowParallel
#           MSnbase
#           CAMERA

################ Load data and define basic variables and functions #############

# Load required packages

library(tools)

library(xcms)

library(CAMERA)

library(BiocParallel)

library(shiny)

library(parallel)

library(LOBtools)


################ User: define locations of data files and database(s) #############

working_dir = "C:/your_wd/goes/here/" # specify working directory
setwd(working_dir) # set working directory to working_dir

# specify directories subordinate to the working directory in which the .mzXML files 
# for xcms can be found; per xcms documentation, use subdirectories within these to 
# divide files according to treatment/primary environmental variable (e.g., station
# number along a cruise transect) and file names to indicate timepoint/secondary 
# environmental variable (e.g., depth)

chosenFileSubset = c("mzXML_ms1_pos_mode/")

# IMPORTANT: ***NO SPACES IN FILE STRINGS*** Make sure that your file names do not 
# include spaces as it interferes with creation of the the file lists.

# check to make sure user has specified at least something in mzXMLdirs

if (!exists("chosenFileSubset")) {
  
  stop("User has not specified any directories containing mzXML files. Specify a value for mzXMLdirs.")
  
}

################## Unchanging Functions and Settings ######################

## Use socket based parallel processing on Windows systems
if (.Platform$OS.type == "unix") {
    register(bpstart(MulticoreParam()))
} else {
    register(bpstart(SnowParam()))
} 

#Start timer to see how long processes take
ptm <- proc.time()

#### Create Helper Functions #####
# doshiny_cent - Launch the centwave UI
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
                   selectInput('mzCenterFunc', 'mzCenterFunction2',selected ="wMeanApex3",c("Intensity Weighted Mean" = "wMean",
                                                                     "Intensity Mean" = "mean",
                                                                     "Peak Apex" = "apex",
                                                                     "Weighted Mean Apex" = "wMeanApex3",
                                                                     "Mean Apex" = "meanApex3")),
                   selectInput('intergrate', 'Intergrate',c( "raw function" = 2,"Mexican Hat function" = 1)),
                   checkboxInput('fitgauss', 'Fit Gauss', value = FALSE, width = 5),
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
                 HTML("<html>
                      <body>
                      
                      <h3>Parameter Info</h2>
                      
                      <p>Here I have put together a guide of what I understand these parameters to really mean and how you should think about setting them. Quoted text is from the XCMS manual. If you have thoughts on this or think something should be added/deleted let Henry know.</p>
                      
                      <h3><p>ppm</p></h3>
                      <p>Defaults to 2.5</p>
                      <p><q><i>numeric(1) defining the maximal tolerated m/z deviation in consecutive scans in parts per million (ppm) for the initial ROI definition.</i></q></p>
                      
                      <p>This is simply the mass tolerence in ppm. For our system this is always 2.5 ppm.</p>
                      
                      <h3><p>min_peakwidth and max_peakwidth</p></h3>
                      <p>Defaults to <q>20</q> and <q>150</q></p>
                      
                      <p><q><i>numeric(2) with the expected approximate peak width in chromatographic space. Given as a range (min, max) in seconds.</i></q></p>
                      
                      <p>This essentially the widest and thinnest your peaks could possibly be. It is set in seconds. The defaults were lowered from a range of 10-60 secs (Patti et al.) to a range of 10-45 sec based on our own Hummel Chromotography. If you are running a new type on chromatography it might be a good idea to see how wide the widest peaks are.</p>
                      
                      <h3><p>snthresh</p></h3>
                      <p>Defaults to 50</p>
                      
                      <p><q><i>numeric(1) defining the signal to noise ratio cutoff.</i></q></p>
                      
                      <p>Peak must be <q>X</q> times the noise to be a peak.</p>
                      
                      <h3><p>mzdiff</p></h3>
                      <p>Defaults to 0.005</p>
                      
                      <p><q><i>numeric(1) representing the minimum difference in m/z dimension for peaks with overlapping retention times; can be negatove to allow overlap.</i></q></p>
                      
                      <p>Essentially saying how different do two peaks at the same retention time need to be in their mass before they are identifed as two different peaks. Shouldnt change too much with good mass accuracy. DNPPE varies across samples by ~0.00015 but for smaller compounds that could be higher.</p>
                      
                      <h3><p>prefilter_k and prefilter_I</p></h3>
                      <p>Defaults to k=3 and I=100</p>
                      
                      <p><q><i>numeric(2): c(k, I) specifying the prefilter step for the first analysis step (ROI detection). Mass traces are only retained if they contain at least k peaks with intensity >= I.</i></q></p>
                      
                      <p>Before the program detects peaks, it detects <q>ROIs</q> which stands for regions of interest. It will then look in these ROIs to find peaks. It will skip a ROI if it doesnt have <q>k</q> peaks with a certian intensity <q>I</q></p>
                      
                      <h3><p>noise</p></h3>
                      <p>Defaults to 500</p>
                      <p><q><i>numeric(1) allowing to set a minimum intensity required for centroids to be considered in the first analysis step (centroids with intensity noise are omitted from ROI detection).</i></q></p>
                      
                      <p>This value define the hieght a peak must be in order to be considered above the noise. We have set it at 500 becuase we can sometimes see useful peaks that are e^4 but anything e^3 is likely noise. This could be set heigher for really rich sample.</p>
                      
                      <h3><p>fitgauss</p></h3>
                      <p>Defaults to True</p>
                      <p><q><i>logical(1) whether or not a Gaussian should be fitted to each peak.</i></q></p>
                      
                      <p>Whether or not you want to smooth your peak shapes with a filter.</p>
                      
                      <h3><p>verbose.columns</p></h3>
                      <p>Defaults to True</p>
                      <p><q><i>logical(1) whether additional peak meta data columns should be returned.</i></q></p>
                      
                      <p>Whether to include Metadata. Keep true.</p>
                      
                      <h3><p>mzCenterFun</p></h3>
                      <p>Defaults to wMeanApex3 </p>
                      <p><q><i>Name of the function to calculate the m/z center of the chromatographic peak. Allowed are: 'wMean': intensity weighted mean of the peak's m/z values, 'mean': mean of the peak's m/z values, 'apex': use the m/z value at the peak apex, 'wMeanApex3': intensity weighted mean of the m/z value at the peak apex and the m/z values left and right of it and 'meanApex3': mean of the m/z value of the peak apex and the m/z values left and right of it.</i></q></p>
                      
                      <h3><p> integrate </p></h3>
                      <p>Defaults to raw function </p>
                      <p><q><i>Integration method. For integrate = 1 peak limits are found through descent on the mexican hat filtered data, for integrate = 2 the descent is done on the real data. The latter method is more accurate but prone to noise, while the former is more robust, but less exact.</p></q></i>
                      
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
        cat(snthresh = input$snthresh, noise = input$noise, ppm= input$ppm, mzdiff = input$mzdiff,
                          prefilter = c(input$prefilter_k, input$prefilter_I), peakwidth = c(input$min_peakwidth,input$max_peakwidth), fitgauss = input$fitgauss,
                          mzCenterFun = input$mzCenterFunc, verboseColumns = input$verbosecol, integrate = input$intergrate)
      })
      observeEvent(input$ending, {
        cwp <- CentWaveParam(snthresh = input$snthresh, noise = input$noise, ppm= input$ppm, mzdiff = input$mzdiff,
                             prefilter = c(input$prefilter_k, input$prefilter_I), peakwidth = c(input$min_peakwidth,input$max_peakwidth), fitgauss = input$fitgauss,
                             mzCenterFun = input$mzCenterFunc, verboseColumns = input$verbosecol, integrate = input$intergrate)
        stopApp(cwp)
      })
    }
                 )
  runApp(app)
  }

# doshiny_group - Launch the grouping UI
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
        tags$p("Grouping Params :", tags$span(id = "valueA", "")),
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

# doshiny_rt - Launch the RT correction UI
doshiny_rt <- function(object) {
  app=shinyApp(
    ui = fluidPage(
      column(
        width = 6,
        numericInput('loess.minfrac', 'Minimum fraction of samples in at least one sample group',"0.9",step = 1),
        numericInput('loess.extra', 'Max Extra Peaks',"1"),
        selectInput('loess.smoothing', 'Smoothing Function',choices = c("loess","linear"),selected = "loess"),
        numericInput('loess.span', 'Degree of Smoothing',"0.2"),
        selectInput('loess.family','Method to be used for loess Smoothing',choices = c("gaussian","symmetric"),selected = "gaussian"),
        selectInput('exclude',"Optional Sample Exclusion",choices = c("None",sampleNames(object)),multiple = TRUE,selected = "None")
      ),
      column(
        width = 6,
        tags$p("Retcor Params :", tags$span(id = "valueA", "")),
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
        
        if (!("None" %in% input$exclude)) {
          excld <-as.character(input$exclude)
          excld_n <- which(sampleNames(object) %in% excld)
          subset <- (1:length(sampleNames(object)))[-excld_n]
        }else{
          subset <- 1:length(sampleNames(object))
        }
        
        pgp <- PeakGroupsParam(minFraction = input$loess.minfrac,
                               extraPeaks = input$loess.extra,
                               smooth = input$loess.smoothing,
                               span = input$loess.span,
                               family = as.character(input$loess.family),
                               subset = subset,
                               subsetAdjust = "average")
        stopApp(pgp)
      })
      session$onSessionEnded(function() {
        stopApp()
      })
    }
    )
  runApp(app)
}

# readinteger: for a given prompt, allows capture of user input as an integer; rejects non-integer input

readinteger = function(prompttext) {
  
  n = readline(prompt=prompttext)
  
  if (!grepl("^[0-9]+$", n)) {
    
    return(readinteger(prompttext))
    
  }
  
  as.integer(n)
  
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

# getFNmatches: returns index(es) of file names in a given file list containing the ID numbers in a match list

getFNmatches = function(filelist,IDnumlist) {
  
  unique(grep(paste(IDnumlist,collapse="|"),filelist, value=FALSE))
  
}

######################## End of Functions Section #################################

######################### Load in mzXML files ######################

# Print Core Feedback

cat(sprintf("PrepOrbit.R running on maximum %d cores for defined dataset",
            ifelse(length(detectCores()), 
                   cores <- detectCores()[1],
                   cores <- 1)))

# load selected subset for processing

mzXMLfiles.raw = list.files(chosenFileSubset, recursive = TRUE, full.names = TRUE)

# verify the ion mode of the data in these files

subset.polarity = getSubsetIonMode(mzXMLfiles.raw)

# provide some feedback to user

print(paste0("Loaded ",length(mzXMLfiles.raw)," mzXML files. These files contain ",subset.polarity,
             " ion mode data. Raw dataset consists of:"))

print(mzXMLfiles.raw)

#####################################################
############ DATA ANAYLSIS STARTS HERE ##############
#####################################################

##### PEAK PICKING

mzXMLfiles <- mzXMLfiles.raw

# read in only msLevel1
rawSpec <- MSnbase::readMSData(mzXMLfiles, centroided=TRUE, mode="onDisk", msLevel = 1)

# read in all levels for MS2 access later 
rawSpec_ms2 <- MSnbase::readMSData(mzXMLfiles, centroided=TRUE, mode="onDisk")

# Launch GUI for picking parameters
cwp<-doshiny_cent()

#find peaks
centWave <- findChromPeaks(rawSpec, param = cwp)

#print info
print(paste0("Peak picking completed"))
print(centWave)

#plot image of detected images
plotChromPeakImage(centWave) 

#print timer
print(proc.time() - ptm)

#### GROUPING

# Next we group identified chromatographic peaks across samples. 
# We use the peak density method [@Smith:2006ic] specifying 
# that a chromatographic peak has to be present in at least a certian number 
# of the samples within each group to be combined to a mz-rt feature.

# Launch GUI for grouping parameters
pdp<-doshiny_group()

pdp@sampleGroups <- rep(1,length(mzXMLfiles))

x_density <- groupChromPeaks(centWave, param = pdp)

#### RETENTION TIME CORRECTION

print(paste0("Doing the loess method RT alignment using the defined settings...."))

pgp <- doshiny_rt(x_density)

rt_adjusted <-adjustRtime(x_density, param = pgp)

## plotting difference between RT adjustment and original results
## Calculate the difference between the adjusted and the raw retention times.
xod <- rt_adjusted
diffRt <- rtime(xod) - rtime(xod, adjusted = FALSE)

## By default, rtime and most other accessor methods return a numeric vector. To
## get the values grouped by sample we have to split this vector by file/sample
diffRt <- split(diffRt, fromFile(xod))

#Plot the boxplot
boxplot(diffRt, main = "alignment results", ylab = "adjusted - raw rt")

#set up some colors for the squid plot
cols <- RColorBrewer::brewer.pal(8,name = "Dark2")
pall <- colorRampPalette(cols)
colors <-pall(length(sampleNames(rt_adjusted)))

#plot the squid plot with and without the legend
plotAdjustedRtime(rt_adjusted,col = colors,peakGroupsPch = 4)
plotAdjustedRtime(rt_adjusted,col = colors,peakGroupsPch = 4)
legend('bottomleft',legend = sampleNames(rt_adjusted),col = colors,pch = 16)

#Second round of "grouping" after RT correction
print(paste0("Performing second peak grouping after application of retcor..."))

#Now that we have RT corrected we can reset our bandwitdh to 5 for tight groups
pdp@bw <- 5

x_2density <- groupChromPeaks(rt_adjusted, param = pdp)

#fill peaks
print(paste0("Filling peaks... NOTE: SERIAL PROCESSING ONLY"))
x_filled <- fillChromPeaks(x_2density, BPPARAM = SnowParam(workers = 7))

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

xset_a = CAMERA::annotate(xset,
                  
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

#Fix for LOBSTAHS XCMS 3 error. LOBSTAHS seems to need a more factor levels to run than there are samples. 
#Here I simply create a factor with double the levels of samples in the "class" slot. Temp fix until we fix dolobscreen()
xset_a@xcmsSet@phenoData$class <- factor(x = rep(1,length(mzXMLfiles)),levels = 1:(length(mzXMLfiles)*2))


