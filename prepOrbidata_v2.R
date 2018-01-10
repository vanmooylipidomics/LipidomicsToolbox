#prepOrbidata_v2.R
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

## Use socket based parallel processing on Windows systems
if (.Platform$OS.type == "unix") {
    register(bpstart(MulticoreParam()))
} else {
    register(bpstart(SnowParam()))
} 

#timer start
ptm <- proc.time()
# ******************************************************************
################ Basic user begin editing here #############
# ******************************************************************

################ User: define locations of data files and database(s) #############

working_dir = "/home/tyronelee/PtH2O2lipids/mzXML" # specify working directory
setwd(working_dir) # set working directory to working_dir

# specify directories subordinate to the working directory in which the .mzXML files for xcms can be found; per xcms documentation, use subdirectories within these to divide files according to treatment/primary environmental variable (e.g., station number along a cruise transect) and file names to indicate timepoint/secondary environmental variable (e.g., depth)

mzXMLdirs = c("Pt_H2O2_mzXML_ms1_pos/","Pt_H2O2_mzXML_ms1_neg/")

# specify which of the directories above you wish to analyze this time through

chosenFileSubset = "Pt_H2O2_mzXML_ms1_pos/"

# specify the ID numbers (i.e., Orbi_xxxx.mzXML) of any files you don't want to push through xcms (e.g., blanks); note that the blanks for the Pt H2O2 dataset (Orbi_0481.mzXML and Orbi_0482.mzXML) have already been removed

excluded.mzXMLfiles = c("0475","0474") # Jamie's notes: specifying removal of Orbi_0475.mzXML and Orbi_0474.mzXML since chromatography was screwy, to the point that weird things started to happen when I used retcor() on them

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
plotcentwave =FALSE

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
print(paste0("Using values of centWave parameters specified by Jamie for peak picking..."))


#read in only msLevel1
rawSpec   <- MSnbase::readMSData(mzXMLfiles, centroided=TRUE, mode="onDisk", msLevel = 1)

#format centwave parameters
cwp <- CentWaveParam(snthresh = centW.snthresh, noise = centW.noise, ppm= centW.ppm, mzdiff = centW.mzdiff, 
                     prefilter = centW.prefilter, peakwidth = c(centW.min_peakwidth,centW.max_peakwidth), fitgauss = centW.fitgauss,
                     mzCenterFun = centW.mzCenterFun, verboseColumns = centW.verbose.columns, integrate = centW.integrate
                     )

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

pdp <- PeakDensityParam(sampleGroups = centWave$sampleNames,
                        bw = density.bw, minFraction = density.minfrac, minSamples = density.minsamp, 
                        binSize = density.mzwid, maxFeatures = density.max)

                        
x_density <- groupChromPeaks(centWave, param = pdp)

## Obiwarp
#print(paste0("Doing the obiwarp alignment using the default settings...."))
#rt_adjusted <- adjustRtime(x_density, param = ObiwarpParam())

## Loess
print(paste0("Doing the loess method RT alignment using the defined settings...."))
PeakGroupsParam(minFraction = lowess.minfrac, extraPeaks = loess.extra, smooth = loess.smoothing,
               span = loess.span, family = loess.family)

rt_adjusted <-adjustRtime(x_density, param = PeakGroupsParam())

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
print(paste0("Filling peaks..."))
#x_filled <- fillChromPeaks(x_2density, BPPARAM = BPPARAM_fillpeaks)

## convert to xset and correct for missing values
xset <- x_2density
xset <-as(xset, "xcmsSet")

## important! you might want to set/adjust the 'sampclass' of the returned xcmSet object before proceeding with the analysis.
## XCMSnExp saves this as "sampleNames" so we will copy from that
xset$class <- centWave$sampleNames

#run fill peaks on xset, use all cores -2 and with memory allocation from multiple of 4
BPPARAM_fillpeaks <- MulticoreParam(min(detectCores()-2,4), progressbar = TRUE)
xset_fill <- do.call(fillPeaks,list(object = xset,BPPARAM = BPPARAM_fillpeaks))



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

xset_a = annotate(xset_fill,
                  
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


