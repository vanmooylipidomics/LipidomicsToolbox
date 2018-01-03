#prepOrbidata_v2.R

library(xcms)

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

# library(Rmpi)

# run two lines below only if IPO hasn't been installed already

# library(devtools)
# install_github("glibiseller/IPO")

library(IPO)

# library(snowfall) # if multicore tasking is desired; seems this line may be unncessary under Windows R

#library(snow) 

## Use socket based parallel processing on Windows systems
if (.Platform$OS.type == "unix") {
    register(bpstart(MulticoreParam()))
} else {
    register(bpstart(SnowParam()))
} 

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

excluded.mzXMLfiles = c("0475","0474") # specifying removal of Orbi_0475.mzXML and Orbi_0474.mzXML since chromatography was screwy, to the point that weird things started to happen when I used retcor() on them

# if planning to use IPO, specify the ID numbers (i.e., Orbi_xxxx.mzXML) of the files you'd like to use for optimization; otherwise, IPO will try to use the entire dataset and you'll probably wind up with errors

# if you aren't planning on running IPO to optimize centWave and/or group/retcor parameters this session, but you have some parameter values from an earlier IPO run saved in a .csv file, you can specify the file paths below. you will still be given the option later to choose which parameters to use.




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


##### PEAK PICKING
print(paste0("Using values of centWave parameters specified by Jamie for peak picking..."))


#read in only msLevel1
rawSpec   <- MSnbase::readMSData(mzXMLfiles, centroided=TRUE, mode="onDisk", msLevel = 1)

#format centwave parameters
cwp <- CentWaveParam(snthresh = 10, 
                     noise = 500, 
                     ppm= 2.5, 
                     mzdiff = 0.005, 
                     prefilter = c(3,7500))

#find peaks
centWave <- findChromPeaks(rawSpec, param = cwp)


print(paste0("Peak picking completed"))

print(centWave)
plotChromPeakImage(centWave) 

print(proc.time() - ptm)

#### GROUPING AND RETCOR

print(paste0("Doing the obiwarp alignment using the default settings...."))

## 
rt_adjusted <- adjustRtime(centWave, param = ObiwarpParam())
xod <- rt_adjusted

## plotting difference between RT adjustment and original results
## Calculate the difference between the adjusted and the raw retention times.
diffRt <- rtime(xod) - rtime(xod, adjusted = FALSE)

## By default, rtime and most other accessor methods return a numeric vector. To
## get the values grouped by sample we have to split this vector by file/sample
diffRt <- split(diffRt, fromFile(xod))

boxplot(diffRt, main = "Obiwarp alignment results", ylab = "adjusted - raw rt") 

#Next we group identified chromatographic peaks across samples. We use the peak density method [@Smith:2006ic] specifying 
#that a chromatographic peak has to be present in at least 2/3 of the samples within each group to be combined to a mz-rt feature.
pdp <- PeakDensityParam(sampleGroups = centWave$sampleNames)

x_density <- groupChromPeaks(rt_adjusted, param = pdp)


#fill peaks
x_filled <- fillChromPeaks(x_density)

#####################################################################################
##### Isotope peak identification, creation of xsAnnotate object using CAMERA #######
#####################################################################################

## convert to xset and correct for missing values
xset <- x_filled
xset <-as(xset, "xcmsSet")

## important! you might want to set/adjust the 'sampclass' of the returned xcmSet object before proceeding with the analysis.
## XCMSnExp saves this as "sampleNames" so we will copy from that
xset$class <- centWave$sampleNames

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

print(proc.time() - ptm)


