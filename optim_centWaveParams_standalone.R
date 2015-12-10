# optim_centWaveParams_standalone.R
#
# Created 11/28/2015 by J.R.C.
#
# Purpose: "Standalone" script to optimize parameters for the xcms command findPeaks.centWave, using the R package "IPO." Currently, the script is written to optimize parameters for peak-picking, alignment, etc., of lipid data from the experiment described in Graff van Creveld et al., 2015, "Early perturbation in mitochondria redox homeostasis in response to environmental stress predicts cell fate in diatoms," ISME Journal 9:385-395. This dataset is used to demonstrate the LOBSTAHS lipidomics pipeline in Collins, J.R., B.R. Edwards, H.F. Fredricks, and B.A.S. Van Mooy, 2015, "Untargeted discovery and identification of oxidative stress biomarkers using a lipidomics pipeline for complex datasets."
#
# This script still functions on its own, but the functionality contained in it has been folded into prepOrbidata.R
#
# IPO is described in Libiseller et al., 2015, "IPO: a tool for automated optimization of XCMS parameters," BMC Bioinformatics 16:118; see https://github.com/glibiseller/IPO/blob/master/vignettes/IPO.Rmd for installation instructions
#
# See https://github.com/vanmooylipidomics/LOBSTAHS for current versions of all pipeline scripts

################ Initial setup and variable definition #############

# load required packages

library(tools) 

library(xcms)

library(CAMERA)

library(rsm)

# run two lines below only if IPO hasn't been installed already

# library(devtools)
# install_github("glibiseller/IPO") 

library(IPO)

library(snowfall) # if multicore tasking is desired

################# User: define locations of data files and database(s) #############

working_dir = "/Users/jrcollins/Dropbox/code/LOBSTAHS/" # specify working directory
setwd(working_dir) # set working directory to working_dir

# specify directories subordinate to the working directory in which the .mzXML files for xcms can be found; per xcms documentation, use subdirectories within these to divide files according to treatment/primary environmental variable (e.g., station number along a cruise transect) and file names to indicate timepoint/secondary environmental variable (e.g., depth)
mzXMLfiles_folder_pos = "Pt_H2O2_mzXML_ms1_pos/" 
mzXMLfiles_folder_neg = "Pt_H2O2_mzXML_ms1_neg/" 

################# Load in mzXML files #############

mzXMLfiles = list.files(mzXMLfiles_folder_pos, recursive = TRUE, full.names = TRUE)

# # exclude any files you don't want to push through xcms (e.g., blanks); note that the blanks for the Pt H2O2 dataset (Orbi_0481.mzXML and Orbi_0482.mzXML) have already been removed
# mzXMLfiles = mzXMLfiles[-c(1,2)]

################# Use IPO to optimize some findPeaks.centWave parameters #############

# will use IPO to optimize settings for method = centWave

# define ranges of parameters to be tested
# if single value is specified for a parameter, or centWave default is used, that parameter will not be optimized

peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')
peakpickingParameters$min_peakwidth <- c(10,20) # centerpoint is 15
peakpickingParameters$max_peakwidth <- c(40,80) # centerpoint is 60
peakpickingParameters$ppm <- 2.5 # want to set this low to avoid peak data insertion errors from centWave; IPO wants to use something like 5.5 ppm if you allow it to "optimize," but this is too high
peakpickingParameters$prefilter <- 3 # a very long optimization routine settled on a value of 2.4
peakpickingParameters$value_of_prefilter <- c(1000,10000)
peakpickingParameters$snthresh <- c(10)
peakpickingParameters$noise <- c(500)

# only going to use 4 0 uM H2O2 treatment files from the dataset for optimization routine
# seems that mzXMLfiles[c(2,5)] throw up the error:
#
#       Error in checkForRemoteErrors(val) : 
#        ... nodes produced an error
#
# so, will use the other four files (mzXMLfiles[c(1,3,4,6)])

# # code to diagnose which files contain scans that are causing the package to give error

# xcmsraw = xcmsRaw(mzXMLfiles[5])
# scan121=getScan(xcmsraw,scan=107)
# scan121[2429,]

resultPeakpicking <- optimizeXcmsSet(files= mzXMLfiles[c(1,3,4,6)], 
                                     params=peakpickingParameters, nSlaves=4, subdir='rsmDirectory')
optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset

################# Export IPO starting value(s) and optimal settings for each parameter to .csv #############

# generate unique timestamp for filename so we don't overwrite any existing output

output_DTG = format(Sys.time(), "%Y-%m-%dT%X%z") # return current time in a good format
output_DTG = gsub(" ", "_", output_DTG) # replace any spaces
output_DTG = gsub(":", "-", output_DTG) # replaces any colons with dashes (Mac compatibility)

# write 3-column table to .csv using write.table()

write.table(cbind(sort(rownames(as.matrix(peakpickingParameters))),
                  as.character(resultPeakpicking$best_settings$parameters[sort(names(resultPeakpicking$best_settings$parameters))]),
                  as.character(peakpickingParameters[sort(rownames(as.matrix(peakpickingParameters)))])),
            file = paste("IPO_centWaveparamfits_",output_DTG,".csv",sep=""),
            col.names = c("centWave_parameter","IPO_optim_value","Starting_value(s)"),
            row.names = FALSE,
            sep=",")