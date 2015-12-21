# Exactive_full_scan_process_ms1+.r
#
# Purpose: Batch convert and process full-scan mass spectrometer data from the Exactive Plus Orbitrap; creates .mzXML files from Thermo .raw files by leveraging the command-line msconvert toolbox
#
# Latest version of script available at https://github.com/vanmooylipidomics/LipidomicsToolbox/
#
# Please direct questions/suggestions/comments to Jamie Collins, james.r.collins@aya.yale.edu, or Helen Fredricks, hfredricks@whoi.edu
#
# Note: This version of the script extracts and processes all (i.e., ms1 and higher) scans
#
################ Revision history #############
#
# Created 8/7/2014 by J.R.C.
# Modified 8/14/2014 by J.R.C.
# Modified 12/12/2014 by J.R.C.
# Mod'fd 1/15/2015 by HFF
#
################ Caveats and prerequisites #############
#
# 1. Unfortunately, this script must be run on a Windows PC since the Thermo .raw conversion tool only runs on Windows
#
# 2. The script presumes you've installed the msconvert command-line tool (part of the ProteoWizard suite)
# and (this is very important) added msconvert to the global path variable so it can be called at the command line.
#
# 3. The script also assumes you've created three directories for the converted files in the master directory where your raw Exactive data is located. We call them:
# r_process_files/two_mode_mzxml_files/ms1, r_process_files/pos_mzxml_files/ms1 and r_process_files/neg_mzxml_files/ms1
# ... but one could easily modify the script for a different directory structure.
#
# 4. This script does check at each step to see if you've already processed the files you're asking it to process. It does this by checking to see if the .mzXML files for the particular sample ID exist already. If for some reason you want to run the script to convert a file again, you'll have to move the existing .mzXML file out of the destination directory first!
#
# 5. For conversion of selected files only, requires a text file "mzxml_convert_list.txt" (see below) with the filenames of your targets
#
# This script contains two code snippets:
#
# 1. The first snippet will batch convert all the .raw files in the Data folder and its subfolders. You probably don't want to do this!
#
# 2. The second snippet will only convert the .raw files you tell it to, using a tab-delimited text file you create that contains the list of .raw filenames

################ Setup and variable definition #############

# load required packages

library(tools)   

# user: define location of the file conversion list, if this is desired

convertlist_file = "mzxml_convert_list.txt"

# user: define location of Thermo .raw files to be converted and where the ms1 .mzXML files should be dumped

data_dir = "/Volumes/Exactive/Data/" # top-level directory (location of .raw files); VML users, change this depending on where you have your data server mapped on your machine

# user: define destination directories where the converted two mode and single ion mode data files will be dumped

twomodemzxmldir_ms1 = paste("/Volumes/Exactive/Data/r_process_files/two_mode_mzxml_files/ms1",sep="")
negmodemzxmldir_ms1 = paste("/Volumes/Exactive/Data/r_process_files/neg_mzxml_files/ms1",sep="")
posmodemzxmldir_ms1 = paste("/Volumes/Exactive/Data/r_process_files/pos_mzxml_files/ms1",sep="")

################# To find and batch convert all .raw files #############
#
## first, look for any previously converted .mzXML files in the target directory already, and compare with list of files to be converted so we don't process files that have already been processed
#
#setwd(twomodemzxmldir_ms1)
#EXISTING_mzXML_FILES = list.files(recursive=TRUE, full.names=FALSE, pattern="\\.mzXML")
#
#setwd(data_dir)
#ALL_RAW_FILES = list.files(recursive=TRUE, full.names=FALSE, pattern="\\.raw")
#POTENTIAL_NEW_FILENAMES = vector(length=length(ALL_RAW_FILES))
#for (i in 1:length(ALL_RAW_FILES))
#{basename = file_path_sans_ext(ALL_RAW_FILES[i])
# POTENTIAL_NEW_FILENAMES[i] = paste(basename,".mzXML",sep="")}
#
## compare the two lists and select only those which haven't been done already
#
#mzXMLs_NOT_YET_CREATED = setdiff(POTENTIAL_NEW_FILENAMES,EXISTING_mzXML_FILES)
#
#RAW_FILES_NOT_YET_CONVERTED = vector(length=length(mzXMLs_NOT_YET_CREATED))
#for (i in 1:length(mzXMLs_NOT_YET_CREATED))
#{basename = file_path_sans_ext(mzXMLs_NOT_YET_CREATED[i])
# RAW_FILES_NOT_YET_CONVERTED[i] = paste(basename,".raw",sep="")}
#
#show(RAW_FILES_NOT_YET_CONVERTED)
#

############### To convert just selected .raw files #############

# read in list of Exactive files to be converted from a text file

setwd(data_dir)

# the text file is called "mzxml_convert_list.txt" and is a tab-delimited text file with the names of the .raw Thermo data files you want to convert
# if using a PC, you should edit this file with Wordpad, not Notepad!
# note that if you want to convert any .raw files that aren't in the top-level Data directory (i.e., they are in a subfolder you've created within Data), you'll need specify the relative file path as well and quotation marks around the entry
# also, you would need to use quotation marks if your filename has any spaces in it
# e.g., "CS TLE ESI/Orbi_1109"

CONVERT_LIST = read.table(convertlist_file)
numfiles = nrow(CONVERT_LIST) # get number of files to be converted
RAW_FILES = matrix(nrow=numfiles, ncol=1) # create a matrix to be populated with filenames + extensions

# append .raw file extension to each file name

for (i in 1:nrow(CONVERT_LIST))
{RAW_FILES[i,1] = paste(CONVERT_LIST[i,1],".raw",sep="")}

# look for any previously converted .mzXML files in the target directory already, and compare with list of files to be converted so we don't process files that have already been processed

setwd(twomodemzxmldir_ms1)
EXISTING_mzXML_FILES = list.files(recursive=TRUE, full.names=FALSE, pattern="\\.mzXML")

POTENTIAL_NEW_FILENAMES = matrix(nrow=numfiles, ncol=1)
for (i in 1:nrow(CONVERT_LIST))
{POTENTIAL_NEW_FILENAMES[i,1] = paste(CONVERT_LIST[i,1],".mzXML",sep="")}

# compare the two lists and select only those which haven't been done already

mzXMLs_NOT_YET_CREATED = setdiff(POTENTIAL_NEW_FILENAMES,EXISTING_mzXML_FILES)

RAW_FILES_NOT_YET_CONVERTED = vector(length=length(mzXMLs_NOT_YET_CREATED))
for (i in 1:length(mzXMLs_NOT_YET_CREATED))
{basename = file_path_sans_ext(mzXMLs_NOT_YET_CREATED[i])
 RAW_FILES_NOT_YET_CONVERTED[i] = paste(basename,".raw",sep="")}

show(RAW_FILES_NOT_YET_CONVERTED)

################ Universally applicable code begins again below #############

################ call msconvert tool to convert files to mzXML #############

# these settings will use peakPicking = 1 to convert data from profile mode to centroid mode, extracting ms1 and higher scans

# change working directory back to top-level data directory

setwd(data_dir)

# run msconvert once to centroid and convert to mzXML

for (i in 1:length(RAW_FILES_NOT_YET_CONVERTED))
{system (paste("msconvert --mzXML --filter \"peakPicking true 1-\" -o ",twomodemzxmldir_ms1," -v \"",RAW_FILES_NOT_YET_CONVERTED[i],"\"",sep=""))}

# look for any previously converted positive mode mzXML files in the target directory already, and compare with list of "full" files to be converted so we don't process files that have already been processed

setwd(posmodemzxmldir_ms1)

EXISTING_pos_mzXML_FILES = list.files(recursive=TRUE, full.names=FALSE, pattern="\\.mzXML")

# compare the two lists and select only those which haven't been done already

NOT_YET_pos_EXTRACTED = setdiff(basename(POTENTIAL_NEW_FILENAMES),EXISTING_pos_mzXML_FILES)

# change working directory back to directory containing "full" .mzXML files

setwd(twomodemzxmldir_ms1)

# run msconvert a second time to extract pos mode data and save

for (i in 1:length(NOT_YET_pos_EXTRACTED))
{system (paste("msconvert --mzXML --filter \"polarity positive\" -o ",posmodemzxmldir_ms1," -v \"",NOT_YET_pos_EXTRACTED[i],"\"",sep=""))}  

# look for any previously converted negative mode mzXML files in the target directory already, and compare with list of "full" files to be converted so we don't process files that have already been processed

setwd(negmodemzxmldir_ms1)

EXISTING_neg_mzXML_FILES = list.files(recursive=TRUE, full.names=FALSE, pattern="\\.mzXML")

# compare the two lists and select only those which haven't been done already

NOT_YET_neg_EXTRACTED = setdiff(basename(POTENTIAL_NEW_FILENAMES),EXISTING_neg_mzXML_FILES)

# change working directory back to directory containing "full" .mzXML files

setwd(twomodemzxmldir_ms1)

# run msconvert a third time to extract neg mode data and save

for (i in 1:length(NOT_YET_neg_EXTRACTED))
{system (paste("msconvert --mzXML --filter \"polarity negative\" -o ",negmodemzxmldir_ms1," -v \"",NOT_YET_neg_EXTRACTED[i],"\"",sep=""))}
