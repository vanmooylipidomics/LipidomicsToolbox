#MSConvert Batch Processing of Thermo.raw mass spectrometry files, for processing with XCMS, CAMERA and LOBSTAHS.
#Jonathan E. Hunter, Woods Hole Oceanographic Institution (https://www.whoi.edu/profile.do?id=jhunter)

#Requires windows and installation of Proteowizard MSConvert (http://proteowizard.sourceforge.net/downloads.shtml)
#Furthermore, MSConvert.exe must be specified in the system PATH (http://www.howtogeek.com/118594/how-to-edit-your-system-path-for-easy-command-line-access/)

#set working directory
setwd("H:/Alina/Microlayer/")

#set raw data directory location
path = ("H:/Alina/Microlayer/")

#prepare msconvert system commands to extract all .raw data files to .mzXML
file.names <- dir(path, pattern =".raw")
msconvcomm <- paste("msconvert ", path, file.names, " --mzXML --filter  \"peakPicking true 1-\" --filter \"scanTime [0,5400]\" -o mzXML_ms1_two_mode -v ", sep="")

#run msconvert system commands
for(i in 1:length(msconvcomm)){
system(msconvcomm[i])
}

#extract positive and negative ms1 data from the .mzXML files
#prepare msconvert system commands
path = "mzXML_ms1_two_mode/"
file.namesmzxml <- dir(path, pattern =".mzXML")
mspolarityseppos <- paste("msconvert ", path, file.namesmzxml, " --mzXML --filter \"polarity positive\" -o mzXML_ms1_pos_mode -v", sep="")
mspolaritysepneg <- paste("msconvert ", path, file.namesmzxml, " --mzXML --filter \"polarity negative\" -o mzXML_ms1_neg_mode -v", sep="")

#run msconvert system commands
for(j in 1:length(mspolarityseppos)){
  system(mspolarityseppos[j])
}

for(k in 1:length(mspolaritysepneg)){
  system(mspolaritysepneg[k])
}

