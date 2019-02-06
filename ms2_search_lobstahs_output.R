# Searching for ms2 via mz/rt window
# 2/1/19

library(xcms)
library(MSnbase)

#load lobstahs output
setwd("C:/Users/TSQ/Desktop/Daniel Lowenstein/Nicole_Coral_Data/Rerun_with_ms2/")
data <- read.csv("Nicole_Rerun_PosMode_LOBSTAHS_screened_peakdata_2019-02-06T10-59-51_AM-0500.csv")

lobs_csv <- data %>% filter(species != "NA")

# reassign centwaves for easy reference
centWave_ms1 <- centWave
centWave <- centWave_ms2

ms1mz <- as.data.frame(precursorMz(centWave))
ms1rt <- as.data.frame(rtime(centWave))

colnames(ms1mz) <- "precursorMz"
colnames(ms1rt) <- "rtime"

#find our detected peaks
Storage <- list()

i <- NULL

for (i in 1:nrow(lobs_csv)) {
  run<-lobs_csv[i,] # pull out one row of LOBSTAHs output

  mz <- run$LOBdbase_mz # pull mz and rt for that one row
  rt <- run$peakgroup_rt
  
  mzrange <- mz*(0.000001*10) # calculate a mass range (in this case, within 10 ppm)
  mzlow <- (mz-mzrange)
  mzhigh <- (mz+mzrange)
  
  rthigh<-rt+30 # and an rt range
  rtlow<-rt-30
  
  # possible ms2 candidates start as a subset of the precursor ions in the mass range
  ms2candid <- subset.data.frame(x = ms1mz,subset = precursorMz>=mzlow & precursorMz<=mzhigh)
  
  # add a column for the retention times by pulling from precursor rts using the rownames
  ms2candid$retention <- ms1rt[rownames(ms2candid),]
  
  # subset that by only taking retention times within the rt window
  ms2matches <- subset.data.frame(ms2candid, subset = retention>=rtlow&retention<=rthigh)
  # create empty row for file names
  ms2matches$file <- rep(0,nrow(ms2matches))
  # add a column with compound name so we can search later
  ms2matches$compound_name <- rep(paste(run$compound_name), nrow(ms2matches))
  
  # if there are ms2 matches within the window, go row by row 
  # and match the fileIdx (file number within the prepOrbidata run)
  # and paste it together with the QE number right before the beginning
  # of the run, so they match up with this run's QE number
  # ****Change QE number!!!!!!!!!!!!!
  if(nrow(ms2matches)>0){
    for (j in 1:nrow(ms2matches)){
      ms2matches[j,"file"] <- paste0("QE00", (6420+as.numeric(centWave@featureData@data[row.names(ms2matches[j,]),"fileIdx"])))
    }
  }
  Storage <- rbind(Storage, ms2matches)
  print(paste0("Searching for matches in ", ms2matches$compound_name[1], " assignment ", i, " of ", nrow(lobs_csv), " total assignments."))
  print(paste0("Matches found:", nrow(ms2matches)))
}
