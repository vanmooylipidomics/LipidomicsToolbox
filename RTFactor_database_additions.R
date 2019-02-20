# Merge RT Factor Dataframes
# 2/19/19

# Match rows from the new ok'ed data to the main Hummel database
# by LOBdbase_mz, and add the new values into new columns



setwd("C:/Users/TSQ/Desktop/Daniel Lowenstein/RT_Factors/")

doc <- read.csv("Updated_Hummel_rtf_dbase.csv")
new <- read.csv("Nicole_RT_Factors.csv")

# going to make it a function...one day. for now, just adding it with a for loop
# add_to_RTF_dbase <- function(Main_RTF_Dbase, New_RTF_mz){
#   if (ncol(New_RTF_mz) != 3){
#     stop("Check your new RTF df. Looks like it doesn't have 3 columns.")
#   }
#   row_number <- which(Main_RTF_Dbase[1] == New_RTF_mz[2])
#   dbase_cols <- length(Main_RTF_Dbase)
#   Main_RTF_Dbase[(dbase_cols + 1)]
#     
# }


# add some new columns
doc$rt_Nicole <- rep("", length(doc$LOBdbase_mz))
doc$DNPPE_rtf_Nicole <- rep("", length(doc$LOBdbase_mz))
doc$class_Nicole <- rep("", length(doc$LOBdbase_mz))
DNPPE_rt <- new$peakgroup_rt[which(grepl("DNPPE", new$compound_name))]

#haven't totally figured out how to use the DNPPE rt that's closer to 930 seconds...
# come back to this later
# if (length(DNPPE_rt) > 1){
#   first_ratio <- 930/DNPPE_rt[1]
#   second_ratio <- 930/DNPPE_rt[2]
#   if first_ratio
# }

for (i in 1:length(new$compound_name)){
  row_number <- which(grepl(paste0("^", new$compound_name[i], "$"), doc$compound_name))
  doc$rt_Nicole[row_number] <- new$peakgroup_rt[i]
  doc$DNPPE_rtf_Nicole[row_number] <- (new$peakgroup_rt[i]/DNPPE_rt)
}
write.csv(doc, "Updated_Hummel_rtf_dbase.csv")

















