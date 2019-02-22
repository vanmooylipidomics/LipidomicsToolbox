#### RT Factor Filtering Function
# Functional Draft 
# Feb. 20, 2019
# input a raw lobset


library(tidyverse)


# where are source files?
setwd("C:/Users/TSQ/Desktop/Daniel Lowenstein/KimT_cleaning/")

# Load data and RT Factor database
# May need to change -X0 or -X1 in source csv
original_data <- read.csv("KimT_LOBSTAHS_screened_peakdata_test.csv")
RT_Factor_Dbase <-read.csv("C:/Users/TSQ/Desktop/Daniel Lowenstein/RT_Factors/Hummel RtF Master Database - rtf_data.csv")

RT_Factor_Sort <- function(original_data, RT_Factor_Dbase){
  
  # make sure peakgroup rt is numeric
  original_data$peakgroup_rt <- as.numeric(original_data$peakgroup_rt)
  
  # Extract correct DNPPE retention time
  DNPPE_RT <- original_data$peakgroup_rt[which(grepl("DNPPE", original_data$compound_name))]
  #DNPPE_RT <- DNPPE_RT[***]  #if there are two DNPPE peaks, index to the correct one
  
  # Add column for DNPPE factor and an empty one for flagging
  original_data <- original_data %>% 
    mutate(DNPPE_Factor = peakgroup_rt/DNPPE_RT, Flag = "None", RTF_Window = NA) %>% 
    filter(species != "NA") 
  
  # isolate major intact polar lipid classes, unoxidized (need to add pigments, etc.)         
  Main_Lipids <- original_data %>%
    filter(degree_oxidation == "0", 
           species == "BLL"| 
           species == "DGCC" |
           species == "DGTS_DGTA" |
           species == "PE" |
           species == "PG" |
           species == "PC" |
           species == "MGDG" |
           species == "DGDG" |
           species == "SQDG"|
           species == "TAG"|
            species == "DAG"|
             species == "FFA"|
             species == "DNPPE")
  
  # isolate oxidized lipids into df
  Ox_Lipids <- original_data %>% 
    filter(degree_oxidation > 0)
  
  # flag known vs unknown by checking whether grepl returns anything in the database
  for (i in 1:length(Main_Lipids$compound_name)){
    which_row <- which(grepl(paste0("^", Main_Lipids$compound_name[i], "$"), RT_Factor_Dbase$compound_name))
  
    if(is.na(RT_Factor_Dbase$Mean_DNPPE_Factor[which_row]) != TRUE ){
      Main_Lipids$Flag[i] = "Known"
      Main_Lipids$RTF_Window[i] = RT_Factor_Dbase$Mean_DNPPE_Factor[which_row]
    }else{
      Main_Lipids$Flag[i] = "Unknown"
    }
    
  }
  
  # separate into two dfs
  Known_RtFs <- Main_Lipids %>% filter(Flag == "Known")
  Unknown_RtFs <- Main_Lipids %>% filter(Flag == "Unknown")
  
  # for each compound in the "Known" df, grab its corresponding row number in RT_Factor_Dbase,
  # then flag it as follows
  # first check if it's in a 10% 
  
  for (i in 1:length(Known_RtFs$compound_name)){
    which_row <- as.numeric(which(grepl(paste0("^", Known_RtFs$compound_name[i], "$"), RT_Factor_Dbase$compound_name)))
    
    
    if(Known_RtFs$DNPPE_Factor[i] < (RT_Factor_Dbase$Mean_DNPPE_Factor[which_row]*1.1) & Known_RtFs$DNPPE_Factor[i] > (RT_Factor_Dbase$Mean_DNPPE_Factor[which_row]*0.9)){
        if(Known_RtFs$DNPPE_Factor[i] < (RT_Factor_Dbase$Mean_DNPPE_Factor[which_row]*1.05) & Known_RtFs$DNPPE_Factor[i] > (RT_Factor_Dbase$Mean_DNPPE_Factor[which_row]*0.95)){
          if(RT_Factor_Dbase$ms2_verified[which_row] == "Yes"){
            Known_RtFs$Flag[i] = "ms2v"
          }else{
            Known_RtFs$Flag[i] = "5%_rtv"
         }
          }else{
          (Known_RtFs$Flag[i] = "10%_rtv")
      }
      }else{
      Known_RtFs$Flag[i] = "Red"
        }
    }
  
  
  
  # combine known and unknown dfs
  Combined <- rbind(Known_RtFs, Unknown_RtFs)
  
  # change levels of colors so they plot in the right order
  Combined$Flag = factor(Combined$Flag, levels = c("Red", "ms2v", "5%_rtv", "10%_rtv", "Unknown"))
  
  Combined <<- Combined
}


RT_Factor_Sort(original_data, RT_Factor_Dbase)

#####################
# From here, giving the option to generate mz vs. rt graphs
# for each of the major classes. 

# add a column for plot labelling by C and DB #
lipidclass <- Combined %>%
  mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+")))

# going through the big nine lipids plus TAGs, DAGs, and FFAs
lipid_classes <- c("DGCC", "DGTS_DGTA", "PC", "PE", "PG", "MGDG", "DGDG", "SQDG", "TAG", "DAG", "FFA")

# and plot each, saving a copy to the working directory
for (i in 1:length(lipid_classes)){
  Lipid <- lipidclass %>% 
    filter(species == paste(lipid_classes[i]))
  
  print(ggplot(Lipid, aes(x = DNPPE_Factor, y = LOBdbase_mz, color =  Flag))+
    geom_point()+
    geom_errorbarh(aes(xmax = RTF_Window*1.1, xmin = RTF_Window*0.9, height = 0.2))+
    geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
    ggtitle(paste0("M/Z vs. RT in ", lipid_classes[i])))
  # ggsave(filename = paste0(lipid_classes[i], "_MZRT_Nicole.tiff"), 
  #        plot = last_plot(),
  #        device = "tiff",
  #        width = 22, height = 17)
  
  }



LOB_lpsolve <- function(LOBpeaklist,choose_class=NULL,use_ms2_RtF) {
  
  library(lpSolve)
  library(ggplot2)
  
  ### Check Inputs ###
  
  if (!class(LOBpeaklist)=="data.frame") {
    
    stop("Input 'LOBpeaklist' is not an 'data.frame' object.\n", 
         "Please use a data.frame generated by 'getLOBpeaklist'.")
    
  }
  
  if (is.null(LOBpeaklist$match_ID)) {
    
    stop("Input data.frame does not contain a 'match_ID' column.", 
         "Please use a data.frame generated by 'getLOBpeaklist'.")
    
  }
  
  if (is.null(LOBpeaklist$compound_name)) {
    
    stop("Input data.frame does not contain a 'compound_name' column.", 
         "Please use a data.frame generated by 'getLOBpeaklist'.")
    
  }
  
  if (is.null(LOBpeaklist$LOBdbase_mz)) {
    
    stop("Input data.frame does not contain a 'LOBdbase_mz' column.", 
         "Please use a data.frame generated by 'getLOBpeaklist'.")
    
  }
  
  if (is.null(LOBpeaklist$peakgroup_rt)) {
    
    stop("Input data.frame does not contain a 'peakgroup_rt' column.", 
         "Please use a data.frame generated by 'getLOBpeaklist'.")
    
  }
  
  if (is.null(LOBpeaklist$FA_total_no_C)) {
    
    stop("Input data.frame does not contain a 'FA_total_no_C' column.", 
         "Please use a data.frame generated by 'getLOBpeaklist'.")
    
  }
  
  if (is.null(LOBpeaklist$FA_total_no_DB)) {
    
    stop("Input data.frame does not contain a 'FA_total_no_DB' column.", 
         "Please use a data.frame generated by 'getLOBpeaklist'.")
    
  }
  
  ### Format our input in a 'run' dataframe
  return <- LOBpeaklist
  LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$degree_oxidation==0),]
  
  if (is.null(choose_class)==FALSE) {
    if(choose_class%in%unique(LOBpeaklist$species)==FALSE){
      stop("Chosen 'choose_class' does not appear in the 'species' column of data.frame.")
    }else{
      LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$species==choose_class),]}
  }else{
    LOBpeaklist <- subset(LOBpeaklist, subset = lipid_class %in% c("IP_DAG","IP_MAG","TAG"))
  }
  # Put what we need in a dataframe 
  PRErun  <- data.frame(LOBpeaklist$match_ID,
                        LOBpeaklist$compound_name,
                        LOBpeaklist$LOBdbase_mz,
                        LOBpeaklist$peakgroup_rt,
                        LOBpeaklist$FA_total_no_C,
                        LOBpeaklist$FA_total_no_DB,
                        LOBpeaklist$species)
  
  #Re-name our column names
  colnames(PRErun) <- c("match_ID",
                        "compound_name",
                        "LOBdbase_mz",
                        "peakgroup_rt",
                        "FA_total_no_C",
                        "FA_total_no_DB",
                        "species")
  
  ### Begin Screening
  
  for (k in 1:length(unique(PRErun$species))) {
    
    run <- PRErun[which(PRErun$species== unique(PRErun$species)[k]),]
    
    #Binary string for each point
    Binary_String <- rep(1, nrow(run)) 
    
    #Empty String we will build a restrictions from
    Empty_String <- rep(0, nrow(run)) 
    
    #Matrix of our Exclusions
    Exclusion_Matrix <- matrix(nrow = 1,ncol = nrow(run))
    
    # Run a loop to find what to exclude for each point
    for (i in 1:nrow(run)) {
      
      #Get our row
      subject <- run[i,]
      
      # Make a table to store our exclusion info
      Exclusion_Table <- run
      Exclusion_Table$Exclude <- rep(FALSE,nrow(run))
      
      #Lets sort the compounds run above and below our point in terms of rt
      lower_rt <- run[which(run$peakgroup_rt < subject$peakgroup_rt),]
      higher_rt <- run[which(run$peakgroup_rt > subject$peakgroup_rt),]
      
      #Now find ones that break the rules for lower and higher and set Exclude to TRUE in the Exclusion_Table
      
      #Exclude the lower
      lower_names <- row.names(lower_rt[lower_rt$FA_total_no_C>=subject$FA_total_no_C & lower_rt$FA_total_no_DB<=subject$FA_total_no_DB,])
      Exclusion_Table[lower_names,"Exclude"] <- TRUE
      #Exclude the higher
      higher_names <- row.names(higher_rt[higher_rt$FA_total_no_C<=subject$FA_total_no_C & higher_rt$FA_total_no_DB>=subject$FA_total_no_DB,])
      Exclusion_Table[higher_names,"Exclude"] <- TRUE
      #Exclude the compounds with the same name
      Exclusion_Table[which(Exclusion_Table$compound_name==subject$compound_name),"Exclude"] <- TRUE
      
      Exclusion_String <- Empty_String
      
      for (j in 1:nrow(run)) {
        if(j!=i){
          if(Exclusion_Table[j,"Exclude"]==TRUE){
            Exclusion_String <- Empty_String
            Exclusion_String[j]<-1
            Exclusion_String[i]<-1
            Exclusion_Matrix <- rbind(Exclusion_Matrix,Exclusion_String)
            rownames(Exclusion_Matrix) <- NULL
            Exclusion_Matrix <- unique(Exclusion_Matrix)
          }
        }
        cat("\r")
        flush.console()
        cat("Writing rules for",as.character(unique(PRErun$species)[k]),"compound number",i,"of",nrow(run),". Number of Rules created:",nrow(Exclusion_Matrix),"...")
      }
    }
    cat(" Done")
    Final_Exclusion_Matrix <- Exclusion_Matrix[-1,]
    
    if(is.null(nrow(Final_Exclusion_Matrix))){
      cat("\nCompound class to small or any lacks noise to screen.")
    }else{
      
      #time to screen
      dir <- rep("<=", nrow(Final_Exclusion_Matrix)) # all constraints '<='
      
      rhs <- rep(1, nrow(Final_Exclusion_Matrix)) # all right hand sides = 1
      
      #all.bin for binary. Set Solution number high to get all solutions.
      
      
      cat("\nApplying lpSolve algorythm...")
      sol <- lpSolve::lp("max", Binary_String, Final_Exclusion_Matrix, dir, rhs,all.bin = TRUE,num.bin.solns = 100) 
      cat(" Done")
      numcols <- nrow(run)
      numsols <- sol$num.bin.solns
      
      solutions <- matrix(head(sol$solution, numcols*numsols), nrow=numsols, byrow=TRUE)
      
      run$Picked <- solutions[1,]
      FINAL <-run[which(run$Picked==1),]
      
      bar <- data.frame(colSums(solutions))
      run$Type<-rep(0,nrow(run))
      run[which(bar==0),"Type"] <- 'No'
      run[which(bar!=nrow(solutions) & bar!=0 ),"Type"] <- 'Maybe'
      run[which(bar==nrow(solutions)),"Type"] <- 'Yes'
      
      print(ggplot(run,aes(x = peakgroup_rt, y = LOBdbase_mz,color=Type)) +
              scale_color_manual(values=c("#e7cd08", "#e70808", "#08e799")) +
              geom_point() +
              geom_text(label=as.character(run$compound_name),nudge_y = 10,size=2,color="black")+
              ggtitle("lpSolve Screened Data") +
              xlab("Peak Group Retention Time (sec)")+
              ylab("Peak Group m/z")
      )
      
      return[return$match_ID %in% run$match_ID,"lpSolve"] <- run$Type
    }
  }
  return(return)
}

LOB_lpsolve(Combined)

