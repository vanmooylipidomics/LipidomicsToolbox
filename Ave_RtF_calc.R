# Make sure we have what we need

library(ggplot2)
library(ggrepel)
library(grid)

## Dont need right now
# Raw Data including ms2
#rawSpec <- rawSpec

# A peak file with peak detection done on it AND MS2!
centWave <- centWave

# A LOBset from both of those
#LOBset <- LOBset
  
# Properly formated list of RtF estimates: Must include columns: "peakgroup_rt","compound_name","peakgroup_mz"
#estimates <- LOBSTAHS_screened_peakdata_2018.09.01T12.29.18.0400
estimates

################################################################

#Define our helper functions 
plot_table_compound <- function(mz,rtlow,rthigh,ppm,xclname) {

  seclow <- (rtlow*60)
  sechigh <- (rthigh*60)
  
  #create our m/z range based on the ppm
  mzrange <- mz*(0.000001*ppm)
  mzlow <- (mz-mzrange)
  mzhigh <- (mz+mzrange)
  
  #make a data frame of our sample names without string
  
  samplenames <- gsub(chosenFileSubset,"",mzXMLfiles)
  
  #get rid of the "/" if there is one
  
  samplenamesnoslash <- gsub("/","",samplenames)
  
  #number each sample in the dataframe
  
  samplenamesframe <- data.frame(samplenamesnoslash,samplenumber =
                                   seq(from=1, to=length(mzXMLfiles)))
  
  #create + extract a lists of peaks that fit our parameters 
  peaks <- chromPeaks(object = centWave, 
                      mz = c(mzlow, mzhigh),
                      rt = c(seclow, sechigh))
  
  if (nrow(peaks)==0) {
    return("No peaks found")
  }
  #turn our matrix into a dataframe
  peaksframe <- as.data.frame(peaks)
  
  #pull out the columns we want
  peaksnumber <- peaksframe[["sample"]]
  peaksmz <- peaksframe[["mz"]]
  peaksrt <- peaksframe[["rt"]]
  peaksintensity <-format(peaksframe[["into"]], scientific = TRUE)
  peaksintensitynon <-format(peaksframe[["into"]])
  
  #Special Time conversions
  rtminutes <- floor(peaksrt/60)
  rtsecs <- format((peaksrt%%60),digits = 0 )
  rtminsec <- paste(rtminutes,"m",rtsecs,"s", sep = " ")
  
  #make them into another dataframe
  samplevalues <- data.frame(name = peaksnumber,
                             mz = peaksmz,
                             rt = rtminsec,
                             rtsecond = peaksrt,
                             intensity = peaksintensity,
                             intensitynon = peaksintensitynon)
  
  #Add the sample names back in
  merged <- merge(samplevalues, samplenamesframe, by.x="name", by.y= "samplenumber")
  
  #Reorder our coulmns so sample name comes seconds
  reordered <- merged[c(1,7,2,3,4,5,6)]
  
  #Added this warning message so we dont break anything if we dont find any peaks
  if(is.na(reordered[1,"name"])== TRUE){
    print("No peaks found in centWave for current settings.")
  }else{
    
    #Make everything a character so we can add a page break in <- I made switch for this but dont like it so i commented it out
    
    #if(input$XYZ == TRUE){
    #ascharacters <- as.data.frame(lapply(reordered, as.character), stringsAsFactors = FALSE)
    
    #Done <- head(do.call(rbind, by(ascharacters, reordered$name, rbind, "")), -1 )
    # }else{
    
    Done <- reordered
    # }
    
    #give our table nice names
    colnames(Done)[1] <- "Sample Number"
    colnames(Done)[2] <- "Sample Name"
    colnames(Done)[3] <- "m/z"
    colnames(Done)[4] <- "Retention Time"
    colnames(Done)[5] <- "Retention Time Sec"
    colnames(Done)[6] <- "Intensity Scientific"
    colnames(Done)[7] <- "Intensity"
  
  write.table(Done, file = paste("~/Desktop/RtF_data/RtF/tables/",xclname,".txt",sep = ""), sep="\t")

}
return(Done)
}
plot_graph <- function(mz,rtlow,rthigh,ppm,xclname){
  
  seclow <- (rtlow*60) 
  sechigh <- (rthigh*60)
  
  #create our m/z range based on the ppm
  mzrange <- mz*(0.000001*ppm)
  mzlow <- (mz-mzrange)
  mzhigh <- (mz+mzrange)
  
  chroms <- chromatogram(object = centWave, mz = c(mzlow,mzhigh),rt = c(seclow, sechigh))
  return(plot(chroms,main=xclname))
}

#Lets find DNPPE
DNPPE <- plot_table_compound(mz = 875.550487,rtlow = 14,rthigh = 17,ppm = 2.5,xclname = "DNPPE")
plot_graph(mz = 875.550487,rtlow = 14,rthigh = 17,ppm = 2.5,xclname = "DNPPE")

ms1mz <- as.data.frame(precursorMz(centWave))
ms1rt <- as.data.frame(rtime(centWave))
colnames(ms1mz) <- "precursorMz"
colnames(ms1rt) <- "rtime"

#find our detected peaks
Storage <- list()

i <- NULL
test<- list()
for (i in 1:nrow(estimates)) {
  run<-estimates[i,]
  high <-run$rt_dan+30
  low <- run$rt_dan-30
  
  high <- high/60
  low <- low/60
  
  name <- as.character(run$compound_name[1])
    
  mz <- run$mz
  rt <- run$rt_dan
    
    mzrange <- mz*(0.000001*5)
    mzlow <- (mz-mzrange)
    mzhigh <- (mz+mzrange)
    
    rthigh<-rt+30
    rtlow<-rt-30
    
    ms2candid <- subset.data.frame(x = ms1mz,subset = precursorMz>=mzlow & precursorMz<=mzhigh)
    
    ms2candid$retention <- ms1rt[rownames(ms2candid),]
    
    ms2matchs <- subset.data.frame(ms2candid, subset = retention>=rtlow&retention<=rthigh)
    
    ms2matchs$file <- rep(0,nrow(ms2matchs))
    
    if(nrow(ms2matchs)>0){
    for (j in 1:nrow(ms2matchs)){
    ms2matchs[j,"file"] <- centWave@featureData@data[row.names(ms2matchs[j,]),"fileIdx"]
    }
      }
  Store <- plot_table_compound(mz = run$mz,rtlow = low,rthigh = high,ppm = 2.5,xclname = name)
  if (class(Store)=="character") {
    Storage[[i]]<- list(Store,ms2matchs)
    names(Storage)[[i]] <- name
  }else{
  Store$estimate_rt <- rep(estimates[i,"rt_dan"],nrow(Store))
  Store$estimate_RtF <- rep(estimates[i,"DNPPE_rtf_dan"],nrow(Store))
  Store$ppm_error <- ((Store$`m/z`-run$mz)/run$mz)*10^6
  Storage[[i]]<- list(Store,ms2matchs)
  names(Storage)[[i]] <- name
  }
}
#see what those chroms look like (off for now kinda buggy)

# i <- NULL
# for (i in 1:nrow(estimates)) {
#   run<-estimates[i,]
#   high <-run$peakgroup_rt+75
#   low <- run$peakgroup_rt-75
#   
#   high <- high/60
#   low <- low/60
#   
#   name <- as.character(run$compound_name[1])
#   
#   setwd("~/Desktop/RtF_data/RtF/plots/")
#   
#   jpeg(filename = name)
#   
#   plot_graph(mz = run$peakgroup_mz,rtlow = low,rthigh = high,ppm = 2.5,xclname = name)
#   
#   dev.off()
#   
# }

#Make some plots
i <- NULL

RtFactors <- data.frame(matrix(nrow = length(Storage),ncol = 4))
colnames(RtFactors) <- c("Name", "RtFmean","Stdev","mass")

for (i in 1:length(Storage)) {
  
  run <- Storage[[i]][[1]]
  
  if (class(run)=="data.frame") {
    if(nrow(run)>=1){
  
  run$DNPPE_name <- DNPPE$`Sample Name`[match(run$`Sample Name`,DNPPE$`Sample Name`)]
  run$DNPPE_rt_sec <- DNPPE$`Retention Time Sec`[match(run$`Sample Name`,DNPPE$`Sample Name`)]
  
  run$RtF <- run$`Retention Time Sec`/run$DNPPE_rt_sec
  
  doops <-run[duplicated(run$`Sample Number`),]
  
  doop_Nums <- doops$`Sample Number`
  
  run$duplicated <- rep(x = FALSE,nrow(run))
  
  run$is_there_ms2 <- rep(x = FALSE,nrow(run))
  
  y<-NULL
  for (y in 1:nrow(run)) {
    if (run[y,"Sample Number"]%in%doop_Nums) {
    run[y,"duplicated"]<- TRUE 
    }
  }
  
  y<-NULL
  for (y in 1:nrow(run)) {
    if (run[y,"Sample Number"]%in%Storage[[i]][[2]]$file) {
      run[y,"is_there_ms2"]<- TRUE 
    }
  }
  
  runGOOD <- run[run$duplicated==FALSE,]
  runms2 <- run[run$is_there_ms2==TRUE,]
  
  RtFactors[i,"Name"]<- names(Storage[i])
  RtFactors[i,"RtFmean"]<- mean(runGOOD$RtF)
  RtFactors[i,"Stdev"]<- sd(runGOOD$RtF)
  RtFactors[i,"mass"]<-mean(runGOOD$`m/z`)
  RtFactors[i,"Number_of_Files"]<-nrow(runGOOD)
  
  plot1<- ggplot(run, aes(x=`Sample Number`, y=`Retention Time Sec`, group=duplicated)) +
    geom_line(linetype="dotted")+
    geom_point(aes(color=duplicated))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "none")+
    ggtitle(names(Storage[i]))+
    geom_hline(yintercept=run[1,"estimate_rt"], linetype="dashed", color = "red")
  
  plot2<- ggplot(run, aes(x=`Sample Number`, y=DNPPE_rt_sec, group=1)) +
    geom_line(linetype="dotted")+
    geom_point()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
  plot3<- ggplot(run, aes(x=`Sample Number`, y=RtF, group=duplicated)) +
    geom_line(linetype="dotted")+
    geom_point(aes(color=is_there_ms2,size=abs(ppm_error),shape=duplicated))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = "bottom")+
    geom_hline(yintercept=run[1,"estimate_RtF"], linetype="dashed", color = "red")+
    geom_text_repel(data = runms2,mapping = aes(label=run[which(run$is_there_ms2==TRUE,),"Sample Name"]),nudge_y = .005,size=3,segment.alpha = .3)
  
  setwd("~/Desktop/RtF_data/RtF/plots/")
  
  jpeg(filename = names(Storage[i]),width = 1000,height = 600)
  
  grid.arrange(plot1,plot2,plot3,ncol=1)
  
  dev.off()
  }}
}

write.csv(x = RtFactors,file = "PE_RtF.csv")

graphms2 <-function(raw,ms2,compound,sample_num){
  list <- ms2[[compound]][2]
  list <- list[[1]]
  list <- list[which(list$file==sample_num),]
  
  result <- vector("list",length(row.names(list[[1]])))
  
  for (i in 1:length(row.names(list[[1]]))) {
    
    sp <-raw[[row.names(list)[i]]]
    
    splabel <-round(sp@mz,digits = 2)
    
    sp@intensity
    
    sp@precursorMz
    
   # jpeg(filename =paste("ms2 mz",round(sp@precursorMz,2),"_RT",round(sp@rt,1),"seconds.jpeg"),
     #    width = 1000,
     #    height = 700)
    
    plots <- print(plot(sp, full=TRUE,centroided=TRUE)+
                     geom_text(aes(label=splabel),hjust=0,nudge_x = 1,check_overlap = TRUE)+
                     ggtitle(label=paste("ms2: Precursor of",round(sp@precursorMz,5),"m/z and RT of",round(sp@rt,1),"seconds")))
    #dev.off()
  }}

graphms2(raw = centWave,
         ms2 = Storage,
         compound = "PE 34:2",sample_num = 35)
