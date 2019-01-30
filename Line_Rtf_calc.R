## Plotting MZ vs RtF clusters

library(ggplot2)

#import our data
RT_Filtering <- read.csv("~/Desktop/RT_Filtering_db_R.csv")
data <- RT_Filtering
#get only rows we have a RtF for
data <- data[which(is.na(data$peakgroup_rt)==FALSE),]

data <- data[data$species=="PE",] #Change this to do a different compound for diff graph

#For Just two plots 
RtF_equation <- function(data,class,plot_text,plot_connections,plot_subplots,strictness){
  
  x <- data$DNPPE_Factor
  y <- data$mz
  ID <- data$frag_ID
  name <- data$compound_name
  CarbNum <- data$FA_total_no_C
  dbNum <- data$FA_total_no_DB
  
  crosslinkframe <- data.frame(x,y,ID,name,CarbNum,dbNum)
  
  # get range of carbons and db
  ###for C
  C_range <- range(CarbNum,na.rm=T)
  C_rangelow <- C_range[[1]]
  C_rangehigh <- C_range[[2]]
  ###for db
  db_range <- range(dbNum,na.rm=T)
  db_rangelow <- db_range[[1]]
  db_rangehigh <- db_range[[2]]
  
  ### make a storage lists for df seperated by C and by db
  C_storage <- list()
  db_storage <-list()
  
  #Loop it Loop it Loop it
  ###seperated by C number
  for (i in C_rangelow:C_rangehigh) {
    result <- data[which(CarbNum==i),]
    C_storage[[(i-C_rangelow)+1]] <- result
  }
  
  ###seperated by db number
  for (i in db_rangelow:db_rangehigh) {
    result <- data[which(dbNum==i),]
    db_storage[[(i-db_rangelow)+1]] <- result
  }
  return <- list()
  return[["db"]] <- db_storage
  return[["C"]] <- C_storage
  
  ### plus 1 db plus 1 carbon delta
  # reference DB creation
  reference <- list()
  
  for (i in 1:nrow(crosslinkframe)){
    #high part
    db_high <- crosslinkframe[i,]$dbNum:30
    Carb_high <- crosslinkframe[i,]$CarbNum:(crosslinkframe[i,]$CarbNum+length(db_high)-1)
    run <- data.frame(carb=Carb_high,db=db_high)
    
    #low part
    db_low <- 0:crosslinkframe[i,]$dbNum
    Carb_low <- (crosslinkframe[i,]$CarbNum-length(db_low)+1):crosslinkframe[i,]$CarbNum
    run2 <- data.frame(carb=Carb_low,db=db_low)
    
    #combo
    reference[[i]] <- rbind(run2,run)
    names(reference)[[i]] <- print(reference[[i]][1,]$carb)
  }
  
  i <- NULL
  names <- unique(names(reference))
  reference_filterd <- list()
  for (i in 1:length(names)) {
    reference_filterd[i]<- reference[names[i]]
    names(reference_filterd)[[i]]<-names[i]
  }
  
  i <- NULL
  for (i in 1:length(reference_filterd)) {
    run<-reference_filterd[[i]]
    run$ref<-rep(i,length(reference_filterd[[i]]))
    if (i==1){
      final<-run
      }else{
        final<-rbind(final,run)}
    }


final <- unique(final)

#time to screen and sort
C_db_storage <- list()
data_ref <- data
for (i in 1:nrow(data_ref)) {
  run <- final[which(data_ref[i,]$FA_total_no_C==final$carb & data_ref[i,]$FA_total_no_DB==final$db),]
  run <- run$ref
  data_ref[i,"REF"] <- run
}
for (i in 1:length(reference_filterd)) {
  result <- data_ref[which(data_ref$REF==i),]
  C_db_storage[[i]] <- result
}
return[["C_db"]] <- C_db_storage

mod_storage <- list()
i <-NULL
  for (i in 1:length(return[["db"]])) {
    if (nrow(return[["db"]][[i]])>strictness) {
   
    #get our data
    x <- return[["db"]][[i]]$DNPPE_Factor
    y <- return[["db"]][[i]]$mz
    ID <- return[["db"]][[i]]$frag_ID
    name <- return[["db"]][[i]]$compound_name

    if (plot_subplots==TRUE) {
      
    print(ggplot(mapping = aes(x=x, y=y))+
      geom_point()+
      geom_text(label=name,size=4)+
      stat_smooth(method = "lm",formula = y ~ x,fullrange = TRUE))
    
    }
      
    mod <- lm(y ~ x)
    smooth_vals = predict(object=mod,newdata = data.frame(x=seq(0,2,0.01)))
    smooth_vals_list <- list(smooth_vals)
    mod_storage[[i]]<-list(mod,smooth_vals)
    
  }
  }
  
  x <- crosslinkframe$x
  y <- crosslinkframe$y
  name <- crosslinkframe$name
  
  print(
    final_plot <- ggplot(mapping = aes(x = x, y = y,xlab= "Rt in Sec",ylab = "mz"))+
      geom_point()+
      geom_text(label=paste(name)))
  
i <- NULL
lineFrame <- NULL
  for (i in 1:length(mod_storage)) {
    if(length(mod_storage[[i]])!=0){
      if (is.null(lineFrame)) {
        lineFrame <- data.frame(mod_storage[[i]][[2]])
        colnames(lineFrame) <- "values"
        lineFrame$type <- rep(i,nrow(lineFrame))
        lineFrame$x_val <- seq(0,2,0.01)
      }else{
        i_frame <- data.frame(mod_storage[[i]][[2]])
        i_frame$type <- rep(i,nrow(i_frame))
        i_frame$x_val <- seq(0,2,0.01)
        colnames(i_frame) <- c("values","type","x_val")
        lineFrame <- rbind(lineFrame,i_frame)
      }
    }
  }
  
print(final_plot + geom_line(data = lineFrame,mapping = aes(x = x_val,y = values, group=type,color="Delta Carbon"))+
                    xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Change in Carbon Line"))

final2<- final_plot + geom_line(data = lineFrame,mapping = aes(x = x_val,y = values, group=type,color="Delta Carbon"))+
  xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Change in Carbon Line")

mod_storage2 <- list()
i <-NULL
for (i in 1:length(return[["C"]])) {
  if (nrow(return[["C"]][[i]])>strictness) {
    
    #get our data
    x <- return[["C"]][[i]]$DNPPE_Factor
    y <- return[["C"]][[i]]$mz
    ID <- return[["C"]][[i]]$frag_ID
    name <- return[["C"]][[i]]$compound_name
    
    if (plot_subplots==TRUE) {
      
      print(ggplot(mapping = aes(x=x, y=y))+
              geom_point()+
              geom_text(label=name,size=4)+
              stat_smooth(method = "lm",formula = y ~ x,fullrange = TRUE))
      
    }
    
    mod <- lm(y ~ x)
    smooth_vals = predict(object=mod,newdata = data.frame(x=seq(0,2,0.01)))
    smooth_vals_list <- list(smooth_vals)
    mod_storage2[[i]]<-list(mod,smooth_vals)
    
  }
}

i <- NULL
lineFrame2 <- NULL
for (i in 1:length(mod_storage2)) {
  if(length(mod_storage2[[i]])!=0){
    if (is.null(lineFrame2)) {
      lineFrame2 <- data.frame(mod_storage2[[i]][[2]])
      colnames(lineFrame2) <- "values"
      lineFrame2$type <- rep(i,nrow(lineFrame2))
      lineFrame2$x_val <- seq(0,2,0.01)
    }else{
      i_frame <- data.frame(mod_storage2[[i]][[2]])
      i_frame$type <- rep(i,nrow(i_frame))
      i_frame$x_val <- seq(0,2,0.01)
      colnames(i_frame) <- c("values","type","x_val")
      lineFrame2 <- rbind(lineFrame2,i_frame)
    }
  }
}
x <- crosslinkframe$x
y <- crosslinkframe$y
name <- crosslinkframe$name

print(final_plot + geom_line(data = lineFrame2,mapping = aes(x = x_val,y = values, group=type,color="Delta DB"))+
        xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Change in DB Line"))

print(final2 + geom_line(data = lineFrame2,mapping = aes(x = x_val,y = values, group=type,color="Delta DB"))+
        xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Combinded Plots"))

final3 <- final2 + geom_line(data = lineFrame2,mapping = aes(x = x_val,y = values, group=type,color="Delta DB"))+
  xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Combinded Plots")

#######
######
####

mod_storage3 <- list()
i <-NULL
x <- crosslinkframe$x
y <- crosslinkframe$y
name <- crosslinkframe$name
for (i in 1:length(return[["C_db"]])) {
  if (nrow(return[["C_db"]][[i]])>strictness) {
    
    #get our data
    x <- return[["C_db"]][[i]]$DNPPE_Factor
    y <- return[["C_db"]][[i]]$mz
    ID <- return[["C_db"]][[i]]$frag_ID
    name <- return[["C_db"]][[i]]$compound_name
    
    if (plot_subplots==TRUE) {
      
      print(ggplot(mapping = aes(x=x, y=y))+
              geom_point()+
              geom_text(label=name,size=4)+
              stat_smooth(method = "lm",formula = y ~ x,fullrange = TRUE))
      
    }
    
    mod <- lm(y ~ x)
    smooth_vals = predict(object=mod,newdata = data.frame(x=seq(0,2,0.01)))
    smooth_vals_list <- list(smooth_vals)
    mod_storage3[[i]]<-list(mod,smooth_vals)
    
  }
}

x <- crosslinkframe$x
y <- crosslinkframe$y
name <- crosslinkframe$name

i <- NULL
lineFrame3 <- NULL
for (i in 1:length(mod_storage3)) {
  if(length(mod_storage3[[i]])!=0){
    if (is.null(lineFrame3)) {
      lineFrame3 <- data.frame(mod_storage3[[i]][[2]])
      colnames(lineFrame3) <- "values"
      lineFrame3$type <- rep(i,nrow(lineFrame3))
      lineFrame3$x_val <- seq(0,2,0.01)
    }else{
      i_frame <- data.frame(mod_storage3[[i]][[2]])
      i_frame$type <- rep(i,nrow(i_frame))
      i_frame$x_val <- seq(0,2,0.01)
      colnames(i_frame) <- c("values","type","x_val")
      lineFrame3 <- rbind(lineFrame3,i_frame)
    }
  }
}

print(final_plot + geom_line(data = lineFrame3,mapping = aes(x = x_val,y = values, group=type,color="Delta Carbon and DB"))+
        xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Change in Carbon and DB Line"))

print(final3 + geom_line(data = lineFrame3,mapping = aes(x = x_val,y = values, group=type,color="Delta Carbon and DB"))+
        xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Change in Carbon and DB Line"))
Doneee <-final3 + geom_line(data = lineFrame3,mapping = aes(x = x_val,y = values, group=type,color="Delta Carbon and DB"))+
  xlim(min(x),max(x))+ylim(min(y),max(y))+ggtitle("Change in Carbon and DB Line")

}

RtF_equation(data = data,class = "NA",plot_text = TRUE,plot_connections = TRUE,plot_subplots = FALSE,strictness = 2)


data$frag_ID <- row.names(data)
data$DNPPE_factor <- data$DNPPE_rtf_dan
