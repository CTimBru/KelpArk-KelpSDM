rm(list=ls())
require(data.table)
require(sf)
require(raster)
require(stars)
require(terra)
require(dplyr)
require(plyr)
require(dismo)
require(randomForest)
require(DescTools)
require(geodata)
require(ggplot2)
require(virtualspecies)
require(geosphere) #distm
require(viridis)
require(stringr)
source("rvar/var.R")

#Set a predictable random number generator seed for reproducibility.
seed.number <- 42
set.seed(seed.number)

#Set working directory
#RVar_wd should be stored in rvar/var.R
setwd(RVar_wd)

#Specific taxa found in the California study area
California_taxa <- c("nereocystis_luetkeana","macrocystis_pyrifera")

#Specific decade of interest starting year, 20#0s
prediction_years <- c(2020,2030,2040,2050,2060,2070,2080,2090)

#Specific future climate SSP of interest
economic_pathways <- c("ssp126","ssp245","ssp585")

#Intake species occurrence data from GBIF into a data table.  Use the variable name taxon in this process.
nereocystis_occurrence_data <- fread("nereocystis_luetkeana_occurrence.csv",sep="\t")
macrocystis_occurrence_data <- fread("macrocystis_pyrifera_occurrence.csv",sep="\t")

#Split into relevant baselines
macrocystis_occurrence_data <- macrocystis_occurrence_data[which(macrocystis_occurrence_data$year>2009)]
nereocystis_occurrence_data <- nereocystis_occurrence_data[which(nereocystis_occurrence_data$year>2009)]

#Spatially thin occurrence data so that no two points are closer than 9.26 km, the resolution of Bio-Oracle data.
nereocystis_coords <- nereocystis_occurrence_data[, c("decimalLongitude", "decimalLatitude")]
macrocystis_coords <- macrocystis_occurrence_data[, c("decimalLongitude", "decimalLatitude")]

#Create distance matrix to calculate distance between every point as a matrix.
nereocystis_dist_matrix <- distm(nereocystis_coords)
macrocystis_dist_matrix <- distm(macrocystis_coords)

#Create vector that denotes if a coordinate is being kept
nereocystis_keep <- rep(TRUE, nrow(nereocystis_dist_matrix))
macrocystis_keep <- rep(TRUE, nrow(macrocystis_dist_matrix))

#Loop through each coordinate, checking distance. If below 9.26km and not 0km (all points are 0km away from themselves)
for (i in 1:(nrow(nereocystis_dist_matrix)-1)) {
  if (nereocystis_keep[i]) {
    close_points <- which(nereocystis_dist_matrix[i, ] < 9260 & nereocystis_dist_matrix[i, ] > 0)
    nereocystis_keep[close_points] <- FALSE
  }
}
for (i in 1:(nrow(macrocystis_dist_matrix)-1)) {
  if (macrocystis_keep[i]) {
    close_points <- which(macrocystis_dist_matrix[i, ] < 9260 & macrocystis_dist_matrix[i, ] > 0)
    macrocystis_keep[close_points] <- FALSE
  }
}

#Eliminate deselected occurrence data due to proximity.
nereocystis_occurrence_data <- nereocystis_occurrence_data[nereocystis_keep, ]
macrocystis_occurrence_data <- macrocystis_occurrence_data[macrocystis_keep, ]

#Spatially thin occurrence data removing duplicated records at same coords.
nereocystis_coords <- nereocystis_occurrence_data[, c("decimalLongitude", "decimalLatitude")]
macrocystis_coords <- macrocystis_occurrence_data[, c("decimalLongitude", "decimalLatitude")]

#Check if any two coords are exactly equal
#nereocystis_equal_coords <- duplicated(nereocystis_coords)
#macrocystis_equal_coords <- duplicated(macrocystis_coords)

#Switch trues to false, we want to keep non-duplicates
#i <- 1
#for(dupe in nereocystis_equal_coords){
  #if(dupe == TRUE){
  #  dupe <- FALSE
  #} else {
  #  dupe <- TRUE
  #}
  #nereocystis_equal_coords[[i]] <- dupe
  #i <- i+1
#}
#i <- 1
#for(dupe in macrocystis_equal_coords){
  #if(dupe == TRUE){
  #  dupe <- FALSE
  #} else {
  #  dupe <- TRUE
  #}
  #macrocystis_equal_coords[[i]] <- dupe
  #i <- i+1
#}


#Remove duplicates
#nereocystis_occurrence_data <- nereocystis_occurrence_data[nereocystis_equal_coords, ]
#macrocystis_occurrence_data <- macrocystis_occurrence_data[macrocystis_equal_coords, ]

#Convert species data to a spatial points object
nereocystis_species_points_2010 <- st_as_sf(nereocystis_occurrence_data, coords=c("decimalLongitude","decimalLatitude"), crs = 4326)
macrocystis_species_points_2010 <- st_as_sf(macrocystis_occurrence_data, coords=c("decimalLongitude","decimalLatitude"), crs = 4326)

#Get lists of all current environmental rasters in tif format.
list_baseline_2010 <- list.files(path="MapLayers",pattern = "^baseline_2010.*\\.tif$")

#Get lists of all 'permanent' environmental rasters in tif format.
list_static <- c(list.files(path="MapLayers",pattern = "^static_.*\\.tif$"),list.files(path="MapLayers",pattern = "^terrain_.*\\.tif$"))

#Build a raster stack of all environmental rasters.
#Rasters are generated using https://github.com/CTimBru/KelpArk-KelpSDM/blob/main/KelpRaster.R
env_layer_2010 <- stack(paste("MapLayers/",list_baseline_2010,sep=""))
env_layer_static <- stack(paste("MapLayers/",list_static,sep=""))

#For each list of files, get the identifying variable name, and rename layer
i <- 1
for(layer in list_baseline_2010){
  split_name <- strsplit(layer,"_")
  model_depth <- split_name[[1]][[length(split_name[[1]])]]
  model_stat <- split_name[[1]][[length(split_name[[1]])-1]]
  model_var <- split_name[[1]][[length(split_name[[1]])-2]]
  list_baseline_2010[[i]] <- paste(model_var,model_stat,model_depth,sep="_")
  i <- i+1
}

#Update raster stack names for current and future environmental raster stacks so they match the environmental layer name
names(env_layer_2010) <- list_baseline_2010
names(env_layer_static) <- list_static

#Filter collinear environmental variables using current environmental raster stack.
#Use a 0.75 Spearman correlation threshold, and 1000 background points, as per https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.10901
#env_retain <- removeCollinearity(env_layer_2010,method='spearman',multicollinearity.cutoff = 0.75, sample.points = TRUE, nb.points=10000, select.variables=TRUE)

#Select subset of layers
#env_layer_2010 <- subset(env_layer_2010, subset=env_retain)

#Build a raster stack of all current environmental rasters using just the filtered layers.
env_layer_2010 <- stack(env_layer_2010,env_layer_static)

#Extract current environmental raster values at species occurrence points
nereocystis_2010_extracted <- as.data.frame(raster::extract(env_layer_2010, nereocystis_species_points_2010))
macrocystis_2010_extracted <- as.data.frame(raster::extract(env_layer_2010, macrocystis_species_points_2010))

#Remove empty rows from the extracted environmental data
nereocystis_2010_extracted <- as.data.frame(nereocystis_2010_extracted[complete.cases(nereocystis_2010_extracted),])
macrocystis_2010_extracted <- as.data.frame(macrocystis_2010_extracted[complete.cases(macrocystis_2010_extracted),])

#Add presence column to extracted environmental data and set its value to 1.
nereocystis_2010_extracted$presence <- 1
macrocystis_2010_extracted$presence <- 1

#Set presence column to be a factor from being numeric for modeling.
nereocystis_2010_extracted$presence <- as.factor(nereocystis_2010_extracted$presence)
macrocystis_2010_extracted$presence <- as.factor(macrocystis_2010_extracted$presence)

#Set California + Oregon/California Baja Nor area boundaries
Pacific <- st_bbox(c(xmin=-124.5,xmax=-114.1,ymin=22.89,ymax=46.25))

#Determine the number of rows in the extracted environmental data
nrow_nereocystis_2010_extracted <- nrow(nereocystis_2010_extracted)
nrow_macrocystis_2010_extracted <- nrow(macrocystis_2010_extracted)

#Generate a set of random points, three times the size of the number of rows in the extracted environmental data, within points_buffer
nereocystis_2010_background_points <- sf::st_sample(Pacific, size=nrow_nereocystis_2010_extracted*3)
macrocystis_2010_background_points <- sf::st_sample(Pacific, size=nrow_macrocystis_2010_extracted*3)

#Convert the single column coordinates in your background points object to standard longitude/latitude columns
nereocystis_2010_background_points <- sf::st_coordinates(nereocystis_2010_background_points)
macrocystis_2010_background_points <- sf::st_coordinates(macrocystis_2010_background_points)

#Convert background points object to a data table
nereocystis_2010_background_points <- as.data.frame(nereocystis_2010_background_points)
macrocystis_2010_background_points <- as.data.frame(macrocystis_2010_background_points)

#Using the longitude (X) and latitude (Y) columns, convert the background points to a sf object.
#Make sure you're using a coordinate reference system (CRS) of 4326.
nereocystis_2010_background_points <- st_as_sf(nereocystis_2010_background_points,coords=c("X","Y"),crs=4326)
macrocystis_2010_background_points <- st_as_sf(macrocystis_2010_background_points,coords=c("X","Y"),crs=4326)

#Extract current environmental raster values at background points
nereocystis_2010_bg_extracted <- as.data.frame(raster::extract(env_layer_2010, nereocystis_2010_background_points))
macrocystis_2010_bg_extracted <- as.data.frame(raster::extract(env_layer_2010, macrocystis_2010_background_points))

#Remove empty rows from the object with environmental data extracted at background points
nereocystis_2010_bg_extracted <- as.data.frame(nereocystis_2010_bg_extracted[complete.cases(nereocystis_2010_bg_extracted),])
macrocystis_2010_bg_extracted <- as.data.frame(macrocystis_2010_bg_extracted[complete.cases(macrocystis_2010_bg_extracted),])

#Add presence column to the object with environmental data extracted at background points.  Set its value to 0.
nereocystis_2010_bg_extracted$presence <- 0
macrocystis_2010_bg_extracted$presence <- 0

#Set presence variable to factor for modeling.
nereocystis_2010_bg_extracted$presence <- as.factor(nereocystis_2010_bg_extracted$presence)
macrocystis_2010_bg_extracted$presence <- as.factor(macrocystis_2010_bg_extracted$presence)

#Model Selection
selectedTaxa <- California_taxa[[2]]

if(selectedTaxa == California_taxa[[1]]){
  presence_extracted <- nereocystis_2010_extracted
  background_extracted <- nereocystis_2010_bg_extracted
  print(selectedTaxa)
} else if(selectedTaxa == California_taxa[[2]]){
  presence_extracted <- macrocystis_2010_extracted 
  background_extracted <- macrocystis_2010_bg_extracted
}

#Create an empty list to store prediction rasters.
raster_predict_list <- c()
#Create an empty list to store future prediction rasters.
future_raster_predict_list <- c()
#Create an empty list to store relative importance outputs.
importance_list <- c()
#Create an empty list to store accuracy outputs.
accuracy_list <- c()
#Create an empty list to store partial plot outputs.
partial_plot_list <- c()
#Create an empty list to store the trained random forest models
rf1_list <- c()
j <- 1
i <- 1
i_max <- 100
for(i in 1:i_max){
  print(paste("beginning run:",i,sep=" "))
  #Create a subset of the presence/background data with the following properties:
  #1. Composed of a randomly selected 80% of rows from env_extracted.
  #2. Composed of rows randomly selected from background_extracted. The number of rows will also be 80% of rows found in env_extracted.
  #3. Merged these two subsets together.
  subset_extracted <- rbind(presence_extracted[sample(nrow(presence_extracted),0.8*nrow(presence_extracted)),],background_extracted[sample(nrow(background_extracted),0.8*nrow(presence_extracted)),])
  
  #Run a random forest model over this data subset.
  #Use the presence column as the model output, and all other columns as inputs.
  rf1 <- suppressWarnings(tuneRF(x=subset_extracted[,!(colnames(subset_extracted) %in% "presence")],y=subset_extracted$presence,stepFactor=1,plot=FALSE,doBest=TRUE))
  
  #Add the random forest to a list to be able to utilize it in future model loops
  rf1_list[[i]] <- rf1
  
  #Make a prediction raster, using current environmental data, from the random forest model and store it as the ith element in raster_predict_list.
  raster_predict_list[[i]] <- dismo::predict(env_layer_2010,rf1,progress='text')
  
  #Plot prediction raster
  plot(raster_predict_list[[i]])
  
  #Store relative importance of variable outputs in the random forest model in a data frame.
  tmp <- as.data.frame(rf1$importance)
  #Set one column to store the variable names from the row names.
  tmp$VariableName <- rownames(tmp)
  #Store this importance data frame in the ith element of importance_list.
  importance_list[[i]] <- tmp
  
  #Calulate the sensitivity of the random forest model from the confusion matrix.
  sensitivity <- rf1$confusion[[1]] / (rf1$confusion[[1]]+rf1$confusion[[2]])
  
  #Calculate the specificity of the random forest model from the confusion matrix.
  specificity <- rf1$confusion[[4]] / (rf1$confusion[[4]]+rf1$confusion[[3]])
  
  #Calculate the true skill statistic TSS to evaluate model accuracy.
  TSS <- sensitivity+specificity-1
  
  #Store TSS results in the ith element of accuracy_list.
  accuracy_list[i] <- TSS
  
  #Loop through each environmental variable and store the partial response outputs in a temporary data frame.
  for(environ_layer in names(env_layer_2010)){
    #Store partial plot chart data in a temporary data frame.
    tmp <- as.data.frame(partialPlot(rf1,subset_extracted[,!(colnames(subset_extracted) %in% "presence")],x.var=c(environ_layer),plot=F))
    #Transform logistic probabilities to regular probabilities.
    tmp$y <- exp(tmp$y) / (1+exp(tmp$y))
    #Rename probability column
    colnames(tmp) <- c(environ_layer,"Detection_Probability")
    #Store partial plot data in a list of data frames.
    partial_plot_list[[j]] <- tmp
    j <- j+1
  }
}

#Stack the list raster_predict_list into a raster brick.
raster_predict <- brick(raster_predict_list)

#Calculate a raster which is the sum of the layers in this raster brick.  This is the prediction frequency raster.
raster_predict <- calc(raster_predict, sum)

#Save this raster output.  Use the variable taxon in naming the file.
writeRaster(raster_predict,paste("decadalPredictions/",selectedTaxa,"_2010_2020.tif",sep=""),overwrite=T)

# Convert the prediction frequency raster to a data frame
raster_df <- as.data.frame(raster(paste("decadalPredictions/",selectedTaxa,"_2010_2020.tif", sep="")), xy = TRUE)

#Rename the third column of this data frame to value
names(raster_df)[3] <- "value"

#Remove elements from this data frame if the value column is less than or equal to the value of 0.5*i_max.
raster_points <- subset(raster_df, value > 0.5*i_max)

#Plot prediction raster
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "lightblue", high = "lightblue", na.value = "grey") +
  guides(fill = "none") +
  geom_point(data = raster_points, aes(x = x, y = y, color = value), size = 1) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Predicted occurrences of\n",gsub("_"," ",selectedTaxa)," (2020)",sep=""),
       x="Longitude degrees East",y = "Latitude degrees North",
       color = paste("Predicted frequency\nout of ",i_max," models",sep=""))+
  theme(legend.position = "bottom")

#Calculate the mean TSS for the models
TSS_Mean <- mean(accuracy_list)
TSS_SD <- sd(accuracy_list)

TSS <- paste("TSS Mean:",TSS_Mean," TSS Standard Deviation:",TSS_SD,sep="")
write(TSS,paste("ModelStatistics/",selectedTaxa,"_TSS.txt",sep=""))

#Convert importance_list to a single data frame.
importance_total <- rbind.fill(importance_list)

#Calculate the mean relative importance for each variable in this variable importance data frame.
importance_total <- aggregate(x=importance_total$MeanDecreaseGini, by = list(importance_total$VariableName), FUN = mean)

#Rename columns.
colnames(importance_total) <- c("VariableName","Importance")

#Convert variable importance to variable rank importance.  Make 1 correspond to the most important variable.
importance_total$Importance <- rank(desc(importance_total$Importance))

#Save rank importance table.  Use the variable name taxon in naming the file.
write.table(importance_total,paste("ModelStatistics/",selectedTaxa,"_2010_2020_Importance_List.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

#Collapse partial plot outputs into single data frame.
partial_plots <- rbind.fill(partial_plot_list)
partial_plots <- as.data.frame(partial_plots)
write.table(partial_plots,paste("ModelStatistics/",selectedTaxa,"_partial_plots.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
partial_plots <- read.table(paste("ModelStatistics/",selectedTaxa,"_partial_plots.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Plot partial dependence heat maps for continuous data.
k <- 1
ggplot(partial_plots, aes(x=!!sym(names(env_layer_2010[[k]])), y='Detection_Probability') )+
  xlab(names(env_layer_2010[[k]]))+ylab("Detection\nProbability")+
  geom_bin2d(bins = 50)+
  scale_fill_continuous(type = "viridis",name=paste("Frequency\n(Out of ",i_max," models)",sep=""))+
  stat_smooth(aes(y = 'Detection_Probability', fill='Detection_Probability'),method="auto",formula=y~x,color="violet",fill="red",n=0.1*sum(!is.na(partial_plots[,names(env_layer_2010[[k]])])))+
  theme_bw(base_size=25)

#Levi
partial_plots <- rbind.fill(partial_plot_list)
partial_plots <- as.data.frame(partial_plots)
write.table(partial_plots,paste("ModelStatistics/",selectedTaxa,"_partial_plots.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
partial_plots <- read.table(paste("ModelStatistics/",selectedTaxa,"_partial_plots.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Plot partial dependence heat maps for continuous data.
k<- 1
for (k in k:length(names(env_layer_2010))){
  ggplot(partial_plots, aes(x=!!sym(names(env_layer_2010[[k]])), y=`Detection_Probability`) )+
    xlab(names(env_layer_2010[[k]]))+ylab("Detection\nProbability")+
    geom_bin2d(bins = 50)+
    scale_fill_continuous(type = "viridis",name=paste("Frequency\n(Out of ",i_max," models)",sep=""))+
    stat_smooth(aes(y = `Detection_Probability`, fill='Detection_Probability'),method="auto",formula=y~x,color="violet",fill="red",n=0.1*sum(!is.na(partial_plots[,names(env_layer_2010[[k]])])))+
    theme_bw(base_size=25)
  ggsave(filename=paste("ModelStatistics/",selectedTaxa,"_",names(env_layer_2010[[k]]),".png",sep=""),plot=last_plot(),height=1044,width=1760,units=c("px"),dpi=100)
}

#Loop Through All Prediction Years & Pathways
for(prediction_year in prediction_years) {
  for(pathway in economic_pathways){
    print(paste("Model Run:",prediction_year,pathway,sep=" "))
    
    #Get list of selected decades/pathways environmental rasters in fit format.
    list_future <- list.files(path="MapLayers",pattern = paste("^",pathway,"_",prediction_year,"_",prediction_year+10,".*\\.tif$",sep=""))
    
    #Build a raster stack of all environmental rasters.
    #Rasters are generated using https://github.com/CTimBru/KelpArk-KelpSDM/blob/main/KelpRaster.R
    env_layer_future <- stack(paste("MapLayers/",list_future,sep=""))
    
    #For each list of files, get the identifying variable name, and rename layer
    i <- 1
    for(layer in list_future){
      split_name <- strsplit(layer,"_")
      model_depth <- split_name[[1]][[length(split_name[[1]])]]
      model_stat <- split_name[[1]][[length(split_name[[1]])-1]]
      model_var <- split_name[[1]][[length(split_name[[1]])-2]]
      list_future[[i]] <- paste(model_var,model_stat,model_depth,sep="_")
      i <- i+1
    }
    
    #Update raster stack names for current and future environmental raster stacks so they match the environmental layer name
    names(env_layer_future) <- list_future
    
    #Select subset of layers
    #env_layer_future <- subset(env_layer_future, subset=env_retain)
    
    #Build a raster stack of all future environmental rasters using just the filtered layers.
    env_layer_future <- stack(env_layer_future,env_layer_static)
    
    j <- 1
    i <- 1
    i_max <- 100
    #Create an empty list to store future prediction rasters.
    future_raster_predict_list <- c()
    for(i in 1:i_max){
      print(paste("beginning run:",i,sep=" "))
      #Make a future prediction raster from the random forest model, and the future environmental raster stack as input.
      #Store it as the ith element in future_raster_predict_list.
      future_raster_predict_list[[i]] <- dismo::predict(env_layer_future,rf1_list[[i]],progress='text')
    }
    #Stack the list raster_predict_list into a raster brick.
    future_raster_predict <- brick(future_raster_predict_list)
    
    #Calculate a raster which is the sum of the layers in this raster brick.  This is the prediction frequency raster.
    future_raster_predict <- calc(future_raster_predict, sum)
    
    #Save this raster output.  Use the variable taxon in naming the file.
    writeRaster(future_raster_predict,paste("decadalPredictions/",selectedTaxa,"_",prediction_year,"_",prediction_year+10,"_",pathway,".tif",sep=""),overwrite=T)
    
  }
  
  #Plot prediction raster
  
}

# Convert the prediction frequency raster for the future kelp model to data frame

#Set the third column name in this data frame to value

#Remove elements from this data frame if the value column is less than or equal to the value of 0.5*i_max.

#Plot prediction raster

#Get geographic range of predicted taxon occurrences
raster_future_df <- as.data.frame(raster(paste("decadalPredictions/",selectedTaxa,"_2010_2020.tif", sep="")), xy = TRUE)

#Rename the third column of this data frame to value
names(raster_future_df)[3] <- "value"

i_max <- 100
#Remove elements from this data frame if the value column is less than or equal to the value of 0.5*i_max.
raster_future_points <- subset(raster_future_df, value > 0.5*i_max)

ggplot() +
  geom_raster(data = raster_future_df, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "lightblue", high = "lightblue", na.value = "grey") +
  guides(fill = "none") +
  geom_point(data = raster_future_points, aes(x = x, y = y, color = value), size = 1) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Predicted occurrences of\n",gsub("_"," ",selectedTaxa)," (2020)",sep=""),
       x="Longitude degrees East",y = "Latitude degrees North",
       color = paste("Predicted frequency\nout of ",i_max," models",sep=""))+
  theme(legend.position = "bottom")

list_predictions <- list.files(path="decadalPredictions",pattern = "\\.tif$")
list <- strsplit(list_predictions, "_")
name <- c()
year1 <- c()
year2 <- c()
ssp <- c()
ssps <- c()
i<-1
for (i in 1:50) {
  name[i] <-list[[i]][1]
  year1[i] <- list[[i]][3]
  year2[i] <- as.numeric(gsub("([0-9]+).*$", "\\1", list[[i]][4]))
  ssp[i] <- list[[i]][5]
  ssps[i] <- substr(ssp[[i]], 1, 6)
}
range_list <- c()
long_min <- c()
long_max <- c()
lat_min <- c()
lat_max <- c()
pixels <- c()
i <- 1
for (i in 1:length(list_predictions)){
  # Load each raster as a data frame
  raster_df <- as.data.frame(raster(paste("decadalPredictions/",list_predictions[[i]],sep="")), xy = TRUE)
  
  #Set Column to Value
  names(raster_df)[3] <- "value"
  
  #Get geographic range of predicted taxon occurrences
  RangeLong <- range(na.omit(raster_df[raster_df$value == i_max,"x"]))
  RangeLat <- range(na.omit(raster_df[raster_df$value == i_max,"y"]))
  
  #Count the number of locations predicted to have suitable habitat.
  RangePixels <- nrow(na.omit(raster_df[raster_df$value == i_max,]))
  long_min[i] <- RangeLong[[1]]
  long_max[i] <- RangeLong[[2]]
  lat_min[i] <- RangeLat[[1]]
  lat_max[i] <- RangeLat[[2]]
  pixels[i] <- RangePixels[[1]]
  
  range_list[i] <- list_predictions[[i]]
  
}
range_df <- data.frame(range_list, name, year1, year2, ssps, long_min, long_max, lat_min, lat_max, pixels)

range_df <- range_df %>% replace(is.na(.), "Current")

#rbind in duplicates of the 'current data' to be the starting point of each of the modls
current_range_df_dupe <- range_df[range_df$ssps == 'Current',]

#rbind in the new rows twice (once for each model)
i<-1
num_rbind <- length(economic_pathways)
for(i in i:num_rbind){
  print(i)
  current_range_df_dupe$ssps <- economic_pathways[i]
  range_df <- rbind(range_df,current_range_df_dupe)
}

range_df <- range_df[range_df$ssps != 'Current',]

write.csv(range_df,paste("ModelStatistics/species_range.csv",sep=""), row.names = FALSE)

#Range Plots Extents
boundaries <- c("North","South","East","West")
latlong <- c("lat_max","lat_min","long_max","long_min")
j<- 1
for (j in j:length(California_taxa)){
  which_species <- strsplit(California_taxa[j], "_")[[1]]
  species_df <- range_df[range_df$name == which_species[1],]
  k<- 1
  for (k in k:length(latlong)){
    print(k)
    ggplot(species_df, aes(x = year2, y = species_df[,latlong[k]], color = ssps)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) + # Add best-fit lines (linear model, no confidence interval)
      labs(title = paste("Maximal",boundaries[k],"Extent vs. Decade by SSP for",str_to_title(which_species[1]),str_to_title(which_species[2])),
           x = "Decade",
           y = "Degrees",
           color = "SSP") # Customize legend title
    ggsave(filename=paste("ModelStatistics/",California_taxa[j],"_",latlong[k],"_extent_shift",".png",sep=""),plot=last_plot(),height=1044,width=1760,units=c("px"),dpi=100)
  }
}

#Get lists of all current environmental rasters in tif format.
list_refugia <- list.files(path="decadalPredictions",pattern=paste("^",selectedTaxa,".*2090_2100.*\\.tif$",sep=""))

#Build a raster stack of all 2090-2100 predictions
refugia_layer <- stack(paste("decadalPredictions/",list_refugia,sep=""))

#Sum the model outputs of the three SSPs
refugia_predict <- calc(refugia_layer, sum)

#Get geographic range of predicted taxon occurrences
refugia_predict_df <- as.data.frame(refugia_predict, xy = TRUE)

#Rename the third column of this data frame to value
names(refugia_predict_df)[3] <- "value"

refugia_predict_df <- refugia_predict_df[refugia_predict_df$value == 300,]

writeRaster(rasterFromXYZ(refugia_predict_df),paste("refugia/",selectedTaxa,"refugia_2090_2100.tif",sep=""), bylayer=FALSE,overwrite=TRUE)


ggplot() +
  geom_raster(data = refugia_predict_df, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "lightblue", high = "lightblue", na.value = "grey") +
  guides(fill = "none") +
  geom_point(data = refugia_predict_df, aes(x = x, y = y, color = value), size = 1) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Predicted occurrences of\n",gsub("_"," ",selectedTaxa)," (2020)",sep=""),
       x="Longitude degrees East",y = "Latitude degrees North",
       color = paste("Predicted frequency\nout of ",i_max," models",sep=""))+
  theme(legend.position = "bottom")

ggsave(filename=paste("ModelStatistics/",selectedTaxa,"_2090_2100_refugia",".png",sep=""),plot=last_plot(),height=1044,width=1760,units=c("px"),dpi=100)

#Load 15 arc-sec bathometry data
bathymetry_clip_layer <- raster("ncTemp/UncroppedTIFs/bathometry_15_arcsec.tif")

values(bathymetry_clip_layer) <- ifelse(values(bathymetry_clip_layer) <= -100 | values(bathymetry_clip_layer) >= 1, 0, 1)

writeRaster(bathymetry_clip_layer,"MapLayers/bathometry_15_arcsec_crop.tif", bylayer=FALSE,overwrite=TRUE)
