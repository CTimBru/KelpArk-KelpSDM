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
source("rvar/var.R")

#Set random number string
set.seed(42)

#Set working directory
#RVar_wd should be stored in rvar/var.R
setwd(RVar_wd)

#Specific taxa found in the California study area
California_taxa <- c("nereocystis_luetkeana","macrocystis_pyrifera")

#Intake species occurrence data from GBIF into a data table.  Use the variable name taxon in this process.
nereocystis_occurrence_data <- fread("nereocystis_luetkeana_occurrence.csv",sep="\t")
macrocystis_occurrence_data <- fread("macrocystis_pyrifera_occurrence.csv",sep="\t")


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
#CHECK IF > ZERO IS A PROBLEM!
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
nereocystis_equal_coords <- duplicated(nereocystis_coords)
macrocystis_equal_coords <- duplicated(macrocystis_coords)

#Remove duplicates
nereocystis_occurrence_data <- nereocystis_occurrence_data[nereocystis_equal_coords, ]
macrocystis_occurrence_data <- macrocystis_occurrence_data[macrocystis_equal_coords, ]

#Convert species data to a spatial points object
nereocystis_species_points <- st_as_sf(nereocystis_occurrence_data, coords=c("decimalLongitude","decimalLatitude"), crs = 4326)
macrocystis_species_points <- st_as_sf(macrocystis_occurrence_data, coords=c("decimalLongitude","decimalLatitude"), crs = 4326)


#Get lists of all current and future environmental rasters in tif format.
list_baseline_2000 <- list.files(path="MapLayers",pattern = "^baseline_2000.*\\.tif$")

#Build a raster stack of all environmental rasters.
#Rasters are generated using https://github.com/CTimBru/KelpArk-KelpSDM/blob/main/KelpRaster.R
env_layer <- c()


#Update raster stack names for current and future environmental raster stacks so they match the environmental raster file names


#Filter collinear environmental variables using current environmental raster stack.
#Use a 0.75 Spearman correlation threshold, and 1000 background points, as per https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.10901

#Build a raster stack of all current environmental rasters using just the filtered layers.

#Build a raster stack of all future environmental rasters using just the filtered layers.

#Update raster stack names for current and future environmental raster stacks so they match the environmental raster file names

#Extract current environmental raster values at species occurrence points

#Remove empty rows from the extracted environmental data

#Add presence column to extracted environmental data and set its value to 1.

#Set presence column to be a factor from being numeric for modeling.

#Set California + Oregon/California Baja Nor area boundaries

#Determine the number of rows in the extracted environmental data

#Generate a set of random points, three times the size of the number of rows in the extracted environmental data, within points_buffer
#These are your background points

#Convert the single column coordinates in your background points object to standard longitude/latitude columns

#Convert background points object to a data table

#Using the longitude (X) and latitude (Y) columns, convert the background points to a sf object.
#Make sure you're using a coordinate reference system (CRS) of 4326.

#Extract current environmental raster values at background points


#Update column names of object with environmental data extracted at background points so the column names match the environmental raster file names


#Remove empty rows from the object with environmental data extracted at background points


#Add presence column to the object with environmental data extracted at background points.  Set its value to 0.


#Set presence variable to factor for modeling.


#Set a predictable random number generator seed for reproducibility.


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
j <- 1
i_max <- 100
for(i in 1:i_max){
  #Create a subset of the presence/background data with the following properties:
  #1. Composed of a randomly selected 80% of rows from env_extracted.
  #2. Composed of rows randomly selected from background_extracted. The number of rows will also be 80% of rows found in env_extracted.
  #3. Merged these two subsets together.
  
  
  #Run a random forest model over this data subset.
  #Use the presence column as the model output, and all other columns as inputs.
  
  
  #Make a prediction raster, using current environmental data, from the random forest model and store it as the ith element in raster_predict_list.
  
  #Plot prediction raster
  
  
  #Make a future prediction raster from the random forest model, and the future environmental raster stack as input.
  #Store it as the ith element in future_raster_predict_list.
  
  #Plot prediction raster
  
  
  #Store relative importance of variable outputs in the random forest model in a data frame.
  
  #Set one column to store the variable names from the row names.
  
  #Store this importance data frame in the ith element of importance_list.
  
  
  #Calulate the sensitivity of the random forest model from the confusion matrix.
  
  #Calculate the specificity of the random forest model from the confusion matrix.
  
  #Calculate the true skill statistic TSS to evaluate model accuracy.
  
  #Store TSS results in the ith element of accuracy_list.
  
  #Loop through each environmental variable and store the partial response outputs in a temporary data frame.
  
  j <- j+1
}

#Stack the list raster_predict_list into a raster brick.

#Calculate a raster which is the sum of the layers in this raster brick.  This is the prediction frequency raster.

#Save this raster output.  Use the variable taxon in naming the file.

# Convert the prediction frequency raster to a data frame

#Rename the third column of this data frame to value

#Remove elements from this data frame if the value column is less than or equal to the value of 0.5*i_max.

#Plot prediction raster

#Get geographic range of predicted taxon occurrences

#Count the number of locations predicted to have suitable habitat.

#Stack the list of future prediction rasters into a raster brick.

#Sum the future prediction rasters in the raster brick into a single prediction frequency raster.

#Save raster output.  Use the variable taxon in naming the file

# Convert the prediction frequency raster for the future kelp model to data frame

#Set the third column name in this data frame to value

#Remove elements from this data frame if the value column is less than or equal to the value of 0.5*i_max.

#Plot prediction raster

#Get geographic range of predicted taxon occurrences

#Calculate the mean TSS for the models

#Convert importance_list to a single data frame.

#Calculate the mean relative importance for each variable in this variable importance data frame.

#Rename columns.

#Convert varible importance to variable rank importance.  Make 1 correspond to the most important variable.

#Save rank importance table.  Use the variable name taxon in naming the file.

#Collapse partial plot outputs into single data frame.

#Plot partial dependence heat maps for continuous data.