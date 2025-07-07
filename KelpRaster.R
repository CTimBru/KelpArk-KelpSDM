rm(list=ls())
require(data.table)
require(raster)
require(sf)
require(terra)
require(devtools)
#devtools::install_github("bio-oracle/biooracler")
require(biooracler)
source("rvar/var.R")

#Set working directory
#RVar_wd should be stored in rvar/var.R
setwd(RVar_wd)

#Check the names of the available layers
#available_layers <- list_layers()
#Surface: Ocean Temperature, Salinity, Sea Water Velocity, Nitrate, Silicate,
#Photosynthetic Available Radiation, Disolved Molecular Oxygen, pH, Iron
#Additional: Topographic/Terrain?
#Benthic Temperature: Ocean Temperature

dataset_id <- "so_baseline_2000_2019_depthmax"

time <- c('2000-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
#Slightly larger than required:
#22.89 - 46.25
latitude <- c(22.80, 46.30)
#-124.5 - -144.1
longitude <- c(-124.60, -114.00)

constraints <- list(time, latitude, longitude)
names(constraints) <- c("time", "latitude", "longitude")

variables <- c("so_min", "so_max", "so_mean")
layers <- download_layers(dataset_id, variables, constraints)

dir <- wd
download_layers(dataset_id, variables, constraints, fmt = "csv", directory = dir)

download_layers(dataset_id, variables, constraints, fmt = "raster", directory = dir)



#Set random number string
set.seed(1)

#Set Pacific area boundaries 46.25N and 22.89N latitude ymin and ymax, and then -124.5W and -114.1W longitude xmin and xmax.
Pacific <- st_bbox(c(xmin=-124.5,xmax=-114.1,ymin=22.89,ymax=46.25))


#Get a list of all environmental map layers with the file extension .nc (NetCDF)
map_layers <- list.files(path="MapLayers",pattern = "\\.nc$")

#Convert all of the NetCDF map layers from Bio-Oracle into clipped rasters in tif format.
for(map_layer in map_layers){
  #Read in marine layer
  marine_layer <- brick(paste("MapLayers/",map_layer,sep=""))
  #Get base file name
  map_layer_name <- gsub("\\.nc$","", map_layer)
  #Convert from NetCDF to tif format.
  writeRaster(marine_layer, paste("MapLayers/",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
  #Read in tiff formatted raster.
  marine_raster <- raster(paste("MapLayers/",map_layer_name,".tif",sep=""))
  #Crop the raster to the Pacific extent.
  pacific_raster <- crop(marine_raster,Pacific)
  #Export the Pacific cropped raster.
  writeRaster(pacific_raster,paste("MapLayers/Pacific_",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
 }