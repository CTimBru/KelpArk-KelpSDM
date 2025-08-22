# Clear all variables from the current environment to ensure a clean start.
rm(list=ls())

# Load required packages
require(data.table)
require(raster)
require(sf)
require(terra)
require(devtools)
#devtools::install_github("bio-oracle/biooracler")
require(biooracler)
source("rvar/var.R")
require(stringr)

#Set random number string
set.seed(1)

#Set working directory
#RVar_wd should be stored in rvar/var.R
setwd(RVar_wd)

#Data Set of variables, which change with respect to time
datasets_temporal <- c('chl_baseline_2000_2018_depthsurf', #Chlorophyll
				'clt_baseline_2000_2020_depthsurf', #Cloud Cover
				'dfe_baseline_2000_2018_depthsurf', #Iron
				'mlotst_baseline_2000_2019_depthsurf', #Mixed Layer Depth
				'no3_baseline_2000_2018_depthsurf', #Nitrate
				'o2_baseline_2000_2018_depthsurf', #Dissolved Molecular Oxygen
				'ph_baseline_2000_2018_depthsurf', #pH
				'phyc_baseline_2000_2020_depthsurf', #Primary Productivity
				'po4_baseline_2000_2018_depthsurf', #Phosphate
				'si_baseline_2000_2018_depthsurf', #Silicate
				'so_baseline_2000_2019_depthsurf', #Salinity
				'sws_baseline_2000_2019_depthsurf', #Sea Water Velocity
				'tas_baseline_2000_2020_depthsurf', #Air Temperature
				'thetao_baseline_2000_2019_depthsurf', #Surface Temperature min, mean, max
				'thetao_baseline_2000_2019_depthmax', #Max Depth Temperature
				#SSP126 Middle Road
				'chl_ssp126_2020_2100_depthsurf', #Chlorophyll
				'clt_ssp126_2020_2100_depthsurf', #Cloud Cover
				'dfe_ssp126_2020_2100_depthsurf', #Iron
				'mlotst_ssp126_2020_2100_depthsurf', #Mixed Layer Depth
				'no3_ssp126_2020_2100_depthsurf', #Nitrate
				'o2_ssp126_2020_2100_depthsurf', #Dissolved Molecular Oxygen
				'ph_ssp126_2020_2100_depthsurf', #pH
				'phyc_ssp126_2020_2100_depthsurf', #Primary Productivity
				'po4_ssp126_2020_2100_depthsurf', #Phosphate
				'si_ssp126_2020_2100_depthsurf', #Silicate
				'so_ssp126_2020_2100_depthsurf', #Salinity
				'sws_ssp126_2020_2100_depthsurf', #Sea Water Velocity
				'tas_ssp126_2020_2100_depthsurf', #Air Temperature
				'thetao_ssp126_2020_2100_depthsurf', #Surface Temperature
				'thetao_ssp126_2020_2100_depthmax', #Max Depth Temperature
				#SSP245 Regional Competition
				'chl_ssp245_2020_2100_depthsurf', #Chlorophyll
				'clt_ssp245_2020_2100_depthsurf', #Cloud Cover
				'dfe_ssp245_2020_2100_depthsurf', #Iron
				'mlotst_ssp245_2020_2100_depthsurf', #Mixed Layer Depth
				'no3_ssp245_2020_2100_depthsurf', #Nitrate
				'o2_ssp245_2020_2100_depthsurf', #Dissolved Molecular Oxygen
				'ph_ssp245_2020_2100_depthsurf', #pH
				'phyc_ssp245_2020_2100_depthsurf', #Primary Productivity
				'po4_ssp245_2020_2100_depthsurf', #Phosphate
				'si_ssp245_2020_2100_depthsurf', #Silicate
				'so_ssp245_2020_2100_depthsurf', #Salinity
				'sws_ssp245_2020_2100_depthsurf', #Sea Water Velocity
				'tas_ssp245_2020_2100_depthsurf', #Air Temperature
				'thetao_ssp245_2020_2100_depthsurf', #Surface Temperature
				'thetao_ssp245_2020_2100_depthmax', #Max Depth Temperature
				#SSP585 Drill bby Drill
				'chl_ssp585_2020_2100_depthsurf', #Chlorophyll
				'clt_ssp585_2020_2100_depthsurf', #Cloud Cover
				'dfe_ssp585_2020_2100_depthsurf', #Iron
				'mlotst_ssp585_2020_2100_depthsurf', #Mixed Layer Depth
				'no3_ssp585_2020_2100_depthsurf', #Nitrate
				'o2_ssp585_2020_2100_depthsurf', #Dissolved Molecular Oxygen
				'ph_ssp585_2020_2100_depthsurf', #pH
				'phyc_ssp585_2020_2100_depthsurf', #Primary Productivity
				'po4_ssp585_2020_2100_depthsurf', #Phosphate
				'si_ssp585_2020_2100_depthsurf', #Silicate
				'so_ssp585_2020_2100_depthsurf', #Salinity
				'sws_ssp585_2020_2100_depthsurf', #Sea Water Velocity
				'tas_ssp585_2020_2100_depthsurf', #Air Temperature
				'thetao_ssp585_2020_2100_depthsurf', #Surface Temperature
				'thetao_ssp585_2020_2100_depthmax' #Max Depth Temperature
				)
				
#Data Set for static data, not calculated for SSP
datasets_static <- c(
        'terrain_characteristics', # Topographic
        'par_mean_baseline_2000_2020_depthsurf', #Photosynthetic Available Radiation
        'kdpar_mean_baseline_2000_2020_depthsurf' #Diffuse Attenuation
				)

#Defines the geographic extent for data download. This bounding box is slightly larger than the final cropping area to ensure edge data is included.
#22.89 - 46.25
latitude <- c(22.80, 46.30)
#-124.5 - -144.1
longitude <- c(-124.60, -114.00)

#One-Off for Photosynthetic Available Radiation
model_year_string <- '2010-01-01T00:00:00Z'
time <- c(model_year_string,model_year_string)
variables <- c("par_mean_mean")
dataset_id <- "par_mean_baseline_2000_2020_depthsurf"

#Creates a named list of constraints (time, lat, lon) to pass to the download function.
constraints <- list(time, latitude, longitude)
names(constraints) <- c("time", "latitude", "longitude")
download_layers(dataset_id, variables, constraints, fmt = "raster", directory = nc_dir)

#One-Off for Diffuse Attenuation
model_year_string <- '2010-01-01T00:00:00Z'
time <- c(model_year_string,model_year_string)
variables <- c("kdpar_mean_mean")
dataset_id <- "kdpar_mean_baseline_2000_2020_depthsurf"
constraints <- list(time, latitude, longitude)
names(constraints) <- c("time", "latitude", "longitude")
download_layers(dataset_id, variables, constraints, fmt = "raster", directory = nc_dir)

#Terrain Characteristics
model_year_string <- '1970-01-01T00:00:00Z'
time <- c(model_year_string,model_year_string)
variables <- c("bathymetry_min")
dataset_id <- "terrain_characteristics"
constraints <- list(time, latitude, longitude)
names(constraints) <- c("time", "latitude", "longitude")
download_layers(dataset_id, variables, constraints, fmt = "raster", directory = nc_dir)


nc_dir <- paste(RVar_wd,"ncTemp",sep="")

#This is the main loop for downloading the time-series data defined in `datasets_temporal`.
layer_names <- c()
i<-1
for(dataset_id in datasets_temporal){
  #Get which model name: baseline, ssp126, ssp245, ssp585
  split_name <- strsplit(dataset_id,"_")
  model_var <- split_name[[1]][[1]]
  model_name <- split_name[[1]][[2]]
  model_year <- split_name[[1]][[3]]
  model_depth <- split_name[[1]][[5]]
# block sets the specific variable names to download.
  if(model_depth == "depthsurf") {
    if(model_var == "thetao"){
      variables <- c(paste(model_var,"_min",sep=""),paste(model_var,"_mean",sep=""),paste(model_var,"_max",sep=""))
    } else {
      variables <- c(paste(model_var,"_mean",sep=""))
    }
  } else {
    variables <- c(paste(model_var,"_mean",sep=""))
  }
  # Sets the final year for the while loop based on the scenario.
  if(model_name == "baseline") {
    model_year_max <- model_year+10
  } else {
    model_year_max <- model_year+70
  }
  
  while(model_year <= model_year_max){
    print(model_year)
    # Creates the time string required by Bio-ORACLE for the specific year.
    model_year_string <- paste(model_year,'-01-01T00:00:00Z',sep='')
    #Assembles the constraints list for the current iteration.
    time <- c(model_year_string,model_year_string)
    constraints <- list(time, latitude, longitude)
    names(constraints) <- c("time", "latitude", "longitude")
    #Creates a descriptive file name for the downloaded layer.
    layer_names[i] <- paste(model_name,"_",model_year,"_",
                            model_year+10,
                            "_",variables[1],"_",model_depth,sep="")
    #rename the downloaded file to the descriptive name created earlier
    download_layers(dataset_id, variables, constraints, fmt = "raster", directory = nc_dir)
    nc_files_table <- file.info(list.files(nc_dir,pattern="*.nc", full.names = T))
    latest_file <- rownames(nc_files_table)[which.max(nc_files_table$ctime)]
    file.rename(latest_file,paste(nc_dir,"/",layer_names[i],".nc",sep=""))
    i <- i+1
    model_year <- model_year + 10
  }
}

#Set Pacific area boundaries 46.25N and 22.89N latitude ymin and ymax, and then -124.5W and -114.1W longitude xmin and xmax.
Pacific <- st_bbox(c(xmin=-124.5,xmax=-114.1,ymin=22.89,ymax=46.25))

#Get a list of all environmental map layers with the file extension .nc (NetCDF)
map_layers <- list.files(path="ncTemp",pattern = "\\.nc$")

#Convert all of the NetCDF map layers from Bio-Oracle into clipped rasters in tif format.
for(map_layer in map_layers){
  #Read in marine layer
  marine_layer <- brick(paste("ncTemp/",map_layer,sep=""))
  #Get base file name
  map_layer_name <- gsub("\\.nc$","", map_layer)
  #Convert from NetCDF to tif format.
  writeRaster(marine_layer, paste("MapLayers/NC_",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
  #Read in tiff formatted raster.
  marine_raster <- raster(paste("MapLayers/NC_",map_layer_name,".tif",sep=""))
  #Crop the raster to the Pacific extent.
  pacific_raster <- crop(marine_raster,Pacific)
  #Export the Pacific cropped raster.
  writeRaster(pacific_raster,paste("MapLayers/",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
}

