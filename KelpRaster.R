rm(list=ls())
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

#Set random number string
set.seed(1)

#Check the names of the available layers
#available_layers <- list_layers()

#Data Set of Current Condition variables, which change with respect to time
datasets_current <- c('chl_baseline_2000_2018_depthsurf', #Chlorophyll
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
				'thetao_baseline_2000_2019_depthmax' #Max Depth Temperature
				)
				
# (Handle separately- Statics)
#'terrain_characteristics' # Topographic
#'par_mean_baseline_2000_2020_depthsurf' #Photosynthetic Available Radiation
#'kdpar_mean_baseline_2000_2020_depthsurf' #Diffuse Attenuation

#Data Set of all Future Condition variables, which change with respect to time, for three SSP
datasets_ssp <- c(
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

time <- c('2000-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
#Slightly larger than required:
#22.89 - 46.25
latitude <- c(22.80, 46.30)
#-124.5 - -144.1
longitude <- c(-124.60, -114.00)

constraints <- list(time, latitude, longitude)
names(constraints) <- c("time", "latitude", "longitude")

nc_dir <- paste(RVar_wd,"ncTemp",sep="")


#Logic Needs adjusting to handle future decades, three separate ssp
layer_names <- c()
i<-1
for(dataset_id in datasets_current){
  variables <- c(paste(str_extract(dataset_id, regex("([^_]+)")),"_mean",sep=""))
  layer_names[i] <- paste(str_extract(time[1],regex("([^-]+)")),"_",str_extract(time[2],regex("([^-]+)")),
                          "_",variables[1],sep="")
  if(layer_names[i] == "2000_2010_thetao_mean"){
    layer_names[i] <- paste(layer_names[i],"_",str_extract(dataset_id, regex("([^_]+)$")),sep="")
  }
  download_layers(dataset_id, variables, constraints, fmt = "raster", directory = nc_dir)
  nc_files_table <- file.info(list.files(nc_dir,pattern="*.nc", full.names = T))
  latest_file <- rownames(nc_files_table)[which.max(nc_files_table$ctime)]
  file.rename(latest_file,paste(nc_dir,"/",layer_names[i],".nc",sep=""))
  i <- i+1
}

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
