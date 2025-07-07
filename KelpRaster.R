require(sf)
require(terra)
require(devtools)
#devtools::install_github("bio-oracle/biooracler")
require(biooracler)
source("rvar/var.R")

#Set random number string
set.seed(1)

#Set working directory
#RVar_wd should be stored in rvar/var.R
#setwd(RVar_wd)
wd <- "/Users/jerelynlee/Documents/GitHub/KelpArk-KelpSDM/"
setwd(wd)

#Check the names of the available layers
#available_layers <- list_layers()
#Benthic Temperature: Ocean Temperature
datasets_current <- c('thetao_baseline_2000_2019_depthsurf', #Surface Temperature
              'so_baseline_2000_2019_depthsurf', #Salinity
              'sws_baseline_2000_2019_depthsurf', #Sea Water Velocity
              'no3_baseline_2000_2018_depthsurf', #Nitrate
              'chl_baseline_2000_2018_depthsurf', #Phosphate
              'si_baseline_2000_2018_depthsurf', #Silicate
              'o2_baseline_2000_2018_depthsurf', #Dissolved Molecular Oxygen
              'dfe_baseline_2000_2018_depthsurf', #Iron
              'par_mean_baseline_2000_2020_depthsurf', #Photosynthetic Available Radiation
              'chl_baseline_2000_2018_depthsurf', #Chlorophyll
              'phyc_baseline_2000_2020_depthsurf', #Primary Productivity
              'ph_baseline_2000_2018_depthsurf', #pH
              'clt_baseline_2000_2020_depthsurf', #Cloud Cover
              'mlotst_baseline_2000_2019_depthsurf', #Mixed Layer Depth
              'tas_baseline_2000_2020_depthsurf', #Air Temperature
              'kdpar_mean_baseline_2000_2020_depthsurf', #Diffuse Attenuation
              'thetao_baseline_2000_2019_depthmax' #Max Depth Temperature
              )
              
#'terrain_characteristics' # Topographic (Hande separately)

datasetsssp <- c('dfe_ssp126_2020_2100_depthmean',
'no3_ssp126_2020_2100_depthmean',
'thetao_ssp126_2020_2100_depthmean',
'ph_ssp126_2020_2100_depthmean',
'po4_ssp126_2020_2100_depthmean',
'phyc_ssp126_2020_2100_depthmean',
'so_ssp126_2020_2100_depthmean',
'sws_ssp126_2020_2100_depthmean',
'si_ssp126_2020_2100_depthmean',
'dfe_ssp245_2020_2100_depthmean')

available_layers <- list_layers()
dataset_id <- "so_baseline_2000_2019_depthmax"

time <- c('2000-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
#Slightly larger than required:
#22.89 - 46.25
latitude <- c(22.80, 46.30)
#-124.5 - -144.1
longitude <- c(-124.60, -114.00)
constraints <- list(time, latitude, longitude)
latitude <- c(22.89, 46.25)
longitude <- c(-124.5, -114.1)

constraints <- list(time, latitude, longitude)
names(constraints) <- c("time", "latitude", "longitude")

nc_dir <- paste(RVar_wd,"ncTemp",sep="")

for(dataset_id in datasets){
  print(dataset_id)
  
  #download_layers(dataset_id, variables, constraints, fmt = "raster", directory = dir)
}
dir <- wd
download_layers(dataset_id, variables, constraints, fmt = "csv", directory = dir)

download_layers(dataset_id, variables, constraints, fmt = "raster", directory = dir)



#Set random number string
set.seed(1)

#Set Pacific area boundaries 46.25N and 22.89N latitude ymin and ymax, and then -124.5W and -114.1W longitude xmin and xmax.
Pacific <- st_bbox(c(xmin=-124.5,xmax=-114.1,ymin=22.89,ymax=46.25))
