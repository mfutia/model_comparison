###----------------------------------------------------------------------------------------------------
# Simulated data for model estimates
# Prepared by MHF and MH
# Description: Generate simulated tracks and estimate regional occupancy using six movement models
###----------------------------------------------------------------------------------------------------


### set working directory to local path
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###----------------------------------------------------------------------------------------------------
### Load packages
###----------------------------------------------------------------------------------------------------
library(sf)
# remotes::install_github('ocean-tracking-network/glatos', build_vignettes = TRUE)
library(glatos)
# devtools::install_github("rossdwyer/VTrack")
library(VTrack)
# remotes::install_github("YuriNiella/RSP") 
library(raster)
library(move)
library(actel)
library(RSP)
library(data.table)
library(tidyverse)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Generate functions for standardizing datasets
###----------------------------------------------------------------------------------------------------

#### ---
# FUNCTION: %!in%
# description: negate based on values in a list
#### ---

`%!in%` <- Negate(`%in%`)

# end of %!in%
#### ---

####---
# FUNCTION: set_region
# description: set region for interpolated positions
#### ---
set_region <- function(data, regions, latitude, longitude){
  data <- data %>% 
    mutate(regions = case_when(is.na(regions) & 
                                 latitude > 44.835 & latitude < 44.970 &
                                 longitude > -73.298 & longitude < -73.226 ~ "Northeast Arm",
                               is.na(regions) & 
                                 latitude > reg_bbox$ymax[reg_bbox$region == "Main_Central"] &
                                 longitude < -73.314 ~ "Main Lake North",
                               is.na(regions) &
                                 latitude > 44.550 & latitude < 44.628 & 
                                 longitude > -73.314 ~ "Malletts Bay",
                               is.na(regions) &
                                 latitude > 44.628 & latitude < 44.922 &
                                 longitude > -73.291 ~ "Northeast Arm",
                               is.na(regions) &
                                 latitude > 44.973 & longitude > -73.227 ~ "Missisquoi Bay",
                               is.na(regions) & 
                                 latitude > reg_bbox$ymax[reg_bbox$region == "Main_South"] ~ "Main Lake Central",
                               is.na(regions) & 
                                 latitude >= reg_bbox$ymin[reg_bbox$region == "Main_South"] ~ "Main Lake South",
                               is.na(regions) & 
                                 latitude <= reg_bbox$ymin[reg_bbox$region == "Main_South"] ~ "South Lake",
                               TRUE ~ as.character(regions))) %>% 
    mutate(area_sqkm = case_when(is.na(area_sqkm) & regions %in% "Missisquoi Bay" ~ 80.459606,
                                 is.na(area_sqkm) & regions %in% "Main Lake North" ~ 351.817017,
                                 is.na(area_sqkm) & regions %in% "Main Lake Central" ~ 201.214557,
                                 is.na(area_sqkm) & regions %in% "Main Lake South" ~ 145.399656,
                                 is.na(area_sqkm) & regions %in% "Northeast Arm" ~ 250.46186,
                                 is.na(area_sqkm) & regions %in% "Malletts Bay" ~ 54.811470,
                                 is.na(area_sqkm) & regions %in% "South Lake" ~ 56.971202,
                                 !is.na(area_sqkm) ~ area_sqkm)) %>%
    mutate(regions = factor(regions,
                            levels = c('Missisquoi Bay','Northeast Arm','Malletts Bay','Main Lake North','Main Lake Central','Main Lake South', 'South Lake')))
}

# end of set_region
#### ---
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Load receiver data
###----------------------------------------------------------------------------------------------------
# load receiver locations
recs <- readRDS("data/OriginalReceiverSummary_2013-2017.rds")

str(recs)
summary(recs)

# correct regions
recs <- recs %>% 
  mutate(regions = as.character(regions)) %>% 
  mutate(regions = case_when(regions %in% c("Gut", "NE Channel", "Inland Sea") ~ "Northeast Arm",
                             regions %in% c("Main North", "Main Central", "Main South", "Malletts", "South Lake") ~ regions)) %>% 
  mutate(regions = factor(regions,
                          levels = c("Northeast Arm", "Malletts", "Main North", "Main Central", "Main South", "South Lake"),
                          labels = c("Northeast Arm", "Malletts Bay", "Main Lake North", "Main Lake Central", "Main Lake South", "South Lake")))

# summarize receiver locations
recs_short <- recs %>% 
  group_by(StationName,regions) %>% 
  reframe(deploy_lon = mean(deploy_long), deploy_lat = mean(deploy_lat))


# convert to sp object
recs_sf <- st_as_sf(x = recs_short, coords = c('deploy_lon', 'deploy_lat'), crs = 4326)

### Summary of receiver deployment
recs_period <- recs %>%
  group_by(StationName, regions) %>%
  reframe(lon = mean(deploy_long),
          lat = round(mean(deploy_lat),4),
          d_date = unique(deploy_date_time),
          r_date = unique(recover_date_time)) %>% 
  unique() %>% 
  arrange(desc(regions),lat) %>% 
  mutate(lat = factor(lat,
                      labels = unique(lat)),
         StationName = factor(StationName,
                              levels = unique(StationName)))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Load lake map
###----------------------------------------------------------------------------------------------------
# setwd(paste(path.shp, "ChamplainOutline/",sep = ""))
lc_outline <- st_read(dsn = "data",
                      layer = "ChamplainOutline")

outline_sf <- lc_outline[,c("GNIS_NAME","geometry")]

# Create transition layer of the lake
lc_trans <- make_transition(outline_sf,
                            res = c(0.0004,0.0004),
                            all_touched = F) 


lc_spat <- as(lc_outline, "Spatial")

# convert outline to raster
ext <- raster::extent(st_bbox(outline_sf))

gridsize <- 0.005
lc_raster <- raster(ext, res=gridsize)
lc_raster <- rasterize(outline_sf, lc_raster)

# load lake region polygon
# setwd(paste0(path.shp, "ChamplainRegions/"))
lc_regions <- st_read(dsn = "data",
                      layer = "ChamplainRegions")

lc_regions_short <- lc_regions %>% 
  mutate(regions = case_when(GNIS_NAME %in% c("Inland_Sea","Gut","NE_Channel") ~ "Northeast Arm",
                             GNIS_NAME %in% "Missisquoi" ~ "Missisquoi Bay",
                             GNIS_NAME %in% "Main_North" ~ "Main Lake North",
                             GNIS_NAME %in% "Main_Central" ~ "Main Lake Central",
                             GNIS_NAME %in% "Main_South" ~ "Main Lake South",
                             GNIS_NAME %in% "Malletts" ~ "Malletts Bay",
                             GNIS_NAME %in% "South_Lake" ~ "South Lake")) %>% 
  group_by(regions) %>% 
  summarize(geometry = sf::st_union(geometry),
            area_sqkm = sum(AREASQKM)) %>% 
  mutate(regions = factor(regions,
                          levels = c("Missisquoi Bay", "Northeast Arm", "Malletts Bay", 
                                     "Main Lake North", "Main Lake Central", "Main Lake South", "South Lake")))

# load lake region bounding boxes
reg_bbox <- as.data.frame(do.call("rbind", lapply(st_geometry(lc_regions), st_bbox)))
reg_bbox$region <- unique(lc_regions$GNIS_NAME)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Create simulated data
###----------------------------------------------------------------------------------------------------
### Create 100 file names for 100 simulated paths
# Five starting locations with 20 replicates per location

# create vector to hold simulated paths
paths <- vector("list", length = 100)

# create metadata table to capture simulation specs
metadata <- data.frame(matrix(ncol = 4, nrow = 100))
colnames(metadata) <- c("sim_id", "theta", "start_region", "velocity")

### Inland Sea simulations
for (i in 1:20) {
  angle = sample(c(5,15,25),1)
  
  paths[[i]] <- crw_in_polygon(polyg = lc_spat,
                               theta = c(0, angle),
                               stepLen = 500,
                               initPos = c(-73.228,44.924),
                               nsteps = 5000,
                               cartesianCRS = 3175,
                               sp_out = T)
  
  metadata$sim_id[i] <- paste0("sim_",i)
  metadata$theta[i] <- angle
  metadata$start_region[i] <- "Inland_Sea"
}

### Main North simulations
for (i in 21:40) {
  angle = sample(c(5,15,25),1)
  
  paths[[i]] <- crw_in_polygon(polyg = lc_spat,
                               theta = c(0, angle),
                               stepLen = 500,
                               initPos = c(-73.352,44.689),
                               nsteps = 5000,
                               cartesianCRS = 3175,
                               sp_out = T)
  
  metadata$sim_id[i] <- paste0("sim_",i)
  metadata$theta[i] <- angle
  metadata$start_region[i] <- "MainLake_North"
}

### Main Central simulations
for (i in 41:60) {
  angle = sample(c(5,15,25),1)
  
  paths[[i]] <- crw_in_polygon(polyg = lc_spat,
                               theta = c(0, angle),
                               stepLen = 500,
                               initPos = c(-73.23214, 44.47108),
                               nsteps = 5000,
                               cartesianCRS = 3175,
                               sp_out = T)
  
  metadata$sim_id[i] <- paste0("sim_",i)
  metadata$theta[i] <- angle
  metadata$start_region[i] <- "MainLake_Central"
}


# Main South simulations
for (i in 61:80) {
  angle = sample(c(5,15,25),1)
  
  paths[[i]] <- crw_in_polygon(polyg = lc_spat,
                               theta = c(0, angle),
                               stepLen = 500,
                               initPos = c(-73.39258, 44.44039),
                               nsteps = 5000,
                               cartesianCRS = 3175,
                               sp_out = T)
  
  metadata$sim_id[i] <- paste0("sim_",i)
  metadata$theta[i] <- angle
  metadata$start_region[i] <- "MainLake_South"
}

# Malletts Bay simulations
for (i in 81:100) {
  angle = sample(c(5,15,25),1)
  
  paths[[i]] <- crw_in_polygon(polyg = lc_spat,
                               theta = c(0, angle),
                               stepLen = 500,
                               initPos = c(-73.23591, 44.60510),
                               nsteps = 5000,
                               cartesianCRS = 3175,
                               sp_out = T)
  
  metadata$sim_id[i] <- paste0("sim_",i)
  metadata$theta[i] <- angle
  metadata$start_region[i] <- "Malletts_Bay"
}

### Convert list of simulated tracks to data.table
# make empty list to hold results
trans <- vector("list", length = 100)

# loop through and calculate transmissions along each path
for(j in 1:100){
  velocity = sample(c(0.1,0.5,0.9),1)
  
  trans[[j]] <- transmit_along_path(paths[[j]], vel = velocity, delayRng = c(120,120), burstDur = 7)
  
  metadata$velocity[j] <- velocity
}

# make empty list to hold results
sim <- vector("list", length = 100)

# Define detection range function (to pass as detRngFun) that returns detection probability for given distance
# assume logistic form of detection range curve where:
#   dm = distance in meters
#   b = intercept and slope
pdrf <- function(dm, b=c(1.8, -1/200)){
  p <- 1/(1+exp(-(b[1]+b[2]*dm)))
  return(p)
}

pdrf(c(100,200,300,400,500,900)) #view detection probs. at some distances

# loop through and calculate transmissions along each path
for(k in 1:100){
  sim[[k]] <- detect_transmissions(trans[[k]], recLoc = recs_sp,
                                   detRngFun = pdrf, show_progress = T)
}

# convert each list of simulated detections into a data.frame
sim_list <- lapply(sim, as.data.frame)

# combine list of each simulated detections data.frame into a single data.table object
sim_df <- rbindlist(sim_list, fill = TRUE, idcol = "virt_fish")

# assign initial regions to each simulation (virt_fish)
sim_df <- sim_df %>% 
  mutate('start_region' = case_when(virt_fish %in% 1:20 ~ "Inland_Sea",
                                    virt_fish %in% 21:40 ~ "MainLake_North",
                                    virt_fish %in% 41:60 ~ "MainLake_Central",
                                    virt_fish %in% 61:80 ~ "MainLake_South",
                                    virt_fish %in% 81:100 ~ "Malletts_Bay",
                                    TRUE ~ 'start_region'))

# create columns for receiver coordinates in simulated dataframe
sim_recs <- data.frame(st_coordinates(sim_df$rec_geometry)) %>% 
  rename("deploy_lon" = X,
         "deploy_lat" = Y)

sim_df <- bind_cols(sim_df, sim_recs)

# create columns for actual coordinates of simulated movements
sim_dets <- data.frame(st_coordinates(sim_df$trns_geometry)) %>% 
  rename("true_lon" = X,
         "true_lat" = Y)

sim_dets_df <- bind_cols(sim_df, sim_dets)

# determine distance between receivers and simulated detections
sim_dets_df$dist <- geosphere::distHaversine(st_coordinates(sim_dets_df$rec_geometry),st_coordinates(sim_dets_df$trns_geometry))


### Structure sim detection data as vemco detection file
# Assign StationNames to detections
sim_dets_full <- left_join(sim_dets_df, recs_short)


### reformat sim_dets_full dataframe
# create initial period to build detection timestamp from
start_time <- as.POSIXct("2022-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# create new columns, set factors, and select final columns
sim_dets_data <- sim_dets_full %>% 
  mutate(detection_timestamp_utc = start_time + time,
         sim_id = paste("sim", virt_fish, sep = "_")) %>% 
  mutate(sim_id = factor(sim_id),
         start_region = factor(start_region)) %>% 
  select(sim_id,start_region,deploy_lon,deploy_lat,true_lon,true_lat,StationName,detection_timestamp_utc,trns_id)


### create new file for complete dataset of all transmission locations
sim_transmissions_list <- lapply(trans, as.data.frame)

sim_transmissions_df <- rbindlist(sim_transmissions_list, fill = TRUE, idcol = "virt_fish")

### create new file for raw positions
sim_pos_list <- lapply(paths, as.data.frame)

sim_pos_df <- rbindlist(sim_pos_list, fill = T, idcol = "virt_sim")
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Summary stats of simulated tracks
###----------------------------------------------------------------------------------------------------
### Duration of each track
sim_transmissions_df <- sim_transmissions_df %>% 
  mutate(sim_id = paste0("sim_",virt_fish))

dur <- sim_transmissions_df %>% 
  mutate(start_region = case_when(virt_fish %in% 1:20 ~ "NEA",
                                  virt_fish %in% 21:40 ~ "MLN",
                                  virt_fish %in% 41:60 ~ "MLC",
                                  virt_fish %in% 61:80 ~ "MLS",
                                  virt_fish %in% 81:100 ~ "MLB")) %>% 
  group_by(virt_fish,start_region) %>% 
  reframe(duration = max(time)/86400)

region_duration <- dur %>% 
  group_by(start_region) %>% 
  reframe(ave_dur = mean(duration),
          sd = sd(duration))


### Determine actual time in regions from simulated transmissions
# create columns for coordinates
true_coords <- data.frame(st_coordinates(sim_transmissions_df$geometry)) %>% 
  rename("true_lon" = X,
         "true_lat" = Y)

sim_comp_full <- bind_cols(sim_transmissions_df, true_coords)

# assign true positions to regions
true_sf <- st_as_sf(sim_comp_full,
                    coords = c("true_lon","true_lat"),
                    crs=4326,
                    remove = F)

# add region to point data
true_occ_sf <- st_join(true_sf,
                       lc_regions_short[c("regions","area_sqkm")],
                       left = T) 

# assign interpolated points with NA region to regions using set_region function
true_occ_sf <- set_region(data = true_occ_sf,
                          regions = true_occ_sf$regions,
                          latitude = true_occ_sf$true_lat,
                          longitude = true_occ_sf$true_lon)

# calculate true regional occupancy
true_occ <- data.frame(true_occ_sf) %>% 
  group_by(virt_fish) %>% 
  mutate(total_positions = as.numeric(length(geometry))) %>% 
  ungroup() %>% 
  group_by(virt_fish,regions,total_positions, .drop = F) %>% 
  reframe(region_positions = as.numeric(length(geometry))) %>% 
  group_by(virt_fish) %>% 
  mutate(total_positions = max(total_positions, na.rm = T),
         percent_region_occ = region_positions*100/total_positions,
         region_positions = if_else(is.na(region_positions),0,region_positions),
         percent_region_occ = if_else(is.na(percent_region_occ),0,percent_region_occ)) %>% 
  ungroup()

summary(true_occ)


### Unique regions used by each track 
# region count
n_regions <- true_occ %>% 
  filter(region_positions > 0) %>% 
  group_by(virt_fish) %>% 
  reframe(n_region = n_distinct(regions))

n_regions %>% 
  mutate(start_region = case_when(virt_fish %in% 1:20 ~ "NEA",
                                  virt_fish %in% 21:40 ~ "MLN",
                                  virt_fish %in% 41:60 ~ "MLC",
                                  virt_fish %in% 61:80 ~ "MLS",
                                  virt_fish %in% 81:100 ~ "MLB")) %>% 
  group_by(start_region) %>% 
  reframe(min_region = min(n_region),
          max_region = max(n_region))


### View regional distribution
true_occ %>% 
  ggplot(aes(x = regions, y = percent_region_occ)) +
  geom_point(shape = 21, fill = "gray40", alpha = 0.6,
             position = position_dodge2(width = 0.4)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.3) +
  labs(x = NULL,
       y = "Distribution (%)") +
  theme_classic()

ave_reg_use <- true_occ %>% 
  group_by(regions) %>%
  reframe(per_use = mean(percent_region_occ),
          sd_use = sd(percent_region_occ))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Generate data frames to use for each model
###----------------------------------------------------------------------------------------------------
### load simulation detection data
sim_dets_data <- readRDS("outputs/Simulations/simulated_detections.rds")

### detection data
# Calculate the elapsed time (seconds) between detections
sim_dets_data <- sim_dets_data %>%
  group_by(sim_id) %>% 
  arrange(detection_timestamp_utc) %>%
  mutate(diff = detection_timestamp_utc - lag(detection_timestamp_utc, default = first(detection_timestamp_utc))) %>% 
  ungroup()

# remove duplicated detections
sim_reduce <- sim_dets_data %>%
  filter(diff > 60 | trns_id %in% 1)

# remove detections at tributary receivers that are not included in analyses
sim_reduce <- sim_reduce %>% 
  filter(StationName %!in% c("Lamoille River Lower", "Otter Creek Lower", "Otter Creek Upper", "Winooski River Upper", "Winooski River Lower"))

str(sim_reduce)

### Assign points to regions
sim_sf <- st_as_sf(sim_reduce,
                   coords = c("deploy_lon","deploy_lat"),
                   crs=crs(lc_regions_short),
                   remove = F)

# add region to point data
sim_sf_regions <- st_join(sim_sf,
                          lc_regions_short[c("regions","area_sqkm")],
                          left = T) 

# convert to data.frame and create columns
sim_base <- data.frame(sim_sf_regions) %>% 
  mutate(transmitter = paste0("trans_",substr(sim_id, 5, nchar(as.character(sim_id)))),
         receiver_sn = paste0(StationName,"_sn"))


### Change formatting to align with glatos package
sim_int <- sim_base %>% 
  select(StationName,receiver_sn,detection_timestamp_utc,deploy_lat,deploy_lon,regions,sim_id, transmitter) %>% 
  rename("animal_id" = sim_id,
         "glatos_array" = StationName,
         "deploy_long" = deploy_lon)


### Change formatting to align with VTrack package
sim_vTrack <- sim_base %>% 
  rename("animal_id" = sim_id) %>% 
  group_by(animal_id) %>% 
  mutate(cap_date = as.Date(min(detection_timestamp_utc)),
         cap_location = start_region) %>% 
  select(animal_id, transmitter, cap_date, cap_location) %>% 
  unique() 

# Reformat VTrack script for current project data
# trace(setupData, edit = T) # corrected code commented out at end of document (see lines 1144-1192)
sim_coa_file <- VTrack::setupData(Tag.Detections = sim_base, 
                                  Tag.Metadata = sim_vTrack, 
                                  Station.Information = recs, 
                                  source = "VEMCO")

# view data
sim_coa_file$Tag.Detections %>% View
sim_coa_file$Tag.Metadata %>% View
sim_coa_file$Station.Information %>% View

# Run COA for 60 minute timesteps
sim_coa_full <- COA(sim_coa_file, timestep = 60) # use "trace(COA, edit = T)" and add one second to "ex" object to avoid parsing failure

# Wrangle df to change column names and remove empty columns
sim_coa <- sim_coa_full %>% 
  ungroup() %>% 
  rename(animal_id = Tag.ID,
         latitude = Latitude.coa,
         longitude = Longitude.coa,
         stations_n = Number.of.Stations,
         detect_n = Number.of.Detections) %>% 
  mutate(animal_id = factor(animal_id)) %>% 
  select("animal_id","detect_n","latitude","longitude","Release.Date","stations_n","TimeStep.coa")


### LOCF
# create event file from GLATOS data
sim_events <- detection_events(sim_int, condense = F)

# convert NA values for timediff to equal difference in time
sim_locf <- sim_events %>% 
  mutate(event2 = if_else(arrive %in% 1, event-1, event),
         time_diff2 = if_else(is.na(time_diff), as.numeric(detection_timestamp_utc - lag(detection_timestamp_utc)), time_diff))


### dBBMM
# remove timestamps with duplicate occurrences
sim_dets_nodups <- sim_dets_data %>% 
  group_by(sim_id) %>% 
  arrange(detection_timestamp_utc) %>% 
  mutate(time_diff = c(NA,diff(detection_timestamp_utc))) %>% 
  filter(time_diff %!in% 0) %>% 
  mutate(sim_id = as.character(sim_id)) %>% 
  arrange(sim_id) %>% 
  ungroup()

str(sim_dets_nodups)

# reduce columns for move object
station_nums <- recs_short %>% 
  dplyr::select(StationName, deploy_lat) %>% 
  arrange(deploy_lat) %>% 
  mutate(station_id = 1:nrow(recs_short),
         deploy_lat = NULL)

sim_dets_short <- sim_dets_nodups %>% 
  left_join(station_nums) %>% 
  mutate(sim_id_num = sub('.*_', '', sim_id),
         start_region = NULL,
         StationName = NULL) %>% 
  mutate(sim_id_num = as.numeric(sim_id_num))  

# create MoveStack object
sim_move <- move(x = sim_dets_short$deploy_lon,
                      y = sim_dets_short$deploy_lat,
                      time = as.POSIXct(sim_dets_short$detection_timestamp_utc, tz = "UTC"),
                      proj = CRS("+proj=longlat +ellps=WGS84"),
                      # data = sim_dets_short,
                      animal = as.character(sim_dets_short$sim_id))

### dBBMM RSP (see analyses/rsp_prep.R for code to format data using RSP)
# load actel dataset
load("data/actel_base_files.R",
     verbose = T)

# run actel preload
actel_data_sim <- preload(biometrics = actel_biometrics,
                          spatial = actel_spatial_sum,
                          deployments = actel_deployments,
                          detections = actel_dets,
                          dot = readLines(spatial_txt_dot),
                          tz = tz)

# veiw summary of data
test_actel <- explore(actel_data_sim, 
                      report=FALSE, 
                      print.releases=FALSE,
                      tz = "UTC")

# create transition layer used to make shortest paths
lake_outline <- loadShape(path = "data",
                          shape = "ChamplainOutline.shp",
                          size = 0.001, # raster pixel size; not sure what the appropriate value is but this resolution looks appropriate and places all receivers in water 
                          spatial = recs_short,
                          coord.x = 'deploy_lon',
                          coord.y = 'deploy_lat',
                          buffer = 0.3, # increase shape file edges due to error in dynBBmm function
                          type = "water")

t_layer <- transitionLayer(lake_outline)

# run rsp function to generate points
sim_rsp <- runRSP(test_actel,
                  t.layer = t_layer,
                  coord.x = 'Longitude',
                  coord.y = 'Latitude',
                  # tags = tags # used to subset tag
                  time.step = 60, # set to 60 min timestep to match other analyses (mins)
                  min.time = 60, # set to match time.step (mins)
                  max.time = 2400, # When the animal is not detected for a period of time longer than the maximum.time argument, a new track is created (hours)
                  er.ad = 8) # increment rate of the position errors for the estimated locations (in meters)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Model 1: Basic residency index (Base)
###----------------------------------------------------------------------------------------------------
### load original data
sim_base <- readRDS("outputs/model_estimates/sim_base.rds") 

### calculate percent of all detections detected at each region for each fish
# Count total detections in each region for each fish and season*year
sim_base_region <- sim_base %>%
  rename("animal_id" = sim_id) %>% 
  group_by(animal_id, regions, .drop = F) %>%
  summarize(region_detect = length(detection_timestamp_utc)) %>%
  ungroup() %>% 
  unique() 

# calculate total number of detections for each virtual id
sim_total_detect <- aggregate(region_detect~animal_id, data = sim_base_region, FUN = sum) %>% 
  rename("total_detect" = region_detect)

# calculate average receiver latitude
sim_rec_lat <- aggregate(deploy_lat~regions, data = sim_base, FUN = mean)

# merge number of detections for each fish in each region with total detections for each fish 
sim_base_region_full <- left_join(sim_base_region, sim_total_detect) %>% 
  left_join(sim_rec_lat)

# calculate seasonal percentage of detections in each region for each fish
sim_base_region_full <- sim_base_region_full %>% 
  mutate(region_percent = round((region_detect*100/total_detect), digits = 1)) %>% 
  mutate(region_percent = if_else(region_percent %in% NaN, 0, region_percent)) %>% 
  rename("latitude" = deploy_lat)

str(sim_base_region_full)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Model 2: Last observation carried forward (LOCF)
###----------------------------------------------------------------------------------------------------
### load original data
sim_locf <- readRDS("outputs/model_estimates/sim_locf.rds")

### calculate percent of time spent in each region for each fish
# calculate total time 
sim_locf_occ <- sim_locf %>%
  group_by(animal_id) %>%
  summarize(total_detect = as.numeric(sum(time_diff2, na.rm = T))) %>%
  ungroup() %>% 
  unique() 

# calculate region time and percent
sim_locf_region <- sim_locf %>% 
  group_by(animal_id, regions, .drop = F) %>% 
  summarize(region_detect = sum(time_diff2)) %>% 
  mutate(region_detect = round(as.numeric(region_detect),2),
         region_detect = if_else(is.na(region_detect),0,region_detect)) %>%
  ungroup() %>% 
  left_join(sim_locf_occ) %>%
  mutate(region_percent = round(as.numeric(region_detect*100/total_detect),2),
         region_percent = if_else(is.na(region_percent),0,region_percent)) %>% 
  ungroup()


# calculate average receiver latitude and join to df
sim_locf_region_full <- left_join(sim_locf_region,
                                 aggregate(deploy_lat~regions, data = sim_locf, FUN = mean)) %>% 
  rename("latitude" = deploy_lat)


str(sim_locf_region_full)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Model 3: Centers of Activity (COA)
###----------------------------------------------------------------------------------------------------
### Data preparation
# Load COA data
sim_coa <- readRDS("outputs/model_estimates/sim_coa.rds")

### Assign COAs to regions
sim_coa_sf <- st_as_sf(sim_coa,
                       coords = c("longitude","latitude"),
                       crs=4326,
                       remove = F)

# add region to point data
sim_coa_sf_regions <- st_join(sim_coa_sf,
                              lc_regions_short[c("regions","area_sqkm")],
                              left = T) 


# assign correct area for NA positions using set_regions function
sim_coa_sf_regions <- set_region(data = sim_coa_sf_regions)

# convert to data.frame
sim_coa_df_regions <- data.frame(sim_coa_sf_regions)

str(sim_coa_df_regions)

summary(sim_coa_df_regions)

### calculate percent of all COAs located in each region for each fish
# Count total COAs in each region for each fish and season*year
sim_coa_region <- sim_coa_df_regions %>%
  group_by(animal_id, regions, .drop = F) %>%
  summarize(region_detect = length(TimeStep.coa)) %>%
  ungroup() %>% 
  unique() 

# calculate total number of detections for each transmitter
sim_total_detect_coa <- aggregate(region_detect~animal_id, data = sim_coa_region, FUN = sum) %>% 
  rename("total_detect" = region_detect)

# calculate average receiver latitude
sim_rec_lat_coa <- aggregate(latitude~regions, data = sim_coa_sf_regions, FUN = mean)

# merge number of detections for each fish in each region with total detections for each fish 
sim_coa_region_full <- left_join(sim_coa_region, sim_total_detect_coa) %>% 
  left_join(sim_rec_lat_coa)

# calculate seasonal percentage of detections in each region for each fish
sim_coa_region_full <- sim_coa_region_full %>% 
  mutate(region_percent = round((region_detect*100/total_detect), digits = 1)) %>% 
  mutate(region_percent = if_else(region_percent %in% NaN, 0, region_percent))

str(sim_coa_region_full)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Model 4: Standardized detections: linear/non-linear interpolated path (Int)
###----------------------------------------------------------------------------------------------------
# Load data
sim_int <- readRDS("outputs/model_estimates/sim_int.rds")


### Interpolate path for each year with 30 min (1800 sec) & 60 min (3600 sec) timestamp
sim_int_60 <- interpolate_path(sim_int,
                               trans = lc_trans$transition,
                               int_time_stamp = 3600,
                               lnl_thresh = 0.999)

# get unique positions by time and simulated track
sim_NLpos_60 <- setDT(sim_int_60)
sim_int_60 <- unique(sim_NLpos_60, by = c("bin_timestamp", "animal_id"))


### Assign interpolated points to regions
# convert to sf object with set crs
sim_int_sf <- st_as_sf(sim_int_60,
                       coords = c("longitude","latitude"), 
                       crs = 4326,
                       remove = F)

# add region to point data
sim_int_sf_regions <- st_join(sim_int_sf,
                              lc_regions_short[c("regions","area_sqkm")],
                              left = T)

# count positions on land 
sim_int_sf_regions %>% 
  filter(is.na(regions)) %>% 
  nrow() # 400 positions on land

# assign interpolated points with NA region to regions using set_region function
sim_int_sf_regions <- set_region(data = sim_int_sf_regions)

# convert to data frame and edit data
sim_int_df_regions <- data.frame(sim_int_sf_regions) %>% 
  mutate(animal_id = factor(animal_id))

str(sim_int_df_regions)

summary(sim_int_df_regions)

### calculate percent of all detections located in each region for each fish
# Count total detections in each region for each fish and season*year
sim_int_region <- sim_int_df_regions %>%
  group_by(animal_id, regions, .drop = F) %>%
  summarize(region_detect = length(bin_timestamp)) %>%
  ungroup() %>% 
  unique() 

# calculate total number of detections for each transmitter
sim_total_detect_int <- aggregate(region_detect~animal_id, data = sim_int_region, FUN = sum) %>% 
  rename("total_detect" = region_detect)

# calculate average receiver latitude
sim_rec_lat_int <- aggregate(latitude~regions, data = sim_int_df_regions, FUN = mean)

# merge number of detections for each fish in each region with total detections for each fish 
sim_int_region_full <- left_join(sim_int_region, sim_total_detect_int) %>% 
  left_join(sim_rec_lat_int)

# calculate percentage of detections in each region for each fish
sim_int_region_full <- sim_int_region_full %>% 
  mutate(region_percent = round((region_detect*100/total_detect), digits = 1)) %>% 
  mutate(region_percent = if_else(region_percent %in% NaN, 0, region_percent))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Model 5: move dynamic Brownian Bridge Movement Model (Move)
###----------------------------------------------------------------------------------------------------
### Load data
sim_move <- readRDS("outputs/model_estimates/sim_move.rds")

### Create underlying raster layer 
ext <- raster::extent(st_bbox(outline_sf))

gridsize <- 0.005
lc_raster <- raster(ext, res=gridsize)
lc_raster <- rasterize(outline_sf, lc_raster)


simProj <- spTransform(sim_move, CRSobj="+proj=utm +zone=18 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

lcProj <- projectRaster(lc_raster, crs = crs(simProj))
lcProj2 <- setExtent(lcProj, extent(lcProj)*1.3)


### calculate max duration of gaps in transmissions for each simulation (seconds)
sim_gaps <- sim_dets_short %>% 
  group_by(sim_id) %>% 
  summarize(max_gap_mins = max(time_diff, na.rm = T)/60, #convert to minutes
            max_gap_hr = max_gap_mins/60) %>% 
  unique()

hist(sim_gaps$max_gap_hr, breaks = seq(0,2500,200))

### Run dBBMM
# create empty list to hold dBBMM output
dBBMM <- list()

# create empty data.frame to hold summary data of excluded track segments
excludeSum <- data.frame(matrix(nrow=length(levels(simProj@trackId)), ncol=6))
names(excludeSum) <- c('sim_ID', 'maxTimeGap', 'propExclude_segs', 'propExclude_time','max_gap_min',"max_gap_hr")

length(levels(simProj@trackId))

for(i in 1:n_distinct(levels(simProj@trackId))){

  # run dBBMM 
  tempind <- simProj[[i]]
  tempID <- as.character(levels(simProj@trackId)[i])
  sim_id <- paste("sim", substr(tempID, 2, nchar(tempID)), sep = "_")
  copyrow <- match(sim_id,sim_gaps$sim_id)
  
  print(i)
  
  temp <- try(brownian.bridge.dyn(object = tempind,
                                  raster = lcProj2,
                                  location.error = 250,
                                  time.step=15), silent = T)
  
  if("try-error" %in% class(temp)){
    j=168 # set minimum gap size to 168 hours
    dbbv <- brownian.motion.variance.dyn(tempind,
                                         location.error=250,
                                         window.size=31,
                                         margin=11)
    
    
    repeat{
      maxTimeGap <- j*60 # convert gap size to minutes
      
      ## In order to get reasonable estimates of where fish when we have confidence in those estimates,
      ## we have to exclude segments where the time lag is too long, because the uncertainty of where the fish
      ## was during that time is too large. Therefore, we exclude segments that have a time lag larger than 168 hours.
      ## The 'dBMvariance' object resulting from the function above contains the slot '@interest' 
      ## in which those segments marked as FALSE won't be included in the calculation of the dBBMM.
      ## Therefore we set all segments with time lag larger than 10080 mins to false
      
      dbbv@interest[timeLag(tempind,"mins")>maxTimeGap] <- FALSE
      propExclude_segs <- length(which(timeLag(tempind,"mins")>maxTimeGap))/length(timeLag(tempind,'mins'))
      propExclude_time <- sum(which(timeLag(tempind,"mins")>maxTimeGap))/sum(timeLag(tempind,'mins'))
      try(temp <- brownian.bridge.dyn(object = dbbv,
                                      raster = lcProj2,
                                      location.error = 250, # taken from example
                                      time.step=15), silent=T)
      if("try-error" %in% class(temp)){j=j-1} else {break} 
    }
    
    print(paste('proportion segments excluded from bbm variance est =', propExclude_segs))
    
    excludeSum[i,1] <- sim_id
    excludeSum[i,2] <- maxTimeGap
    excludeSum[i,3] <- propExclude_segs
    excludeSum[i,4] <- propExclude_time
    excludeSum[i,5] <- sim_gaps[copyrow, 'max_gap_mins']
    excludeSum[i,6] <- sim_gaps[copyrow, 'max_gap_hr']
    
  } else {
    maxTimeGap = 168*60 # time gap in minutes
    propExclude_segs <- length(which(timeLag(tempind,"mins")>maxTimeGap))/length(timeLag(tempind,'mins')) 
    print(paste('proportion segments excluded from bbm variance est =', propExclude_segs))
    
    excludeSum[i,1] <- sim_id
    excludeSum[i,2] <- NA
    excludeSum[i,3] <- 0
    excludeSum[i,4] <- 0
    excludeSum[i,5] <- sim_gaps[copyrow, 'max_gap_mins']
    excludeSum[i,6] <- sim_gaps[copyrow, 'max_gap_hr']
    
  }
  
  dBBMM[[i]] <- temp
  
  if(i == 1){
    coords <- data.frame(coordinates(temp))
    values <- values(temp)
    sim_ID <- as.character(rep(sim_id,length(temp)))
  } else{
    coords <- bind_rows(coords, data.frame(coordinates(temp)))
    values <- c(values, values(temp))
    sim_ID <- c(sim_ID,rep(sim_id,length(temp)))
  }
}

### Combine model coords with probability values and simulation ID 
dBBMM_full <- bind_cols(coords,values,sim_ID) %>% 
  rename("cell_prob" = `...3`,
         "sim_id" = `...4`)
table(dBBMM_full$sim_id)

# set as sf object with UTM projection
dBBMM_full_sf <- st_as_sf(dBBMM_full,
                          coords = 1:2,
                          crs = crs(simProj))

# convert UTM projection to lat lon
dBBMM_latlon <- st_transform(dBBMM_full_sf,
                             crs = crs(lc_regions_short))

# calculate Utilization Distribution (UD)
full_ud <- lapply(dBBMM, FUN = getVolumeUD)


### Assign positions to regions and calculate regional occupancy
# add region to point data
dBBMM_sf_regions <- st_join(dBBMM_latlon,
                            lc_regions_short[c("regions","Area_sqkm")],
                            left = T) 

dBBMM_df_regions <- bind_cols(data.frame(dBBMM_sf_regions),
                              st_coordinates(dBBMM_sf_regions)) %>% 
  mutate(geometry = NULL) %>% 
  rename("longitude" = X,
         "latitude" = Y)

### Calculate regional occupancy
# identify land cells with probability of occurrence > 0
land <- dBBMM_df_regions %>% 
  filter(is.na(regions) &
           cell_prob > 0)
table(land$sim_id)

nrow(land)/nrow(dBBMM_df_regions) # 10.4%

# calculate new probability for lake cells only
dBBMM_lake <- dBBMM_df_regions %>% 
  filter(!is.na(regions))

# calculate total probability for each simulation
dBBMM_lake <- dBBMM_lake %>% 
  group_by(sim_id) %>% 
  mutate(total_occ = sum(cell_prob)) %>% 
  ungroup()

# calculate percent occupancy for each region by simulation  
sim_region_move_full <- data.frame(dBBMM_lake) %>% 
  group_by(sim_id,regions) %>% 
  summarize(region_detect = sum(cell_prob),
            total_detect = mean(total_occ),
            latitude = mean(latitude)) %>% 
  mutate(region_percent = region_detect*100/total_detect) %>% 
  rename("animal_id" = sim_id) %>% 
  unique() %>% 
  ungroup() #%>% 
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Model 6: RSP dynamic Brownian Bridge Movement Model (RSP)
###----------------------------------------------------------------------------------------------------
# load data
sim_rsp <- readRDS(file = "outputs/model_estimates/sim_rsp.rds")

# run dBBMM
rsp_dBBMM <- dynBBMM(input = sim_rsp,
                     base.raster = lake_outline,
                     UTM = "18",
                     # tags = MLB,
                     verbose = T,
                     window.size = 31,
                     margin = 11)

# create list of start regions
regions <- sim_rsp$spatial$array.order$all[1:5]

# create object to hold probability values
rsp_spat <- matrix(nrow = 317262,
                     ncol = 1) 

# extract probability values for each tag by start region
for (i in 1:n_distinct(sim_rsp$bio$Release.site)) {
  object_df = raster::rasterToPoints(rsp_dBBMM$group.rasters[[i]]) %>% 
    data.frame()
  
  rsp_spat = bind_cols(rsp_spat,object_df)
}

# remove duplicated columns and rename df column
rsp_dBBMM_short <- rsp_spat[!duplicated(lapply(rsp_spat, summary))] %>% 
  rename("x" = `x...2`,
         "y" = `y...3`) %>% 
  mutate(`...1` = NULL)

# convert to long format with new "Tag" column for identifying individuals
rsp_dBBMM_long <- melt(setDT(rsp_dBBMM_short),
                   id.vars = c("x","y"),
                   variable.name = "Tag") %>% 
  data.frame() %>% 
  mutate(Tag = str_extract(string = Tag,
                           pattern = "sim.\\d+"))  

# convert to sf object
rsp_dBBMM_sf <- st_as_sf(rsp_dBBMM_long,
                     coords = c(x = "x", y = "y"),
                     crs = "+proj=utm +zone=18 +datum=WGS84 +ellps=WGS84")  

# crop to lake outline
lake_outline_utm <- projectRaster(lake_outline,
                                  crs = crs(rsp_dBBMM$dbbmm$MLB))

rsp_dBBMM_crop <- st_crop(rsp_dBBMM_sf,
                          lake_outline_utm)

# change format to lat long
rsp_dBBMM_latlon <- st_transform(rsp_dBBMM_crop,
                                 crs = crs(lc_regions_short))

# add region to point data
rsp_regions <- st_join(rsp_dBBMM_latlon, 
                       lc_regions_short[c("regions","area_sqkm")],
                       left = T) 

# convert to df 
rsp_regions_full <- rsp_regions %>%
  data.frame()

# idenify number of positions on land
land <- rsp_regions_full %>% 
  filter(is.na(regions) & value %!in% c(0,1))# 23,886 land

nrow(land)*100/nrow(rsp_regions_full) # 0.07%

# remove points on land
rsp_regions_df <-  rsp_regions %>%
  filter(!is.na(regions)) %>% 
  unique()


### Restructure to match format of other files
rsp_regions_ref <- rsp_regions_df %>% 
  rename("animal_id" = Tag) %>% 
  mutate(animal_id = gsub(pattern = ".", # reformat animal id
                          replacement = "_",
                          fixed = T,
                          animal_id))

### Regional use by probability
# calculate total time 
sim_rsp_occ.1 <- rsp_regions_ref %>%
  mutate(prob = 1-value, # take inverse value for space use to get ~probability
         prob = if_else(prob == 1, 0, prob)) %>% # remove cells that had values of 0
  group_by(animal_id) %>%
  reframe(total_prob = sum(prob, na.rm = T)) %>%
  unique()

# calculate region time and percent
sim_rsp_region_full <- rsp_regions_ref %>% 
  mutate(prob = 1-value, # take inverse value for space use to get ~probability
         prob = if_else(prob == 1, 0, prob)) %>%  # remove cells that had values of 0
  group_by(animal_id, regions, .drop = F) %>% 
  reframe(region_detect = sum(prob)) %>% 
  mutate(region_detect = round(as.numeric(region_detect),2)) %>%
  left_join(sim_rsp_occ.1) %>%
  mutate(region_percent = round(as.numeric(region_detect*100/total_prob),2),
         region_percent = if_else(is.na(region_percent),0,region_percent))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Code for updating VTrack loading function (replace code from VTrack::setupData with the code below)
###----------------------------------------------------------------------------------------------------
function (Tag.Detections, Tag.Metadata, Station.Information,
          source = NULL, tzone = "UTC", crs = NULL)
{
  detection_timestamp <- transmitter_id <- station_name <- receiver_name <- latitude <- longitude <- NULL
  sensor_value <- sensor_unit <- Date.and.Time..UTC. <- Transmitter <- Station.Name <- Receiver <- Latitude <- NULL
  Longitude <- Sensor.Value <- Sensor.Unit <- tag_id <- scientific_name <- common_name <- tag_project_name <- NULL
  release_latitude <- release_longitude <- ReleaseDate <- tag_expected_life_time_days <- tag_status <- sex <- NULL
  measurement <- installation_name <- project_name <- deploymentdatetime_timestamp <- recoverydatetime_timestamp <- NULL
  station_latitude <- station_longitude <- status <- NULL
  if (is.null(source))
    stop("Can't recognize the source of your tag detection data.\n'source' should be either 'IMOS' or 'VEMCO'")
  if (source %in% "IMOS") {
    Tag.Detections = as_tibble(Tag.Detections) %>% transmute(Date.Time = lubridate::ymd_hms(detection_timestamp,
                                                                                            tz = tzone), Transmitter = transmitter_id, Station.Name = station_name,
                                                             Receiver = receiver_name, Latitude = latitude, Longitude = longitude,
                                                             Sensor.Value = sensor_value, Sensor.Unit = sensor_unit)
  }
  if (source %in% "VEMCO") {
    Tag.Detections = as_tibble(Tag.Detections) %>% transmute(Date.Time = lubridate::ymd_hms(detection_timestamp_utc,
                                                                                            tz = "UTC"), Transmitter = transmitter, Station.Name = StationName,
                                                             Receiver = receiver_sn, Latitude = deploy_lat, Longitude = deploy_lon,
                                                             Sensor.Value = NA, Sensor.Unit = NA, Receiver.Group = regions)
  }
  object <- structure(list(Tag.Detections = Tag.Detections,
                           Tag.Metadata = as_tibble(Tag.Metadata) %>% transmute(Tag.ID = animal_id,
                                                                                Transmitter = transmitter, Sci.Name = NA, Common.Name = NA,
                                                                                Tag.Project = NA, Release.Latitude = NA, Release.Longitude = NA,
                                                                                Release.Date = lubridate::as_date(cap_date), Tag.Life = NA,
                                                                                Tag.Status = NA, Sex = NA, Bio = NA, Clip = NA),
                           Station.Information = as_tibble(Station.Information) %>%
                             transmute(Station.Name = StationName, Receiver = receiver_sn,
                                       Installation = NA, Receiver.Group = regions,
                                       Deployment.Date = lubridate::as_date(deploy_date_time),
                                       Recovery.Date = lubridate::as_date(recover_date_time),
                                       Station.Latitude = deploy_lat, Station.Longitude = deploy_long,
                                       Receiver.Status = NA)), class = "ATT")
  if (inherits(crs, "CRS")) {
    attr(object, "CRS") <- crs
  }
  else {
    message("Geographic projection for detection positions not recognised, reverting to WGS84 global coordinate reference system")
    attr(object, "CRS") <- CRS("+init=epsg:4326")
  }
  return(object)
}
###----------------------------------------------------------------------------------------------------