###----------------------------------------------------------------------------------------------------
# Project: dynmaic Brownian Bridge using RSP
# Prepared by MHF and MH
# Description: calculate dBBMM for simulated data using functions from actel and RSP packages
###----------------------------------------------------------------------------------------------------

### set working directory to local path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###----------------------------------------------------------------------------------------------------
### Load packages
###----------------------------------------------------------------------------------------------------
library(raster)
library(terra)
library(sf)
library(lemon)
library(actel) 
# remotes::install_github("YuriNiella/RSP", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(RSP) 
library(magrittr)
library(data.table)
library(tidyverse)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Generate functions for standardizing datasets ####
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
### Load lake map ####
###----------------------------------------------------------------------------------------------------
lc_outline <- st_read(dsn = "data/Shapefiles",
                      layer = "ChamplainOutline")

st_crs(lc_outline) <- 4326

lc_spat <- as(lc_outline, "Spatial")



# load lake region polygon
lc_regions <- st_read(dsn = "data/Shapefiles",
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
### load receiver data
###----------------------------------------------------------------------------------------------------
recs <- readRDS("data/OriginalReceiverSummary_2013-2017.rds")

# correct regions
recs <- recs %>% 
  mutate(regions = as.character(regions),
         regions = if_else(regions %in% c("Gut", "NE Channel", "Inland Sea"),
                           "Northeast Arm",
                           regions),
         regions = factor(regions,
                          levels = c("Northeast Arm", "Malletts", "Main North", "Main Central", "Main South", "South Lake"),
                          labels = c("Northeast_Arm", "Malletts_Bay", "Main_Lake_North", "Main_Lake_Central", "Main_Lake_South", "South_Lake")))


# summarize receiver locations
recs_actel <- recs %>% 
  mutate(receiver_sn = as.numeric(receiver_sn)) %>% 
  group_by(StationName,regions) %>% 
  summarize(deploy_lon = mean(deploy_long),
            deploy_lat = mean(deploy_lat),
            receiver_sn = min(receiver_sn)) %>% 
  mutate("Array" = case_when(regions %in% "Northeast_Arm" ~ "NEA",
                             regions %in% "Malletts_Bay" ~ "MLB",
                             regions %in% "Main_Lake_North" ~ "MLN",
                             regions %in% "Main_Lake_Central" ~ "MLC",
                             regions %in% "Main_Lake_South" ~ "MLS",
                             regions %in% "South_Lake" ~ "SOL")) %>%
  ungroup()
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### load simulation data
###----------------------------------------------------------------------------------------------------
# all simulated transmissions
sim_trans <- readRDS("outputs/Simulations/simulated_trans_regions.rds")

# simulated detected transmissions
sim_dets <- readRDS(paste0(path.sim, "simulated_detections.rds"))

# Calculate the elapsed time (seconds) between detections
sim_dets <- sim_dets %>%
  group_by(sim_id) %>% 
  arrange(detection_timestamp_utc) %>%
  mutate(diff = detection_timestamp_utc - lag(detection_timestamp_utc, default = first(detection_timestamp_utc))) %>% 
  ungroup()

# remove duplicated detections
sim_reduce <- sim_dets %>%
  filter(diff > 60 | trns_id %in% 1)

# remove detections at tributary receivers
sim_reduce <- sim_reduce %>% 
  filter(StationName %!in% c("Lamoille River Lower", "Otter Creek Lower", "Otter Creek Upper", "Winooski River Upper", "Winooski River Lower"))

# Convert to sf object
sim_sf <- st_as_sf(sim_reduce,
                   coords = c("deploy_lon","deploy_lat"),
                   crs = crs(lc_regions_short),
                   remove = F)

# Assign points to regions
sim_sf_regions <- st_join(sim_sf,
                          lc_regions_short[c("regions","area_sqkm")],
                          left = T)

# add convert back to df and correct region names
sim_df_regions <- sim_sf_regions %>% 
  mutate(geometry = NULL,
         start_region = str_replace(start_region, "MainLake", "Main_Lake"), # correct start regions
         start_region = str_replace(start_region, "Inland_Sea", "Northeast_Arm"),
         Group = case_when(start_region %in% "Northeast_Arm" ~ "NEA", # create Array column
                           start_region %in% "Malletts_Bay" ~ "MLB",
                           start_region %in% "Main_Lake_North" ~ "MLN",
                           start_region %in% "Main_Lake_Central" ~ "MLC",
                           start_region %in% "Main_Lake_South" ~ "MLS",
                           start_region %in% "South_Lake" ~ "SOL"),
         regions = str_replace_all(regions, " ", "_"),
         Array = case_when(regions %in% "Northeast_Arm" ~ "NEA", # create Array column
                           regions %in% "Malletts_Bay" ~ "MLB",
                           regions %in% "Main_Lake_North" ~ "MLN",
                           regions %in% "Main_Lake_Central" ~ "MLC",
                           regions %in% "Main_Lake_South" ~ "MLS",
                           regions %in% "South_Lake" ~ "SOL")) %>% 
  data.frame()

str(sim_df_regions)


### create object for actel data
# create transmitter list of 001-100
vals <- seq(1:100) %>% 
  str_pad(width = 3, pad = "0") %>% 
  paste0("sim-",.)

# create general object
sim_actel <- sim_df_regions %>% 
  mutate(Transmitter = paste0("sim-", # add transmitter column
                              str_pad(str_extract(sim_id, "\\d+"),
                                      width = 3,
                                      pad = "0")), 
         sim_id = factor(Transmitter,
                         levels = vals),
         Serial.nr = "sim",
         Signal = as.integer(str_extract(Transmitter,"\\d+")),
         Group = factor(Group,
                        levels = c("NEA","MLB","MLN","MLC","MLS","SOL")),
         Release.date = as.POSIXct("2022-01-01 00:00:00.1",
                                   tz = "UTC")) %>% 
  left_join(recs_actel[,c("StationName","receiver_sn")]) %>% 
  rename("glatos_array" = StationName,
         "animal_id" = sim_id,
         "Release.site" = start_region,
         "deploy_long" = deploy_lon)

str(sim_actel)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Create required dataframes for actel preload
###----------------------------------------------------------------------------------------------------
### Create biometrics dataframe from tag data ---
actel_biometrics <- sim_actel %>% 
  # filter(sim_id %in% "sim_027") %>%  # select specific simulations
  arrange(animal_id) %>% 
  select(Release.date,Serial.nr,Signal,Group,Release.site,Transmitter) %>% 
  unique()


lubridate::tz(actel_biometrics$Release.date)
str(actel_biometrics)


### Create receiver deployment dataframe ---
# add hh:mm:ss to deploy and recover times
actel_deployments <- recs_actel %>% 
  mutate(Start = as.POSIXct(min(sim_trans$time)-1,
                            tz = "UTC",
                            origin = "2022-01-01 00:00:00"),
         Stop = as.POSIXct(max(sim_trans$time)+1,
                           tz = "UTC",
                           origin = "2022-01-01 00:00:00")) %>% 
  rename("Station.name" = StationName,
         "Receiver" = receiver_sn,
         "deploy_long" = deploy_lon) %>%
  arrange(Receiver, Start) %>% 
  data.frame()


### Create the detection dataframe ---
actel_dets <- sim_actel %>% 
  rename("Receiver" = receiver_sn,
         "Timestamp" = detection_timestamp_utc) %>% 
  mutate(CodeSpace = Serial.nr,
         Receiver = as.integer(Receiver)) %>% 
  # animal_id = as.character(animal_id),
  # glatos_array = as.character(glatos_array)),
  select(Timestamp,Receiver,CodeSpace,Signal,Transmitter,animal_id,Array) %>% 
  data.frame()

str(actel_dets)


### Create the spatial dataframe ---
# receiver data
actel_receivers <- recs_actel %>% 
  rename("Station.name" = StationName,
         "Latitude" = deploy_lat,
         "Longitude" = deploy_lon) %>% 
  mutate(Type = "Hydrophone") %>%
  select(Station.name, Latitude, Longitude, Array, Type) %>%
  distinct(Station.name, Latitude, Longitude, Array, Type)


### prepare and style entries for tag releases 
actel_tag_releases <- sim_actel %>%
  rename("Station.name" = Release.site) %>%
  mutate(Latitude = case_when(Station.name %in% "Northeast_Arm" ~ 44.924,
                              Station.name %in% "Main_Lake_North" ~ 44.689,
                              Station.name %in% "Main_Lake_Central" ~ 44.47108,
                              Station.name %in% "Main_Lake_South" ~ 44.44039,
                              Station.name %in% "Malletts_Bay" ~ 44.60510),
         Longitude = case_when(Station.name %in% "Northeast_Arm" ~ -73.228,
                               Station.name %in% "Main_Lake_North" ~ -73.352,
                               Station.name %in% "Main_Lake_Central" ~ -73.23214,
                               Station.name %in% "Main_Lake_South" ~ -73.39258,
                               Station.name %in% "Malletts_Bay" ~ -73.23591),
         Type = "Release") %>% 
  filter(!is.na(Latitude)) %>% 
  select(Station.name, Latitude, Longitude, Array, Type) %>% 
  distinct()

### Combine Releases and Receivers 
actel_spatial <- bind_rows(actel_receivers,actel_tag_releases)

# group by station name and take the mean lat and lon of each station deployment history.
actel_spatial_sum <- actel_spatial %>% 
  group_by(Station.name, Type) %>%
  summarize(Latitude = mean(Latitude),
            Longitude = mean(Longitude),
            Array = first(Array)) %>% 
  mutate(Range = 250)


# load text document for connectivity among regions
spatial_txt_dot <- paste0(path.mvt,'RecLocs_MC/actel_spatial_sim.txt')


### Create actel object from all files
# set timezone object
tz <- "UTC"


### Save actel files and compile with actel function
# # save actel files
# save(actel_biometrics,actel_spatial_sum,actel_deployments,actel_dets,spatial_txt_dot,tz,
#      file = paste0(path.sim,"actel_base_files_Nov2023.R"))

load(paste0(path.sim,"actel_base_files_Nov2023.R"),
     verbose = T)

# run actel preload
actel_data_sim <- preload(biometrics = actel_biometrics,
                          spatial = actel_spatial_sum,
                          deployments = actel_deployments,
                          detections = actel_dets,
                          dot = readLines(spatial_txt_dot),
                          tz = tz)

lubridate::tz(actel_biometrics$Release.date)
lubridate::tz(actel_deployments$Start)
lubridate::tz(actel_deployments$Stop)
lubridate::tz(actel_dets$Timestamp)
lubridate::tz(actel_dets$cap_date)

# # save file
# write_rds(actel_data_sim, file = paste0(path.sim, "actel_sims_Nov2023.rds"))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Prep data for RSP
###----------------------------------------------------------------------------------------------------
### create RSP object
# # load actel dataset - DOESN'T WORK - NEED TO RECREATE FILE EACH SESSION IN PREVIOUS CHUNK 
# actel_data_sim <- readRDS(file = paste0(path.sim, "actel_sims_Nov2023.rds"))

# veiw summary of data
test_actel <- explore(actel_data_sim, 
                      report=F, 
                      print.releases=F,
                      tz = "UTC")

str(test_actel$bio)

# create transition layer used to make shortest paths
setwd(path.mvt)
lake_outline <- loadShape(shape = "ChamplainOutline.shp",
                          size = 0.001, # raster pixel size; not sure what the appropriate value is but this resolution looks appropriate and places all receivers in water 
                          spatial = recs_actel,
                          coord.x = 'deploy_lon',
                          coord.y = 'deploy_lat',
                          buffer = 0.3, # increase shape file edges due to error in dynBBmm function
                          type = "water")

t_layer <- transitionLayer(lake_outline)

# visualize receivers over transition layer
plotRaster(input = test_actel, 
           base.raster = lake_outline_utm, 
           coord.x = "Longitude", 
           coord.y = "Latitude",
           size = 1)

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


# plot tracks
plotTracks(sim_rsp,
           base.raster = lake_outline,
           type = "lines",
           tag = "sim-47")


# extract detections to a new data.frame
rsp_detections <- bind_rows(sim_rsp$detections)

# extract track summary to new data.frame
rsp_tracks <- bind_rows(sim_rsp$tracks)

# # save full rsp file
# write_rds(sim_rsp, file = paste0(path.sim, "rsp_sims_Dec2023.rds"))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Run dynamic Brownian Bridge Movement Models
###---------------------------------------------------------------------------------------------------
# load RSP data
sim_rsp <- readRDS(file = paste0(path.sim, "rsp_sims_Dec2023.rds"))


# ### Test run dBBMM
# # create list of transmitters
# tags <- levels(sim_rsp$detections$`sim-1`$Transmitter)
# tags
# 
# tags_test <- tags[c(1,20,40,60,81)]
# 
# # test run dBBMM for single fish (tag 27)
# test_dBBMM <- dynBBMM(input = sim_rsp,
#                       base.raster = lake_outline,
#                       tags = "sim-47",
#                       UTM = "18",
#                       verbose = T)
# 
# ### assess dBBMM output
# # check valid tracks
# test_dBBMM$valid.tracks
# 
# str(test_dBBMM$valid.tracks)
# 
# 
# # visualize dBBMM
# plotContours(input = test_dBBMM, tag = "sim-81")



### Run dBBMM for all tags
# create vectors with tags by start region
tags <- actel_biometrics[,c("Transmitter","Group")]

MLB <- tags %>% 
  filter(Group %in% "MLB") %>% 
  mutate(Transmitter = paste0(substr(Transmitter, 1, 4),
                              str_sub(Transmitter, -2)),
         Transmitter = if_else(Transmitter %in% "sim-00", "sim-100", Transmitter)) %>% 
  .[["Transmitter"]]

MLN <- tags %>% 
  filter(Group %in% "MLN") %>% 
  mutate(Transmitter = paste0(substr(Transmitter, 1, 4),
                              str_sub(Transmitter, -2))) %>% 
  .[["Transmitter"]]

MLC <- tags %>% 
  filter(Group %in% "MLC") %>% 
  mutate(Transmitter = paste0(substr(Transmitter, 1, 4),
                              str_sub(Transmitter, -2))) %>% 
  .[["Transmitter"]]

# nMLC <- tags %>% 
#   filter(Group %!in% "MLC") %>% 
#   mutate(Transmitter = case_when(Transmitter %in% paste0("sim-",sprintf("%03d",1:9)) ~ paste0(substr(Transmitter, 1, 4), str_sub(Transmitter, -1)),
#                                  Transmitter %in% paste0("sim-",sprintf("%0d",10:99)) ~ paste0(substr(Transmitter, 1, 4), str_sub(Transmitter, -2)),
#                                  Transmitter %in% "sim-100" ~ "sim-100")) %>% 
#   .[["Transmitter"]]
# 

MLS <- tags %>% 
  filter(Group %in% "MLS") %>% 
  mutate(Transmitter = paste0(substr(Transmitter, 1, 4),
                              str_sub(Transmitter, -2))) %>% 
  .[["Transmitter"]]


NEA <- tags %>% 
  filter(Group %in% "NEA") %>% 
  mutate(Transmitter = if_else(Transmitter %in% paste0("sim-",sprintf("%03d",1:9)),
                               paste0(substr(Transmitter, 1, 4), str_sub(Transmitter, -1)),
                               paste0(substr(Transmitter, 1, 4), str_sub(Transmitter, -2)))) %>% 
  .[["Transmitter"]]

# Run dynBBMM for each region
rsp_dBBMM_NEA <- dynBBMM(input = sim_rsp,
                         base.raster = lake_outline,
                         UTM = "18",
                         tags = NEA,
                         verbose = T) # Time spent: 02:09:53

rsp_dBBMM_MLN <- dynBBMM(input = sim_rsp,
                         base.raster = lake_outline,
                         UTM = "18",
                         tags = MLN,
                         verbose = T) # Time spent: 01:19:02

# rsp_dBBMM_MLC.47 <- dynBBMM(input = sim_rsp,
#                             base.raster = lake_outline,
#                             UTM = "18",
#                             tags = "sim-47",
#                             verbose = T) # Time spent: 01:37:55

rsp_dBBMM_MLC <- dynBBMM(input = sim_rsp,
                         base.raster = lake_outline,
                         UTM = "18",
                         tags = MLC,
                         verbose = T) # Time spent: Time spent: 03:11:59

rsp_dBBMM_MLS <- dynBBMM(input = sim_rsp,
                         base.raster = lake_outline,
                         UTM = "18",
                         tags = MLS,
                         verbose = T) # Time spent: 00:42:46

rsp_dBBMM_MLB <- dynBBMM(input = sim_rsp,
                         base.raster = lake_outline,
                         UTM = "18",
                         tags = MLB,
                         verbose = T) # Time spent: 00:18:12

# save(rsp_dBBMM_NEA, rsp_dBBMM_MLN, rsp_dBBMM_MLC, rsp_dBBMM_MLS, rsp_dBBMM_MLB,
#      file = paste0(path.sim,"RSP_dBBMM_sims/all_sims_Dec2023.RData"))


rsp_dBBMM <- dynBBMM(input = sim_rsp,
                     base.raster = lake_outline,
                     UTM = "18",
                     # tags = MLB,
                     verbose = T,
                     window.size = 31,
                     margin = 11) # Time spent: 

# write_rds(rsp_dBBMM,
#           file = paste0(path.sim,"RSP_dBBMM_sims/dBBMM_full_Mar2024.rds"))
# rsp_dBBMM <- readRDS(paste0(path.sim,"RSP_dBBMM_sims/dBBMM_full_Dec2023.rds"))

# run time
# run_time <- c("00:18:10","03:11:50","01:18:37","00:42:43","02:09:45") 
run_time <- c("00:24:25","03:19:01","01:26:25","00:48:58","02:16:28") 
run_time <- chron::chron(times = run_time)
sum(run_time) # 08:15:17

# visualize dBBMM
plotContours(input = rsp_dBBMM,
             tag = "sim-27")


# ### assess dBBMM output
# # check valid tracks
# rsp_dBBMM_all$valid.tracks
# 
# # visualize dBBMM
# plotContours(input = rsp_dBBMM_all, 
#              tag = rsp_dBBMM_all$valid.tracks$Tag[1],
#              breaks = c(0.999999999,.95,.5))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Calculate regional use
###----------------------------------------------------------------------------------------------------
# ### load data
# sim_rsp <- readRDS(file = paste0(path.sim, "rsp_sims_Dec2023.rds"))
# 
# rsp_dBBMM <- readRDS(file = paste0(path.sim,"RSP_dBBMM_sims/dBBMM_full_Dec2023.rds"))

### 
# create list of start regions
regions <- sim_rsp$spatial$array.order$all[1:5]

# create object to hold probability values
dBBMM_spat <- matrix(nrow = 317262, # length(rsp_dBBMM$group.rasters$NEA) #rsp_dBBMM
                     ncol = 1) 

# extract probability values for each tag by start region
for (i in 1:n_distinct(sim_rsp$bio$Release.site)) {
  # object <- mask(, #rsp_dBBMM_all
  #                mask = lake_outline_utm)
  object_df = raster::rasterToPoints(rsp_dBBMM$group.rasters[[i]]) %>% 
    data.frame()
  
  dBBMM_spat = bind_cols(dBBMM_spat,object_df)
}

# remove duplicated columns and rename df column
dBBMM_final <- dBBMM_spat[!duplicated(lapply(dBBMM_spat, summary))] %>% 
  rename("x" = `x...2`,
         "y" = `y...3`) %>% 
  mutate(`...1` = NULL)

# convert to long format with new "Tag" column for identifying individuals
dBBMM_long <- melt(setDT(dBBMM_final),
                   id.vars = c("x","y"),
                   variable.name = "Tag") %>% 
  data.frame() %>% 
  # filter(value >= 0.0001 & value < 0.99) %>% 
  mutate(Tag = str_extract(string = Tag,
                           pattern = "sim.\\d+"))  

# convert to sf object
dBBMM_sf <- st_as_sf(dBBMM_long,
                     coords = c(x = "x", y = "y"),
                     crs = "+proj=utm +zone=18 +datum=WGS84 +ellps=WGS84")  

# crop to lake outline
lake_outline_utm <- projectRaster(lake_outline,
                                  crs = crs(rsp_dBBMM$dbbmm$MLB))

dBBMM_crop <- st_crop(dBBMM_sf,
                      lake_outline_utm)

# change crs to lat long
dBBMM_latlon <- st_transform(dBBMM_crop,
                             crs = crs(lc_regions_short))

# add region to point data
dBBMM_regions <- st_join(dBBMM_latlon, 
                         lc_regions_short[c("regions","area_sqkm")],
                         left = T) 

# convert to df 
dBBMM_regions_full <-  dBBMM_regions %>%
  data.frame()

# idenify number of positions on land
land <- dBBMM_regions_full %>% 
  filter(is.na(regions) & value %!in% c(0,1))# 23,886 land

nrow(land)*100/nrow(dBBMM_regions_full) # 0.07%

# remove points on land
dBBMM_regions_df <-  dBBMM_regions %>%
  filter(!is.na(regions)) %>% 
  unique()

# write_rds(dBBMM_regions_df,
#           file = paste0(path.sim,"RSP_regions_Mar2024.rds"))
# 
# write_rds(dBBMM_regions_full,
#           file = paste0(path.sim,"RSP_full_regions_Mar2024.rds"))


# dBBMM_regions_df <- readRDS(paste0(path.sim,"RSP_regions_Mar2024.rds"))
# dBBMM_regions_full <- readRDS(paste0(path.sim,"RSP_full_regions_Mar2024.rds"))

# dist_df <- raster::rasterToPoints(test_dBBMM$dbbmm$MLN) %>% 
#   data.frame() 
# 
# # reformat to long orientation
# mln_long <- melt(setDT(dBBMM_spat),
#                  id.vars = c("x","y"),
#                  variable.name = "Tag") %>% 
#   data.frame() %>% 
#   filter(value >= 0.1 & value < 0.99) %>% 
#   mutate(Tag = str_extract(string = Tag,
#                            pattern = "sim.\\d+"))
# 
# # convert to sf object
# mlb_sf <- st_as_sf(mln_long,
#                    coords = c(x = "x", y = "y"),
#                    crs = "+proj=utm +zone=18 +datum=WGS84 +ellps=WGS84")
# 
# # change crs to lat long
# mlb_sf <- st_transform(mlb_sf,
#                        crs = crs(lc_regions_short))
# 
# # add region to point data
# sim_mlb_sf_regions <- st_join(mlb_sf,
#                               lc_regions_short[c("regions","area_sqkm")],
#                               left = T)
# 
# sim_df_regions <- sim_mlb_sf_regions %>% 
#   data.frame() %>% 
#   unique()
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Calculate regional distribution
###----------------------------------------------------------------------------------------------------
### load regional distribution data
# # full data
dBBMM_regions_full <- readRDS(paste0(path.sim,"RSP_full_regions_Mar2024.rds"))

# Lake only data
dBBMM_regions_df <- readRDS(paste0(path.sim,"RSP_regions_Mar2024.rds"))

### Regional use by probability
# average probability by individual and region

# calculate total time 
sim_rsp_occ.1 <- dBBMM_regions_df %>%
  mutate(prob = 1-value, # take inverse value for space use to get ~probability
         prob = if_else(prob == 1, 0, prob)) %>% # remove cells that had values of 0
  group_by(Tag) %>%
  reframe(total_prob = sum(prob, na.rm = T)) %>%
  unique()

# calculate region time and percent
sim_rsp_region.1 <- dBBMM_regions_df %>% 
  mutate(prob = 1-value, # take inverse value for space use to get ~probability
         prob = if_else(prob == 1, 0, prob)) %>% # remove cells that had values of 0
  group_by(Tag, regions, .drop = F) %>% 
  reframe(region_detect = sum(prob)) %>% 
  mutate(region_detect = round(as.numeric(region_detect),2)) %>%
  left_join(sim_rsp_occ.1) %>%
  mutate(region_percent = round(as.numeric(region_detect*100/total_prob),2),
         region_percent = if_else(is.na(region_percent),0,region_percent))


### Regional use by count
# calculate total time based on lake only cells
sim_rsp_occ <- dBBMM_regions_df %>%
  filter(value %!in% c(0,1)) %>% 
  group_by(Tag) %>%
  reframe(total_detect = n_distinct(geometry)) %>%
  unique() 

sim_rsp_region <- dBBMM_regions_df %>% 
  filter(value %!in% c(0,1)) %>% 
  group_by(Tag, regions, .drop = F) %>% 
  reframe(region_detect = n_distinct(geometry)) %>% 
  mutate(region_detect = if_else(is.na(region_detect),0,region_detect)) %>%
  left_join(sim_rsp_occ) %>%
  mutate(region_percent = round(as.numeric(region_detect*100/total_detect),2),
         region_percent = if_else(is.na(region_percent),0,region_percent)) %>% 
  unique()
# aggregate(reg_prop~regions, data = dBBMM_reg_count, FUN = n_distinct)
###----------------------------------------------------------------------------------------------------


