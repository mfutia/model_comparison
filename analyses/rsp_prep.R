###----------------------------------------------------------------------------------------------------
# Data preparation for actel and RSP packages
# Prepared by MHF
# Description: reformat simulated data for functions from actel and RSP packages
###----------------------------------------------------------------------------------------------------

### set working directory to local path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###----------------------------------------------------------------------------------------------------
### Load packages
###----------------------------------------------------------------------------------------------------
library(sf)
library(actel) 
# remotes::install_github("YuriNiella/RSP", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(RSP) 
library(magrittr)
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
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Load lake map
###----------------------------------------------------------------------------------------------------
lc_outline <- st_read(dsn = "data",
                      layer = "ChamplainOutline")

st_crs(lc_outline) <- 4326

lc_spat <- as(lc_outline, "Spatial")



# load lake region polygon
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
### load receiver data
###----------------------------------------------------------------------------------------------------
recs <- readRDS("data/OriginalReceiverSummary_2013-2017.rds")

# change region names to condensed options
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
sim_trans <- readRDS("outputs/Simulations/complete_simulated_transmissions.rds")

# simulated detected transmissions
sim_dets <- readRDS("outputs/Simulations/simulated_detections.rds")

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

# convert back to df and change region names
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
### Create biometrics dataframe from tag data
actel_biometrics <- sim_actel %>% 
  arrange(animal_id) %>% 
  select(Release.date,Serial.nr,Signal,Group,Release.site,Transmitter) %>% 
  unique()


lubridate::tz(actel_biometrics$Release.date)
str(actel_biometrics)


### Create receiver deployment dataframe
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


### Create the detection dataframe
actel_dets <- sim_actel %>% 
  rename("Receiver" = receiver_sn,
         "Timestamp" = detection_timestamp_utc) %>% 
  mutate(CodeSpace = Serial.nr,
         Receiver = as.integer(Receiver)) %>% 
  select(Timestamp,Receiver,CodeSpace,Signal,Transmitter,animal_id,Array) %>% 
  data.frame()

str(actel_dets)


### Create the spatial dataframe
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
spatial_txt_dot <- "data/actel_spatial_sim.txt"
###----------------------------------------------------------------------------------------------------


