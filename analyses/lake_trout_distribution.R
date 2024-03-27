#### ---
# Project: Lake trout distribution
# Start Date: 2023-12-21
# Description: Case study of lake trout distribution for model comparison paper using non-linear interpolation
# mhf
#### ---

###----------------------------------------------------------------------------------------------------
### Load packages
###----------------------------------------------------------------------------------------------------
library(raster)
library(sf)
library(lemon) # facet_rep_wrap()
library(glmmTMB) # mixed effects models
library(glatos) # interpolate_path(); non-linear interpolation
library(emmeans)
library(MuMIn)
library(magrittr)
library(data.table)
library(tidyverse)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Identify path roots
###----------------------------------------------------------------------------------------------------
path.drive <- "C:/Users/mttft/OneDrive - University of Vermont/"
path.root <- paste0(path.drive, "Dissertation_Final/Telemetry/") # Telemetry folder
path.mvt <- paste0(path.root, "ModelComparison/") # Model comparison folder
path.merge <- paste0(path.mvt, "MergedDetections_MC/") # Merged detection data folder for this project (2013-2017) for base, COA, INT, and dBBMM
path.sim <- paste0(path.mvt, "SimulatedData_MC/") # Simulation data: original and for each model
path.result <- paste0(path.mvt, "Results_MC/") # folder for model output (occupancy)
path.figs <- paste0(path.mvt, "Figures_MC/") # figures
path.data <- paste0(path.root, "Detections_2013-2019/") # original detection data (2013-2019)
path.shp <- paste0(path.drive,"ArcGIS Pro 2.2/Shapefiles/LakeChamplain/") # shapefiles for lake data
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


#### ---
# FUNCTION: detect_range
# description: generate sf object for circle polygons with set center and radius (km)
#### ---

detect_range <- function(station, latitude, longitude, radius, nPoints = 100){
  # data: the data frame of receivers with ID
  # radius: radius measured in kilometer
  #
  meanLat <- mean(latitude)
  # length per longitude changes with lattitude, so need correction
  radiusLon <- radius /111 / cos(meanLat/57.3) 
  radiusLat <- radius / 111
  circleDF <- data.frame(ID = rep(station, each = nPoints))
  angle <- seq(0,2*pi,length.out = nPoints)
  
  circleDF$lon <- unlist(lapply(longitude, function(x) x + radiusLon * cos(angle)))
  circleDF$lat <- unlist(lapply(latitude, function(x) x + radiusLat * sin(angle)))
  return(circleDF)
}

# end of detect_range
#### ---


#### ---
# FUNCTION: seasons
# description: generate columns for julian date, season, and season_year
#### ---

seasons <- function(data, timestep){
  data <- data %>% 
    mutate(date_day = yday(timestep),
           year_detect = lubridate::year(timestep)) %>% 
    mutate(season = case_when(year_detect %!in% 2016 & date_day > 334 | year_detect %!in% 2016 & date_day <= 90 ~ "Winter",
                              year_detect %in% 2016 & date_day > 335 | year_detect %in% 2016 & date_day <= 91 ~ "Winter",
                              year_detect %in% 2016 & date_day > 91 & date_day <= 152 ~ "Spring",
                              year_detect %!in% 2016 & date_day > 90 & date_day <= 151 ~ "Spring",
                              year_detect %in% 2016 & date_day > 152 & date_day <= 274 ~ "Summer",
                              year_detect %!in% 2016 & date_day > 151 & date_day <= 273 ~ "Summer",
                              year_detect %in% 2016 & date_day > 274 & date_day <= 335 ~ "Fall",
                              year_detect %!in% 2016 & date_day > 273 & date_day <= 334 ~ "Fall",
                              TRUE ~ 'season'),
           year_detect = as.numeric(year_detect)) %>% 
    mutate(season = factor(season, levels = c("Winter", "Spring", "Summer", "Fall"))) %>% 
    mutate(season_year = case_when(season %in% c("Spring","Summer","Fall") ~ paste(season, substr(timestep, 1, 4), sep = " "),
                                   season %in% "Winter" & timestep < as.POSIXct("2014-06-01 00:00:00") ~ "Winter 2013-2014",
                                   season %in% "Winter" & timestep < as.POSIXct("2015-06-01 00:00:00") ~ "Winter 2014-2015",
                                   season %in% "Winter" & timestep < as.POSIXct("2016-06-01 00:00:00") ~ "Winter 2015-2016",
                                   season %in% "Winter" & timestep < as.POSIXct("2017-06-01 00:00:00") ~ "Winter 2016-2017",
                                   TRUE ~ "season_year")) %>% 
    mutate(season_year = factor(season_year,
                                levels = c("Fall 2013","Winter 2013-2014","Spring 2014","Summer 2014",
                                           "Fall 2014","Winter 2014-2015","Spring 2015","Summer 2015",
                                           "Fall 2015","Winter 2015-2016","Spring 2016","Summer 2016",
                                           "Fall 2016","Winter 2016-2017","Spring 2017","Summer 2017",
                                           "Fall 2017"),
                                ordered = T),
           year_group = if_else(season %in% "Winter" & date_day <= 91, year_detect - 1, year_detect)) %>% 
    mutate(year_group = factor(year_group))
}

# end of seasons
#### ---


####---
# FUNCTION: set_region
# description: set region for interpolated positions
#### ---

set_region <- function(data){
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
                                 latitude > reg_bbox$ymax[reg_bbox$region == "Main_South"] ~ "Main Lake Central",
                               is.na(regions) & 
                                 latitude >= reg_bbox$ymin[reg_bbox$region == "Main_South"] ~ "Main Lake South",
                               is.na(regions) & 
                                 latitude <= reg_bbox$ymin[reg_bbox$region == "Main_South"] ~ "South Lake",
                               !is.na(regions) ~ regions)) %>% 
    mutate(Area_sqkm = case_when(is.na(Area_sqkm) & regions %in% "Main Lake North" ~ 351.817017,
                                 is.na(Area_sqkm) & regions %in% "Main Lake Central" ~ 201.214557,
                                 is.na(Area_sqkm) & regions %in% "Main Lake South" ~ 145.399656,
                                 is.na(Area_sqkm) & regions %in% "Northeast Arm" ~ 245.8366,
                                 is.na(Area_sqkm) & regions %in% "Malletts Bay" ~ 54.811470,
                                 is.na(Area_sqkm) & regions %in% "South Lake" ~ 56.971202,
                                 !is.na(Area_sqkm) ~ Area_sqkm)) %>% 
    mutate(regions = factor(regions,
                            levels = c('Northeast Arm','Malletts Bay','Main Lake North','Main Lake Central','Main Lake South', 'South Lake')))
}

# end of set_region
#### ---


#### ---
# FUNCTION: region_barplot
# description: stacked barplot of regional use for each fish
#### ---

region_barplot <- function(data, tag_name, plot_title){
  data %>% 
    ggplot(aes(tag_name, region_percent, fill = regions)) +
    geom_col(position = "stack")+
    facet_rep_wrap(~season_year, ncol = 4)+ 
    scale_fill_manual(values = c("#543005","#8c510a","#bf812d","#dfc27d","#80cdc1","#35978f","#01665e","#003c30")) +
    ggtitle(plot_title)+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.title.x = element_blank(),
          # legend.position = "none",
          text = element_text(size = 14))
}

# end of region_barplot
#### ---


#### ---
# FUNCTION: region_boxplot
# description: boxplot for regional use by season & year
#### ---

region_boxplot <- function(data, plot_title){
  data %>% 
    ggplot(aes(regions, region_percent)) +
    geom_boxplot(outlier.shape = NA)+
    geom_point(shape = 1, position = position_jitter(width = 0.2))+
    labs(title = plot_title, x = "Region", y = "Percent occurrence")+
    facet_rep_wrap(~season, scales = "fixed", repeat.tick.labels = F) + # look into removing x axis text
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    theme_classic() +
    theme(panel.spacing.x = unit(0.5, "lines"),
          panel.spacing.y = unit(-3.5, "lines"),
          axis.text.x = element_text(angle = 45, hjust=1))
}

# end of region_boxplot
#### ---
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Load receiver data
###----------------------------------------------------------------------------------------------------
recs_full <- read_rds(paste0(path.root,"ReceiverLocations/TotalReceiverSummary_2013-2022_Apr2023.rds"))
str(recs_full)

# Wrangle data to remove recent deployments and assign regions
recs <- recs_full %>% 
  filter(deploy_date_time <= "2017-01-01" &
           region %!in% c("Winooski River","Lamoille River","Otter Creek","Missisquoi River") & # remove receivers in tributaries
           StationName %!in% "BurlingtonBay_RangeTest") %>% 
  rename("regions" = region) %>% 
  mutate(receiver_sn = as.character(receiver_sn),
         StationName = as.character(StationName),
         regions = factor(regions,
                          levels = c('Northeast Channel','Gut','Inland Sea','Main Lake North','Main Lake Central','Main Lake South','Malletts', 'South Lake'),
                          labels = c('NE Channel','Gut','Inland Sea','Main North','Main Central','Main South','Malletts', 'South Lake'))) %>% 
  mutate(receiver_sn = factor(receiver_sn),
         StationName = factor(StationName))

# Corrections to recs file
recs$recover_date_time[recs$StationName == "Gut" & recs$recover_date_time == "2017-05-05"] <- "2017-06-23" # change recover time due to later detection and unknown recover date in catos log


# saveRDS(recs_og, file = paste0(path.mvt,"OriginalReceiverSummary_2013-2017_May2023.rds"))

# Unique receivers and station names
recs_summary <- recs %>% 
  select(StationName,receiver_sn,regions) %>% 
  unique()
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Load fish tagging data
###----------------------------------------------------------------------------------------------------
### Load biometrics data
fish_full <- read.csv(paste0(path.root,'SurgeryLog.csv'))

# Wrangle dataset
fish <- fish_full %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%Y")) %>%
  filter(Date < as.Date("2019-01-01") & Species %in% "Lake Trout" & !is.na(Vemco.ID)) %>%
  select(Date, transmitter, Vemco.ID, Cap.Site, Sex, Length, Fin.Clip) %>%
  rename("cap_date" = Date,
         "animal_id" = Vemco.ID,
         "cap_site" = Cap.Site,
         "sex" = Sex,
         "length" = Length,
         "fin_clip" = Fin.Clip) %>%
  mutate(fin_clip = factor(fin_clip),
         transmitter = factor(transmitter),
         animal_id = factor(animal_id),
         cap_site = factor(cap_site),
         sex = factor(sex))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Load lake map
###----------------------------------------------------------------------------------------------------
# load as sf object
lc_outline <- st_read(dsn = paste0(path.shp, "ChamplainOutline/"),
                      layer = "ChamplainOutline")

outline_sf <- lc_outline[,c("GNIS_NAME","geometry")]

# Create transition layer of the lake
lc_trans <- make_transition(outline_sf,
                            res = c(0.001,0.001),
                            all_touched = F) # pixels must be at least 50% covered by polygon to be coded as water


# load lake regions polygon
lc_regions <- st_read(dsn = paste0(path.shp, "ChamplainRegions/"),
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
            Area_sqkm = sum(AREASQKM)) %>% 
  mutate(regions = factor(regions,
                          levels = c("Missisquoi Bay", "Northeast Arm", "Malletts Bay", 
                                     "Main Lake North", "Main Lake Central", "Main Lake South", "South Lake")))

reg_bbox <- read.csv(paste0(path.shp,"ChamplainRegions/RegionsBBox.csv"))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Generate interpolated data - can skip and load data at start of "Evaluate lake trout distribution"
###----------------------------------------------------------------------------------------------------
# ### Load data formatted for GLATOS (from CompareDistributionModels_2013-2017)
# LT_glatos <- readRDS(file = paste0(path.merge, "LT_glatos_2013-2017_Mar2023.rds"))
# 
# remove tagging season_year detections (minimum of 21 days removed)
# lt_glatos_cut <- LT_glatos %>%
#   mutate(cut = case_when(cap_date < "2014-01-01" & season_year %in% "Fall 2013" ~ 0,
#                          cap_date > "2014-01-01" & season_year %in% "Fall 2014" ~ 0,
#                          .default = 1)) %>% 
#   filter(cut %in% 1)

# # remove dead fish (based on abacus plot)
# lt_glatos_live <- lt_glatos_cut %>%
#   filter(animal_id %!in% animal_id[c(3,16,28,53)]) %>%
#   filter(animal_id %!in% animal_id[72] & glatos_array %!in% "Willsboro")
# 
# ### Interpolate path for each year with 60 min (3600 sec) timestamp
# glatos_60 <- interpolate_path(lt_glatos_live,
#                               trans = lc_trans$transition,
#                               int_time_stamp = 3600,
#                               lnl_thresh = 0.999)

# # save interpolation data
# write_rds(int_60,
#           file = paste0(path.merge,"raw_interpolated_data_2014-2017_live_Jan2024.rds"))
#
# # Load interpolation file
# glatos_60 <- readRDS(file = paste0(path.merge,"raw_interpolated_data_2014-2017_live_Jan2024.rds"))
# 

# NLpos_60 <- setDT(glatos_60)
# int_60 <- unique(NLpos_60, by = c("bin_timestamp", "animal_id"))
# 

# ### Assign interpolated points to regions
# # convert to sf object
# int_sf <- st_as_sf(int_60,
#                    coords = c("longitude","latitude"),
#                    crs=4326,
#                    remove = F)
# 
# # add regions to point data
# int_sf_regions <- st_join(int_sf,
#                           lc_regions_short#,
#                           # join = st_nearest_feature)
# )
# 
# # number of positions on land
# land_int <- int_sf_regions %>%
#   filter(is.na(regions)) %>%
#   nrow() # 4,674
# 
# # percent of positions on land
# land_int/nrow(int_sf_regions)*100 # 0.38%
# 
# 
# # assign interpolated points with NA region to regions using set_region function
# int_sf_regions <- set_region(int_sf_regions)
# 
# 
# # add season and julian date columns using seasons function
# int_sf_seasons <- seasons(data = int_sf_regions,
#                           timestep = int_sf_regions$bin_timestamp)
# 
# # # add weighting factor as number of days per season
# # sea_wt <- int_sf_seasons %>%
# #   group_by(animal_id,season_year) %>%
# #   reframe(sea_wt = n_distinct(date_day))
# #
# # int_sf_seasons <- int_sf_seasons %>%
# #   left_join(sea_wt)
# 
# # add fish data and convert to df
# int_df_full <- left_join(int_sf_seasons, fish) %>%
#   mutate(date_detect = date(bin_timestamp)) %>% # add column for date detected
#   data.frame()
# 
# # remove dates with less than 24 positions
# int_full_day <- int_df_full %>%
#   mutate(date_detect = date(bin_timestamp)) %>%
#   group_by(animal_id, date_detect) %>%
#   reframe(n_pos = n_distinct(bin_timestamp)) %>%
#   unique() %>%
#   filter(n_pos == 24)
# 
# int_df_final <- int_full_day %>%
#   left_join(int_df_full)
# 
# # ### plot interpolated paths for each fish
# # for (i in unique(int_df_final$animal_id)) {
# #
# #   plot_lkt <- int_df_final %>%
# #     filter(animal_id %in% i) %>%
# #     mutate(long = unlist(map(geometry,1)),
# #          lat = unlist(map(geometry,2)),
# #          time = as.numeric(bin_timestamp)) %>%
# #     ggplot() +
# #     geom_point(aes(long, lat, color = time), size = 0.5) +
# #     geom_point(data = recs, aes(x = deploy_long, y = deploy_lat), color = "red", inherit.aes = F) +
# #     geom_sf(data = lc_outline, inherit.aes = F, fill = NA) +
# #     # scale_color_continuous(trans = "reverse") +
# #     labs(main = paste0("LKT: ", i)) +
# #     theme_void()
# #
# #   ggsave(plot = plot_lkt,
# #          filename = paste0(path.mvt,"LKT_Data_MC/Int_Paths/Track_LKT_",i,".svg"),
# #          height = 7,
# #          width = 5)
# #
# # }
# 
# 
# # save data
# write_rds(int_df_final,file = paste0(path.merge,"interpolated_data_daily_2014-2017_live_Jan2024.rds"))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
###  Regional and seasonal distribution by Int model  - can skip
###----------------------------------------------------------------------------------------------------
### Load lake trout positions
int_df_final <- readRDS(file = paste0(path.merge,"interpolated_data_daily_2014-2017_live_Jan2024.rds"))

### convert grouping variables to factors
# create list of all unique dates in dataframe to set date levels
all_periods_int <- sort(unique(int_df_final$date_detect)) 

# date detect and region to factors
int_df_final <- int_df_final %>% 
  mutate(date_detect = factor(date_detect,
                              levels = all_periods_int,
                              ordered = T),
         regions = factor(regions,
                          levels = c("Missisquoi Bay", "Northeast Arm", "Malletts Bay",
                                     "Main Lake North", "Main Lake Central", "Main Lake South", "South Lake")))

### Calculate daily region use
basin_glatos_day <- int_df_final %>%
  group_by(animal_id, regions, date_detect, .drop = F) %>% 
  reframe(region_detect = length(bin_timestamp)) %>% 
  unique() %>% 
  left_join(unique(int_df_final[,c("date_detect","season_year")])) %>% 
  filter(animal_id %in% unique(int_df_final$animal_id))

# add year_group back in to data frame
basin_glatos_day <- left_join(basin_glatos_day, unique(int_df_final[,c("season_year","year_group")])) %>% 
  na.omit() # remove periods with no detections (doesn't change anything)

# calculate total number of detections for each transmitter by season and year
total_detect_glatos_day <- basin_glatos_day %>% 
  group_by(animal_id, date_detect) %>% 
  reframe(day_detect = sum(region_detect)) %>% 
  unique()

# calculate average receiver latitude
rec_lat_glatos <- aggregate(latitude~regions, data = int_df_final, FUN = mean)

# merge number of detections for each fish in each region with total detections for each fish 
basin_glatos_full_day <- left_join(basin_glatos_day, total_detect_glatos_day) %>% 
  left_join(rec_lat_glatos) %>% 
  left_join(fish[,c("animal_id","cap_site","sex")]) # add cap site

# calculate seasonal percentage of detections in each basin for each fish
basin_glatos_full_day <- basin_glatos_full_day %>% 
  mutate(region_percent = round((region_detect*100/day_detect), digits = 1), # replace monthly_detect with day_detect
         region_percent = if_else(region_percent %in% NaN, 0, region_percent), # convert NaN values to 0
         # region_prop = if_else(region_percent %in% 0,
         #                       0.000001,
         #                       (region_percent-0.001)/100),
         region_prop = region_percent/100,
         region_beta = (region_prop*(nrow(basin_glatos_full_day)-1)+0.5)/nrow(basin_glatos_full_day))


### remove periods when fish were not active
# create factor list of all dates
all_dates <- sort(unique(int_df_final$date_detect))

# identify active periods for each fish (season/years between first and last detection)
active_days_int <- int_df_final %>% 
  group_by(animal_id) %>% 
  reframe(first = min(date_detect), # replace mm_yy with season_year
          last = max(date_detect), # replace mm_yy with season_year
          all = list(all_dates[all_dates >= first & all_dates <= last])) %>% 
  unique()

# merge active periods with basin occupancy
basin_glatos_full_day <- left_join(basin_glatos_full_day, active_days_int)

# remove periods when fish were not active (before first detection; after last detection)
basin_glatos_final_day <- basin_glatos_full_day %>% 
  group_by(animal_id) %>% 
  filter(date_detect %in% unlist(all)) %>% # remove periods outside of detection window
  ungroup() %>% 
  filter(region_percent %!in% NaN) # remove NaN values for incomplete seasons (first and last), shouldn't change anything

### create season column 
basin_glatos_final_day <- basin_glatos_final_day %>% 
  mutate(season = case_when(grepl("Spring", season_year) ~ "Spring",
                            grepl("Summer", season_year) ~ "Summer",
                            grepl("Fall", season_year) ~ "Fall",
                            grepl("Winter", season_year) ~ "Winter")) %>% 
  mutate(season = factor(season, 
                         levels = c("Winter", "Spring", "Summer", "Fall"))) # set order


# save final dataframe of monthly regional detections
# write_rds(basin_glatos_final_day,
#           file = paste0(path.merge, "interpolated_daily_region_use_2014-2017_live_Mar2024.rds"))
# basin_glatos_final_day <- read_rds(file = paste0(path.merge, "interpolated_daily_region_use_2014-2017_live_Mar2024.rds"))

### create new summarized dataframe
glatos_percent <- basin_glatos_final_day %>% 
  select(region_percent, regions, season) %>%
  group_by(regions, season) %>%
  summarise_all(list(mean, sd)) %>% 
  rename("mean_glatos" = fn1,
         "sd_glatos" = fn2)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Evaluate lake trout distribution
###----------------------------------------------------------------------------------------------------
### Summary of LKT distribution
# load data 
lt_glatos_live <- readRDS(file = paste0(path.merge, "LT_glatos_2013-2017_Mar2023.rds")) %>% 
  filter(detection_timestamp_utc > cap_date+15) %>% 
  filter(animal_id %!in% animal_id[c(3,16,28,53)]) %>% 
  filter(animal_id %!in% animal_id[72] & glatos_array %!in% "Willsboro")

basin_glatos_final_day <- readRDS(file = paste0(path.merge,"interpolated_daily_region_use_2014-2017_live_Mar2024.rds"))


# Percent of fish in each region based on detection locations and interpolated positions
lt_glatos_live %>% 
  mutate(regions = case_when(regions %in% c("NE Channel","Gut","Inland Sea") ~ "Northeast Arm",
                             regions %in% "Malletts" ~ "Malletts Bay",
                             regions %in% "Main North" ~ "Main Lake North",
                             regions %in% "Main Central" ~ "Main Lake Central",
                             regions %in% "Main South" ~ "Main Lake South",
                             regions %in% "South Lake" ~ "South Lake")) %>% 
  aggregate(animal_id~regions, FUN = n_distinct)

basin_glatos_final_day %>%
  filter(region_percent > 0) %>%
  group_by(regions) %>% 
  reframe(n_fish = n_distinct(animal_id)) %>% 
  mutate(total_fish = n_distinct(basin_glatos_final_day$animal_id),
         per_fish = n_fish*100/total_fish)


### Generate glmm for lake trout distribution
# visualize data
hist(basin_glatos_final_day$region_beta)

summary(basin_glatos_final_day$region_beta)

summary(basin_glatos_final_day)
str(basin_glatos_final_day)


# create glmm using proportional data
dist_mod <- glmmTMB(region_beta~regions*cap_site+sex+(1|animal_id),
                    data = basin_glatos_final_day,
                    family = "beta_family",
                    na.action = "na.fail")

# dist_mod_per <- glmmTMB(region_percent~regions*cap_site+(1|animal_id),
#                     data = basin_glatos_final,
#                     family = "tweedie",
#                     na.action = "na.fail")


# view model output
summary(dist_mod) # multiple significant differences among regions

# dredge
mod_sel_dist <- dredge(dist_mod,
                       trace = T)

best_meta_dist <- subset(mod_sel_dist, delta<2)
best_meta_dist <- subset(best_meta_dist, !is.na(delta))
rowind_dist <- which(best_meta_dist$df==min(best_meta_dist$df))
best_meta_dist <- get.models(best_meta_dist, subset=rowind_dist)[[1]]

mod_sum_dist <- summary(best_meta_dist) 
mod_sum_dist

# model output as df
mod_dist_df <- data.frame(mod_sum_dist$coefficients$cond) %>% 
  round(., 4) %>% 
  rownames_to_column(var = "Parameter") 

mod_dist_df


### visualize differences among regions and cap locations by season
# calculate seasonal average of region use for each fish
lkt_region_season <- basin_glatos_final_day %>% 
  group_by(animal_id,regions,season_year,season) %>% 
  reframe(region_percent = mean(region_percent),
          cap_site = unique(cap_site))

# make function to calculate 5th and 95th percentiles
box_plot_aes <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# make plot
reg_plot <- lkt_region_season %>% 
  ggplot(aes(x = regions,y = region_percent))+
  # geom_point(shape = 21, alpha = 0.6,
  #            aes(fill = cap_site),
  #            position = position_jitterdodge()) +
  stat_summary(fun.data=box_plot_aes, geom="boxplot",
               outlier.shape = NA, alpha = 0.3, # **** CHANGE WHISKERS TO 0.05 and 0.95 percentiles
               aes(fill = cap_site),
               position = position_dodge()) +
  # geom_boxplot(outlier.shape = NA, alpha = 0.3, 
  #              aes(fill = cap_site)) +
  scale_fill_manual(values = c("darkgreen","purple")) +
  facet_rep_wrap(~season) +
  labs(x = NULL) +
  theme_classic() +
  theme(legend.position = "none")

reg_plot

# ggsave(plot = reg_plot,
#        filename = paste0(path.figs,"lkt_int_seasonal_region_use_Jan2024.svg"),
#        height = 7.5,
#        width = 7.5)


### Average basin use as a function of region and release site
lkt_reg_use <- basin_glatos_final_day %>% 
  group_by(cap_site,regions) %>% 
  reframe(ave_use = round(mean(region_percent),1),
          sd_use = round(sd(region_percent),1))

lkt_reg_use


# averages based on predict function
newdat <- basin_glatos_final_day %>%
  group_by(animal_id, cap_site, regions, year_group, .drop = F) %>% 
  reframe(`cap_site:regions` = paste(regions, cap_site, sep = " ")) %>% 
  unique()

est_use <- predict(best_meta_dist, 
                   newdat,
                   se.fit = T)

est_use_df <- data.frame(matrix(unlist(est_use),
                                nrow = nrow(newdat),
                                byrow = F),
                         stringsAsFactors = F)

newdata <- bind_cols(newdat,est_use_df) %>% 
  mutate(ll_logit = X1-X2,
         ul_logit = X1+X2,
         ave = exp(X1)/(1+exp(X1))*100,
         se_ll = exp(ll_logit)/(1+exp(ll_logit))*100,
         se_ul = exp(ul_logit)/(1+exp(ul_logit))*100) %>% 
  unique()


sum_predict <- newdata %>% 
  group_by(regions, cap_site) %>% 
  reframe(ave_pred = mean(ave),
          sd_pred = sd(ave))

sum_predict

### emmeans to test pairwise comparisons
mod_emm_dist <- emmeans(best_meta_dist, specs = c("regions","cap_site"))
mod_p_dist <- pwpm(mod_emm_dist)

mod_p_dist_df <- matrix(mod_p_dist,
                        ncol = 14,
                        dimnames = list(rownames(mod_p_dist),rownames(mod_p_dist))) %>% 
  data.frame()


# write.csv(mod_p_dist_df,
#           file = paste0(path.result,"lkt_int_daily_eval_emmeans_live_Mar2024.csv"))


### Number of individuals from North and South capture sites with no use of region
# number of fish tagged at each site
stock_n <- fish %>% 
  filter(animal_id %in% basin_glatos_final_day$animal_id) %>% 
  group_by(cap_site) %>% 
  reframe(n_fish = n_distinct(animal_id))

# percent of north and south fish that didn't use alternative site
no_use <- lkt_region_season %>% 
  filter(regions %in% c("Main Lake North", "Main Lake South")) %>% 
  group_by(animal_id,regions,cap_site) %>% 
  reframe(region_percent = mean(region_percent)) %>% 
  filter(region_percent == 0) %>% 
  group_by(cap_site, regions) %>% 
  reframe(n_zero = n_distinct(animal_id)#,
          # fish_id = c(unique(animal_id))) %>% 
  ) %>% 
  left_join(stock_n) %>% 
  mutate(per_zero = n_zero*100/n_fish)

no_use
###----------------------------------------------------------------------------------------------------
