###----------------------------------------------------------------------------------------------------
# Lake trout distribution
# Prepared by MHF
# Description: Case study of lake trout regional occupancy using linear/non-linear interpolation
###----------------------------------------------------------------------------------------------------


### set working directory to local path
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###----------------------------------------------------------------------------------------------------
### Load packages
###----------------------------------------------------------------------------------------------------
library(sf)
library(mgcv)
library(glatos)
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
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Load data
###----------------------------------------------------------------------------------------------------
### Load receiver data
recs_full <- read_rds("data/OriginalReceiverSummary_2013-2017.rds")
str(recs_full)


### Load fish tagging data
surgery_log <- read.csv("data/surgery_log.csv")

# Wrangle dataset
fish_data <- surgery_log %>%
  mutate(cap_date = as.Date(cap_date, format = "%m/%d/%Y")) %>%
  select(cap_date, transmitter, animal_id, cap_site, sex, length) %>%
  mutate(transmitter = factor(transmitter),
         animal_id = factor(animal_id),
         cap_site = factor(cap_site),
         sex = factor(sex))


### Load lake map
# load as sf object
lc_outline <- st_read(dsn = "data",
                      layer = "ChamplainOutline")

outline_sf <- lc_outline[,c("GNIS_NAME","geometry")]

# Create transition layer of the lake
lc_trans <- make_transition(outline_sf,
                            res = c(0.001,0.001),
                            all_touched = F) # pixels must be at least 50% covered by polygon to be coded as water


# load lake regions polygon
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
            Area_sqkm = sum(AREASQKM)) %>% 
  mutate(regions = factor(regions,
                          levels = c("Missisquoi Bay", "Northeast Arm", "Malletts Bay", 
                                     "Main Lake North", "Main Lake Central", "Main Lake South", "South Lake")))

reg_bbox <- read.csv("data/RegionsBBox.csv")
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Generate interpolated data
###----------------------------------------------------------------------------------------------------
### Load detection data and reformat for glatos package
# load data
lkt_detects <- readRDS(file = "data/lkt_detections_2013-2017.rds")

# change variable names
lkt_detects <- lkt_detects %>% 
  rename(glatos_array = StationName)

# remove tagging season_year detections (minimum of 21 days removed)
lkt_detects_cut <- lkt_detects %>%
  mutate(cut = case_when(cap_date < "2014-01-01" & season_year %in% "Fall 2013" ~ 0,
                         cap_date > "2014-01-01" & season_year %in% "Fall 2014" ~ 0,
                         .default = 1)) %>%
  filter(cut %in% 1)

# remove dead fish 
lkt_detects_live <- lkt_detects_cut %>%
  filter(animal_id %!in% c(24322,24335,24347,24381)) %>%
  filter(animal_id %!in% 26786 | glatos_array %!in% "Willsboro")

### Interpolate path for each year with 60 min (3600 sec) timestamp
lkt_int_60 <- interpolate_path(lkt_detects_live,
                              trans = lc_trans$transition,
                              int_time_stamp = 3600,
                              lnl_thresh = 0.999)

lkt_int_60 <- setDT(lkt_int_60)
lkt_int_60 <- unique(lkt_int_60, by = c("bin_timestamp", "animal_id"))


### Assign interpolated points to regions
# convert to sf object
lkt_int_sf <- st_as_sf(lkt_int_60,
                       coords = c("longitude","latitude"),
                       crs=4326,
                       remove = F)

# add regions to point data
lkt_int_sf_regions <- st_join(lkt_int_sf,
                              lc_regions_short)

# number of positions on land
land_int <- lkt_int_sf_regions %>%
  filter(is.na(regions)) %>%
  nrow()

# percent of positions on land
land_int/nrow(lkt_int_sf_regions)*100

# assign interpolated points with NA region to regions using set_region function
lkt_int_sf_regions <- set_region(lkt_int_sf_regions)


# add season and julian date columns using seasons function
lkt_int_sf_seasons <- seasons(data = lkt_int_sf_regions,
                              timestep = lkt_int_sf_regions$bin_timestamp)

# add fish data and convert to df
lkt_int_df_full <- left_join(lkt_int_sf_seasons, fish_data) %>%
  mutate(date_detect = date(bin_timestamp)) %>% # add column for date detected
  data.frame()

# remove dates with less than 24 positions (first and last dates detected)
lkt_int_full_day <- lkt_int_df_full %>%
  mutate(date_detect = date(bin_timestamp)) %>%
  group_by(animal_id, date_detect) %>%
  reframe(n_pos = n_distinct(bin_timestamp)) %>%
  unique() %>%
  filter(n_pos == 24)

lkt_int_final <- lkt_int_full_day %>%
  left_join(lkt_int_df_full)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
###  Regional and seasonal distribution
###----------------------------------------------------------------------------------------------------
### convert grouping variables to factors
# create list of all unique dates in dataframe to set date levels
all_periods_int <- sort(unique(lkt_int_final$date_detect)) 

# date detect and region to factors
lkt_int_final <- lkt_int_final %>% 
  mutate(date_detect = factor(date_detect,
                              levels = all_periods_int,
                              ordered = T),
         regions = factor(regions,
                          levels = c("Missisquoi Bay", "Northeast Arm", "Malletts Bay",
                                     "Main Lake North", "Main Lake Central", "Main Lake South", "South Lake")))

### Calculate daily region use
region_int_day <- lkt_int_final %>%
  group_by(animal_id, regions, date_detect, .drop = F) %>% 
  reframe(region_detect = length(bin_timestamp)) %>% 
  unique() %>% 
  left_join(unique(lkt_int_final[,c("animal_id","cap_site","sex","date_detect","season_year")])) %>% # add "animal_id","cap_site","sex" to this line
  filter(animal_id %in% unique(lkt_int_final$animal_id)) %>% 
  mutate(day_detect = 24) ##### add this line

# calculate daily percentage of detections in each basin for each fish
region_int_day <- region_int_day %>%  
  mutate(region_percent = round((region_detect*100/day_detect), digits = 1), 
         region_percent = if_else(region_percent %in% NaN, 0, region_percent), # convert NaN values to 0
         region_prop = region_percent/100,
         region_beta = (region_prop*(nrow(region_int_day)-1)+0.5)/nrow(region_int_day))


### Remove periods when fish were not active
# create factor list of all dates
all_dates <- sort(unique(lkt_int_final$date_detect))

# identify active periods for each fish (season/years between first and last detection)
active_days_int <- lkt_int_final %>% 
  group_by(animal_id) %>% 
  reframe(first = min(date_detect),
          last = max(date_detect),
          all = list(all_dates[all_dates >= first & all_dates <= last])) %>% 
  unique()

# merge active periods with region occupancy
region_int_day <- left_join(region_int_day, active_days_int)

# remove periods when fish were not active (before first detection; after last detection)
region_int_day <- region_int_day %>% 
  group_by(animal_id) %>% 
  filter(date_detect %in% unlist(all)) %>% # remove periods outside of detection window
  ungroup() %>% 
  filter(region_percent %!in% NaN) # remove NaN values for incomplete seasons (first and last), shouldn't change anything

### Create season column 
region_int_day <- region_int_day %>% 
  mutate(season = case_when(grepl("Spring", season_year) ~ "Spring",
                            grepl("Summer", season_year) ~ "Summer",
                            grepl("Fall", season_year) ~ "Fall",
                            grepl("Winter", season_year) ~ "Winter")) %>% 
  mutate(season = factor(season, 
                         levels = c("Winter", "Spring", "Summer", "Fall")))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Evaluate lake trout distribution
###----------------------------------------------------------------------------------------------------
### Summary of LKT distribution
# Percent of fish in each region based on detection locations and interpolated positions
lkt_detects_live %>% 
  mutate(regions = case_when(regions %in% c("NE Channel","Gut","Inland Sea") ~ "Northeast Arm",
                             regions %in% "Malletts" ~ "Malletts Bay",
                             regions %in% "Main North" ~ "Main Lake North",
                             regions %in% "Main Central" ~ "Main Lake Central",
                             regions %in% "Main South" ~ "Main Lake South",
                             regions %in% "South Lake" ~ "South Lake")) %>% 
  aggregate(animal_id~regions, FUN = n_distinct)

region_int_day %>%
  filter(region_percent > 0) %>%
  group_by(regions) %>% 
  reframe(n_fish = n_distinct(animal_id)) %>% 
  mutate(total_fish = n_distinct(region_int_day$animal_id),
         per_fish = n_fish*100/total_fish)


### Generate mixed multinomial logistic regression (Mlogit) for lake trout distribution
# remove regions with no detections
region_short <- region_int_day %>% 
  filter(regions != "Missisquoi Bay") %>% 
  mutate(regions = factor(regions, 
                          levels = c("Northeast Arm","Malletts Bay","Main Lake North","Main Lake Central","Main Lake South","South Lake")))


# run global Mlogit model
lkt_mlogit_full <- gam(list(region_beta~cap_site*regions + sex + s(animal_id, bs='re')),
                       family = multinom(K = 1),
                       data = region_short)

# evaluate significance of model parameters
summary(lkt_mlogit_full) # sex insignificant

# run Mlogit model excluding sex
lkt_mlogit_nsex <- gam(list(region_beta~cap_site*regions + s(animal_id, bs='re')),
                       family = multinom(K = 1),
                       data = region_short)

# compare AIC
lkt_mlogit_full$aic # -296004.2
lkt_mlogit_nsex$aic # -296006.2 * lowest AIC and 1 less factor

# identify the animal ID associated with the smallest effect to use for calculating marginal means
cffs <- lkt_mlogit_nsex$coefficients  
min_lkt <- cffs[grep("animal_id", names(cffs))]
lkt_min <- min_lkt[which(abs(min_lkt) == min(abs(min_lkt)))]

# create object with all combinations of fixed effects
regions <- levels(region_short$regions)
cap_site <- levels(region_short$cap_site)
animal_id <- levels(region_short$animal_id)[6]

newdat <- expand.grid(regions, cap_site, animal_id)
names(newdat) <- c("regions", "cap_site", "animal_id")

# calculate marginal means for each fixed effect
preds <- predict.gam(lkt_mlogit_nsex,
                     newdata = newdat, 
                     type = "response", 
                     se.fit = T)

# add marginal means to the new data object
newdat <- newdat %>% 
  mutate(est_prop = preds$fit[,2],
         se_prop = preds$se.fit[,2],
         CI_lo = est_prop - se_prop*1.96,
         CI_up = est_prop + se_prop*1.96)

# visualize estimated marginal means with 95% CI
predict_plot <- newdat %>% 
  ggplot(aes(x = regions, y = est_prop)) +
  geom_errorbar(aes(ymin = CI_lo, ymax = CI_up, color = cap_site), position = position_dodge(0.9), width = 0.25) +
  geom_point(aes(fill = cap_site, shape = cap_site), position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#B3D1B3","#E3BCFB")) +
  scale_color_manual(values = c("black","black")) +
  scale_shape_manual(values = 21:22) +
  theme_classic()

predict_plot


### Visualize differences among regions and cap locations by season
# calculate seasonal average of region use for each fish
lkt_region_season <- region_int_day %>% 
  group_by(animal_id,regions,season_year,season) %>% 
  reframe(region_percent = mean(region_percent),
          cap_site = unique(cap_site))


### Average region use as a function of region and release site
lkt_reg_use <- region_int_day %>% 
  group_by(cap_site,regions) %>% 
  reframe(ave_use = round(mean(region_percent),1),
          sd_use = round(sd(region_percent),1))

lkt_reg_use


### Number of individuals from North and South capture sites with no use of region
# number of fish tagged at each site
stock_n <- fish_data %>% 
  filter(animal_id %in% region_int_day$animal_id) %>% 
  group_by(cap_site) %>% 
  reframe(n_fish = n_distinct(animal_id))

# percent of north and south fish that didn't use alternative site
no_use <- lkt_region_season %>% 
  filter(regions %in% c("Main Lake North", "Main Lake South")) %>% 
  group_by(animal_id,regions,cap_site) %>% 
  reframe(region_percent = mean(region_percent)) %>% 
  filter(region_percent == 0) %>% 
  group_by(cap_site, regions) %>% 
  reframe(n_zero = n_distinct(animal_id)) %>% 
  left_join(stock_n) %>% 
  mutate(per_zero = n_zero*100/n_fish)
###----------------------------------------------------------------------------------------------------
