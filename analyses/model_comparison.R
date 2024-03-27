###----------------------------------------------------------------------------------------------------
# Movement model assessment and comparison
# Prepared by MHF and MH
# Description: Evaluate model accuracy and compare among movement models
###----------------------------------------------------------------------------------------------------

### set working directory to local path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###----------------------------------------------------------------------------------------------------
### Load packages
###----------------------------------------------------------------------------------------------------
library(sf)
library(raster)
library(lemon)
library(glmmTMB)
library(MuMIn)
library(emmeans)
library(magrittr)
library(data.table)
library(tidyverse)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### load data
###----------------------------------------------------------------------------------------------------
### load model occupancy estimates
sim_base_region_full <- read_rds(file = "outputs/model_estimates/sim_base_region_full.rds")
sim_locf_region_full <- read_rds(file = "outputs/model_estimates/sim_locf_region_full.rds")
sim_coa_region_full <- read_rds(file = "outputs/model_estimates/sim_coa_region_full.rds")
sim_int_region_full <- read_rds(file = "outputs/model_estimates/sim_int_region_full.rds")
sim_move_region_full <- read_rds(file = "outputs/model_estimates/sim_move_region_full.rds")
sim_rsp_region_full <- read_rds(file = "outputs/model_estimates/sim_rsp_region_full.rds")

# load true occupancy data
true_use <- readRDS(file = "outputs/simulations/simulated_true_data.rds")

# load metadata
metadata <- readRDS(paste0(path.sim, "simulations_metadata.rds")) %>% 
  rename("animal_id" = sim_id)


### merge dataframes into one
# create column for each model
sim_base_region_full$model <- "Base"
sim_locf_region_full$model <- "LOCF"
sim_coa_region_full$model <- "COA"
sim_int_region_full$model <- "Int"
sim_move_region_full$model <- "Move"
sim_rsp_region_full$model <- "RSP"

# combine model output
model_est <- bind_rows(sim_base_region_full,
                       sim_locf_region_full,
                       sim_coa_region_full,
                       sim_int_region_full,
                       sim_move_region_full,
                       sim_rsp_region_full)

# add true region use and calculate total and regional simulation duration
model_full <- true_use %>% 
  select(animal_id,regions,region_positions,total_positions) %>% 
  right_join(model_est) %>% 
  mutate(region_dur = region_positions*127, # calculate duration by region for each simulation (127 = transmission delay (120 sec) + transmission duration (7 sec))   
         total_dur = total_positions*127, # calculate total duration of simulation in seconds
         region_wt = region_dur/total_dur, # calculate proportion of time spent in each region
         abs_error_per = abs(region_percent-(region_wt*100)), # calculate absolute model error for each sim/region ***** THINK ABOUT THIS ****
         act_error_per = region_percent-(region_wt*100),
         wt_abs_error = abs_error_per*region_wt, # calculate occupancy error for each sim/region/model
         wt_act_error = act_error_per*region_wt) # calculate actual weighted error

str(model_full)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Criteria 1: Weighted mean occupancy error (MOE) 
###----------------------------------------------------------------------------------------------------
### calculate MOE for each simulation per model
# calculate indvidual error for each simulation/region
moe_total <- model_full %>% 
  group_by(animal_id,model) %>% 
  reframe(mean_wt_abs_err = sum(wt_abs_error)) %>% 
  mutate(model = factor(model,
                        levels = c("Base", "LOCF","COA", "Int", "Move", "RSP")))


# plot error among models
moe_plot <- moe_total %>% 
  ggplot(aes(x = model, y = mean_wt_abs_err)) +
  geom_point(shape = 21, fill = "gray40", alpha = 0.6,
             position = position_dodge2(width = 0.4)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  # ylim(c(-0.01,0.15)) +
  labs(x = NULL,
       y = "Mean Weighted Occupancy Error (%)") +
  theme_classic()

moe_plot


### Evaluate MOE among models
# reformat model data and add metadata
mod_moe_perf <- model_full %>% 
  group_by(animal_id,model) %>% 
  summarize(mean_wt_abs_err = sum(wt_abs_error)) %>% 
  mutate(prop_mean_wt_err = mean_wt_abs_err/100,
         model = factor(model,
                        levels = c("Base", "LOCF", "COA", "Int", "Move", "RSP"))) %>% 
  left_join(metadata) %>% 
  mutate(start_region = factor(start_region,
                               levels = c("MainLake_North","MainLake_Central","MainLake_South","Inland_Sea","Malletts_Bay"),
                               labels = c("Main Lake North","Main Lake Central","Main Lake South","Northeast Arm","Malletts Bay")),
         prop_mean_wt_err2 = if_else(prop_mean_wt_err == 0, 1e-10, prop_mean_wt_err),
         animal_id = factor(animal_id))


# MOE parameter assessment
summary(mod_moe_perf)

# create glmm
mod_perf <- glmmTMB(mean_wt_abs_err ~ start_region+velocity+theta+model+(1|animal_id),
                    family = "tweedie",
                    data = mod_moe_perf,
                    na.action = "na.fail")

# view model output
summary(mod_perf) # only model significant 

# model selection based on dAIC and number of model parameters
mod_sel_moe <- dredge(mod_perf, trace = T)

best_meta_moe <- subset(mod_sel_moe, delta<2)
best_meta_moe <- subset(best_meta_moe, !is.na(delta))
rowind_moe <- which(best_meta_moe$df==min(best_meta_moe$df))
best_meta_moe <- get.models(best_meta_moe, subset=rowind_moe)[[1]]

# model summary
summary(best_meta_moe)

# emmeans to test pairwise comparisons of model parameters
mod_emm_moe <- emmeans(best_meta_moe, specs = c("model"))
mod_p_moe <- pwpm(mod_emm_moe)

mod_p_moe_df <- matrix(mod_p_moe,
                       ncol = 6,
                       dimnames = list(rownames(mod_p_moe),rownames(mod_p_moe))) %>%
  data.frame()

mod_p_moe_df
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Criteria 2 & 3: Regional occupancy error (ROE and aROE)
###----------------------------------------------------------------------------------------------------
### Unique regions used by each track
# calculate average use for true distribution and modeled distribution
ave_reg_use <- true_use %>% 
  group_by(regions) %>%
  reframe(per_use = mean(percent_region_use),
          sd_use = sd(percent_region_use),
          med_use = median(percent_region_use))

est_reg <- model_est %>% 
  group_by(regions,model,.drop = F) %>% 
  reframe(ave_use = mean(region_percent),
          sd_use = sd(region_percent))

# calculate average regional difference (actual and absolute) in use for each method
ave_diff <- est_reg %>% 
  left_join(ave_reg_use[,c("regions","per_use")]) %>% 
  mutate(diff = ave_use - per_use,
         abs_diff = abs(diff))

# plot estimated regional distribution for each method (include true median)
est_plot <- model_est %>% 
  mutate(model = factor(model,
                        levels = c("Base", "LOCF", "COA", "Int", "Move", "RSP"))) %>% 
  ggplot(aes(x = regions, y = region_percent)) +
  geom_point(shape = 21, fill = "gray40", alpha = 0.6,
             position = position_dodge2(width = 0.3)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.3)+
  geom_point(data = ave_reg_use, aes(x = regions, y = med_use),
             shape = 24, fill = "red",
             inherit.aes = F) +
  facet_rep_wrap(~model,
                 ncol = 2) +
  theme_classic()

est_plot


### Absolute regional occupancy error (aROE) associated with each model
# table for mean +/- sd of aROE
ave_err <- model_full %>% 
  group_by(regions,model) %>% 
  reframe(abs_mean = round(mean(abs_error_per),1),
          abs_sd = round(sd(abs_error_per),1),
          act_mean = round(mean(act_error_per),1),
          act_sd = round(sd(act_error_per),1)) %>% 
  mutate(model = factor(model,
                        levels = c("Base", "LOCF", "COA", "Int", "Move", "RSP"))) %>% 
  arrange(regions, model)

# calculate difference between estimated use and true use
full_diff <- true_use %>% 
  rename("true_percent" = percent_region_use) %>% 
  left_join(model_est) %>% 
  mutate(diff = region_percent - true_percent,
         model = factor(model,
                        levels = c("Base", "LOCF", "COA", "Int", "Move", "RSP"))) %>%
  select(animal_id,model,regions,true_percent,region_percent,diff)


### map filling region by error
map_error <- lc_regions_short %>% 
  left_join(ave_err[,c("regions", "model", "act_mean")]) %>% # act_mean | abs_mean
  ggplot() +
  geom_sf(aes(fill = act_mean)) + # act_mean | abs_mean
  facet_wrap(~model,
             nrow = 1)+
  scico::scale_fill_scico(palette = "vik", midpoint = 0)+
  theme_void() 

map_error


###  Evaluate model performance
# ROE and aROE data
mod_reg_perf <- model_full %>%
  left_join(metadata) %>%
  mutate(start_region = factor(start_region,
                               levels = c("MainLake_North","MainLake_Central","MainLake_South","Inland_Sea","Malletts_Bay"),
                               labels = c("Main Lake North","Main Lake Central","Main Lake South","Northeast Arm","Malletts Bay")),
         model = factor(model,
                        levels = c("Base", "LOCF", "COA", "Int", "Move", "RSP")))


### Regional performance base on actual error 
mod_perf_act <- glmmTMB(act_error_per ~ factor(velocity)+factor(theta)+start_region+regions*model+(1|animal_id),
                        data = mod_reg_perf,
                        na.action = "na.fail")

# view model output
summary(mod_perf_act) # multiple significant differences among regions

# model selection based on dAIC and number of model parameters
mod_sel_act <- dredge(mod_perf_act,
                      trace = T)

best_meta_act <- subset(mod_sel_act, delta<2)
best_meta_act <- subset(best_meta_act, !is.na(delta))
rowind_act <- which(best_meta_act$df==min(best_meta_act$df))
best_meta_act <- get.models(best_meta_act, subset=rowind_act)[[1]]

mod_sum_act <- summary(best_meta_act)
mod_sum_act

# model output as df
mod_act_df <- data.frame(mod_sum_act$coefficients$cond) %>%
  round(., 4) %>%
  rownames_to_column(var = "Parameter")


# emmeans to test pairwise comparisons
mod_emm_act <- emmeans(best_meta_act, specs = c("regions","model"))
mod_p_act <- pwpm(mod_emm_act)

mod_p_act_df <- matrix(mod_p_act,
                       ncol = 42,
                       dimnames = list(rownames(mod_p_act),rownames(mod_p_act))) %>%
  data.frame()


### Regional performance based on absolute error
mod_perf_abs <- glmmTMB(abs_error_per ~ factor(velocity)+factor(theta)+start_region+regions*model+(1|animal_id),
                        family = tweedie,
                        data = mod_reg_perf,
                        na.action = "na.fail")

# view model output
summary(mod_perf_abs) 

# model selection based on dAIC and number of model parameters
mod_sel_abs <- dredge(mod_perf_abs,
                      trace = T)

best_meta_abs <- subset(mod_sel_abs, delta<2)
best_meta_abs <- subset(best_meta_abs, !is.na(delta))
rowind_abs <- which(best_meta_abs$df==min(best_meta_abs$df))
best_meta_abs <- get.models(best_meta_abs, subset=rowind_abs)[[1]]

mod_sum_abs <- summary(best_meta_abs) 
mod_sum_abs

# model output as df
mod_abs_df <- data.frame(mod_sum_abs$coefficients$cond) %>% 
  round(., 4) %>% 
  rownames_to_column(var = "Parameter") 

mod_abs_df

# emmeans to test pairwise comparisons **** REDO EMMEANS WITH TRUE ERROR??? *****
mod_emm_abs <- emmeans(best_meta_abs, specs = c("regions","model"))
mod_p_abs <- pwpm(mod_emm_abs)

mod_p_abs_df <- matrix(mod_p_abs,
                       ncol = 42,
                       dimnames = list(rownames(mod_p_abs),rownames(mod_p_abs))) %>% 
  data.frame()
###----------------------------------------------------------------------------------------------------

##### Start here
###----------------------------------------------------------------------------------------------------
### Daily region assignment accuracy
###----------------------------------------------------------------------------------------------------
### Determine actually daily zone assignment
# load simulated transmissions
true_sf_regions <- readRDS(paste0(path.sim,"simulated_trans_regions.rds"))

# convert to df
true_df <- true_sf_regions %>% 
  mutate(geometry = NULL) %>% 
  data.frame()

# quantify daily region transmissions
true_daily_regions <- true_df %>% 
  mutate(day = ceiling(time/86400)) %>%
  group_by(virt_fish,regions,day) %>% 
  summarize(points = n_distinct(time)) %>% 
  ungroup() %>% 
  unique()

# select region with greatest number of points for each day/simulation
max_daily_regions <- true_daily_regions %>% 
  group_by(virt_fish,day) %>% 
  mutate(max_pt = max(points),
         is_max = if_else(points %in% max_pt, 1, 0),
         sim_id = paste0("sim_",virt_fish),
         sim_id = factor(sim_id,
                         levels = unique(sim_id))) %>% 
  ungroup() %>% 
  group_by(sim_id, day) %>% 
  filter(is_max %in% 1) %>% 
  reframe(regions = regions[1]) %>% 
  rename("true_region" = regions) %>% 
  unique()

### Determine daily region assignment for each simulation/model
# base model  ### --- CAN SKIP THIS MODEL ??? ---
sim_locf <- readRDS(paste0(path.sim,"sim_locf_detections.rds"))

# daily assignment
daily_dist_base <- sim_locf %>%
  mutate(day = yday(detection_timestamp_utc)) %>%
  rename("sim_id" = animal_id) %>%
  group_by(sim_id,regions,day) %>%
  reframe(points = sum(time_diff2)) %>%
  group_by(sim_id, day) %>% 
  filter(points %in% max(points)) %>% 
  reframe(points = max(points),
          regions = regions[1]) %>% 
  unique()

# # Identify duplicated dates by sim id
# t <- daily_dist_base %>% select(sim_id,day) %>% arrange(sim_id,day) %>% duplicated() %>% data.frame()

# select region with greatest number of points for each day/simulation
max_daily_base <- daily_dist_base %>%
  # group_by(sim_id,day) %>%
  mutate(# max_pt = max(points),
    # is_max = if_else(points %in% max_pt, 1, 0),
    model = "base") %>%
  # filter(is_max %in% 1) %>%
  # ungroup() %>%
  select(sim_id,regions,day,model) %>%
  unique()


# COA
load(paste0(path.sim, "simulated_COA_data.Rdata"), verbose = T)

# determine regional duration
coa_reform <- sim_coa_df_regions %>% 
  rename("detection_timestamp_utc" = TimeStep.coa,
         "deploy_lat" = latitude,
         "deploy_long" = longitude) %>% 
  mutate(glatos_array = regions)

# calculate time of detection events
coa_events <- detection_events(coa_reform,
                               condense = F)

# convert NA values for timediff to equal difference in time
coa_locf <- coa_events %>% 
  mutate(event2 = if_else(arrive %in% 1, event-1, event)) %>% 
  group_by(animal_id,event2) %>% 
  mutate(time_diff = duration(time_diff, units = "seconds"),
         time_diff2 = if_else(is.na(time_diff), 
                              detection_timestamp_utc - lag(detection_timestamp_utc), 
                              time_diff),
         time_diff_fin = if_else(time_diff2 > 86400,
                                 period_to_seconds(hms(format(as.POSIXct(detection_timestamp_utc),
                                                              format = "%H:%M:%S"))),
                                 as.numeric(time_diff2)))

# daily assignment
daily_dist_coa <- coa_locf %>% 
  mutate(day = yday(detection_timestamp_utc)) %>%
  rename("sim_id" = animal_id) %>% 
  group_by(sim_id,regions,day) %>% 
  reframe(points = sum(time_diff2, na.rm = T)) %>% 
  group_by(sim_id, day) %>% 
  filter(points %in% max(points)) %>% 
  reframe(points = max(points),
          regions = regions[1]) %>% 
  unique()

# # fill missing dates based on last observed region - skip because we don't want to fill undetected days
# daily_coa_full <- max_daily_regions[,c("sim_id","day")] %>% 
#   left_join(daily_dist_coa) %>% 
#   arrange(sim_id,day) %>% 
#   mutate(regions = zoo::na.locf(regions),
#          points = zoo::na.locf(points))


# select region with greatest number of points for each day/simulation
max_daily_coa <- daily_dist_coa %>% 
  # group_by(sim_id,day) %>% 
  mutate(# max_pt = max(points),
    # is_max = if_else(points %in% max_pt, 1, 0),
    model = "coa") %>% 
  # filter(is_max %in% 1) %>% 
  # ungroup() %>% 
  select(sim_id,regions,day,model) %>% 
  unique()


# LOCF
sim_locf <- readRDS(paste0(path.sim,"sim_locf_detections.rds")) 

# daily assignment
daily_dist_locf <- sim_locf %>% 
  mutate(day = yday(detection_timestamp_utc)) %>%
  rename("sim_id" = animal_id) %>% 
  group_by(sim_id,regions,day) %>% 
  reframe(points = sum(time_diff2,na.rm = T)) %>% 
  group_by(sim_id, day) %>% 
  filter(points %in% max(points)) %>% 
  reframe(points = max(points),
          regions = regions[1]) %>% 
  unique()

# fill missing dates based on last observed region
daily_locf_full <- max_daily_regions[,c("sim_id","day")] %>% 
  left_join(daily_dist_locf) %>% 
  arrange(sim_id,day) %>% 
  mutate(regions = zoo::na.locf(regions),
         points = zoo::na.locf(points))

# select region with greatest number of points for each day/simulation
max_daily_locf <- daily_locf_full %>% 
  # group_by(sim_id,day) %>% 
  mutate(# max_pt = max(points),
    # is_max = if_else(points %in% max_pt, 1, 0),
    model = "locf") %>% 
  # filter(is_max %in% 1) %>% 
  # ungroup() %>% 
  select(sim_id,regions,day,model) %>% 
  unique()


# nlINT
load(file = paste0(path.sim, "simulated_GLATOS_data.Rdata"), verbose = T)

# unique days by individual
ud <- sim_int_df_full %>% 
  mutate(day = yday(bin_timestamp)) %>% 
  group_by(animal_id) %>% 
  reframe(day_m = max(day),
          day_n = n_distinct(day)) %>% 
  unique()

# View(ud)

# daily assignment
daily_dist_int <- sim_int_df_full %>% 
  mutate(day = yday(bin_timestamp)) %>%
  rename("sim_id" = animal_id) %>% 
  group_by(sim_id,regions,day) %>% 
  reframe(points = n_distinct(bin_timestamp)) %>% 
  group_by(sim_id, day) %>% 
  filter(points %in% max(points)) %>% 
  reframe(points = max(points),
          regions = regions[1]) %>% 
  unique()

# fill missing dates based on last observed region
daily_int_full <- max_daily_regions[,c("sim_id","day")] %>% 
  left_join(daily_dist_int) %>% 
  arrange(sim_id,day) #%>% 
# mutate(regions = zoo::na.locf(regions),
#        points = zoo::na.locf(points))

# select region with greatest number of points for each day/simulation
max_daily_int <- daily_dist_int %>% 
  # group_by(sim_id,day) %>% 
  mutate(# max_pt = max(points),
    # is_max = if_else(points %in% max_pt, 1, 0),
    model = "int") %>% 
  # filter(is_max %in% 1) %>% 
  # ungroup() %>% 
  select(sim_id,regions,day,model) %>% 
  unique()


### calculate daily accuracy
# Merge daily distributions from all models and true daily assignment
daily_assn <- bind_rows(max_daily_base,
                        max_daily_coa,
                        max_daily_locf,
                        max_daily_int)

# Add true daily region assignment and calculate accuracy
daily_accy <- left_join(max_daily_regions, daily_assn) %>% 
  mutate(correct = if_else(true_region == regions, 1, 0)) %>% 
  group_by(sim_id,model) %>% 
  summarize(day_n = max(day),
            total_corr = sum(correct)) %>% 
  unique() %>% 
  mutate(per_corr = total_corr*100/day_n,
         model = factor(model,
                        levels = c("base","locf","coa","int"),
                        labels = c("Base","LOCF","COA","Int")),
         day_error = 100-per_corr) %>% 
  ungroup()


### Plot accuracy and stats
# plot error among models
day_plot <- daily_accy %>% 
  ggplot(aes(x = model, y = day_error)) +
  geom_point(shape = 21, fill = "gray40", alpha = 0.6,
             position = position_dodge2(width = 0.4)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.3) +
  ylim(c(0,100)) +
  labs(x = NULL,
       y = "Daily region assignment error (%)") +
  theme_classic()

day_plot

ggsave(plot = day_plot,
       filename = paste0(path.figs,"dayerror_4mod_Jan2024.svg"))

# summary of daily accuracy and stats comparison
daily_accy %>% 
  group_by(model) %>% 
  reframe(ave_cor = mean(day_error),
          sd_cor = sd(day_error))

shapiro.test(daily_accy$day_error)
plot(lm(day_error~model, data = daily_accy))

friedman.test(y = daily_accy$day_error, group = daily_accy$model, block = daily_accy$sim_id)

pairwise.wilcox.test(daily_accy$day_error, daily_accy$model,
                     p.adjust.method = "bonf")

### summary of daily assignments by model
# average and percent of total days assigned
sim_dur <- max_daily_regions %>% 
  group_by(sim_id) %>% 
  reframe(dur = max(day))

daily_assn %>%
  left_join(sim_dur) %>% 
  group_by(model,sim_id) %>% 
  mutate(days = n_distinct(day)) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  reframe(ave_days = mean(days),
          per_days = mean(days*100/dur),
          per_days_min = min(days*100/dur),
          per_days_max = max(days*100/dur)) %>% 
  unique()
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Total percent occurrence by region
###----------------------------------------------------------------------------------------------------
### Determine actually percent of simulations that used each region
# load simulated transmissions
true_sf_regions <- readRDS(paste0(path.sim,"simulated_trans_regions.rds"))

# convert to df
true_df <- true_sf_regions %>% 
  mutate(geometry = NULL) %>% 
  data.frame()

# calculate percent of fish that occupied each region
true_reg_oc <- true_df %>% 
  group_by(regions) %>% 
  reframe(n_sim = n_distinct(virt_fish)) %>% 
  unique()


### Determine modeled number of simulations that used each region
# load regional model estimates
model_est <- readRDS(file = paste0(path.sim, "model_regional_estimates_Feb2024.rds"))

# calculate unique simulations by region*model
model_reg_oc <- model_est %>% 
  filter(region_detect > 0) %>% 
  group_by(model, regions, .drop = F) %>% 
  reframe(n_sim_mod = n_distinct(animal_id)) %>% 
  unique()


### Calculate error between true and modeled number occupied
# join true and model data
region_occ <- model_reg_oc %>% 
  full_join(true_reg_oc) %>% 
  mutate(oc_error = n_sim_mod - n_sim,
         model = factor(model,
                        levels = c("Base","LOCF","COA","Int","Move","RSP"))) %>% 
  arrange(model)

# write.csv(region_occ,
#           file = paste0(path.result,"sims_per_region_Feb2024.csv"),
#           row.names = F)

# # plot error
# ggplot(region_occ, aes(x = model, y = oc_error)) +
#   geom_point(shape = 21, alpha = 0.8, size = 3,
#              aes(fill = regions)) +
#   scale_fill_manual(values = c('#8c510a','#f6e8c3','#bf812d','#80cdc1','#01665e','#35978f','#c7eae5'))+
#   theme_classic()

###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Evaluate error for non-zero occurrence ### SKIP
###----------------------------------------------------------------------------------------------------
#-- Error is likely largely diminished by 
#-- accurate lack of presence for individuals
#-- that don't show up in a region (lots of 
#-- accurate 0s). This section assess model
#-- accuracy when excluding true zeros to 
#-- compare error for missed detections

### Load full model output
model_full <- readRDS(paste0(path.sim, "model_error_full_Feb2024.rds"))

### Remove region/sim combos with no region occurrence
dets_fail <- model_full %>% 
  filter(total_region > 0) %>% 
  mutate(model = factor(model,
                        levels = c("Base","COA","LOCF","Int","Move","RSP")))


fail_err_plot <- dets_fail %>% 
  ggplot(aes(x = regions, y = occ_error)) +
  geom_point(fill = "gray40", shape = 21, alpha = 0.6, size = 1,
             position = position_dodge2(width = 0.4)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  facet_rep_wrap(~model,
                 ncol = 2) +
  labs(x = NULL,
       y = "Model Error (%)") +
  theme_classic()

fail_err_plot

ggsave(plot = fail_err_plot,
       filename = paste0(path.figs,"fail_mod_error_region.svg"),
       width = 7.5,
       height = 6)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Summary stats of simulated tracks
###----------------------------------------------------------------------------------------------------
### Duration of each track
sim_full_df <- readRDS(paste0(path.sim,"complete_simulated_transmissions.rds")) %>% 
  mutate(sim_id = paste0("sim_",virt_fish))

dur <- sim_full_df %>% 
  mutate(s_region = case_when(virt_fish %in% 1:20 ~ "NEA",
                              virt_fish %in% 21:40 ~ "MLN",
                              virt_fish %in% 41:60 ~ "MLC",
                              virt_fish %in% 61:80 ~ "MLS",
                              virt_fish %in% 81:100 ~ "MLB")) %>% 
  group_by(virt_fish,s_region) %>% 
  reframe(duration = max(time)/86400)

reg_dur <- dur %>% 
  group_by(s_region) %>% 
  reframe(ave_dur = mean(duration),
          sd = sd(duration))

### Unique regions used by each track
# region count
true_use <- readRDS(file = paste0(path.sim,"simulated_true_data.rds"))

n_regs <- true_use %>% 
  filter(total_region > 0) %>% 
  group_by(virt_fish) %>% 
  reframe(n_region = n_distinct(regions))

n_regs %>% 
  mutate(s_region = case_when(virt_fish %in% 1:20 ~ "NEA",
                              virt_fish %in% 21:40 ~ "MLN",
                              virt_fish %in% 41:60 ~ "MLC",
                              virt_fish %in% 61:80 ~ "MLS",
                              virt_fish %in% 81:100 ~ "MLB")) %>% 
  group_by(s_region) %>% 
  reframe(min_region = min(n_region),
          max_region = max(n_region))

# regional distribution
reg_use <- true_use %>% 
  # group_by(regions) %>% 
  # reframe(per_use = mean(percent_region),
  #         sd_use = sd(percent_region)) %>% 
  ggplot(aes(x = regions, y = percent_region)) +
  geom_point(shape = 21, fill = "gray40", alpha = 0.6,
             position = position_dodge2(width = 0.4)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.3) +
  labs(x = NULL,
       y = "Distribution (%)") +
  theme_classic()

reg_use

ggsave(plot = reg_use,
       filename = paste0(path.figs,"true_regional_dist.svg"),
       width = 7.5, height = 2.25)

ave_reg_use <- true_use %>% 
  group_by(regions) %>%
  reframe(per_use = mean(percent_region),
          sd_use = sd(percent_region))

### Number/percent of detected transmissions
dets <- readRDS(paste0(path.sim,"simulated_detections.rds"))

sim_full_df <- readRDS(paste0(path.sim,"complete_simulated_transmissions.rds")) %>% 
  mutate(sim_id = paste0("sim_",virt_fish)) %>% 
  group_by(sim_id) %>% 
  reframe(sim_full_df = length(time))

n_dets <- dets %>% 
  group_by(sim_id,start_region) %>% 
  left_join(sim_full_df) %>% 
  reframe(sim_full_df = length(detection_timestamp_utc),
          sim_full_df = unique(sim_full_df),
          p_dets = n_dets*100/sim_full_df)

n_dets %>% 
  group_by(start_region) %>% 
  reframe(ave_dets = mean(p_dets),
          sd_dets = sd(p_dets))

### Summary of metadata
metadata <- readRDS(paste0(path.sim, "simulations_metadata.rds"))

metadata %>% 
  group_by(start_region) %>% 
  reframe(ave_theta = mean(theta),
          sd_theta = sd(theta),
          ave_vel = mean(velocity),
          sd_vel = sd(velocity))
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Receiver cover
###----------------------------------------------------------------------------------------------------
### Lake region detection probability
lc_regions <- st_read(dsn = ".",
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

recs_loc <- readRDS(paste0(path.mvt,"OriginalReceiverSummary_2013-2017_May2023.rds")) %>% 
  select(StationName,deploy_lat,deploy_long) %>% 
  group_by(StationName) %>% 
  reframe(lat = first(deploy_lat),
          lon = first(deploy_long))


recs_sf <- st_as_sf(recs_loc,
                    coords = c('lon', 'lat'),
                    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

reg_recs <- st_join(recs_sf,lc_regions_short) 

reg_sum <- lc_regions_short %>% 
  data.frame() %>% 
  select(regions, area_sqkm) %>% 
  left_join(reg_recs) %>% 
  group_by(regions) %>% 
  reframe(n_rec = if_else(is.na(StationName), 0,
                          n_distinct(StationName)),
          cover = n_rec*(pi*(0.250^2)),
          rec_rate = n_rec/area_sqkm,
          per_cover = cover*100/area_sqkm) %>% 
  unique()

# write.csv(reg_sum,
#           paste0(path.result,"regional_cover.csv"),
#           row.names = F)

mean(reg_sum$per_cover)


### Relationship between coverage and error
# load model and true distribution data
model_est <- readRDS( file = paste0(path.sim, "model_regional_estimates.rds"))
true_use <- readRDS(file = paste0(path.sim,"simulated_true_data.rds"))

# calculate difference between estimated use and true use
full_diff <- true_use %>% 
  rename("true_percent" = percent_region) %>% 
  mutate(animal_id = paste0("sim_",virt_fish)) %>% 
  left_join(model_est) %>% 
  mutate(diff = region_percent - true_percent,
         model = factor(model,
                        levels = c("Base", "COA", "LOCF", "Int", "dBBMM", "RSP"),
                        labels = c("Base", "COA", "LOCF", "Int", "move", "RSP"))) %>%
  select(animal_id,model,regions,true_percent,region_percent,diff)

# summary table of difference
full_diff %>% 
  group_by(regions,model) %>% 
  reframe(mean_dff = mean(diff),
          sd_diff = sd(diff)) %>% 
  print(n = 42)

# merge with receiver cover
full_diff_rec <- full_diff %>% 
  left_join(reg_sum)

rec_cov_mod <- glmmTMB(diff~per_cover+model+(1|animal_id),
                       data = full_diff_rec)

summary(rec_cov_mod)

# plot error for all simulation*model*region combos
full_diff_plot <- full_diff_rec %>% 
  ggplot(aes(x = per_cover, y = diff))+ 
  geom_smooth(method = "lm", se = F, alpha = 0.8,
              aes(color = model)) +
  geom_hline(yintercept = 0, color = "gray75", linetype = "dashed") +
  geom_point(size = 2.5, alpha = 0.6,
             aes(fill = model, shape = model)) +
  scale_fill_manual(values = c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
  scale_shape_manual(values = c(8,21:25)) +
  scale_color_manual(values = c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
  labs(x = "Receiver Coverage (%)", y = "Average Model Error (%)") +
  theme_classic()

full_diff_plot



### Unique regions used by each track
# load model and true distribution data
model_est <- readRDS( file = paste0(path.sim, "model_regional_estimates_Feb2024.rds"))
true_use <- readRDS(file = paste0(path.sim,"simulated_true_data.rds"))

# calculate average use for true distribution and modeled distribution
ave_reg_use <- true_use %>% 
  group_by(regions) %>%
  reframe(per_use = mean(percent_region),
          sd_use = sd(percent_region),
          med_use = median(percent_region))

est_reg <- model_est %>% 
  group_by(regions,model,.drop = F) %>% 
  reframe(ave_use = mean(region_percent),
          sd_use = sd(region_percent))

# calculate average regional difference in use for each method
ave_diff <- est_reg %>% 
  left_join(ave_reg_use[,c("regions","per_use")]) %>% 
  mutate(diff = ave_use - per_use)

# merge with receiver cover
ave_diff_rec <- ave_diff %>% 
  left_join(reg_sum)

# plot difference as a function of receiver cover
rec_diff <- ave_diff_rec %>% 
  mutate(model = factor(model,
                        levels = c("Base","COA","LOCF","Int","dBBMM","RSP"),
                        labels = c("Base","COA","LOCF","Int","move","RSP"))) %>% 
  ggplot(aes(x = per_cover, y = diff))+
  # geom_smooth(method = "lm", se = F, alpha = 0.8,
  #             aes(color = model)) +
  geom_hline(yintercept = 0, color = "gray75", linetype = "dashed") +
  geom_point(aes(fill = model, shape = model),
             size = 2.5, alpha = 0.6)+
  scale_fill_manual(values = c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
  scale_shape_manual(values = c(8,21:25)) +
  # scale_color_manual(values = c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
  labs(x = "Receiver Coverage (%)", y = "Average Model Error (%)") +
  theme_classic()

rec_diff

ggsave(plot = rec_diff,
       filename = paste0(path.figs,"error_by_rec_cover.svg"),
       scale = 1.25)

aggregate(diff~per_cover+model, data = ave_diff_rec, FUN = mean)
###----------------------------------------------------------------------------------------------------
