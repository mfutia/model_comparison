###----------------------------------------------------------------------------------------------------
# Movement model assessment and comparison
# Prepared by MHF and MH
# Description: Evaluate model accuracy and compare among movement models
###----------------------------------------------------------------------------------------------------

### set working directory to local path
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


###----------------------------------------------------------------------------------------------------
### Load packages
###----------------------------------------------------------------------------------------------------
library(raster)
library(magrittr)
library(data.table)

library(sf)
library(lemon)
library(emmeans)
library(MuMIn)
library(glmmTMB)
library(tidyverse)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### load data
###----------------------------------------------------------------------------------------------------
### load lake map data
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


### load model occupancy estimates
sim_base_region_full <- readRDS(file = "outputs/model_estimates/sim_base_region_full.rds")
sim_locf_region_full <- readRDS(file = "outputs/model_estimates/sim_locf_region_full.rds")
sim_coa_region_full <- readRDS(file = "outputs/model_estimates/sim_coa_region_full.rds")
sim_int_region_full <- readRDS(file = "outputs/model_estimates/sim_int_region_full.rds")
sim_move_region_full <- readRDS(file = "outputs/model_estimates/sim_move_region_full.rds")
sim_rsp_region_full <- readRDS(file = "outputs/model_estimates/sim_rsp_region_full.rds")

# load true occupancy data
true_occ <- readRDS(file = "outputs/Simulations/simulated_true_data.rds")

# load metadata
metadata <- readRDS(file = "outputs/Simulations/simulations_metadata.rds") %>% 
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

# select required columns
model_est_final <- model_est %>%
  select(animal_id, regions, region_percent, model)

# add true region use and calculate total and regional simulation duration
model_full <- true_occ %>% 
  select(virt_fish,regions,region_positions,total_positions) %>% 
  mutate(animal_id = paste0("sim_", virt_fish)) %>% # create animal_id column for merging tables
  right_join(model_est_final) %>% 
  mutate(region_dur = region_positions*127, # calculate duration by region for each simulation (127 = transmission delay (120 sec) + transmission duration (7 sec))   
         total_dur = total_positions*127, # calculate total duration of simulation in seconds
         region_wt = region_dur/total_dur, # calculate proportion of time spent in each region
         abs_error_per = abs(region_percent-(region_wt*100)), # calculate absolute model error for each sim/region 
         act_error_per = region_percent-(region_wt*100),
         wt_abs_error = abs_error_per*region_wt, # calculate occupancy error for each sim/region/model
         wt_act_error = act_error_per*region_wt) # calculate actual weighted error

str(model_full)
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Criterion 1: Weighted mean occupancy error (MOE) 
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
ave_reg_use <- true_occ %>% 
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
full_diff <- true_occ %>% 
  rename("true_percent" = percent_region_use) %>% 
  left_join(model_est) %>% 
  mutate(diff = region_percent - true_percent,
         model = factor(model,
                        levels = c("Base", "LOCF", "COA", "Int", "Move", "RSP"))) %>%
  select(animal_id,model,regions,true_percent,region_percent,diff)


### map filling region by error
map_error <- lc_regions_short %>% 
  left_join(ave_err[,c("regions", "model", "act_mean")]) %>%
  ggplot() +
  geom_sf(aes(fill = act_mean)) +
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

# emmeans to test pairwise comparisons
mod_emm_abs <- emmeans(best_meta_abs, specs = c("regions","model"))
mod_p_abs <- pwpm(mod_emm_abs)

mod_p_abs_df <- matrix(mod_p_abs,
                       ncol = 42,
                       dimnames = list(rownames(mod_p_abs),rownames(mod_p_abs))) %>% 
  data.frame()
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Criterion 4: Daily region assignment error
###----------------------------------------------------------------------------------------------------
### Determine actually daily zone assignment
# load simulated transmissions
true_occ_sf <- readRDS("outputs/simulations/simulated_occ_sf.rds")

# convert to df
true_occ_df <- true_occ_sf %>% 
  mutate(geometry = NULL) %>% 
  data.frame()

# quantify daily region transmissions
true_daily_regions <- true_occ_df %>% 
  mutate(day_n = ceiling(time/86400)) %>%
  group_by(virt_fish,regions,day_n) %>% 
  summarize(points = n_distinct(time)) %>% 
  ungroup() %>% 
  unique()

# select region with greatest number of points for each day/simulation
max_daily_regions <- true_daily_regions %>% 
  group_by(virt_fish,day_n) %>% 
  mutate(max_pt = max(points),
         is_max = if_else(points %in% max_pt, 1, 0),
         sim_id = paste0("sim_",virt_fish),
         sim_id = factor(sim_id,
                         levels = unique(sim_id))) %>% 
  ungroup() %>% 
  group_by(sim_id, day_n) %>% 
  filter(is_max %in% 1) %>% 
  reframe(regions = regions[1]) %>% 
  rename(true_region = regions) %>% 
  unique()

### Determine daily region assignment for each simulation/model
# Base model
sim_base <- readRDS("outputs/model_estimates/sim_base.rds") 

sim_base_daily <- sim_base %>% 
  mutate(day_n = yday(detection_timestamp_utc)) %>% 
  group_by(day_n, regions, sim_id) %>% 
  reframe(n_dets = n_distinct(detection_timestamp_utc)) %>% 
  group_by(day_n, sim_id) %>% 
  filter(n_dets == max(n_dets)) %>% 
  ungroup() %>% 
  mutate(model = "base") %>% 
  select(sim_id,regions,day_n,model) %>%
  unique()

# LOCF model
sim_locf <- readRDS("outputs/model_estimates/sim_locf.rds")

sim_locf_daily <- sim_locf %>% 
  mutate(day_n = yday(detection_timestamp_utc)) %>%
  rename("sim_id" = animal_id) %>% 
  group_by(sim_id, regions, day_n) %>% 
  reframe(duration = sum(time_diff2,na.rm = T)) %>% 
  group_by(day_n, sim_id) %>% 
  filter(duration == max(duration)) %>% 
  ungroup() %>% 
  select(sim_id,regions,day_n) %>%
  unique() %>% 
  right_join(max_daily_regions[,c("sim_id","day_n")]) %>% # add all dates
  arrange(sim_id,day_n) %>% 
  mutate(regions = zoo::na.locf(regions), # fill region for added dates based on last observed region
         model = "locf") 

# COA model
sim_coa_df_regions <- readRDS("outputs/model_estimates/sim_coa_regions.rds")

sim_coa_daily <- sim_coa_df_regions %>% 
  rename(sim_id = animal_id) %>% 
  mutate(day_n = yday(TimeStep.coa)) %>% 
  group_by(day_n, regions, sim_id) %>% 
  reframe(n_dets = n_distinct(TimeStep.coa)) %>% 
  group_by(day_n, sim_id) %>% 
  filter(n_dets == max(n_dets)) %>% 
  ungroup() %>% 
  mutate(model = "coa") %>% 
  select(sim_id,regions,day_n,model) %>%
  unique()

# Int Model
sim_int_df_regions <- readRDS("outputs/model_estimates/sim_int_regions.rds")

sim_int_daily <- sim_int_df_regions %>% 
  rename(sim_id = animal_id) %>% 
  mutate(day_n = yday(bin_timestamp)) %>% 
  group_by(day_n, regions, sim_id) %>% 
  reframe(n_dets = n_distinct(bin_timestamp)) %>% 
  group_by(day_n, sim_id) %>% 
  filter(n_dets == max(n_dets)) %>% 
  ungroup() %>% 
  mutate(model = "int") %>% 
  select(sim_id,regions,day_n,model) %>%
  unique()


### calculate daily accuracy
# Merge daily distributions from all models and true daily assignment
daily_assn <- bind_rows(sim_base_daily,
                        sim_locf_daily,
                        sim_coa_daily,
                        sim_int_daily)

# Add true daily region assignment and calculate accuracy
daily_accy <- left_join(max_daily_regions, daily_assn) %>% 
  mutate(correct = if_else(true_region == regions, 1, 0)) %>% 
  group_by(sim_id, model) %>% 
  summarize(day_n = max(day_n),
            total_corr = sum(correct, na.rm = T)) %>% 
  unique() %>% 
  mutate(per_corr = total_corr*100/day_n,
         day_error = 100-per_corr) %>% 
  ungroup() %>% 
  mutate(model = factor(model,
                        levels = c("base", "locf", "coa", "int")))


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

# ggsave(plot = day_plot,
#        filename = paste0(path.figs,"dayerror_4mod_Mar2024.svg"))

# summary of daily accuracy and stats comparison
daily_accy %>% 
  group_by(model) %>% 
  reframe(ave_err = mean(day_error),
          sd_err = sd(day_error))

friedman.test(y = daily_accy$day_error, group = daily_accy$model, block = daily_accy$sim_id)

pairwise.wilcox.test(daily_accy$day_error, daily_accy$model,
                     p.adjust.method = "bonf")

### summary of daily assignments by model
# average and percent of total days assigned
sim_dur <- max_daily_regions %>% 
  group_by(sim_id) %>% 
  reframe(dur = max(day_n))

daily_assn %>%
  left_join(sim_dur) %>% 
  group_by(model,sim_id) %>% 
  mutate(days = n_distinct(day_n)) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  reframe(ave_days = mean(days),
          per_days = mean(days*100/dur),
          per_days_min = min(days*100/dur),
          per_days_max = max(days*100/dur)) %>% 
  unique()
###----------------------------------------------------------------------------------------------------


###----------------------------------------------------------------------------------------------------
### Criterion 5: Within-region occurrence error
###----------------------------------------------------------------------------------------------------
### Determine  number of simulations that used each region
# calculate actual percent of fish that occupied each region
true_reg_occ <- true_occ_df %>% 
  group_by(regions) %>% 
  reframe(n_sim = n_distinct(virt_fish)) %>% 
  unique()

# calculate unique simulations by region*model
model_reg_occ <- model_est %>% 
  filter(region_detect > 0) %>% 
  group_by(model, regions, .drop = F) %>% 
  reframe(n_sim_mod = n_distinct(animal_id)) %>% 
  unique()


### Calculate error between true and modeled number occupied
# join true and model data
region_occ <- model_reg_occ %>% 
  full_join(true_reg_occ) %>% 
  mutate(occ_error = n_sim_mod - n_sim,
         model = factor(model,
                        levels = c("Base","LOCF","COA","Int","Move","RSP"))) %>% 
  arrange(model)

summary(region_occ)
###----------------------------------------------------------------------------------------------------
