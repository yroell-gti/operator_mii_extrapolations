# -----------------------------------------------------------------------------
# Script Name   : EST_createsyntheticdata.R
# Author        : Yannik Roell
# Last Modified : 2025-06-01
#
# Description   :
#   Create synthetic data from real data by adding noise
# -----------------------------------------------------------------------------

library(dplyr)
library(lubridate)

### datasets loaded in
# load production data
production = read.csv("../../../Data/Production/production_ghgrp.csv")

# load aerial datasets for both flights, renaming and selecting necessary columns
aerialscan1 = read.csv("../../../Data/Aerial/facility_first_scan.csv") %>%
  rename(pad_id = asset_number) %>%
  select(-job_id, -facility_name)
aerialscan2 = read.csv("../../../Data/Aerial/facility_second_scan.csv") %>%
  rename(pad_id = asset_number) %>%
  select(-job_id, -facility_name)

# merge aerial datasets and filter for sites with emission data in both months
aerial = aerialscan1 %>%
  full_join(aerialscan2, by = "pad_id", suffix = c("_scan1", "_scan2")) %>%
  filter(!is.na(sum_of_emissions_first_scan) & !is.na(sum_of_emissions_second_scan))

# load and preprocess PSN A data
psna = read.csv("../../../Data/PSNA/psna_data.txt") %>%
  rename(pad_id = SHORT_NAME,
         emission = EMISSION,
         methane = METHANE) %>%
  mutate(timestamp = dmy_hms(DATE)) %>%
  select(pad_id, timestamp, emission, methane) %>%
  arrange(pad_id, timestamp, desc(emission)) %>%
  distinct(pad_id, timestamp, .keep_all = TRUE)   # remove duplicate entries

# load and preprocess PSN B data
psnb = read.csv("../../../Data/PSNB/psnb_data.csv") %>%
  rename(pad_id = site_alias,
         emission = rate_kg_per_hour,
         volume = volume_kg) %>%
  mutate(timestamp = ymd_hms(time)) %>%
  select(pad_id, timestamp, emission, volume)

### generate synthetic site list
# define number of sites per category
aerial_count = 100
psna_count = 50
psnb_count = 50

# load emission source data created while working on milestone
list_compared = read.csv("../../Data/all_emission_sources.csv")

# select sites for aerial
aerial100 = list_compared %>%
  filter(is.na(psna_nan0assumpt_ch4tonnes) & is.na(psnb_ch4tonnes) & !is.na(ghgrp_ch4tonnes),
         pad_id %in% unique(production$pad_id)) %>%
  filter(pad_id %in% aerial$pad_id) %>%
  slice_sample(n = aerial_count, replace = FALSE) %>%
  mutate(source = "aerial") %>%
  select(pad_id, source)

# select sites for PSN A
psna50 = list_compared %>%
  filter(!is.na(psna_nan0assumpt_ch4tonnes) & is.na(psnb_ch4tonnes) & !is.na(ghgrp_ch4tonnes),
         pad_id %in% unique(production$pad_id)) %>%
  filter(pad_id %in% aerial$pad_id) %>%
  slice_sample(n = psna_count, replace = FALSE) %>%
  mutate(source = "psna") %>%
  select(pad_id, source)

# select sites for PSN B
psnb50 = list_compared %>%
  filter(is.na(psna_nan0assumpt_ch4tonnes) & !is.na(psnb_ch4tonnes) & !is.na(ghgrp_ch4tonnes),
         pad_id %in% unique(production$pad_id)) %>%
  filter(pad_id %in% aerial$pad_id) %>%
  slice_sample(n = psnb_count, replace = FALSE) %>%
  mutate(source = "psnb") %>%
  select(pad_id, source)

# combine all selected sites into one list
site_list = rbind(aerial100, psna50, psnb50)

### getting real data
# extract real production for selected sites
real_production = production %>%
  rename(ghgrp = GHGRP23_CH4_TONNES_CO2E_PER_YEAR) %>%
  filter(pad_id %in% unique(site_list$pad_id)) %>%
  select(pad_id, ghgrp, starts_with("total"), starts_with("age"))

# extract real emission data for selected sites
real_aerial = aerial %>% filter(pad_id %in% unique(site_list$pad_id))
real_psna = psna %>% filter(pad_id %in% unique(site_list$pad_id))
real_psnb = psnb %>% filter(pad_id %in% unique(site_list$pad_id))
  
### generate synthetic reduction values per site
# define reduction parameters
rnorm_mean = 0.45
rnorm_sd = 0.025
monthly_noise = 0.05
beta = 3
red_factor_offset_value = 0.35
bandwidth = 8

# compute reduction factors
red_factor = rnorm(length(unique(site_list$pad_id)), mean = rnorm_mean, sd = rnorm_sd)
alpha = red_factor * beta / (1 - red_factor)

# assign reduction factors to site list
site_list$red_factor = red_factor + red_factor_offset_value
site_list$beta = beta
site_list$alpha = alpha

### generate synthetic data for production, aerial, PSN A, and PSN B
# adjust production emissions based on reduction factor
syn_production = merge(real_production, site_list, by = "pad_id") %>%
  select(-beta, -alpha, -total_prod_2023) %>%
  rowwise() %>%
  mutate(ghgrp = ghgrp * red_factor,
         across(starts_with("total"), ~ . * (red_factor + runif(1, -monthly_noise, monthly_noise)))) %>%
  ungroup() %>%
  mutate(across(starts_with("total"), ~ replace(., is.na(.), 0)),
         age2023_min = round(age2023_min),
         age2023_max = round(age2023_max),
         age2023_mean = round(age2023_mean)) %>%
  mutate(total_prod_2023 = rowSums(select(., starts_with("total"))))

# adjust aerial emissions based on reduction factor
syn_aerial = merge(real_aerial, site_list, by = "pad_id") %>%
  mutate(sum_of_emissions_first_scan = sum_of_emissions_first_scan * (red_factor + runif(nrow(real_aerial), -monthly_noise, monthly_noise)),
         sum_of_emissions_second_scan = sum_of_emissions_second_scan * (red_factor + runif(nrow(real_aerial), -monthly_noise, monthly_noise))) %>%
  select(-beta, -alpha)

# process PSN A synthetic emissions using beta distribution and smoothing
syn_psna = merge(real_psna, site_list, by = "pad_id") %>%
  arrange(pad_id, timestamp)
syn_psna$smwn = NA
syn_psna$new_emission = NA

psna_uniquesites = site_list %>%
  filter(source == "psna")

for (i in seq_along(psna_uniquesites$pad_id)) {
  site_id = psna_uniquesites$pad_id[i]
  site_rows = syn_psna$pad_id == site_id
  
  n = sum(site_rows)
  wn = rbeta(n, shape1 = syn_psna$alpha[i], shape2 = beta) + red_factor_offset_value
  smwn = ksmooth(1:n, wn, bandwidth = bandwidth)
  
  syn_psna$smwn[site_rows] = smwn$y
  syn_psna$new_emission[site_rows] = syn_psna$smwn[site_rows] * syn_psna$emission[site_rows]
}

syn_psna = syn_psna %>%
  select(-beta, -alpha, -emission, -methane) %>%
  rename(emission = new_emission) %>%
  select(pad_id, timestamp, emission, red_factor, smwn, source)

# process PSN B synthetic emissions using beta distribution and smoothing
syn_psnb = merge(real_psnb, site_list, by = "pad_id") %>%
  arrange(pad_id, timestamp)
syn_psnb$smwn = NA
syn_psnb$new_emission = NA

psnb_uniquesites = site_list %>%
  filter(source == "psnb")

for (i in seq_along(psnb_uniquesites$pad_id)) {
  site_id = psnb_uniquesites$pad_id[i]
  site_rows = syn_psnb$pad_id == site_id
  
  n = sum(site_rows)
  wn = rbeta(n, shape1 = syn_psnb$alpha[i], shape2 = beta) + red_factor_offset_value
  smwn = ksmooth(1:n, wn, bandwidth = bandwidth)
  
  if (length(smwn$y) == n) {
    syn_psnb$smwn[site_rows] = smwn$y
    syn_psnb$new_emission[site_rows] = syn_psnb$smwn[site_rows] * syn_psnb$emission[site_rows]
  } else {
    # a few sites only had a few entries so smoothing function didn't work the same; 
    # only multiplying by reduction factor and not extra noise from smoothing function
    syn_psnb$new_emission[site_rows] = syn_psnb$red_factor[site_rows] * syn_psnb$emission[site_rows] 
  }
}

syn_psnb = syn_psnb %>%
  select(-beta, -alpha, -emission, -volume) %>%
  rename(emission = new_emission) %>%
  select(pad_id, timestamp, emission, red_factor, smwn, source)
  
### determine install date for continous monitors
cm_installdates = rbind(syn_psna, syn_psnb) %>%
  group_by(pad_id, source) %>%
  summarize(first_date = as.Date(min(timestamp))) %>%
  ungroup() %>%
  mutate(random_days = if_else(source == "psnb", 
                               sample(0:90, n(), replace = TRUE), 
                               0),
         install_date = first_date - days(random_days),
         install_date = if_else(install_date < as.Date("DATE"),   # actual date has been replaced
                                as.Date("DATE"),   # actual date has been replaced
                                install_date)) %>%
  select(pad_id, install_date)
site_list = site_list %>%
  left_join(cm_installdates, by = "pad_id")

### save real and synthetic data into rds
saveRDS(list(real_production = real_production, real_aerial = real_aerial,
             real_psna = real_psna, real_psnb = real_psnb,
             syn_production = syn_production, syn_aerial = syn_aerial,
             syn_psna = syn_psna, syn_psnb = syn_psnb,
             site_list = site_list %>% select(-beta, -alpha)), "../AnalyzeSyntheticData/synthetic_data.rds")
