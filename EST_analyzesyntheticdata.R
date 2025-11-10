# -----------------------------------------------------------------------------
# Script Name   : EST_analyzesyntheticdata.R
# Author        : Yannik Roell
# Created Date  : 2025-06-01
#
# Description   :
#   Analyze synthetic data to get emissions for various extrapolation methods
# -----------------------------------------------------------------------------

##### Load and process data #####

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(purrr)

# load in synthetic data
site_list = read.csv("site_list.csv") %>%
  mutate(
    install_date_clean = gsub('=\\"|\\"', '', install_date),
    install_date_parsed = as.Date(paste0("2023-", install_date_clean), format = "%Y-%m-%d"),
    install_month = month(install_date_parsed)
  ) %>%
  select(pad_id, source, install_date_parsed, install_month) %>%
  rename(install_date = install_date_parsed)
syn_production_raw = read.csv("syn_production.csv")
syn_aerial_raw = read.csv("syn_aerial.csv")
syn_psna_raw = read.csv("syn_psna.csv")
syn_psnb_raw = read.csv("syn_psnb.csv")

# process synthetic production data
ghgrp_to_tonnes_conversion = 28   # https://www.epa.gov/system/files/documents/2025-01/ghg-emission-factors-hub-2025.pdf
prod_to_boe_conversion = 6   # https://www.eia.gov/petroleum/wells/pdf/WDR2024_Full%20Report.pdf
prod_to_tonnes_conversion = 0.0192    # https://www.eei.org/-/media/Project/EEI/Documents/Issues-and-Policy/NGSI_MethaneIntensityProtocol.pd
syn_production = syn_production_raw %>%
  mutate(
    ghgrp = ghgrp / ghgrp_to_tonnes_conversion,   # convert ghgrp to tonnes
    prod_strata = as.factor(
      case_when(
        # categorize production levels
        total_prod / prod_to_boe_conversion / 365 == 0 |
          is.na(total_prod) ~ "Inactive",
        total_prod / prod_to_boe_conversion / 365 > 0 &
          total_prod / prod_to_boe_conversion / 365 < 15 ~ "Marginal",
        total_prod / prod_to_boe_conversion / 365 >= 15 &
          total_prod / prod_to_boe_conversion / 365 < 300 ~ "Standard",
        total_prod / prod_to_boe_conversion / 365 >= 300 ~ "High"
      )
    )
  )

# process synthetic aerial data
syn_aerial = syn_aerial_raw %>%
  mutate(row_id = row_number()) %>%
  mutate(
    sum_of_emissions_first_scan = if_else(
      row_id %in% sample(row_id[sum_of_emissions_first_scan == 0], 46),
      0.922,
      sum_of_emissions_first_scan
    )
  ) %>%
  mutate(
    sum_of_emissions_second_scan = if_else(
      row_id %in% sample(row_id[sum_of_emissions_second_scan == 0], 26),
      0.922,
      sum_of_emissions_second_scan
    )
  ) %>%
  mutate(
    sum_of_emissions_first_scan = sum_of_emissions_first_scan / 1000,   # convert scans to tonnes
    sum_of_emissions_second_scan = sum_of_emissions_second_scan / 1000   # convert scans to tonnes
  ) 

# process synthetic PSN A data
syn_psna = syn_psna_raw %>%
  mutate(timestamp_clean = gsub('=\\"|\\"', '', timestamp),
    timestamp_parsed = ymd_hms(paste0("2023-", timestamp_clean)),
    month = month(timestamp_parsed)) %>%
  mutate(emission = if_else(emission < 0, 0, emission)) %>%   # change negative readings to 0
  mutate(emission = emission / 4 / 1000) %>%   # convert readings to tonnes
  select(pad_id, source, timestamp_parsed, month, emission) %>%
  rename(timestamp = timestamp_parsed)
  
# define monthly time intervals in 15-minute increments
monthly_intervals = data.frame(month = c(1:12),
                               days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
monthly_intervals$intervals_15min = monthly_intervals$days * 24 * 4
monthly_intervals = monthly_intervals %>%
  mutate(month = ym(paste0("2023-",month)))

# process synthetic PSN B data
syn_psnb = syn_psnb_raw %>%
  mutate(timestamp_clean = gsub('=\\"|\\"', '', timestamp),
         timestamp_parsed = ymd_hms(paste0("2023-", timestamp_clean)),
         month = month(timestamp_parsed)) %>%
  mutate(emission = if_else(emission < 0, 0, emission)) %>%
  select(pad_id, source, timestamp_parsed, month, emission) %>%
  rename(timestamp = timestamp_parsed)

generate_pad_timestamps = function(pad_id, install_date) {
  current_month = floor_date(install_date, "month")
  final_timestamps = c()
  
  while (format(current_month, "%Y-%m") %in% format(monthly_intervals$month, "%Y-%m")) {
    # Lookup month data
    month_info = monthly_intervals %>% filter(month == current_month)
    if (nrow(month_info) == 0) break
    
    # Start and end times
    if (current_month == floor_date(install_date, "month")) {
      start_time = as.POSIXct(install_date, tz = "UTC")
    } else {
      start_time = as.POSIXct(current_month, tz = "UTC")
    }
    
    end_time = as.POSIXct(ceiling_date(current_month + months(1), "month") - minutes(15), tz = "UTC")
    
    # Generate sequence
    timestamps = seq(from = start_time, to = end_time, by = "15 min")
    final_timestamps = c(final_timestamps, timestamps)
    
    current_month = current_month %m+% months(1)
  }
  
  tibble(
    pad_id = pad_id,
    timestamp = as.POSIXct(final_timestamps, origin = "1970-01-01", tz = "UTC")
  )
}

psnb_site_list = site_list %>% filter(source == "psnb")

# Generate for all PSN B pad_ids
full_time_series = pmap_dfr(
  list(psnb_site_list$pad_id, psnb_site_list$install_date),
  generate_pad_timestamps
)

final_emissions <- full_time_series %>%
  left_join(syn_psnb, by = c("pad_id", "timestamp")) %>%
  mutate(emission = replace_na(emission, 0),) %>%
  arrange(pad_id, timestamp, desc(emission)) %>%
  distinct(pad_id, timestamp, .keep_all = TRUE)

# Set 67.1% of rows to NA
zero_rows <- which(final_emissions$emission == 0)
n_zero <- nrow(final_emissions)
n_na <- floor(0.671 * sum(final_emissions$emission == 0))
na_indices <- sample(zero_rows, n_na)
final_emissions$emission[na_indices] <- NA_real_

syn_psnb = final_emissions %>%
  mutate(month = month(timestamp)) %>%
  mutate(emission = emission / 4 / 1000)   # convert readings to tonnes

# define monthly time intervals in 15-minute increments
monthly_intervals = data.frame(month = c(1:12),
                               days = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
monthly_intervals$intervals_15min = monthly_intervals$days * 24 * 4

##### Data end #####

##### Temporal assumption #####

### AERIAL
aerial_calculations = syn_aerial %>%
  left_join(syn_production, by = "pad_id") %>%
  mutate(
    first_half_producing_months = rowSums(across(
      matches("Jan|Feb|Mar|Apr|May|Jun"),
      ~ !is.na(.) & . != 0
    )),
    second_half_producing_months = rowSums(across(
      matches("Jul|Aug|Sep|Oct|Nov|Dec"),
      ~ !is.na(.) & . != 0
    )),
    
    # extrapolate emissions using 6 month assumption
    aerial_emissions_first_6monthassumpt = sum_of_emissions_first_scan * 6 * 730,
    aerial_emissions_second_6monthassumpt = sum_of_emissions_second_scan * 6 * 730,
    
    # extrapolate emissions based on production months
    aerial_emissions_first_monthconsider = sum_of_emissions_first_scan * first_half_producing_months * 730,
    aerial_emissions_second_monthconsider = sum_of_emissions_second_scan * second_half_producing_months * 730,
    
    # total emissions for each method
    aerial_emissions_total_6monthassumpt = aerial_emissions_first_6monthassumpt + aerial_emissions_second_6monthassumpt,
    aerial_emissions_total_monthconsider = aerial_emissions_first_monthconsider + aerial_emissions_second_monthconsider
  ) %>%
  mutate(
    # Reference productions from flyover months
    feb_prod = total_gas_prod_mcf_vol_Feb,
    aug_prod = total_gas_prod_mcf_vol_Aug,

    # Avoid divide-by-zero errors
    feb_prod_safe = if_else(is.na(feb_prod) | feb_prod == 0, NA_real_, feb_prod),
    aug_prod_safe = if_else(is.na(aug_prod) | aug_prod == 0, NA_real_, aug_prod),

    # Emission per unit production in Feb and Aug
    feb_emission_per_unit = sum_of_emissions_first_scan / feb_prod_safe,
    aug_emission_per_unit = sum_of_emissions_second_scan / aug_prod_safe
  ) %>%

  rowwise() %>%
  mutate(
    # Scale emissions by production per month
    aerial_emissions_first_prodscaled = sum(c_across(
      c(total_gas_prod_mcf_vol_Jan,
        total_gas_prod_mcf_vol_Feb,
        total_gas_prod_mcf_vol_Mar,
        total_gas_prod_mcf_vol_Apr,
        total_gas_prod_mcf_vol_May,
        total_gas_prod_mcf_vol_Jun)) * feb_emission_per_unit, na.rm = TRUE) * 730,
    aerial_emissions_second_prodscaled = sum(c_across(
      c(total_gas_prod_mcf_vol_Jul,
        total_gas_prod_mcf_vol_Aug,
        total_gas_prod_mcf_vol_Sep,
        total_gas_prod_mcf_vol_Oct,
        total_gas_prod_mcf_vol_Nov,
        total_gas_prod_mcf_vol_Dec)) * aug_emission_per_unit, na.rm = TRUE) * 730,

    aerial_emissions_total_prodscaled = aerial_emissions_first_prodscaled + aerial_emissions_second_prodscaled
  ) %>%
  ungroup()

aerial_bysource = aerial_calculations %>%
  group_by(source.x) %>%
  summarize(aerial_6month_bysource = sum(aerial_emissions_total_6monthassumpt),
            aerial_monthconsider_bysource = sum(aerial_emissions_total_monthconsider),
            aerial_prodscaled_bysource = sum(aerial_emissions_total_prodscaled))


aerial_prodscaled = syn_aerial %>%
  left_join(syn_production, by = "pad_id") %>%
  summarize(sum_feb_emission = sum(sum_of_emissions_first_scan),
            sum_aug_emission = sum(sum_of_emissions_second_scan),
            sum_feb_production = sum(total_gas_prod_mcf_vol_Feb),
            sum_aug_production = sum(total_gas_prod_mcf_vol_Aug)) %>%
  mutate(
    feb_emission_per_unit = sum_feb_emission / sum_feb_production,
    aug_emission_per_unit = sum_aug_emission / sum_aug_production
  )

# Now apply these values across the full production data
total_prod_scaled = syn_production %>%
  summarize(
    total_first_half_prod = sum(c_across(matches("Jan|Feb|Mar|Apr|May|Jun")), na.rm = TRUE),
    total_second_half_prod = sum(c_across(matches("Jul|Aug|Sep|Oct|Nov|Dec")), na.rm = TRUE)
  ) %>%
  bind_cols(aerial_prodscaled) %>%
  mutate(
    aerial_emissions_first_prodscaled_agg = total_first_half_prod * feb_emission_per_unit * 730,
    aerial_emissions_second_prodscaled_agg = total_second_half_prod * aug_emission_per_unit * 730,
    aerial_emissions_total_prodscaled_agg = aerial_emissions_first_prodscaled_agg + aerial_emissions_second_prodscaled_agg
  )

### PSN A
# aggregate PSN A by month
psna_monthly = syn_psna %>%
  group_by(pad_id, month) %>%
  summarize(
    total = sum(emission, na.rm = TRUE),
    intervals = n(),
    missing_intervals = sum(is.na(emission))
  ) %>%
  ungroup() %>%
  inner_join(monthly_intervals, by = "month") %>%
  left_join(site_list, by = "pad_id") %>%
  mutate(expected_intervals = (days - day(install_date) + 1) * 24 * 4)

# account for missing months
psna_distinctmonths = syn_psna %>%
  group_by(pad_id) %>%
  summarize(cm_months = n_distinct(month)) %>%
  ungroup() %>%
  left_join(syn_production, by = "pad_id") %>%
  mutate(
    prod_months = rowSums(!is.na(across(
      starts_with("total_gas")
    )) & across(starts_with("total_gas")) != 0),
    total_months = if_else(cm_months < prod_months, prod_months, cm_months),
    unmonitored_months = if_else(cm_months < prod_months, prod_months - cm_months, 0)
  ) %>%
  select(pad_id,
         cm_months,
         prod_months,
         total_months,
         unmonitored_months)

# naive scaling
psna_emission_ns = sum(psna_monthly$total)
mm_avg_emission_ns = psna_emission_ns / (
  sum(psna_distinctmonths$total_months) - sum(psna_distinctmonths$unmonitored_months)
)
mm_u_emission_ns = mm_avg_emission_ns * sum(psna_distinctmonths$unmonitored_months)
psna_emission_ns = psna_emission_ns + mm_u_emission_ns

# observed scaling for all
psna_avg_emission = sum(psna_monthly$total) / (sum(psna_monthly$intervals) - sum(psna_monthly$missing_intervals)) # get average emissions per interval observed
u_emission_os = psna_avg_emission * sum(psna_monthly$missing_intervals) # multiply average emission by intervals not monitored
t_emission_os = sum(psna_monthly$total) + u_emission_os # add monitored and unmonitored emissions for total
mm_avg_emission_os = t_emission_os / (
  sum(psna_distinctmonths$total_months) - sum(psna_distinctmonths$unmonitored_months)
) # get average emissions per observed months
mm_u_emission_os = mm_avg_emission_os * sum(psna_distinctmonths$unmonitored_months) # get emissions missed during unmonitored months
psna_emission_os = t_emission_os + mm_u_emission_os # add total emissions of monitored months with total emissions of unmonitored months

# observed scaling by site
psna_emission_site = psna_monthly %>%
  mutate(
    psna_avg_emission_site = total / (intervals - missing_intervals),
    u_emission_site = psna_avg_emission_site * missing_intervals,
    t_emission_site = total + if_else(is.na(u_emission_site), 0, u_emission_site)
  )

psna_emission_site_allmonths = psna_monthly %>%
  mutate(
    psna_avg_emission_site = total / (intervals - missing_intervals),
    u_emission_site = psna_avg_emission_site * missing_intervals,
    t_emission_site = total + if_else(is.na(u_emission_site), 0, u_emission_site)
  ) %>%
  group_by(pad_id) %>%
  summarize(psna_monitored_emission_site = sum(t_emission_site)) %>%
  ungroup() %>%
  left_join(psna_distinctmonths, by = "pad_id") %>%
  mutate(
    mm_avg_emission_site = psna_monitored_emission_site / (total_months - unmonitored_months),
    mm_u_emission_site = mm_avg_emission_site * unmonitored_months,
    psna_emission_site = psna_monitored_emission_site + mm_u_emission_site
  )

# linear scaling
fill_na_with_average <- function(data, leading_w_0 = TRUE) {
  if (leading_w_0) {
    # Handle leading NAs: Assume they are 0
    if (is.na(data[1])) {
      data[1] = 0
    }
    if (is.na(data[length(data)])) {
      data[length(data)] = 0
    }
  }
  
  filled_data <- data  # Copy the data
  
  for (i in seq_along(data)) {
    if (is.na(data[i])) {
      # Look back for the most recent non-NA value
      j = i - 1
      while (j > 0 && is.na(data[j])) {
        j = j - 1
      }
      
      # Look forward for the next non-NA value
      k = i + 1
      while (k <= length(data) && is.na(data[k])) {
        k = k + 1
      }
      
      # Fill the NA with the average of the previous and next non-NA values
      if (j > 0 && k <= length(data)) {
        filled_data[i] = mean(c(data[j], data[k]), na.rm = TRUE)
      } else if (j > 0) {
        # If no next non-NA value exists, use the previous value
        filled_data[i] = data[j]
      } else if (k <= length(data)) {
        # If no previous non-NA value exists, use the next value
        filled_data[i] = data[k]
      }
    }
  }
  
  return(filled_data)
}

emission_normal = c()
emission_filled_leadingT = c()
emission_filled_leadingF = c()
n = 1
for (id in unique(syn_psna$pad_id)) {
  print(paste(id, n, sep = " - "))
  subset_data = syn_psna[syn_psna$pad_id == id, ]
  subset_data$leadingT = fill_na_with_average(subset_data$emission, leading_w_0 = TRUE)
  subset_data$leadingF = fill_na_with_average(subset_data$emission, leading_w_0 = FALSE)
  emission_normal[n] = sum(subset_data$emission, na.rm = TRUE)
  emission_filled_leadingT[n] = sum(subset_data$leadingT, na.rm = TRUE)
  emission_filled_leadingF[n] = sum(subset_data$leadingF, na.rm = TRUE)
  
  n = n + 1
}
sum(emission_normal)
linear_leadingw0_T = sum(emission_filled_leadingT)
linear_leadingw0_F = sum(emission_filled_leadingF)

linear_avg_mm = linear_leadingw0_F / (
  sum(psna_distinctmonths$total_months) - sum(psna_distinctmonths$unmonitored_months)
)
linear_mm_u = linear_avg_mm * sum(psna_distinctmonths$unmonitored_months)
linear_emissions = linear_leadingw0_F + linear_mm_u

linear_sites = data.frame(
  pad_id = unique(syn_psna$pad_id),
  linear_emission = emission_filled_leadingF
) %>%
  left_join(psna_distinctmonths, by = "pad_id") %>%
  mutate(
    mm_avg_emission_site = linear_emission / (total_months - unmonitored_months),
    mm_u_emission_site = mm_avg_emission_site * unmonitored_months,
    linear_emission_site = linear_emission + mm_u_emission_site
  )

### PSN B
# aggregate PSN B by month
psnb_monthly = syn_psnb %>%
  group_by(pad_id, month) %>%
  summarize(
    total = sum(emission, na.rm = TRUE),
    intervals = n(),
    missing_intervals = sum(is.na(emission))
  ) %>%
  ungroup() %>%
  inner_join(monthly_intervals, by = "month") %>%
  left_join(site_list, by = "pad_id") %>%
  mutate(expected_intervals = (days - day(install_date) + 1) * 24 * 4)

# account for missing months
psnb_distinctmonths = syn_psnb %>%
  group_by(pad_id) %>%
  summarize(cm_months = n_distinct(month)) %>%
  ungroup() %>%
  left_join(syn_production, by = "pad_id") %>%
  mutate(
    prod_months = rowSums(!is.na(across(
      starts_with("total_gas")
    )) & across(starts_with("total_gas")) != 0),
    total_months = if_else(cm_months < prod_months, prod_months, cm_months),
    unmonitored_months = if_else(cm_months < prod_months, prod_months - cm_months, 0)
  ) %>%
  select(pad_id,
         cm_months,
         prod_months,
         total_months,
         unmonitored_months)

# naive scaling
psnb_emission_ns = sum(psnb_monthly$total)
mm_avg_emission_ns = psnb_emission_ns / (
  sum(psnb_distinctmonths$total_months) - sum(psnb_distinctmonths$unmonitored_months)
)
mm_u_emission_ns = mm_avg_emission_ns * sum(psnb_distinctmonths$unmonitored_months)
psnb_emission_ns = psnb_emission_ns + mm_u_emission_ns

# observed scaling for all
psnb_avg_emission = sum(psnb_monthly$total) / (sum(psnb_monthly$intervals) - sum(psnb_monthly$missing_intervals)) # get average emissions per interval observed
u_emission_os = psnb_avg_emission * sum(psnb_monthly$missing_intervals) # multiply average emission by intervals not monitored
psnb_t_emission_os = sum(psnb_monthly$total) + u_emission_os # add monitored and unmonitored emissions for total
mm_avg_emission_os = psnb_t_emission_os / (
  sum(psnb_distinctmonths$total_months) - sum(psnb_distinctmonths$unmonitored_months)
) # get average emissions per observed months
mm_u_emission_os = mm_avg_emission_os * sum(psnb_distinctmonths$unmonitored_months) # get emissions missed during unmonitored months
psnb_emission_os = psnb_t_emission_os + mm_u_emission_os # add total emissions of monitored months with total emissions of unmonitored months

# observed scaling by site
psnb_emission_site = psnb_monthly %>%
  mutate(
    psnb_avg_emission_site = total / (intervals - missing_intervals),
    u_emission_site = psnb_avg_emission_site * missing_intervals,
    psnb_emission_site = total + if_else(is.na(u_emission_site), 0, u_emission_site)
  )

psnb_emission_site_allmonths = psnb_monthly %>%
  mutate(
    psnb_avg_emission_site = total / (intervals - missing_intervals),
    u_emission_site = psnb_avg_emission_site * missing_intervals,
    psnb_emission_site = total + if_else(is.na(u_emission_site), 0, u_emission_site)
  ) %>%
  group_by(pad_id) %>%
  summarize(psnb_monitored_emission_site = sum(psnb_emission_site)) %>%
  ungroup() %>%
  left_join(psnb_distinctmonths, by = "pad_id") %>%
  mutate(
    mm_avg_emission_site = psnb_monitored_emission_site / (total_months - unmonitored_months),
    mm_u_emission_site = mm_avg_emission_site * unmonitored_months,
    psnb_emission_site = psnb_monitored_emission_site + mm_u_emission_site
  )

psnb_emission_normal = c()
psnb_emission_filled_leadingT = c()
psnb_emission_filled_leadingF = c()
n = 1
for (id in unique(syn_psnb$pad_id)) {
  print(paste(id, n, sep = " - "))
  subset_data = syn_psnb[syn_psnb$pad_id == id, ]
  subset_data$leadingT = fill_na_with_average(subset_data$emission, leading_w_0 = TRUE)
  subset_data$leadingF = fill_na_with_average(subset_data$emission, leading_w_0 = FALSE)
  psnb_emission_normal[n] = sum(subset_data$emission, na.rm = TRUE)
  psnb_emission_filled_leadingT[n] = sum(subset_data$leadingT, na.rm = TRUE)
  psnb_emission_filled_leadingF[n] = sum(subset_data$leadingF, na.rm = TRUE)
  
  n = n + 1
}
sum(psnb_emission_normal)
psnb_linear_leadingw0_T = sum(psnb_emission_filled_leadingT)
psnb_linear_leadingw0_F = sum(psnb_emission_filled_leadingF)

psnb_linear_avg_mm = psnb_linear_leadingw0_F / (
  sum(psnb_distinctmonths$total_months) - sum(psnb_distinctmonths$unmonitored_months)
)
psnb_linear_mm_u = psnb_linear_avg_mm * sum(psnb_distinctmonths$unmonitored_months)
psnb_linear_emissions = psnb_linear_leadingw0_F + psnb_linear_mm_u

psnb_linear_sites = data.frame(
  pad_id = unique(syn_psnb$pad_id),
  linear_emission = psnb_emission_filled_leadingF
) %>%
  left_join(psnb_distinctmonths, by = "pad_id") %>%
  mutate(
    mm_avg_emission_site = linear_emission / (total_months - unmonitored_months),
    mm_u_emission_site = mm_avg_emission_site * unmonitored_months,
    linear_emission_site = linear_emission + mm_u_emission_site
  )

### temporal results by method
aerial_temporal_emissions = data.frame(
  source = "aerial",
  extrapolation = "temporal",
  method = c("6monthassumption", "prodmonthconsider", "prodscaling"),
  emissions = c(
    sum(aerial_calculations$aerial_emissions_total_6monthassumpt),
    sum(aerial_calculations$aerial_emissions_total_monthconsider),
    sum(total_prod_scaled$aerial_emissions_total_prodscaled_agg)
  )
)

psna_temporal_emissions = data.frame(
  source = "psna",
  extrapolation = "temporal",
  method = c(
    "naive",
    "naive_allmonths",
    "observedall",
    "observedall_allmonths",
    "observedsites",
    "observedsites_allmonths",
    "linear_leadingw0_F",
    "linear_leadingw0_F_allmonths"
  ),
  emissions = c(
    sum(psna_monthly$total),
    psna_emission_ns,
    t_emission_os,
    psna_emission_os,
    sum(psna_emission_site$t_emission_site),
    sum(psna_emission_site_allmonths$psna_emission_site),
    linear_leadingw0_F,
    sum(linear_sites$linear_emission_site)
  )
)

psnb_temporal_emissions = data.frame(
  source = "psnb",
  extrapolation = "temporal",
  method = c(
    "naive",
    "naive_allmonths",
    "observedall",
    "observedall_allmonths",
    "observedsites",
    "observedsites_allmonths",
    "linear_leadingw0_F",
    "linear_leadingw0_F_allmonths"
  ),
  emissions = c(
    sum(psnb_monthly$total),
    psnb_emission_ns,
    psnb_emission_os,
    psnb_emission_os,
    sum(psnb_emission_site$psnb_emission_site),
    sum(psnb_emission_site_allmonths$psnb_emission_site),
    psnb_linear_leadingw0_F,
    sum(psnb_linear_sites$linear_emission_site)
  )
)

temporal_emissions = rbind(aerial_temporal_emissions,
                           psna_temporal_emissions,
                           psnb_temporal_emissions)

##### Temporal end #####

##### Spatial assumption #####

cm_temporal_emissions_observedall_allmonths = psna_emission_os + psnb_emission_os
cm_temporal_emissions_observedsites_allmonths = sum(psna_emission_site_allmonths$psna_emission_site) + sum(psnb_emission_site_allmonths$psnb_emission_site)
cm_temporal_emissions_linear_allmonths = sum(linear_sites$linear_emission_site) + sum(psnb_linear_sites$linear_emission_site)

# multiply by percent of unseen sites
cm_naive_spatial_observedall = cm_temporal_emissions_observedall_allmonths * nrow(site_list) / sum(site_list$source == "aerial")
cm_naive_spatial_observedsites = cm_temporal_emissions_observedsites_allmonths * nrow(site_list) / sum(site_list$source == "aerial")
cm_naive_spatial_linear = cm_temporal_emissions_linear_allmonths * nrow(site_list) / sum(site_list$source == "aerial")

# multiply by production strata
psna_temporal_emissions = psna_emission_site_allmonths %>%
  select(pad_id, psna_emission_site) %>%
  rename(cm_emission = psna_emission_site)
psnb_temporal_emissions = psnb_emission_site_allmonths %>%
  select(pad_id, psnb_emission_site) %>%
  rename(cm_emission = psnb_emission_site)
cm_prodstrata = rbind(psna_temporal_emissions, psnb_temporal_emissions) %>%
  left_join(syn_production, by = "pad_id") %>%
  group_by(prod_strata) %>%
  summarize(avg_yearly_emission = mean(cm_emission)) %>%
  ungroup()
cm_strata_spatial = syn_production %>%
  filter(source == "aerial") %>%
  group_by(prod_strata) %>%
  summarize(site_count = n()) %>%
  ungroup() %>%
  full_join(cm_prodstrata, by = "prod_strata") %>%
  mutate(total_yearly_emission = site_count * avg_yearly_emission) %>%
  summarize(
    strata_emission_observedsites = sum(total_yearly_emission, na.rm = TRUE) + cm_temporal_emissions_observedsites_allmonths,
  )

psna_temporal_emissions_linear = linear_sites %>%
  select(pad_id, linear_emission_site) %>%
  rename(cm_emission = linear_emission_site)
psnb_temporal_emissions_linear = psnb_linear_sites %>%
  select(pad_id, linear_emission_site) %>%
  rename(cm_emission = linear_emission_site)
cm_prodstrata_linear = rbind(psna_temporal_emissions_linear, psnb_temporal_emissions_linear) %>%
  left_join(syn_production, by = "pad_id") %>%
  group_by(prod_strata) %>%
  summarize(avg_yearly_emission = mean(cm_emission)) %>%
  ungroup()
cm_strata_spatial_linear = syn_production %>%
  filter(source == "aerial") %>%
  group_by(prod_strata) %>%
  summarize(site_count = n()) %>%
  ungroup() %>%
  full_join(cm_prodstrata_linear, by = "prod_strata") %>%
  mutate(total_yearly_emission = site_count * avg_yearly_emission) %>%
  summarize(
    strata_emission_linear = sum(total_yearly_emission, na.rm = TRUE) + cm_temporal_emissions_linear_allmonths
  )

# model using age, production, and ghgrp
cm_modeldata = rbind(psna_temporal_emissions, psnb_temporal_emissions) %>%
  left_join(syn_production, by = "pad_id") %>%
  mutate(log_emission = log10(cm_emission))
unmonitored_modeldata = syn_production %>%
  filter(source == "aerial")
emission_model = lm(log_emission ~ total_prod + ghgrp + age, data = cm_modeldata)
summary(emission_model)
u_modeled_emission = sum(10^predict(emission_model, newdata = unmonitored_modeldata))
modeled_emissions = cm_temporal_emissions_observedsites_allmonths + u_modeled_emission

# model using age, production, and ghgrp
cm_modeldata_l = rbind(linear_sites %>% select(pad_id, linear_emission_site) %>% rename(cm_emission = linear_emission_site), 
                       psnb_linear_sites %>% select(pad_id, linear_emission_site) %>% rename(cm_emission = linear_emission_site)) %>%
  left_join(syn_production, by = "pad_id") %>%
  mutate(log_emission = log10(cm_emission))
unmonitored_modeldata = syn_production %>%
  filter(source == "aerial")
emission_model = lm(log_emission ~ total_prod + ghgrp + age, data = cm_modeldata_l)
summary(emission_model)
u_modeled_emission = sum(10^predict(emission_model, newdata = unmonitored_modeldata))
modeled_emissions_l = cm_temporal_emissions_linear_allmonths + u_modeled_emission

### spatial results by source
spatial_emissions = data.frame(
  extrapolation = "spatial",
  method = c(
    "naive_times2_observedall",
    "naive_times2_observedsites",
    "naive_times2_linear",
    "production_strata_observedsites",
    "production_strata_linear",
    "modeling_unseen",
    "modeling_linear"
  ),
  emissions = c(
    cm_naive_spatial_observedall,
    cm_naive_spatial_observedsites,
    cm_naive_spatial_linear,
    cm_strata_spatial$strata_emission_observedsites,
    cm_strata_spatial_linear$strata_emission_linear,
    modeled_emissions,
    modeled_emissions_l
  )
)

##### Spatial end #####