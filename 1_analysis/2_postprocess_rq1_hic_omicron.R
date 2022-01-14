name <- "rq1_hic_abmodel_omicron_SD1"

# join the runs and link to parameters
scenarios <- read_csv(paste0("scenarios_", name, ".csv"), show_col_types = FALSE)
df <- list.files(path = paste0("raw_outputs/output_", name, "/"), pattern = ".rds")
df <- map(paste0("raw_outputs/output_", name, "/", df), readRDS)
df <- do.call(rbind,df)

df <- left_join(df, scenarios, by = "scenario") 

vfr_vec <- unique(df$vfr)[order(unique(df$vfr))]
min_vfr <- vfr_vec[1]
central_vfr <- vfr_vec[2]
max_vfr <- vfr_vec[3]

# name the options
df <- df %>%
  mutate(strategy_name = 
           if_else(vaccine_doses == 2 & age_groups_covered == 15, "10y+ 2 doses, no booster",
                   if_else(vaccine_doses == 3 & age_groups_covered == 15 & age_groups_covered_d3 == 5, "10y+ 2 doses, booster 60y+",
                           if_else(vaccine_doses == 3 & age_groups_covered == 15 & age_groups_covered_d3 == 9, "10y+ 2 doses, booster 40y+",
                                   if_else(vaccine_doses == 3 & age_groups_covered == 15 & age_groups_covered_d3 == 13, "10y+ 2 doses, booster 20y+",
                                           if_else(vaccine_doses == 3 & age_groups_covered == 15 & age_groups_covered_d3 == 15, "10y+ 2 doses, booster 10y+", "NA")))))) %>%
  mutate(rollout_rate = if_else(vacc_per_week == 0.05, "Default", if_else(vacc_per_week < 0.05, "Slower rollout", "None"))) %>%
  mutate(dose_3_timing = if_else(t_d3 == 180, "6 months (default)", if_else(t_d3 == 90, "3 months", if_else(t_d3 == 360, "12 months", "NA")))) %>%
  mutate(severity = factor(hosp_scal_omicron, levels = c(0.4, 0.5, 0.6), labels = c("60%", "50% (default)", "40%"))) %>%
  mutate(vfr_lab = case_when(vfr==min_vfr ~ "Optimistic scenario",
                             vfr==central_vfr ~ "Central scenario",
                             vfr==max_vfr ~ "Pessimistic scenario"))

m <- unique(df$strategy_name)
m
df <- df %>%
  mutate(strategy_name = factor(strategy_name, levels = c(m[1], m[4], m[5], m[2], m[3]), ordered = TRUE))

# summarise totals over repetitions
df <- df %>%
  group_by(income_group, target_pop, max_coverage, age_groups_covered, age_groups_covered_d3, vaccine_doses, vacc_start, strategy_name, ab_model_infection, vacc_per_week, t_d3, dose_3_timing, rollout_rate, max_Rt, max_Rt_omicron, vfr, vfr_time1, vfr_time2, mu_ab_infection, hosp_scal_omicron, ICU_scal_omicron, severity, vfr_lab) %>%
  mutate(deaths_med = median(deaths),
         deaths_lower = quantile(deaths, 0.025),
         deaths_upper = quantile(deaths, 0.975),
         hosp_med = median(cum_hosp),
         hosp_lower = quantile(cum_hosp, 0.025),
         hosp_upper = quantile(cum_hosp, 0.975),
         icu_med = median(cum_ICU),
         icu_lower = quantile(cum_ICU, 0.025),
         icu_upper = quantile(cum_ICU, 0.975),
         inc_med = median(inc),
         inc_lower = quantile(inc, 0.025),
         inc_upper = quantile(inc, 0.975),
         prop_R_med = median(prop_R),
         total_doses_med = median(total_doses)) %>%
  ungroup() 

df_summarise_totals <- df %>%
  select(-deaths, -cum_hosp, -cum_hosp_all, -cum_ICU, -cum_ICU_all, -inc, -total_doses, -prop_R, -cols, -repetition, - scenario) %>%
  unique()

m3 <- unique(df_summarise_totals$strategy_name)
m3
df_summarise_totals <- df_summarise_totals %>%
  mutate(strategy_name = factor(strategy_name, levels = c(m3[1], m3[5], m3[2], m3[3], m3[4]), ordered = TRUE))

# summarise temporal dynamics over repetitions
df_summarise <- df %>%
  unnest(cols) %>%
  select(-c(deaths, cum_hosp, cum_hosp_all, cum_ICU, cum_ICU_all, prop_R, inc)) %>%
  group_by(timestep, income_group, target_pop, max_coverage, age_groups_covered, age_groups_covered_d3, vaccine_doses, vacc_start, deaths_med, deaths_lower, deaths_upper, inc_med, inc_lower, inc_upper, prop_R_med, total_doses_med, strategy_name, vacc_per_week, t_d3, dose_3_timing, rollout_rate, max_Rt, max_Rt_omicron, vfr, vfr_time1, vfr_time2, mu_ab_infection, hosp_scal_omicron, ICU_scal_omicron, severity, vfr_lab) %>%
  summarise(deaths_t = median(D_count),
            deaths_tmin = quantile(D_count, 0.025),
            deaths_tmax = quantile(D_count, 0.975),
            hosp_t = median(hosp),
            hosp_tmin = quantile(hosp, 0.025),
            hosp_tmax = quantile(hosp, 0.975),
            ICU_t = median(ICU),
            ICU_tmin = quantile(ICU, 0.025),
            ICU_tmax = quantile(ICU, 0.975),
            inc_t = median(incidence),
            inc_tmin = quantile(incidence, 0.025),
            inc_tmax = quantile(incidence, 0.975),
            vax_ab_med = median(vax_ab_mean),
            vax_ab_lower = median(vax_ab_lower),
            vax_ab_upper = median(vax_ab_upper),
            nat_ab_med = median(nat_ab_mean),
            nat_ab_lower = median(nat_ab_lower),
            nat_ab_upper = median(nat_ab_upper),
            nat_med = median(nat_mean),
            nat_lower = median(nat_lower),
            nat_upper = median(nat_upper),
            sp_med = median(sp_mean),
            sp_lower = median(sp_lower),
            sp_upper = median(sp_upper),
            vaccines_t = median(X1_count + X2_count * 2 + X3_count * 3),
            dose1_t = median(X1_count),
            dose2_t = median(X2_count),
            dose3_t = median(X3_count),
            prop_R = median(round(sp_mean * 100,2)),
            vax_ab = median(vax_ab_mean),
            nat_ab = median(nat_ab_mean),
            nat = median(nat_mean),
            Rt = median(Rt),
            .groups = 'drop') %>%
  unique() %>%
  mutate(date = timestep + as.Date("2020-02-01"))

# tidy up period before vacc introduction
df0_pre_vacc <- df_summarise %>%
  filter(date < as.Date("2021-01-01")) %>%
  filter(strategy_name == "10y+ 2 doses, no booster") %>%
  mutate(strategy_name = "Pre-vaccine introduction")

df0_post_vacc <- df_summarise %>%
  filter(date >= as.Date("2021-01-01"))

df_summarise <- rbind(df0_pre_vacc, df0_post_vacc)

saveRDS(df_summarise, paste0("processed_outputs/df_summarise_", name, ".rds"))
saveRDS(df_summarise_totals, paste0("processed_outputs/df_summarise_totals_", name, ".rds"))

