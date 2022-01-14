name <- "rq2_lmic_abmodel_omicron_SD1"

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

df <- df %>%
  mutate(vaccine_doses = factor(vaccine_doses, levels = c(2,3), labels = c("Expand 2 doses", "2 doses + booster"))) %>%
  mutate(age_groups_covered_label = factor(age_groups_covered, levels = c(5,9), labels = c("60+ years vaccinated before booster introduced", "40+ years vaccinated before booster introduced"))) %>%
  mutate(rollout_rate = if_else(vacc_per_week == 0.02, "Default", if_else(vacc_per_week == 0.015, "Slower rollout", "None"))) %>%
  mutate(dose_3_timing = if_else(t_d3 == 180, "6 months (default)", if_else(t_d3 == 90, "3 months", if_else(t_d3 == 360, "12 months", "NA"))))%>%
  mutate(severity = factor(hosp_scal_omicron, levels = c(0.4, 0.5, 0.6), labels = c("60%", "50% (default)", "40%"))) %>%
  mutate(vfr_lab = case_when(vfr==min_vfr ~ "Optimistic scenario",
                             vfr==central_vfr ~ "Central scenario",
                             vfr==max_vfr ~ "Pessimistic scenario"))

# summarise totals over repetitions
df <- df %>%
  group_by(income_group, target_pop, max_coverage, age_groups_covered, age_groups_covered_label, vaccine_doses, vacc_start, rollout_rate, vacc_per_week, max_Rt, max_Rt_omicron, t_d3, dose_3_timing, vfr, vfr_time1, vfr_time2, mu_ab_infection, hosp_scal_omicron, ICU_scal_omicron, severity, vfr_lab) %>%
  mutate(deaths_med = median(deaths),
         deaths_lower = quantile(deaths, 0.025),
         deaths_upper = quantile(deaths, 0.975),
         hosp_med = median(cum_hosp),
         hosp_lower = quantile(cum_hosp, 0.025),
         hosp_upper = quantile(cum_hosp, 0.975),
         hosp_all_med = median(cum_hosp_all),
         hosp_all_lower = quantile(cum_hosp_all, 0.025),
         hosp_all_upper = quantile(cum_hosp_all, 0.975),
         icu_med = median(cum_ICU),
         icu_lower = quantile(cum_ICU, 0.025),
         icu_upper = quantile(cum_ICU, 0.975),
         icu_all_med = median(cum_ICU_all),
         icu_all_lower = quantile(cum_ICU_all, 0.025),
         icu_all_upper = quantile(cum_ICU_all, 0.975),
         inc_med = median(inc),
         inc_lower = quantile(inc, 0.025),
         inc_upper = quantile(inc, 0.975),
         prop_R_med = median(prop_R),
         total_doses_med = median(total_doses)) %>%
  ungroup()  

df_summarise_totals <- df %>%
  select(-deaths, -cum_hosp, -cum_hosp_all, -cum_ICU, -cum_ICU_all, -inc, -total_doses, -prop_R, -cols, -repetition, - scenario) %>%
  unique()

# summarise temporal dynamics over repetitions
df_summarise <- df %>%
  unnest(cols) %>%
  select(-c(deaths, prop_R, inc)) %>%
  group_by(timestep, income_group, target_pop, max_coverage, age_groups_covered,age_groups_covered_label, vaccine_doses, vacc_start, max_Rt, vfr, max_Rt_omicron, rollout_rate, t_d3, vacc_per_week, deaths_med, deaths_lower, deaths_upper, inc_med, inc_lower, inc_upper, prop_R_med, total_doses_med, vfr_time1, vfr_time2, mu_ab_infection, hosp_scal_omicron, ICU_scal_omicron, severity, vfr_lab) %>%
  summarise(deaths_t = median(D_count),
            deaths_tmin = quantile(D_count, 0.025),
            deaths_tmax = quantile(D_count, 0.975),
            hosp_t = median(hosp),
            hosp_tmin = quantile(hosp, 0.025),
            hosp_tmax = quantile(hosp, 0.975),
            hosp_all_t = median(hosp_all),
            hosp_all_tmin = quantile(hosp_all, 0.025),
            hosp_all_tmax = quantile(hosp_all, 0.975),
            ICU_t = median(ICU),
            ICU_tmin = quantile(ICU, 0.025),
            ICU_tmax = quantile(ICU, 0.975),
            ICU_all_t = median(ICU_all),
            ICU_all_tmin = quantile(ICU_all, 0.025),
            ICU_all_tmax = quantile(ICU_all, 0.975),
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
            prop_R = round(median(R_count)/target_pop * 100,2),
            Rt = median(Rt),
            .groups = 'drop') %>%
  unique() %>%
  mutate(date = timestep + as.Date("2020-02-01"))

# tidy up period before vacc introduction
df0_pre_vacc <- df_summarise %>%
  filter(vacc_per_week == 0.02) %>%
  filter(date < as.Date("2021-04-01")) %>%
  filter(vaccine_doses == "Expand 2 doses") %>%
  mutate(vaccine_doses = "Pre-vaccine introduction")

df0_post_vacc <- df_summarise %>%
  filter(date >= as.Date("2021-04-01"))

df_summarise <- rbind(df0_pre_vacc, df0_post_vacc)

saveRDS(df, paste0("processed_outputs/df_", name, ".rds"))
saveRDS(df_summarise, paste0("processed_outputs/df_summarise_", name, ".rds"))
saveRDS(df_summarise_totals, paste0("processed_outputs/df_summarise_totals_", name, ".rds"))
