
name <- "rq2_lmic_abmodel_omicron_SD1"
ages_covered <- 9
dose_3_t <- 180
prior_Rt <- 4
central_Rt <- 4
central_VFR <- 3.9

df_summarise <- readRDS(paste0("processed_outputs/df_summarise_", name, ".rds")) %>%
  filter(age_groups_covered == ages_covered,
         t_d3 == dose_3_t) %>%
  mutate(vaccine_doses = factor(vaccine_doses, levels = c("Expand 2 doses", "2 doses + booster", "Pre-vaccine introduction"), ordered = TRUE)) %>%
  mutate(vfr_lab = factor(vfr_lab, levels=c("Optimistic scenario", "Central scenario", "Pessimistic scenario")))

df_summarise_totals <- readRDS(paste0("processed_outputs/df_summarise_totals_", name, ".rds")) %>%
  filter(age_groups_covered == ages_covered)%>%
  mutate(vaccine_doses = factor(vaccine_doses, levels = c("Expand 2 doses", "2 doses + booster"), ordered = TRUE)) %>%
  mutate(vfr_lab = factor(vfr_lab, levels=c("Optimistic scenario", "Central scenario", "Pessimistic scenario")))

#######################################################################
### calculate differences in peaks between delta wave and omicron waves
#######################################################################

df2 <- df_summarise %>%
  filter(vacc_per_week == 0.02, t_d3== 180, vaccine_doses=="Expand 2 doses") %>%
  select(timestep, vfr_lab, date,vaccine_doses,mu_ab_infection, hosp_scal_omicron, inc_med, inc_t, inc_tmin, inc_tmax, 
         hosp_t, hosp_tmin, hosp_tmax, deaths_t, target_pop) %>%
  mutate(variant = ifelse(date <= "2021-11-15", "delta","omicron"))


peak <- df2 %>%
  group_by(variant,vfr_lab) %>%
  summarise(hosp = max(hosp_t), deaths = max(deaths_t)) %>%
  pivot_wider(names_from = variant, values_from = c(hosp,deaths)) %>%
  mutate(percentage_reduction_hosp = (1-hosp_omicron/hosp_delta)*100,
         percentage_reduction_deaths = (1-deaths_omicron/deaths_delta)*100)

peak



#################################################
# plot total doses over time
ggplot(data = df_summarise, aes(x = as.Date(date), y = vaccines_t/target_pop, col = vaccine_doses)) +
  geom_line() +
  facet_wrap(~rollout_rate)

df_summarise %>%
  select(rollout_rate, vaccines_t, target_pop, t_d3, vaccine_doses, date, severity) %>%
  group_by(rollout_rate, target_pop, t_d3, vaccine_doses, severity) %>%
  summarise(vaccines_t = max(vaccines_t)/max(target_pop))

#cumsum(rev(pop$n/1e6))*0.8

#################################################
# blue-green doses barplot
df_doses <- df_summarise %>%
  filter(vaccine_doses != "Pre-vaccine introduction",
         t_d3 == 180,
         max_Rt_omicron == central_Rt,
         vfr == central_VFR,
         severity == "50% (default)") %>%
  rename("Dose 1" = "dose1_t", "Dose 2" = "dose2_t", "Booster" = "dose3_t") %>%
  filter(rollout_rate == "Default") %>%
  pivot_longer(cols = c("Dose 1", "Dose 2", "Booster"), names_to = "dose") %>%
  mutate(dose = factor(dose, levels = c("Dose 1", "Dose 2", "Booster"), ordered = TRUE))

df1_doses_month <- df_doses %>%
  # filter to last date of each month
  mutate(year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         date = lubridate::floor_date(date, "month")) %>%
  group_by(income_group, target_pop, age_groups_covered, vaccine_doses, dose, year, month) %>% 
  mutate(max_day = max(day)) %>%
  ungroup() %>%
  filter(day == max_day)

plot_doses <- ggplot(data = df1_doses_month, aes(x = as.Date(date), y = value/target_pop*100, fill = dose)) +
  geom_bar(stat = "identity") +
  facet_grid( ~ vaccine_doses) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        axis.text.x = element_text(angle = 335, vjust = 0.3, hjust=0.2)) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Time", y = "Vaccinated (%)", fill = "Dose number")
plot_doses

df1 <- df_summarise

deaths_omicron <- ggplot(data = filter(df1, vacc_per_week == 0.02), aes(x = as.Date(date), y = deaths_t/target_pop * 1e6, col = vaccine_doses)) +
  geom_ribbon(aes(ymin = deaths_tmin/target_pop * 1e6, ymax = deaths_tmax/target_pop * 1e6, fill = vaccine_doses), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c(col6, col8, "black")) +
  scale_fill_manual(values = c(col6, col8, "grey")) +
  labs(x = "Time", y = "Daily deaths per million", col = "Dose scenario", fill = "Dose scenario")

deaths_omicron

hosp_omicron <- ggplot(data = filter(df1, vacc_per_week == 0.02), aes(x = as.Date(date), y = hosp_t/target_pop * 1e6, col = vaccine_doses)) +
  geom_ribbon(aes(ymin = hosp_tmin/target_pop * 1e6, ymax = hosp_tmax/target_pop * 1e6, fill = vaccine_doses), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c(col6, col8, "black")) +
  scale_fill_manual(values = c(col6, col8, "grey")) +
  labs(x = "Time", y = "Daily hospital admission per million", col = "Dose scenario", fill = "Dose scenario")

hosp_omicron


df1_nat <- df1 %>%
  filter(max_Rt == prior_Rt,
         t_d3 == 180,
         vacc_per_week == 0.02,
         timestep < max(df_summarise$timestep)) %>%
  mutate(titre = nat_med, nat_label = "All") %>%
  select(date, titre, vfr_lab, nat_label, severity, vaccine_doses)

df1_vax_ab <- df1 %>%
  filter(vacc_per_week == 0.02,
         max_Rt == prior_Rt,
         t_d3 == 180,
         timestep < max(df_summarise$timestep)) %>%
  mutate(titre = vax_ab_med, nat_label = "Vaccination") %>%
  select(date, titre, vfr_lab, nat_label, severity, vaccine_doses)

df1_nat_ab <- df1 %>%
  filter(vacc_per_week == 0.02,
         max_Rt == prior_Rt,
         t_d3 == 180,
         timestep < max(df_summarise$timestep)) %>%
  mutate(titre = nat_ab_med, nat_label = "Infection-induced") %>%
  select(date, titre, vfr_lab, nat_label, severity, vaccine_doses)

df1_nat <- rbind(df1_nat, df1_vax_ab, df1_nat_ab)

nat_omicron <- ggplot(data = filter(df1_nat, vaccine_doses != "Expand 2 doses"), aes(x = as.Date(date), y = titre, col = nat_label)) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_viridis_d(option = "A", begin = 0.2, end = 0.8) +
  labs(x = "Time", y = "Mean neutralizing antibody titre", col = "Immunity Type", fill = "Immunity Type")

nat_omicron

#########################
# summary barchart
df_om <- df_summarise_totals %>%
  select(deaths_med, hosp_med, inc_med, vfr_lab, vaccine_doses, rollout_rate, dose_3_timing, max_Rt_omicron, target_pop) #%>%
  #mutate(scenario = paste0("VFR = ", vfr)) %>%
 # filter(max_Rt_omicron == central_Rt)

df_barchart <- df_om %>%
  filter(!(dose_3_timing == "12 months" & rollout_rate == "Slower rollout")) %>%
  mutate(sensitivity_scenario = if_else(dose_3_timing == "6 months (default)" & rollout_rate == "Default", "Default", if_else(dose_3_timing == "6 months (default)" & rollout_rate == "Slower rollout", "Slower rollout", if_else(dose_3_timing == "12 months" & rollout_rate == "Default", "12 months to booster", "NA")))) %>%
  filter(sensitivity_scenario != "NA") %>%
  mutate(sensitivity_scenario = factor(sensitivity_scenario, levels = c("Default", "Slower rollout", "12 months to booster")))

x <- df_barchart %>%
  filter(sensitivity_scenario == "Default") %>%
  select(-sensitivity_scenario, -rollout_rate, -dose_3_timing) 
x
### 9 age groups
(1516  -1379)/1516
(5262-4755)/5262

### 5 age groups
(1900-1810)/1900
(5664-5438)/5664


# barplot summary of deaths
p_deaths_summary <- ggplot(data = df_barchart, aes(x = vfr_lab, y = deaths_med/target_pop * 1e6, fill = vaccine_doses)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
  labs(x = "Scenario", y = "Total deaths per million", col = "Dose scenario", fill = "Dose scenario") +
  facet_wrap(~sensitivity_scenario) +
  scale_fill_manual(values = c(col6, col8)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)

p_deaths_summary

# barplot summary of incidence
p_inc_summary <- ggplot(data = df_barchart, aes(x = vfr_lab, y = inc_med/target_pop * 1e6, fill = vaccine_doses)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
  labs(x = "Scenario", y = "Total incidence per million", col = "Dose scenario", fill = "Dose scenario") +
  facet_wrap(~sensitivity_scenario) +
  scale_fill_manual(values = c(col6, col8)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)

p_inc_summary
ggsave(paste0("plots/fig3_inc_age_groups_covered_", ages_covered, ".png"),p_inc_summary, height = 4, width = 11)


#############################################################
#########################
# summary roll-out rate
df_om_2 <- df_summarise_totals %>%
  select(deaths_med, inc_med, vfr, vaccine_doses, vacc_per_week, t_d3, max_Rt_omicron, target_pop) %>%
  filter(vfr == central_VFR,
         t_d3 == 180,
         ) %>%
  mutate(max_Rt_omicron = paste0("Rt = ", max_Rt_omicron))


# combine plots
library(patchwork)
layout <- "
AA
BB
CC
DD
"
combined <- plot_doses + hosp_omicron + deaths_omicron+ p_deaths_summary  + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect") + 
  plot_layout(ncol = 2, nrow = 3, design = layout, widths = c(2,1))

combined
ggsave(paste0("plots/fig3_age_groups_covered_", ages_covered, ".png"),combined, height = 11, width = 11)

nat_omicron
ggsave("plots/nat_omicron_rq2.png",nat_omicron, height = 3.5, width = 11)

