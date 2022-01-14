prior_Rt <- 3.5
central_Rt <- 3.5
central_vfr <- 3.9
# plot deaths - omicron scenario
name <- "rq1_hic_abmodel_omicron_SD1"
df_summarise_om <- readRDS(paste0("processed_outputs/df_summarise_", name, ".rds"))
df_summarise_totals_om <- readRDS(paste0("processed_outputs/df_summarise_totals_", name, ".rds"))

fit1 <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE"
vacc_params <- readRDS(paste0("data/param_list_",fit1,".rds")) %>%
  mutate(std10 = 0.44) 

#############################################
# numbers for text
x <- df_summarise_totals_om %>% 
  filter(vacc_per_week == 0.05,
         t_d3 == 180,
         max_Rt_omicron == central_Rt,
         vfr == central_vfr,
         strategy_name %in% c("10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+")) %>%
  select(target_pop, strategy_name, severity, vfr_lab, deaths_med, hosp_med, icu_med, inc_med)
x

(x[which(x$strategy_name == "10y+ 2 doses, no booster"),]$deaths_med - x[which(x$strategy_name == "10y+ 2 doses, booster 60y+"),]$deaths_med) / x[which(x$strategy_name == "10y+ 2 doses, no booster"),]$deaths_med

(x[which(x$strategy_name == "10y+ 2 doses, no booster"),]$hosp_med - x[which(x$strategy_name == "10y+ 2 doses, booster 60y+"),]$hosp_med) / x[which(x$strategy_name == "10y+ 2 doses, no booster"),]$hosp_med

################################
m <- unique(df_summarise_om$strategy_name)
m

df_summarise_om <- df_summarise_om %>%
  mutate(strategy_name = factor(strategy_name, levels = m, ordered = TRUE)) %>%
  mutate(vfr_lab = factor(vfr_lab, levels=c("Optimistic scenario", "Central scenario", "Pessimistic scenario")))

df_summarise_totals_om <- df_summarise_totals_om %>%
  mutate(strategy_name = factor(strategy_name, levels = m, ordered = TRUE)) %>%
  mutate(vfr_lab = factor(vfr_lab, levels=c("Optimistic scenario", "Central scenario", "Pessimistic scenario")))

df_summarise_om <- df_summarise_om %>%
  mutate(dose_3_timing = factor(dose_3_timing, levels = c("6 months (default)", "3 months")))

df_summarise_totals_om <- df_summarise_totals_om %>%
  mutate(dose_3_timing = factor(dose_3_timing, levels = c("6 months (default)", "3 months")))

p_deaths_omicron <- ggplot(data = filter(df_summarise_om,
                                 strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"),
                                 t_d3 == 180)
                   , aes(x = as.Date(date), y = deaths_t/target_pop * 1e6, col = strategy_name)) +
  geom_ribbon(aes(ymin =deaths_tmin/target_pop * 1e6, ymax = deaths_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c("grey20", col_set_3)) +
  scale_fill_manual(values = c("grey20", col_set_3)) +
  labs(x = "Time", y = "Daily deaths per million", col = "Dose scenario", fill = "Dose scenario")

p_deaths_omicron

p_ICU_omicron <- ggplot(data = filter(df_summarise_om,
                                         strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"),
                                         t_d3 == 180)
                           , aes(x = as.Date(date), y = ICU_t/target_pop * 1e6, col = strategy_name)) +
  geom_ribbon(aes(ymin =ICU_tmin/target_pop * 1e6, ymax = ICU_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c("grey20", col_set_3)) +
  scale_fill_manual(values = c("grey20", col_set_3)) +
  labs(x = "Time", y = "Daily ICU admissions per million", col = "Dose scenario", fill = "Dose scenario")

p_ICU_omicron

p_hosp_omicron <- ggplot(data = filter(df_summarise_om,
                                      strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"),
                                      t_d3 == 180)
                        , aes(x = as.Date(date), y = hosp_t/target_pop * 1e6, col = strategy_name)) +
  geom_ribbon(aes(ymin =hosp_tmin/target_pop * 1e6, ymax = hosp_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c("grey20", col_set_3)) +
  scale_fill_manual(values = c("grey20", col_set_3)) +
  labs(x = "Time", y = "Daily hospital admissions per million", col = "Dose scenario", fill = "Dose scenario")

p_hosp_omicron


p_inc_omicron <- ggplot(data = filter(df_summarise_om,
                                         strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"),
                                         t_d3 == 180)
                                      , aes(x = as.Date(date), y = (inc_t/target_pop*1e6), col = strategy_name)) +
   geom_ribbon(aes(ymin =inc_tmin/target_pop * 1e6, ymax = inc_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c("grey20", col_set_3)) +
  scale_fill_manual(values = c("grey20", col_set_3)) +
#  scale_y_continuous(trans = "log10") +
  labs(x = "Time", y = "Daily incidence per million", col = "Dose scenario", fill = "Dose scenario")

p_inc_omicron

p_sp_omicron <- ggplot(data = filter(df_summarise_om,
                                      strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"),
                                      t_d3 == 180)
                                       , aes(x = as.Date(date), y = sp_med, col = strategy_name)) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none") +
  scale_color_manual(values = c("grey20", col_set_3)) +
  scale_fill_manual(values = c("grey20", col_set_3)) +
  #scale_y_continuous(trans = "log10") +
  labs(x = "Time", y = "NAT > 50% threshold", col = "Dose scenario", fill = "Dose scenario")
p_sp_omicron

p_nat_omicron <- ggplot(data = filter(df_summarise_om,
                                     strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"),
                                     t_d3 == 180)
                       , aes(x = as.Date(date), y = nat_med, col = strategy_name)) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_manual(values = c("grey20", col_set_3)) +
  scale_fill_manual(values = c("grey20", col_set_3)) +
  #scale_y_continuous(trans = "log10") +
  labs(x = "Time", y = "Mean neutralizing antibody titre", col = "Dose scenario", fill = "Dose scenario")
p_nat_omicron

df1_nat <- df_summarise_om %>%
  filter(strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, booster 10y+"),
         max_Rt == prior_Rt,
         t_d3 == 180) %>%
  mutate(titre = nat_med, nat_label = "All") %>%
  select(date, titre, vfr_lab, nat_label)

df1_vax_ab <- df_summarise_om %>%
  filter(strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, booster 10y+"),
         max_Rt == prior_Rt,
         t_d3 == 180) %>%
  mutate(titre = vax_ab_med, nat_label = "Vaccination") %>%
  select(date, titre, vfr_lab, nat_label)

df1_nat_ab <- df_summarise_om %>%
  filter(strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, booster 10y+"),
         max_Rt == prior_Rt,
         t_d3 == 180) %>%
  mutate(titre = nat_ab_med, nat_label = "Infection-induced") %>%
  select(date, titre, vfr_lab, nat_label)

df1_nat <- rbind(df1_nat, df1_vax_ab, df1_nat_ab)

nat_omicron <- ggplot(data = df1_nat, aes(x = as.Date(date), y = titre, col = nat_label)) +
  geom_line() +
  facet_wrap(~vfr_lab, nrow = 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_viridis_d(option = "A", begin = 0.2, end = 0.8) +
  labs(x = "Time", y = "Mean neutralizing antibody titre", col = "Immunity type", fill = "Immunity Type")

nat_omicron

############################################################################
#changed to include all scenarios for the text 

y1 <- df_summarise_om %>%
  filter(date > as.Date("2021-01-01")) %>%
  filter(vacc_per_week == 0.05,
         t_d3 == 180,
         strategy_name %in% c("Pre-vaccine introduction", "10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"),
         ) %>%
  group_by(target_pop, strategy_name, vfr, max_Rt_omicron, severity) %>%
  summarise(deaths = max(deaths_t),
            hosp = max(hosp_t),
            inf = max(inc_t))

## ranges 

y2 <- y1 %>%
      group_by(strategy_name) %>%
      summarise (max_deaths = max(deaths), min_deaths=min(deaths), max_hosp=max(hosp), min_hosp=min(hosp))
y2


#########################
# summary barchart
df_barchart <- df_summarise_totals_om %>%
  filter(t_d3 == 180) %>%
  select(deaths_med, icu_med, hosp_med, inc_med, vfr, vfr_lab, strategy_name, max_Rt_omicron, target_pop, severity, t_d3) 

df_barchart_90 <- df_summarise_totals_om %>%
  filter(t_d3 == 90) %>%
  filter(vfr == central_vfr) %>%
  select(deaths_med, icu_med, hosp_med, inc_med, vfr, vfr_lab, strategy_name, max_Rt_omicron, target_pop, severity, t_d3) %>%
  mutate(vfr_lab = paste0(vfr_lab, "\n3 month boost")) 

df_barchart <- rbind(df_barchart, df_barchart_90) %>%
  mutate(vfr_lab = case_when(vfr_lab == "Optimistic scenario" ~ "Optimistic\nscenario",
                             vfr_lab == "Central scenario" ~ "Central\nscenario",
                             vfr_lab == "Central scenario\n3 month boost" ~ "Central scenario\n3 month boost",
                             vfr_lab == "Pessimistic scenario" ~ "Pessimistic\nscenario"))  %>%
  mutate(vfr_lab = factor(vfr_lab, levels=c("Optimistic\nscenario", "Central\nscenario", "Central scenario\n3 month boost","Pessimistic\nscenario")))
  

# barplot summary of deaths
p_deaths_summary <- ggplot(data = df_barchart,  aes(x = vfr_lab, y = deaths_med/target_pop * 1e6, fill = strategy_name)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = col_set) +
  labs(x = "Scenario", y = "Total deaths per million", col = "Dose scenario", fill = "Dose scenario") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)
p_deaths_summary

# barplot summary of hospitalisations
p_hosp_summary <- ggplot(data = df_barchart , aes(x = vfr_lab, y = hosp_med/target_pop * 1e6, fill = strategy_name)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = col_set) +
  labs(x = "Scenario", y = "Total hospitalisations per million", col = "Dose scenario", fill = "Dose scenario") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)

p_hosp_summary


# barplot summary of incidence
p_inc_summary_180 <- ggplot(data = df_barchart, aes(x = vfr_lab, y = inc_med/target_pop * 1e6, fill = strategy_name)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = col_set) +
  labs(x = "Scenario", y = "Total incidence per million", col = "Scenario", fill = "Dose scenario") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)

p_inc_summary_180

# plot total doses over time
fig_doses_time <- ggplot(data = filter(df_summarise_om,
                                       strategy_name != "Pre-vaccine introduction",
                                       vfr == central_vfr
                                      ), aes(x = as.Date(date), y  = vaccines_t/target_pop, col = strategy_name, linetype = dose_3_timing)) +
  geom_line(size = 1) +
  lims(x = c(as.Date("2021-01-01"), as.Date("2022-06-30"))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line() ) + #,
      #  legend.text.align = 0) +
  theme(legend.position = "none") +
  labs(x = "Time", y = "Cumulative doses per person", col = "Dose scenario", linetype = "Booster dose timing") +
  scale_color_manual(values = col_set)

fig_doses_time
ggsave("plots/rq1_fig_doses_time.png", fig_doses_time, height = 3.5, width = 7)



############################
library(patchwork)
layout <- "
AB
CC
DD
EE
"
combined <-  fig_doses_time + p_hosp_summary + p_hosp_omicron + p_deaths_omicron + p_sp_omicron + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect") + 
  plot_layout(ncol = 2, nrow = 3, design = layout)

combined
ggsave("plots/Figure 2.png",combined, height = 11, width = 12)

layout2 <- "
AA
BB
"
inc_combined <- p_inc_omicron +p_inc_summary_180 + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect") + 
  plot_layout(ncol = 2, nrow = 2, design = layout2)
inc_combined
ggsave("plots/fig2_incidence.png",inc_combined, height = 6, width = 11)

nat_omicron
ggsave("plots/nat_omicron_rq1.png",nat_omicron, height = 3.5, width = 11)

