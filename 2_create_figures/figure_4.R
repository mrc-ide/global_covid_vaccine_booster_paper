central_Rt <- 3.0
central_VFR <- 3.9

name <- "rq3_hic_abmodel_omicron_SD1"

df_summarise <- readRDS(paste0("processed_outputs/df_summarise_", name, ".rds"))
df_summarise_totals <- readRDS(paste0("processed_outputs/df_summarise_totals_", name, ".rds"))


m <- unique(df_summarise$strategy_name)
m <- c(m[1], m[3], m[2], m[4])
df_summarise <- df_summarise %>%
  mutate(strategy_name = factor(strategy_name, levels = m, ordered = TRUE))
df_summarise_totals <- df_summarise_totals %>%
  mutate(strategy_name = factor(strategy_name, levels = m, ordered = TRUE))

df1_omicron <- filter(df_summarise, vacc_per_week == 0.05) %>%
  filter(strategy_name %in% c("10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+")) %>%
 # filter(vfr == central_VFR) %>%
  mutate(Rt_lift_t = factor(Rt_lift_t, levels = c("Sept '21 lift", "Nov '21 lift", "April '22 lift", "Slow April '22 lift"))) %>%
  filter(Rt_lift_t == "Sept '21 lift" | (Rt_lift_t == "Nov '21 lift" & (strategy_name %in% c("10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"))) | (Rt_lift_t == "April '22 lift" & strategy_name == "10y+ 2 doses, booster 10y+") | (Rt_lift_t == "Slow April '22 lift" & strategy_name == "10y+ 2 doses, booster 10y+")) %>%
  mutate(vfr_lab = factor(vfr_lab, levels=c("Optimistic scenario", "Central scenario", "Pessimistic scenario")))

df2_omicron <- filter(df_summarise_totals, vacc_per_week == 0.05) %>%
  filter(strategy_name %in% c("10y+ 2 doses, no booster", "10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+")) %>%
  filter(vfr == central_VFR) %>%
  mutate(Rt_lift_t = factor(Rt_lift_t, levels = c("Sept '21 lift", "Nov '21 lift", "April '22 lift", "Slow April '22 lift"))) %>%
  filter(Rt_lift_t == "Sept '21 lift" | (Rt_lift_t == "Nov '21 lift" & (strategy_name %in% c("10y+ 2 doses, booster 60y+", "10y+ 2 doses, booster 10y+"))) | (Rt_lift_t == "April '22 lift" & strategy_name == "10y+ 2 doses, booster 10y+") | (Rt_lift_t == "Slow April '22 lift" & strategy_name == "10y+ 2 doses, booster 10y+"))

# plot outputs: deaths
deaths_omicron <- ggplot(data = filter(df1_omicron, t_d3 == 180, max_Rt_omicron == central_Rt), aes(x = as.Date(date), y = deaths_t/target_pop * 1e6, col = strategy_name)) +
  geom_ribbon(aes(ymin = deaths_tmin/target_pop * 1e6, ymax = deaths_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~ Rt_lift_t, nrow = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none") +
  scale_color_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Time", y = "Daily deaths per million", col = "Dose scenario", fill = "Dose scenario")

deaths_omicron

hosp_omicron <- ggplot(data = filter(df1_omicron, t_d3 == 180, max_Rt_omicron == central_Rt), aes(x = as.Date(date), y = hosp_t/target_pop * 1e6, col = strategy_name)) +
  geom_ribbon(aes(ymin = hosp_tmin/target_pop * 1e6, ymax = hosp_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~ Rt_lift_t, nrow = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none") +
  scale_color_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Time", y = "Daily hospital admissions per million", col = "Dose scenario", fill = "Dose scenario")
hosp_omicron

# plot of all scenarios for the Supplementary

hosp_omicron_all <- ggplot(data = filter(df1_omicron, t_d3 == 180), aes(x = as.Date(date), y = hosp_t/target_pop * 1e6, col = strategy_name)) +
  geom_ribbon(aes(ymin = hosp_tmin/target_pop * 1e6, ymax = hosp_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_grid( Rt_lift_t ~ vfr_lab) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "bottom") +
  scale_color_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Time", y = "Daily hospital admissions per million", col = "Dose scenario", fill = "Dose scenario")
hosp_omicron_all

ggsave("plots/rq3_hic_allhosp.png", hosp_omicron_all, height = 10, width = 10)


sp_omicron <- ggplot(data = filter(df1_omicron, t_d3 == 180, max_Rt_omicron == central_Rt), aes(x = as.Date(date), y = sp_med, col = strategy_name)) +
  geom_line() +
  facet_wrap(~ Rt_lift_t, nrow = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none") +
  scale_color_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Time", y = "NAT > 50% threshold", col = "Dose scenario", fill = "Dose scenario")
sp_omicron

df3_omicron <- df2_omicron %>%
  mutate(Rt_lift_t = str_sub(Rt_lift_t,1,nchar(as.character(Rt_lift_t))-5)) %>%
  mutate(Rt_lift_t = case_when(Rt_lift_t == "April '22" ~ "Apr '22",
                               Rt_lift_t== "Sept '21" ~ "Sep '21",
                               Rt_lift_t == "Slow April '22" ~ "Apr '22 (S)",
                               Rt_lift_t == "Nov '21" ~ "Nov '21")) %>%
  mutate(Rt_lift_t =factor(Rt_lift_t, c("Sep '21", "Nov '21", "Apr '22", "Apr '22 (S)")))

deaths_summary_omicron <- ggplot(data = filter(df3_omicron, t_d3 == 180, max_Rt_omicron == central_Rt), aes(x = Rt_lift_t, y = deaths_med / target_pop * 1e6, fill = strategy_name)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Rt lift scenario", y = "Total deaths per million", fill = "Dose scenario")

deaths_summary_omicron

hosp_summary_omicron <- ggplot(data = filter(df3_omicron, t_d3 == 180, max_Rt_omicron == central_Rt), aes(x = Rt_lift_t, y = hosp_med / target_pop * 1e6, fill = strategy_name)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Rt lift scenario", y = "Total hospitalisations per million", fill = "Dose scenario")

hosp_summary_omicron

###########################

# plot outputs: infections
infections_omicron <- ggplot(data = filter(df1_omicron, t_d3 == 180, max_Rt_omicron == central_Rt), aes(x = as.Date(date), y = inc_t/target_pop * 1e6, col = strategy_name)) +
  geom_ribbon(aes(ymin = inc_tmin/target_pop * 1e6, ymax = inc_tmax/target_pop * 1e6, fill = strategy_name), alpha = 0.5, col = NA) +
  geom_line() +
  facet_wrap(~ Rt_lift_t, nrow = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none") +
  scale_color_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Time", y = "Daily incidence per million", col = "Dose scenario", fill = "Dose scenario")

infections_omicron

infections_summary_omicron <- ggplot(data = filter(df3_omicron, t_d3 == 180, max_Rt_omicron == central_Rt), aes(x = Rt_lift_t, y = inc_med / target_pop * 1e6, fill = strategy_name)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_fill_viridis_d(option = "E", begin = 0.2, end = 0.9) +
  labs(x = "Rt lift scenario", y = "Total incidence per million", fill = "Dose scenario")

infections_summary_omicron


df1_nat <- df1_omicron %>%
  mutate(titre = nat_med, nat_label = "All") %>%
  select(date, titre, Rt_lift_t, nat_label, strategy_name)

df1_vax_ab <- df1_omicron %>%
  mutate(titre = vax_ab_med, nat_label = "Vaccination") %>%
  select(date, titre, Rt_lift_t, nat_label, strategy_name)

df1_nat_ab <- df1_omicron %>%
  mutate(titre = nat_ab_med, nat_label = "Infection-induced") %>%
  select(date, titre, Rt_lift_t, nat_label, strategy_name)

df1_nat <- rbind(df1_nat, df1_vax_ab, df1_nat_ab) %>%
  filter(strategy_name == "10y+ 2 doses, booster 10y+")



nat_omicron <- ggplot(data = df1_nat, aes(x = as.Date(date), y = titre, col = nat_label)) +
  geom_line() +
  facet_wrap(~ Rt_lift_t, nrow = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) +
  scale_color_viridis_d(option = "A", begin = 0.2, end = 0.8) +
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8) +
  labs(x = "Time", y = "Mean neutralising antibody titre", col = "Dose scenario", fill = "Dose scenario")
nat_omicron

###########################
library(patchwork)
layout <- "
ABC
ABD
ABE
"
combined <- hosp_omicron + sp_omicron + infections_summary_omicron + hosp_summary_omicron + deaths_summary_omicron + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", design = layout, widths = c(1.2,1.2,1)) 
combined

ggsave("plots/fig4.png", combined, height = 8, width = 12)

ggsave("plots/nat_omicron_rq3.png", nat_omicron, height = 8, width = 6)

