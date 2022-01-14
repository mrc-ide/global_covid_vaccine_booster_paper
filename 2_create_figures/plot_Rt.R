central_VFR <- 3.9

name <- "rq1_hic_abmodel_omicron_SD1"
df1o <- readRDS(paste0("processed_outputs/df_summarise_", name, ".rds")) %>%
  filter(vaccine_doses == 2,
         t_d3 == 180,
         vfr == central_VFR) %>%
  select(date, Rt, max_Rt_omicron)

name <- "rq2_lmic_abmodel_omicron_SD1"
df2o <- readRDS(paste0("processed_outputs/df_summarise_", name, ".rds")) %>%
  filter(vaccine_doses != "2 doses + booster",
         t_d3 == 180,
         vfr == central_VFR) %>%
  select(date, Rt, max_Rt_omicron) 

name <- "rq3_hic_abmodel_omicron_SD1"
df3o <- readRDS(paste0("processed_outputs/df_summarise_", name, ".rds")) %>%
  filter(vaccine_doses == 3,
         t_d3 == 180,
         strategy_name == "10y+ 2 doses, booster 10y+") %>%
         select(date, Rt, max_Rt_omicron, Rt_lift_t) %>%
         mutate(Rt_lift_t = factor(Rt_lift_t, levels = c("Sept '21 lift", "Nov '21 lift", "April '22 lift", "Slow April '22 lift"), ordered = TRUE))

p_Rt1 <- ggplot(data = df1o, aes(x = as.Date(date), y = Rt, linetype = factor(max_Rt_omicron))) +
  geom_line(size = 0.8) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none",
        plot.title = element_text(size=11)) +
  labs(x = "Time", y = "Rt", col = "", title = "Substantial prior transmission, HIC", linetype = "")

p_Rt1

p_Rt2 <- ggplot(data = df2o, aes(x = as.Date(date), y = Rt, linetype = factor(max_Rt_omicron))) +
  geom_line(size = 0.8) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "none",
        plot.title = element_text(size=11)) +
  labs(x = "Time", y = "Rt", col = "", title = "Substantial prior transmission, LMIC", linetype = "")

p_Rt2

p_Rt3 <- ggplot(data = df3o, aes(x = as.Date(date), y = Rt, col = Rt_lift_t)) +
  geom_line(size = 0.8) +
  facet_wrap(~Rt_lift_t, nrow = 4) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0,
        legend.position = "right",
        plot.title = element_text(size=11)) +
  labs(x = "Time", y = "Rt", col = "", title = "Minimal prior transmission, HIC", linetype = "")

p_Rt3

library(patchwork)
layout <- "
AC
BC
"
combined <- p_Rt1 + p_Rt2 + p_Rt3 +
  plot_annotation(tag_levels = "A") + 
  #plot_layout(guides = "collect") + 
  plot_layout(ncol = 2, design = layout)

combined
ggsave("plots/p_Rt.png", combined, height = 5, width = 8)
