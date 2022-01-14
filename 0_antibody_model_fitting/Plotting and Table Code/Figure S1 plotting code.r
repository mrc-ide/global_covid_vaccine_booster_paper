library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(emdbook)
library(drjacoby)

setwd(dirname(getActiveDocumentContext()$path)) 
data <- "Imperial"
#data <- "PHE"

##################################
##### LOAD THE PARAMETERS 
##################################

if(data=="Imperial") {
  load("../Model Fits/imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_mcmc_chain.Rdata")  
  
}

if(data=="PHE") {
  load("../Model Fits/phe_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_mcmc_chain.Rdata")  
}

draws <- sample_chains(mcmc, 5000)

draws_transform <- draws %>%
  select(-sample, -AZ_ns_off ) %>%
  mutate(d2_AZ_log10 = log10(32/59) - fold_red_AZ,
         d2_PF_log10 = log10(223/94) - fold_red_PF,
         d2_PF = 10^(d2_PF_log10),
         d2_AZ = 10^(d2_AZ_log10),
         d1_AZ = 10^(d2_AZ_log10 + d1_AZ),
         d1_PF = 10^(d2_PF_log10 + d1_PF),
         d3 = 10^(d2_PF_log10 + bst_PF),
         ab50 = 10^(d2_PF_log10 + ni50),
         ab50_s = 10^(d2_PF_log10 + ns50), 
         ab50_d = 10^(d2_PF_log10 + nd50),
         om_red = 10^(om_red),
         fold_red_PF = 10^(fold_red_PF),
         fold_red_AZ = 10^(fold_red_AZ)) %>%
  select(-ni50, -ns50, -nd50, -d2_AZ_log10, -d2_PF_log10, -bst_AZ, -bst_PF) 

posterior_median_transform <-draws_transform %>%
  summarise( 
    across(where(is.numeric), median)
  )%>%
  mutate(measure = "median")


n<-lseq(0.01,100,length.out=100)


#len<-length(mcmc$output$k)
#sam<-mcmc$parameters$samples

#ni50v        <- 10^mcmc$output$ni50[seq(len-sam,len,100)]
#ns50v        <- 10^mcmc$output$ns50[seq(len-sam,len,100)]
#kv           <- mcmc$output$k[seq(len-sam,len,100)]

ni50v        <- draws_transform$ab50
ns50v        <- draws_transform$ab50_s
nd50v        <- draws_transform$ab50_d
kv           <- draws_transform$k

ru <-NULL

for (i in (1:length(kv))){
  
  ni50         <- ni50v[i] # titre relative to convalescent required to provide 50% protection from infection, on linear scale
  ns50         <- ns50v[i]
  nd50         <- nd50v[i]
  k            <- kv[i] # shape parameter of efficacy curve
  
  ef_infection <- 1 / (1 + exp(-k * (log10(n) - log10(ni50))))
  ef_severe    <- 1 / (1 + exp(-k * (log10(n) - log10(ns50))))
  ef_death    <- 1 / (1 + exp(-k * (log10(n) - log10(nd50))))
  
  r <- data.frame(n=n, run=rep(i,length(n)), ef_infection=ef_infection, ef_severe=ef_severe, ef_death=ef_death)  
  ru<-rbind(ru,r)

}

ru <- ru %>%
  pivot_longer(cols = c("ef_infection", "ef_severe", "ef_death"), names_to = "type") %>%
  mutate(type = factor(type, levels = c("ef_infection", "ef_severe", "ef_death")))

ru_summary <- ru %>%
  group_by(type, n) %>%
  summarise(median = median(value),
            upper = quantile(value, 0.975),
            lower = quantile(value, 0.025),)


ni50      <- posterior_median_transform$ab50
ns50      <- posterior_median_transform$ab50_s
nd50      <- posterior_median_transform$ab50_d
k         <- posterior_median_transform$k


ef_infection <- 1 / (1 + exp(-k * (log10(n) - log10(ni50))))
ef_severe <- 1 / (1 + exp(-k * (log10(n) - log10(ns50))))
r1 <- data.frame(n,ef_infection,ef_severe, ef_death)

g1 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_infection"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkblue") +
  geom_line(data = r1, aes(x = n, y = 100*ef_infection), size = 1,col = "darkblue") +  
  labs(x = "NAT (fold of convalescent)", y = "vaccine efficacy mild (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))

g2 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_severe"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkgray") +
  geom_line(data = r1, aes(x = n, y = 100*ef_severe), size = 1,col = "black") +  
  labs(x = "NAT  (fold of convalescent)", y = "vaccine efficacy severe (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))

g3 <- ggplot(data = r1) +
  geom_ribbon(data = filter(ru_summary, type == "ef_death"), aes(x = n, ymin = 100*lower, ymax = 100*upper), alpha = 0.2, fill = "darkred") +
  geom_line(data = r1, aes(x = n, y = 100*ef_death), size = 1,col = "darkred") +  
  labs(x = "NAT  (fold of convalescent)", y = "vaccine efficacy death (%)") +
  scale_x_log10(limits = c(0.01,100),breaks = c(0.01,0.1,1,10,100)) +
  scale_y_continuous(limits = c(0,100)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1))

combined <- g1 + g2  + g3 + plot_annotation(tag_levels = "A") 
combined
if(data=="Imperial") {
  ggsave("../Figures/FigureS1.png",combined, height = 5, width = 15)  
}
if(data=="PHE") {
  ggsave("../Figures/FigureS1_phe.png",combined, height = 5, width = 15)  
}
