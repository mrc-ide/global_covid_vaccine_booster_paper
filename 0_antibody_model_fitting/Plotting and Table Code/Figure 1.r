library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(drjacoby)

setwd(dirname(getActiveDocumentContext()$path)) 

source("R/create_profile.R")

data <- "Imperial"
#data <- "PHE"

##################################
##### LOAD THE PARAMETERS AND DATA
##################################

if(data=="Imperial") {
  load("../Model Fits/imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_mcmc_chain.Rdata")  
  all_data <- read.csv("../Model Fits/Imperial_VE.csv")
}

if(data=="PHE") {
  load("../Model Fits/phe_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_mcmc_chain.Rdata")  
  all_data <- read.csv("../Model Fits/PHE_Data_Long_Omicron_add_2.csv")
}

chain <- sample_chains(mcmc, 10000)

posterior_median <- chain %>%
  summarise( 
  across(where(is.numeric), median)
  )

name <- c("AZ","PF")
  
log10_d2_AZ <- log10(32/59)- posterior_median$fold_red_AZ 
log10_d2_PF <- log10(223/94) - posterior_median$fold_red_PF

d1_AZ       <- 10^( log10_d2_AZ + posterior_median$d1_AZ)
d1_PF       <- 10^( log10_d2_PF + posterior_median$d1_PF)
fold_red_AZ <- 10^( posterior_median$fold_red_AZ) 
fold_red_PF <- 10^( posterior_median$fold_red_PF)

d3_PF       <- 10^(log10_d2_PF + posterior_median$bst_PF)

ab_50_PF       <- 10^( log10_d2_PF + posterior_median$ni50) 
ab_50_severe_PF <- 10^( log10_d2_PF + posterior_median$ns50)
ab_50_death_PF  <- 10^( log10_d2_PF + posterior_median$nd50)

ab_50_AZ       <- 10^( log10_d2_PF + posterior_median$ni50 ) 
ab_50_severe_AZ <- 10^( log10_d2_PF + posterior_median$ns50 ) 
ab_50_death_AZ  <- 10^( log10_d2_PF + posterior_median$nd50 ) 

ab_50 <-c(ab_50_AZ, ab_50_PF)
ab_50_severe <-c(ab_50_severe_AZ, ab_50_severe_PF)
ab_50_death <-c(ab_50_death_AZ, ab_50_death_PF)

k           <- posterior_median$k
hl_s        <- posterior_median$hl_s
hl_l        <- posterior_median$hl_l
period_s    <- posterior_median$period_s
#t_period_l  <- posterior_median$period_l

# fixed parameter
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 scale

om_red <- 10^(posterior_median$om_red)
vfr <- c(1,om_red)

mu_ab_d1 <- c(d1_AZ,d1_PF)
mu_ab_d2 <- c(32/59, 223/94)/c(fold_red_AZ,fold_red_PF)
mu_ab_d3 <- c(d3_PF,d3_PF)   
dose_3_fold_increase <- mu_ab_d3/mu_ab_d2

# transforms
dr_s <- -log(2)/hl_s  # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

# Timing of doses
max_t <- 365*2 # number of days to model
t <- 0:max_t
t_d2 <- 84 # timing of second dose relative to first
t_d3 <- 180 # timing of third dose relative to second dose

param_list <- data.frame(name,mu_ab_d1,mu_ab_d2, mu_ab_d3, t_d2, t_d3, dr_s, dr_l,period_s,ab_50,ab_50_severe) 
  
## manipulate data for plotting
## remove data points not used in fitting

all_data <- all_data %>%
            mutate(t_orig=floor((t_max+t_min)/2)) %>%
            filter(t_orig>=14) %>%
            mutate(t_fit= case_when( dose==1 ~ t_orig-21,
                                     dose>1 ~ t_orig-14)) %>%
            filter((dose==1 & t_fit<t_d2)|(dose==2 & t_fit<t_d3)|dose==3) %>%
            mutate(t= case_when( dose==1 ~ t_fit,
                                 dose==2 ~ t_d2+t_fit,
                                 dose==3 ~ t_d2+t_d3+t_fit )) %>%
            mutate(type=case_when( endpoint==1 ~ "Efficacy",
                                   endpoint==2 ~ "Efficacy_Severe",
                                   endpoint==3 ~ "Death")) %>%
            filter(endpoint==1 | endpoint == 2) 

# initialise other parameters

nt <- NULL
r1_summary <- NULL
summary_stats <- NULL
plots <- NULL
plotlist <- list()

##########################
###### OUTPUT OPTIONS
##########################

output_plots <- 1    # boolean for generating plots
num_ind <- 100     # number of individuals to simulate - keep to 100 for plots, 1000 or more for summary statistics

#####################################################################################
### loop through the vaccines to calculate the profiles of NAT and efficacy over time
#####################################################################################

m <- 2 ### set to 1 for delta, 2 for omicron


for (j in 1:2){

  r2 <- NULL
  r2_summary <- NULL
  sub2 <- NULL
  
          mu_ab_d1 <- param_list$mu_ab_d1[j]/vfr[m]
          mu_ab_d2 <- param_list$mu_ab_d2[j]/vfr[m]
          mu_ab_d3 <- param_list$mu_ab_d3[j]/vfr[m]
          dr_s <- param_list$dr_s[j]
          dr_l <- param_list$dr_l[j]
          period_s <- param_list$period_s[j]
          t_d2 <- param_list$t_d2[j]
          t_d3 <- param_list$t_d3[j]
          t_d4 <- param_list$t_d4[j]
          ab_50 <- param_list$ab_50[j]
          ab_50_severe <- param_list$ab_50_severe[j]

          ### generate array of simulations and process into plots
          r1 <- NULL
          for (i in 1:num_ind){
            out <- draw(mu_ab_d1 = mu_ab_d1, mu_ab_d2 = mu_ab_d2, mu_ab_d3 = mu_ab_d3, std10 = std10, ab_50 = ab_50, ab_50_severe = ab_50_severe, 
                dr_s = dr_s, dr_l = dr_l, period_s = period_s, t = t, t_d2 = t_d2, t_d3, k = k)
            sub <- data.frame(t = t, run = rep(i, length(t)), Titre = out$titre, Efficacy = out$VE, Efficacy_Severe = out$VE_severe)
            r1 <- rbind(r1, sub)
          }
          out0 <- draw(mu_ab_d1 = mu_ab_d1, mu_ab_d2 = mu_ab_d2, mu_ab_d3 = mu_ab_d3, std10 = 0, ab_50 = ab_50, ab_50_severe = ab_50_severe, 
              dr_s = dr_s, dr_l = dr_l, period_s = period_s, t = t, t_d2 = t_d2, t_d3, k = k)
          sub0<- data.frame(t = t, run = 0, variant=m, Titre = out0$titre, Efficacy = out0$VE, Efficacy_Severe = out0$VE_severe)

          ## cross-check code
          #  look <- sub0 %>%
          #  filter(t==t_d2+t_d3)
        
          sub0 <- sub0 %>%
              mutate(Efficacy = Efficacy * 100,
              Efficacy_Severe = Efficacy_Severe * 100) %>%
              pivot_longer(cols = c("Titre", "Efficacy", "Efficacy_Severe"), names_to = "type") %>%
              mutate(type = factor(type, levels = c("Titre", "Efficacy", "Efficacy_Severe"))) 

          look <- sub0 %>% 
             filter(t==20 )
          
          r1 <- r1 %>%
              mutate(variant=m) %>%
              mutate(Efficacy = Efficacy * 100,
              Efficacy_Severe = Efficacy_Severe * 100) %>%
              pivot_longer(cols = c("Titre", "Efficacy", "Efficacy_Severe"), names_to = "type") %>%
              mutate(type = factor(type, levels = c("Titre", "Efficacy", "Efficacy_Severe")))
  
          r1_summary <- r1 %>%
              group_by(type, t) %>%
              summarise(median = median(value),
                  upper = quantile(value, 0.975),
                  lower = quantile(value, 0.025),
                  ) %>%
              mutate(name=name[j]) %>%
              mutate(vaccine = name[j]) %>%
              mutate(variant=m)
          
          sub2 <- rbind(sub2, sub0)
          r2 <- rbind(r2,r1)
          r2_summary <- rbind(r2_summary, r1_summary)
  
  if(j==1){
    uk_data <- all_data %>%
    filter(vaccine=="AZ") 
  }
  if(j==2) {
    uk_data <- all_data %>%
    filter(vaccine=="PF" | vaccine=="Pfizer")
  }
  if(j==3) {
    uk_data <- all_data %>%
    filter(vaccine=="MD")
  }
    
  r2_summary <- r2_summary %>%
    left_join(uk_data,by=c("t","type","variant"))
  
  r2_summary <- r2_summary %>%
    left_join(sub2,by=c("t","type","variant")) %>%
    # dont allow data CIs to be <0
    mutate(L95 = if_else(L95 < 0, 0, L95))

  summary_stats <-rbind(summary_stats,r2_summary)
  
   ###################################
   ##### generate plots
   ###################################
 
  if(output_plots == 1){
      
      g1 <- ggplot(data = filter(r2, type == "Titre", variant==m)) +
      geom_line(aes(x = t, y = value, group = run), col = "grey") +
      geom_line(data = filter(r2_summary, type == "Titre", variant==m), aes(x = t, y = value), size = 1,col = "black") +  
      geom_ribbon(data = filter(r2_summary, type == "Titre", variant==m), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "darkgreen") +
      labs(x = "time (days)", y = "NAT") +
      scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
      scale_y_log10(limits = c(1e-3,1e2)) +
      theme_bw() +
      theme(strip.background = element_rect(fill = NA, color = "white"),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_text(angle=60, hjust = 1))
    
      g2 <- ggplot(data = filter(r2, type == "Efficacy", variant==m), na.rm=FALSE) +
      geom_line(aes(x = t, y = value, group = run), col = "grey") +
      geom_line(data = filter(r2_summary, type == "Efficacy", variant==m), aes(x = t, y = value), size = 1,col = "black") +  
      geom_ribbon(data = filter(r2_summary, type == "Efficacy", variant==m), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "darkblue") +
      geom_point(data = filter(r2_summary, type == "Efficacy", variant==m),aes(x = t, y=VE), color = "blue", size = 1.5) +
      geom_errorbar(data = filter(r2_summary, type == "Efficacy"),aes(x = t, ymin=L95, ymax=U95), color = "blue", width=15) +
      lims(y = c(0, 100)) +
      labs(x = "time (days)", y = "vaccine efficacy mild (%)") +
      scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
      theme_bw() +
      theme(strip.background = element_rect(fill = NA, color = "white"),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_text(angle=60, hjust = 1))
    
      g3 <- ggplot(data = filter(r2, type == "Efficacy_Severe", variant==m)) +
      geom_line(aes(x = t, y = value, group = run), col = "grey") +
      geom_line(data = filter(r2_summary, type == "Efficacy_Severe", variant==m), aes(x = t, y = median), size = 1) +
      geom_ribbon(data = filter(r2_summary, type == "Efficacy_Severe", variant==m), aes(x = t, ymin = lower, ymax = upper), alpha = 0.4, fill = "darkgrey") +
      geom_point(data = filter(r2_summary, type == "Efficacy_Severe", variant==m),aes(x = t, y = VE), color = "black", size = 1.5) +
      geom_errorbar(data = filter(r2_summary, type == "Efficacy_Severe"),aes(x = t, ymin=L95, ymax=U95), color = "black", width=15) +
      lims(y = c(0, 100)) +
      labs(x = "time (days)", y = "vaccine efficacy severe (%)") +
      scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
      theme_bw() +
      theme(strip.background = element_rect(fill = NA, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(),
            axis.text.x=element_text(angle=60, hjust = 1))
     
      plots[[j*3-2]] <- g1
      plots[[j*3-1]] <- g2
      plots[[j*3]] <- g3
      
  }
}

## plots for the 3 vaccines, natural infection, and natural infection followed by Pfizer ##

plot_AZ_profiles_delta <- (plots[[1]] | plots[[2]] | plots[[3]])  
plot_PF_profiles_delta <- (plots[[4]] | plots[[5]] | plots[[6]]) 

plot_AZ_profiles_omicron <- (plots[[1]] | plots[[2]] | plots[[3]])  
plot_PF_profiles_omicron <- (plots[[4]] | plots[[5]] | plots[[6]]) 


text = paste("AZ primary PF boost: Delta")
cap1<- ggplot() + annotate("text", x = 0, y = 0, size=4, label = text) + theme_void()+ plot_layout(tag_level = 'new')

text = paste("AZ primary PF boost: Omicron")
cap2<- ggplot() + annotate("text", x = 0, y = 0, size=4, label = text) + theme_void()+ plot_layout(tag_level = 'new')

text = paste("PF primary PF boost: Delta")
cap3<- ggplot() + annotate("text", x = 0, y = 0, size=4, label = text) + theme_void()+ plot_layout(tag_level = 'new')

text = paste("PF primary PF boost: Omicron")
cap4<- ggplot() + annotate("text", x = 0, y = 0, size=4, label = text) + theme_void()+ plot_layout(tag_level = 'new')

combined <- cap1 / plot_AZ_profiles_delta / cap2 / plot_AZ_profiles_omicron /
              cap3 / plot_PF_profiles_delta / cap4 / plot_PF_profiles_omicron +
           plot_annotation(tag_levels = "A") +   plot_layout(heights = c(2,5,2,5,2,5,2,5))
combined

if(data=="Imperial") {
  ggsave("../Figures/Figure 1 Imperial.png", combined, height = 13, width = 10)
}

if(data=="PHE") {
  ggsave("../Figures/Figure 1 PHE.png", combined, height = 13, width = 10)
}

