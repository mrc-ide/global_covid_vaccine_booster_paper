library(tidyverse)
library(dplyr)
library(purrr)
library(drjacoby)

##################################
##### LOAD THE PARAMETERS AND DATA
##################################

if(data=="Imperial") {
    name <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE"
}

###################################
### parameters from MCMC chain  
###################################

load(paste0("model fits/",name,"_mcmc_chain.Rdata"))

chain <- sample_chains(mcmc, 10000)

draws_transform <- chain %>%
  select(-sample, -AZ_ns_off ) %>%
  mutate(d2_AZ_log10 = log10(32/59) - fold_red_AZ,
         d2_PF_log10 = log10(223/94) - fold_red_PF,
         d1_AZ = 10^(d2_AZ_log10 + d1_AZ),
         d1_PF = 10^(d2_PF_log10 + d1_PF),
         fold_red_PF = 10^(fold_red_PF),
         fold_red_AZ = 10^(fold_red_AZ),
         d3_PF = 10^(d2_PF_log10 + bst_PF),
         ab50 = 10^(d2_PF_log10 + ni50),
         ab50_d = 10^(d2_PF_log10 + nd50),
         # note we are going to use the death efficacy estimates for our hospitalisations efficacy
         ab50_s = ab50_d,
         om_red = 10^(om_red),
         d2_PF = 10^(d2_PF_log10),
         d2_AZ = 10^(d2_AZ_log10),
         dose_3_fold_increase_PF = d3_PF/d2_PF,
         d3_AZ = d2_AZ * dose_3_fold_increase_PF
  ) %>%
  select(-ni50, -ns50, -nd50, -d2_AZ_log10, -d2_PF_log10, -bst_AZ, -bst_PF, -dose_3_fold_increase_PF)

posterior_median <-draws_transform %>%
  summarise( 
    across(where(is.numeric), median))

posterior_upper <- draws_transform %>%
  summarise( 
    across(where(is.numeric), quantile, 0.975))

posterior_lower <- draws_transform %>%
  summarise( 
    across(where(is.numeric), quantile, 0.025))
 

#####################################################
### calculate alternative decay rate vector
### use the three-parameter biphasic exponential decay
######################################################

hl_s <- posterior_median$hl_s
hl_l <- posterior_median$hl_l
period_s <- posterior_median$period_s

max_t     <- 730
t         <- 0:(max_t-1) #vaccinated on day 0
dr_s      <- -log(2)/hl_s # Corresponding decay rate in days for half life above
dr_l      <- -log(2)/hl_l

# simple biphasic decay implemented as sum of decaying exponentials

denom=log(exp(dr_l*period_s)+exp(dr_s*period_s))
cum_dr_vec=log(exp(dr_s*t+dr_l*period_s)+exp(dr_l*t+dr_s*period_s))-denom
dr_vec=c(0,diff(cum_dr_vec,1))

D1 <- dr_vec
D2 <- dr_vec

# for dose three have half the period_s to account for more long-lived cells

period_s <- period_s/2

denom=log(exp(dr_l*period_s)+exp(dr_s*period_s))
cum_dr_vec=log(exp(dr_s*t+dr_l*period_s)+exp(dr_l*t+dr_s*period_s))-denom
dr_vec=c(0,diff(cum_dr_vec,1))

D3 <- dr_vec

# natural immunity is the same as D3

N <- D3

# output to data frame for reading in

dr_vec_new <- data.frame(t,D1,D2,D3,N)

saveRDS(dr_vec_new, paste0("data/dr_vec_", name, "_SD1.rds"))

#################################
### save parameter list for runs
#################################

pm <- posterior_median
vaccine <- c("AZ", "PF")
vfr <- c(1,round(posterior_median$om_red, 1), round(posterior_lower$om_red,1),round(posterior_upper$om_red,1))

mu_ab_d1 <- c(pm$d1_AZ, pm$d1_PF)
mu_ab_d2 <- c(pm$d2_AZ, pm$d2_PF)
mu_ab_d3 <- c(pm$d3_AZ, pm$d3_PF)
k <- pm$k
hl_s <- pm$hl_s
hl_l <- pm$hl_l
period_s <- pm$period_s
period_l <- pm$period_l
dose_3_fold_increase <- mu_ab_d3/mu_ab_d2
ab_50 <- pm$ab50
ab_50_severe <- pm$ab50_d

param_list <- data.frame(vaccine,mu_ab_d1,mu_ab_d2, k, dose_3_fold_increase, hl_s, hl_l, period_s, period_l, ab_50, ab_50_severe) 

param_list_out <- 
  # Create input options
  expand_grid(
    vfr = vfr,
    vaccine = c("AZ", "PF")) %>%
  # Join with MCMC samples
  left_join(param_list, by = "vaccine") 

saveRDS(param_list_out, paste0("data/param_list_", name, ".rds"))

