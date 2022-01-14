library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(drjacoby)
library(matrixStats)

setwd(dirname(getActiveDocumentContext()$path))

# Read in parameters from MCMC chain 

load("../../Model Fitting/imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_mcmc_chain.Rdata")  

#############################################
# get the params we are using


chain <- mcmc$output %>%
  filter(phase == "sampling") %>%
  select(-c(chain, phase, iteration, logprior, loglikelihood)) 

#### pull out parameter medians - just need bounds on om_red ##

draws <- sample_chains(mcmc, 10000)

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

posterior_upper <- draws_transform %>%
  summarise( 
    across(where(is.numeric), quantile, 0.975)
  )%>%
  mutate(measure = "upper")

posterior_lower <- draws_transform %>%
  summarise( 
    across(where(is.numeric), quantile, 0.025)
  ) %>%
  mutate(measure = "lower")

params_est <- posterior_median_transform %>%
  rbind(posterior_upper) %>%
  rbind(posterior_lower) %>%
  pivot_longer(cols = c(d1_AZ, d1_PF, d2_PF, d2_AZ,d3,fold_red_AZ, fold_red_PF, om_red,ab50, ab50_s, ab50_d, k, hl_s, hl_l, period_s, period_l)) %>%
  pivot_wider(names_from = measure, values_from = value)

om_red_est <- params_est %>%
       filter(name=="om_red")

posterior_median <- chain %>%
  summarise( 
    across(where(is.numeric), median)
  )

log10_d2_PF <- log10(223/94)  - posterior_median$fold_red_PF

fold_red <- 10^posterior_median$fold_red_PF

ab_50       <- log10_d2_PF + posterior_median$ni50 
ab_50_s <- log10_d2_PF + posterior_median$ns50
ab_50_d  <- log10_d2_PF + posterior_median$nd50

k           <- posterior_median$k
hl_s        <- posterior_median$hl_s
hl_l        <- posterior_median$hl_l
period_s    <- posterior_median$period_s
period_l  <- posterior_median$period_l

om_red <- posterior_median$om_red


dr_s <- -log(2)/hl_s  # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

lg10 <- log(10)

## assume natural immunity is similar to dose 3 and generates a slower decay by halving the short period

period_s <- period_s/2

#############################################################################################################
### Calculate the mean efficacy over 1st year for a range of NAT boosts  
#############################################################################################################

max = 5.0
min = -2.0
stepsize = 0.01
array_size <- (max-min/stepsize)
titre <- c(rep(0,array_size))  
mean_eff_infection <- c(rep(0,array_size))  
mean_eff_severe <- c(rep(0,array_size))  
mean_eff_death <- c(rep(0,array_size))  

j <- 0  
for(i in seq (from = min, to = max, by=stepsize)) {

  j <- j + 1
  titre[j] <- i
  
  t <- 1:365
  nt <- c(rep(i,length(t)))
  
## simple biphasic decay  
  
  denom=log10(exp(dr_l*period_s)+exp(dr_s*period_s))
  cum_dr_vec=log10(exp(dr_s*t+dr_l*period_s)+exp(dr_l*t+dr_s*period_s))-denom
  dr_vec=c(0,diff(cum_dr_vec,1))*lg10
  
  nt <- nt + cum_dr_vec
  
  # relate titre to efficacy over time - using log-10 parameters
  ef_infection <- 1/(1+exp(-k*(nt -ab_50)))
  ef_severe <- 1/(1+exp(-k*(nt -ab_50_s)))
  ef_death <- 1/(1+exp(-k*(nt -ab_50_d)))
  
  mean_eff_infection[j] <- mean(ef_infection)
  mean_eff_severe[j] <- mean(ef_severe)
  mean_eff_death[j] <- mean(ef_death)
}

plot1 <- ggplot(data=NULL, aes(x=titre, y=mean_eff_infection) ) +
  geom_line()
plot1

plot2 <- ggplot(data=NULL, aes(x=titre, y=mean_eff_severe) ) +
  geom_line()
plot2

## match 90% over 1 year - NAT of 1.7 

l <- match(0.90, round(mean_eff_infection,2))
mean_eff_infection[l]
NAT <- 10^(titre[l])
NAT

mean_eff_severe[l]

RR_ratio = (1-mean_eff_severe[l])/(1-mean_eff_infection[l])
RR_ratio


# omicron reduction

new_t <- round(titre[l]-log10(om_red_est$median),2)
new_t 
i <-match(round(new_t,2),round(titre,2))
titre[i]
mean_eff_infection[i] # 60% protection from re-infection 

mean_eff_severe[i] #90% protection from hospitalisation

RR_ratio = (1-mean_eff_severe[i])/(1-mean_eff_infection[i])
RR_ratio

#upper 95%

new_t <- titre[l]-log10(om_red_est$lower)
new_t
i <- match(round(new_t,2),round(titre,2))
titre[i]
mean_eff_infection[i] # 68% protection from re-infection 
mean_eff_severe[i] #93% protection from hospitalisation


#lower 95%

new_t <- titre[l]-log10(om_red_est$upper)
new_t
i <- match(round(new_t,1),round(titre,2))
titre[i]
mean_eff_infection[i] # 50% protection from re-infection 
mean_eff_severe[i] #85% protection from hospitalisation


#################################################################
####  titre of 1.7 gives mean_eff_infection of 90% over 1st year
#################################################################



