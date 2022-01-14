
get_representative_country <- function(income_group){
  case_when(income_group == "HIC" ~ "Malta",
            income_group == "UMIC" ~ "Grenada",
            income_group == "LMIC" ~ "Nicaragua",
            income_group == "LIC" ~ "Madagascar")
}


get_capacity <- function(country, income_group, pop, hs_constraints){
  hc <- squire::get_healthcare_capacity(country = country)
  
  # Unconstrained healthcare
  if(hs_constraints == "Absent"){
    hc$hosp_beds <- 1000000
    hc$ICU_beds <- 1000000
  }
  
  if(hs_constraints == "Present"){
    if(income_group %in% c("HIC", "UMIC")){
      hc$hosp_beds <- 1000000
      hc$ICU_beds <- 1000000
    }
    if(income_group %in% c("LMIC", "LIC")){
      hc$ICU_beds <- 0
    }
  }
  
  hc$hosp_beds <- round(hc$hosp_beds * sum(pop) / 1000)
  hc$ICU_beds <- round(hc$ICU_beds * sum(pop) / 1000)
  
  return(hc)
}

get_prob_non_severe_death_treatment <- function(income_group, hs_constraints){
  psdt <- squire:::probs$prob_non_severe_death_treatment
  
  if(income_group  == "LIC" & hs_constraints == "Present"){
    psdt <- c(rep(0.25, 16), 0.5804312)
  }
  return(psdt)
}

get_vaccine_pars <- function(
  vaccine,
  mu_ab_d1,
  mu_ab_d2,
  vaccine_doses,
  dose_3_fold_increase,
  ab_50,
  ab_50_severe,
  std10 = 0.44,
  k,
  t_d2 = 28,
  t_d3,
  hl_s,
  hl_l,
  period_s,
  period_l
){
  mu_ab_list <- data.frame(name = vaccine,
                           mu_ab_d1 = mu_ab_d1,
                           mu_ab_d2 = mu_ab_d2) %>%
    mutate(mu_ab_d3 = mu_ab_d2 * dose_3_fold_increase)
  
  ab_parameters <- safir::get_vaccine_ab_titre_parameters(
    vaccine = vaccine, max_dose = vaccine_doses, correlated = TRUE,
    hl_s = hl_s, hl_l = hl_l, period_s = period_s, t_period_l = period_l,
    ab_50 = ab_50, ab_50_severe = ab_50_severe, std10 = std10, k = k,
    mu_ab_list = mu_ab_list
  )
  ab_parameters$max_ab <- 5 # max titre on natural log scale
  return(ab_parameters) 
}
