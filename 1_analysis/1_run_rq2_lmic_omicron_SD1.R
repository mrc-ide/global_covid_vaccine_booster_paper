name <- "rq2_lmic_abmodel_omicron_SD1"

fit1 <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE"
fit <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_SD1"

R_profile <- read_csv("data/category_2_Rt.csv") 
saveRDS(R_profile,"data/R_profile_lmic.rds")

#### Get vaccine parameters  ##############################################
vaccine <- "Oxford-AstraZeneca"

vacc_names <- data.frame(vaccine = c("Pfizer", "Oxford-AstraZeneca"), vacc = c("PF", "AZ"))
vaccine_set <- vaccine
vacc_params <- readRDS(paste0("data/param_list_",fit1,".rds")) %>%
  rename(vacc = vaccine) %>%
  left_join(vacc_names, by = "vacc") %>%
  filter(vaccine == vaccine_set) %>%
  mutate(std10 = 0.44) %>%
  filter(vfr > 1) %>%
  select(-c(vacc))

#### Set up other simulation parameters  ##############################################

target_pop <- 1e6
income_group <- "LMIC"
hs_constraints <- "Present"
dt <- 0.25
repetition <- 1:20
vacc_start <- "4/1/2021"
vaccine_doses <- c(2,3)
max_coverage <- 0.8
age_groups_covered <- c(5, 9)
seeding_cases <- 10
vacc_per_week <- c(0.02, 0.015)
ab_model_infection <- TRUE
strategy <- "same_doses"
t_d3 <- 180
max_Rt <- 4
#max_Rt_omicron <- 4.25 #c(4.25,4.5,4.75)
vfr_time1 <- "11/27/2021"
vfr_time2 <- "12/31/2021"
vfr <- unique(vacc_params$vfr)
#hosp_scal_omicron <- c(0.4,0.5,0.6)
ICU_scal_omicron <-  0.42
mu_ab_infection <- 1.7 

central_vfr <- vfr[1]
corr_params <- data.frame(vfr = vfr[order(vfr)], max_Rt_omicron = c(3.75, 4, 4.25), hosp_scal_omicron = c(0.3,0.5,0.6))

#### Create scenarios ##########################################################

scenarios <- expand_grid(fit = fit,
                         income_group = income_group,
                         target_pop = target_pop,
                         hs_constraints = hs_constraints,
                         vaccine_doses = vaccine_doses,
                         vaccine = vaccine,
                         max_coverage = max_coverage,
                         age_groups_covered = age_groups_covered,
                         vacc_start = vacc_start,
                         dt = dt,
                         repetition = repetition,
                         seeding_cases = seeding_cases,
                         vacc_per_week = vacc_per_week,
                         ab_model_infection = ab_model_infection,
                         t_d3 = t_d3,
                         max_Rt = max_Rt,
                         #max_Rt_omicron = max_Rt_omicron,
                         vfr = vfr,
                         vfr_time1 = vfr_time1,
                         vfr_time2 = vfr_time2,
                         mu_ab_infection = mu_ab_infection,
               #          hosp_scal_omicron = hosp_scal_omicron,
                         ICU_scal_omicron = ICU_scal_omicron) %>%
  filter(!(vfr != central_vfr & vacc_per_week != 0.02 & t_d3 != 180)) %>%
#  filter(!(vfr != central_vfr & hosp_scal_omicron != 0.5)) %>%
#  filter(!(vacc_per_week != 0.02 & hosp_scal_omicron != 0.5)) %>%
  unique()

scenarios$scenario <- 1:nrow(scenarios)
scenarios$name <- name
scenarios$strategy <- strategy

scenarios <- left_join(scenarios, vacc_params, by = c("vaccine", "vfr")) %>%
  left_join(corr_params, by = "vfr")

nrow(scenarios)

write_csv(scenarios, paste0("scenarios_", name, ".csv"))

## test on PC

source("R/run_function_abmodel_omicron_lmic.R")
source("R/utils.R")
source("R/vaccine_strategy.R")
plan(multicore, workers = 4)
system.time({out <- future_pmap(scenarios, run_scenario, .progress = TRUE)})



#### Run the model on cluster ###############################################
# Load functions
sources <- c("R/run_function_abmodel_omicron_lmic.R", "R/utils.R", "R/vaccine_strategy.R")
src <- conan::conan_sources(c("mrc-ide/safir", "mrc-ide/squire", "mrc-ide/nimue"))
ctx <- context::context_save("context",
                             sources = sources,
                             packages = c("tibble", "dplyr", "tidyr", "countrycode", "safir", "nimue", "squire", "data.table"),
                             package_sources = src)

config <- didehpc::didehpc_config(use_rrq = FALSE, use_workers = FALSE, cluster="fi--didemrchnb")
#config <- didehpc::didehpc_config(use_rrq = FALSE, use_workers = FALSE, cluster="fi--dideclusthn")

# Create the queue
run <- didehpc::queue_didehpc(ctx, config = config)
# Summary of all available clusters
# run$cluster_load(nodes = FALSE)
# Run
runs <- run$enqueue_bulk(scenarios, run_scenario, do_call = TRUE, progress = TRUE)
runs$status()
