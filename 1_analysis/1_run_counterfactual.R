name <- "counterfactual"
fit1 <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE"

R_profile <- read_csv("data/category_1_Rt.csv") 
saveRDS(R_profile,"data/R_profile.rds")

#### Get vaccine parameters  ##############################################
vaccine <- "Pfizer"
t_d3_set <- 180

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
income_group <- "HIC"
hs_constraints <- "Present"
dt <- 0.25
repetition <- 1:20
vacc_start <- "1/1/2021"
vaccine_doses <- 2
max_coverage <- 0
age_groups_covered <- 15
age_groups_covered_d3 <- 5
seeding_cases <- 10
vacc_per_week <- 0.05
ab_model_infection <- TRUE
strategy <- "realistic"
t_d3 <- t_d3_set
max_Rt <- 3.5
max_Rt_omicron <- 3.5
vfr_time1 <- "11/27/2021"
vfr_time2 <- "12/31/2021"
vfr <- unique(vacc_params$vfr)[1]
hosp_scal_omicron <- 0.5
ICU_scal_omicron <-  0.42
mu_ab_infection <- 1.7 

#### Create scenarios ##########################################################

scenarios <- expand_grid(fit = fit,
                         income_group = income_group,
                         target_pop = target_pop,
                         hs_constraints = hs_constraints,
                         vaccine_doses = vaccine_doses,
                         vaccine = vaccine,
                         max_coverage = max_coverage,
                         age_groups_covered = age_groups_covered,
                         age_groups_covered_d3 = age_groups_covered_d3,
                         vacc_start = vacc_start,
                         dt = dt,
                         repetition = repetition,
                         seeding_cases = seeding_cases,
                         vacc_per_week = vacc_per_week,
                         ab_model_infection = ab_model_infection,
                         t_d3 = t_d3,
                         max_Rt = max_Rt,
                         max_Rt_omicron = max_Rt_omicron,
                         vfr = vfr,
                         vfr_time1 = vfr_time1,
                         vfr_time2 = vfr_time2,
                         mu_ab_infection = mu_ab_infection,
                         hosp_scal_omicron = hosp_scal_omicron,
                         ICU_scal_omicron = ICU_scal_omicron) %>%
  filter((vaccine_doses == 2 & age_groups_covered_d3 == 5 ) | (vaccine_doses == 3) ) %>%
  unique()

scenarios$scenario <- 1:nrow(scenarios)
scenarios$name <- name
scenarios$strategy <- strategy

scenarios <- left_join(scenarios, vacc_params, by = c("vaccine", "vfr"))

nrow(scenarios)

write_csv(scenarios, paste0("scenarios_", name, ".csv"))

#### Run the model on cluster ###############################################
# Load functions
sources <- c("R/run_function_abmodel_omicron.R", "R/utils.R", "R/vaccine_strategy.R")
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

