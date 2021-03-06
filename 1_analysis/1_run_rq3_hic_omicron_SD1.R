name <- "rq3_hic_abmodel_omicron_SD1"

fit1 <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE"
fit <- "imp_v2_20211219_AZPD2=FALSE_SB=TRUE_NewDecay=TRUE_SD1"


#### Get vaccine parameters  ##############################################
vaccine <- "Pfizer"

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
vaccine_doses <- c(2,3)
max_coverage <- 0.9
age_groups_covered <- 15
age_groups_covered_d3 <- c(5, 15)
seeding_cases <- 10
vacc_per_week <- 0.05
ab_model_infection <- TRUE
strategy <- "realistic"
t_d3 <- 180
max_Rt <- 3
vfr_time1 <- "11/27/2021"
vfr_time2 <- "12/31/2021"
vfr <- unique(vacc_params$vfr)
ICU_scal_omicron <-  0.42
mu_ab_infection <- 1.7 

central_vfr <- vfr[1]
corr_params <- data.frame(vfr = vfr[order(vfr)], max_Rt_omicron = c(2.75, 3, 3.25), hosp_scal_omicron = c(0.3,0.5,0.6))

# how many days to reach 80% vaccinated with 2 doses, at 5% vacc per dose per week?
days <- floor(0.8 / (vacc_per_week/7))*2
d1 <- as.Date(x = (as.Date("1/1/2021", format = "%m/%d/%Y") + days), format = "%m/%d/%Y")

R0_t3_in <- c("8/1/2021", "10/1/2021", "3/1/2022", "3/2/2022")


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
                         R0_t3_in = R0_t3_in,
                         max_Rt = max_Rt,
                         #max_Rt_omicron = max_Rt_omicron,
                         vfr = vfr,
                         vfr_time1 = vfr_time1,
                         vfr_time2 = vfr_time2,
                         mu_ab_infection = mu_ab_infection,
                        # hosp_scal_omicron = hosp_scal_omicron,
                         ICU_scal_omicron = ICU_scal_omicron
                         )  %>%
  filter((vaccine_doses == 2 & age_groups_covered_d3 == 15 ) | (vaccine_doses == 3) ) %>%
  filter(R0_t3_in == "8/1/2021" | (R0_t3_in == "10/1/2021" & (age_groups_covered_d3 %in% c(5,15))) | (R0_t3_in == "3/1/2022" & age_groups_covered_d3 == 15) | (R0_t3_in == "3/2/2022" & age_groups_covered_d3 == 15)) %>%
#  filter(vfr == central_vfr) %>%
  unique()

scenarios$scenario <- 1:nrow(scenarios)
scenarios$name <- name
scenarios$strategy <- strategy

scenarios <- left_join(scenarios, vacc_params, by = c("vaccine", "vfr")) %>%
  left_join(corr_params, by = "vfr")


nrow(scenarios)
write_csv(scenarios, paste0("scenarios_", name, ".csv"))

## test on PC
source("R/run_function_abmodel_omicron_zerocovid.R")
source("R/utils.R")
source("R/vaccine_strategy.R")
plan(multicore, workers = 4)
system.time({out <- future_pmap(scenarios, run_scenario, .progress = TRUE)})

#### Run the model on cluster ###############################################
# Load functions
sources <- c("R/run_function_abmodel_omicron_zerocovid.R", "R/utils.R", "R/vaccine_strategy.R")
src <- conan::conan_sources(c("mrc-ide/safir", "mrc-ide/squire", "mrc-ide/nimue"))
ctx <- context::context_save("context",
                             sources = sources,
                             packages = c("tibble", "dplyr", "tidyr", "countrycode", "safir", "nimue", "squire", "data.table"),
                             package_sources = src)

config <- didehpc::didehpc_config(use_rrq = FALSE, use_workers = FALSE, cluster="fi--didemrchnb")

# Create the queue
run <- didehpc::queue_didehpc(ctx, config = config)

# Run
runs <- run$enqueue_bulk(scenarios, run_scenario, do_call = TRUE, progress = TRUE)
runs$status()

