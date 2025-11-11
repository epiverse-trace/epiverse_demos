# ============================================================================== #
# SETUP AND DATA PREPARATION
# ============================================================================== #

# This script aims to reconstruct the dynamics of SARS-CoV-2 infections in the UK
# using daily death data. It follows the methodology outlined in the Epiverse-TRACE
# howto guide: https://epiverse-trace.github.io/howto/analyses/reconstruct_transmission/estimate_infections.html

# Load required packages
library(incidence2) # for uk covid daily deaths
library(EpiNow2) # to estimate time-varying reproduction number
library(epiparameter) # to access delay distributions
library(dplyr) # to format input and outputs
library(ggplot2) # to generate plots

# Set number of cores
withr::local_options(list(mc.cores = 4))

# Extract data on UK COVID deaths and format for EpiNow2
incidence_data <- incidence2::covidregionaldataUK %>% 
  # preprocess missing values
  tidyr::replace_na(list(deaths_new = 0)) %>%
  # compute the daily incidence
  incidence2::incidence(
    date_index = "date",
    counts = "deaths_new",
    count_values_to = "confirm",
    date_names_to = "date",
    complete_dates = TRUE
  ) %>%
  dplyr::select(-count_variable) %>% 
  # Focus on early 2020 period and sort by ascending date
  dplyr::filter(date<"2020-07-01" & date>="2020-03-01") %>% 
  # convert to tibble format for simpler data frame output
  dplyr::as_tibble()

# Preview data
incidence_data

# Define parameters
# Extract infection-to-death distribution (from Aloon et al)
incubation_period_in <-
  epiparameter::epiparameter_db(
    disease = "covid",
    epi_name = "incubation",
    single_epiparameter = TRUE
  )

# Summarise distribution and type
print(incubation_period_in)

# Get parameters and format for EpiNow2 using LogNormal input
incubation_params <- epiparameter::get_parameters(incubation_period_in)

# Find the upper 99.9% range by the interval
incubation_max <- round(quantile(incubation_period_in,0.999))

incubation_period <- EpiNow2::LogNormal(
  meanlog = incubation_params[["meanlog"]], 
  sdlog = incubation_params[["sdlog"]], 
  max = incubation_max
)

## Set onset to death period (from Linton et al)
onset_to_death_period_in <-
  epiparameter::epiparameter_db(
    disease = "covid",
    epi_name = "onset to death",
    single_epiparameter = TRUE
  )

# Summarise distribution and type
print(onset_to_death_period_in)

# Get parameters and format for EpiNow2 using LogNormal input
onset_to_death_params <- epiparameter::get_parameters(onset_to_death_period_in)

# Find the upper 99.9% range by the interval
onset_to_death_max <- round(quantile(onset_to_death_period_in,0.999))

onset_to_death_period <- LogNormal(
  meanlog = onset_to_death_params[["meanlog"]], 
  sdlog = onset_to_death_params[["sdlog"]], 
  max = onset_to_death_max
)

## Combine infection-to-onset and onset-to-death
infection_to_death <- incubation_period + onset_to_death_period

# Plot underlying delay distributions
# plot(infection_to_death)

# Extract serial interval distribution distribution (from Yang et al)
serial_interval_in <-
  epiparameter::epiparameter_db(
    disease = "covid",
    epi_name = "serial",
    single_epiparameter = TRUE
  )

# Summarise distribution and type
print(serial_interval_in)

# Discretise serial interval for input into EpiNow2
serial_int_discrete <- epiparameter::discretise(serial_interval_in)

# Find the upper 99.9% range by the interval
serial_int_discrete_max <- quantile(serial_int_discrete,0.999)

# Get parameters
serial_params <- epiparameter::get_parameters(serial_int_discrete)

# Define parameters using LogNormal input
serial_interval_covid <- LogNormal(
  meanlog = serial_params[["meanlog"]],
  sdlog = serial_params[["sdlog"]],
  max = serial_int_discrete_max
)
# Run infection estimation model
epinow_estimates <- epinow(
  data = incidence_data, # time series data
  # assume generation time = serial interval
  generation_time = generation_time_opts(serial_interval_covid),
  # delay from infection-to-death
  delays = delay_opts(infection_to_death),
  # no Rt estimation
  rt = NULL,
  # change default Gaussian Process priors
  gp = gp_opts(alpha = Normal(0, 0.05))
)

# Extract infection estimates from the model output
infection_estimates <- epinow_estimates$estimates$summarised %>% 
  dplyr::filter(variable=="infections")

# Plot output
epinow_estimates$plots$infections +
  geom_vline(aes(xintercept = as.Date("2020-03-16")), linetype = 3) +
  geom_text(aes(x = as.Date("2020-03-16"), 
                y = 3000,
                label = "Non-essential contact advice"),
            hjust = 0) +
  geom_vline(aes(xintercept = as.Date("2020-03-23")), linetype = 3) +
  geom_text(aes(x = as.Date("2020-03-23"), 
                y = 2500,
                label = "Stay-at-home order (i.e. lockdown)"),
            hjust = 0) +
  labs(
    title = "Estimated dynamics of SARS-CoV-2 infections
    among those with subsequent fatal outcomes in the UK,
    reconstructed using data on reported deaths.",
    subtitle = "Dashed lines show dates of
    UK non-essential contact advice (16 Mar)
    and lockdown (23 Mar)."
  )

#' Further exploration

#' - In EpiNow2::epinow(), explore the following priors on the Gaussian process and observe what changes in the output:
#'    - alpha: gp = gp_opts(alpha = Normal(0, 0.05))
#'    - length scale: gp = gp_opts(ls = LogNormal(mean = 14, sd = 7))
