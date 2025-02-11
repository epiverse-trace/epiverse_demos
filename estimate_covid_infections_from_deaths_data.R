# ============================================================================== #
# SETUP AND DATA PREPARATION
# ============================================================================== #

# This script aims to reconstruct the dynamics of SARS-CoV-2 infections in the UK
# using daily death data. It follows the methodology outlined in the Epiverse-TRACE
# howto guide: https://epiverse-trace.github.io/howto/analyses/reconstruct_transmission/estimate_infections.html

# Load necessary packages for analysis
library(incidence2)   # For UK COVID daily deaths data
library(EpiNow2)      # To estimate time-varying reproduction number
library(epiparameter) # To access delay distributions
library(dplyr)        # For data manipulation
library(ggplot2)      # For plotting
library(parallel)     # For parallel processing
library(withr)        # For setting local options

# Set the number of cores for faster processing
withr::local_options(list(mc.cores = parallel::detectCores() - 1))

# Extract and preprocess data on UK COVID deaths for EpiNow2
uk_covid_deaths <- incidence2::covidregionaldataUK %>%
  tidyr::replace_na(list(deaths_new = 0)) %>%
  incidence2::incidence(
    date_index = "date",
    counts = "deaths_new",
    count_values_to = "confirm",
    date_names_to = "date",
    complete_dates = TRUE
  ) %>%
  dplyr::select(-count_variable) %>%
  dplyr::filter(date < "2020-07-01" & date >= "2020-03-01") %>%
  dplyr::as_tibble()

# Display the preprocessed incidence data
uk_covid_deaths

# ============================================================================== #
# DEFINE EPIDEMIOLOGICAL PARAMETERS AND DISTRIBUTIONS
# ============================================================================== #

# Extract distribution the incubation period for COVID-19
covid_incubation_dist <- epiparameter::epiparameter_db(
  disease = "covid",
  epi_name = "incubation",
  single_epiparameter = TRUE
)

# Display the incubation period distribution information
covid_incubation_dist

# Extract parameters for EpiNow2 using LogNormal distribution
incubation_params <- epiparameter::get_parameters(covid_incubation_dist)
incubation_max_days <- round(quantile(covid_incubation_dist, 0.999))  # Upper 99.9% range needed for EpiNow2

# Create a LogNormal object for the incubation period
incubation_lognormal <- EpiNow2::LogNormal(
  meanlog = incubation_params[["meanlog"]],
  sdlog = incubation_params[["sdlog"]],
  max = incubation_max_days
)

# Get the onset-to-death distribution
onset_to_death_dist <- epiparameter::epiparameter_db(
  disease = "covid",
  epi_name = "onset to death",
  single_epiparameter = TRUE
)

# Display the onset-to-death distribution information
onset_to_death_dist

# Extract parameters for EpiNow2 using LogNormal distribution
onset_death_params <- epiparameter::get_parameters(onset_to_death_dist)
onset_death_max_days <- round(quantile(onset_to_death_dist, 0.999))  # Upper 99.9% range

# Create an EpiNow2 LogNormal object for the onset-to-death distribution
onset_death_lognormal <- EpiNow2::LogNormal(
  meanlog = onset_death_params[["meanlog"]],
  sdlog = onset_death_params[["sdlog"]],
  max = onset_death_max_days
)

# Use EpiNow2 to convolve the infection-to-onset and onset-to-death distributions
infection_to_death_dist <- incubation_lognormal + onset_death_lognormal

# Extract the serial interval distribution
serial_interval_dist <- epiparameter::epiparameter_db(
  disease = "covid",
  epi_name = "serial",
  single_epiparameter = TRUE
)

# Display the serial interval distribution information
serial_interval_dist

# Find the upper 99.9% range for the serial interval
serial_interval_max_days <- round(quantile(serial_interval_dist, 0.999))

# Extract parameters for the serial interval
serial_interval_params <- epiparameter::get_parameters(serial_interval_dist)

# Create a LogNormal object for the serial interval
serial_interval_lognormal <- EpiNow2::LogNormal(
  meanlog = serial_interval_params[["meanlog"]],
  sdlog = serial_interval_params[["sdlog"]],
  max = serial_interval_max_days
)

# ============================================================================== #
# ESTIMATE INFECTIONS AND VISUALIZE RESULTS
# ============================================================================== #

# Estimate infections using the non-mechanistic model
infection_estimates <- EpiNow2::epinow(
  data = uk_covid_deaths,  # Time series data
  generation_time = EpiNow2::generation_time_opts(serial_interval_lognormal),  # Generation time
  delays = EpiNow2::delay_opts(infection_to_death_dist),  # Delay from infection to death
  rt = NULL  # Rt is not estimated
)

# Extract and filter infection estimates from the model output
estimated_infections <- infection_estimates$estimates$summarised %>%
  dplyr::filter(variable == "infections")

# Annotate the infections plot with key interventions in the UK
infection_estimates$plots$infections +
  geom_vline(aes(xintercept = as.Date("2020-03-16")), linetype = 3) +
  geom_text(
    aes(
      x = as.Date("2020-03-16"),
      y = 3000,
      label = "Non-essential contact advice"
    ),
    hjust = 0
  ) +
  geom_vline(aes(xintercept = as.Date("2020-03-23")), linetype = 3) +
  geom_text(
    aes(
      x = as.Date("2020-03-23"),
      y = 2500,
      label = "Stay-at-home order (i.e. lockdown)"
    ),
    hjust = 0
  ) +
  labs(
    title = "Estimated Dynamics of SARS-CoV-2 Infections in the UK",
    subtitle = "Reconstructed using data on reported deaths. Dashed lines indicate key intervention dates."
  )
