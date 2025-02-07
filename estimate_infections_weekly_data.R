# ============================================================================== #
# SETUP AND DATA PREPARATION
# ============================================================================== #

# This script aims to estimate Rt from weekly reported SARS-CoV-2 infections
# in the UK using EpiNow2 and EpiEstim. The EpiEstim example follows the
# methodology outlined in the EpiEstim vignette in
# https://mrc-ide.github.io/EpiEstim/articles/EpiEstim_aggregated_data.html.

# Load necessary packages for analysis
library(EpiNow2)      # To estimate time-varying reproduction number
library(EpiEstim)     # To estimate time-varying reproduction number
library(epiparameter) # To extract epidemiological parameters
library(patchwork)    # For plot composition
library(data.table)   # For data manipulation
library(parallel)     # For parallel processing
library(withr)        # For setting local options
library(dplyr)

# Set the number of cores for faster processing
local_options(list(mc.cores = parallel::detectCores() - 1))

# Use the example data for confirmed cases from EpiNow2
reported_cases <- EpiNow2::example_confirmed
reported_cases_weekly <- data.table::copy(reported_cases)

# Aggregate the daily cases to weekly cases (sum of daily cases)
reported_cases_weekly[, confirm := frollsum(confirm, 7)]
reported_cases_weekly <- reported_cases_weekly[seq(7, nrow(reported_cases_weekly), 7)]

# ============================================================================== #
# DEFINE EPIDEMIOLOGICAL PARAMETERS AND DISTRIBUTIONS
# ============================================================================== #

# Extract distribution the incubation period for COVID-19
covid_incubation_dist <- epiparameter_db(
    disease = "covid",
    epi_name = "incubation",
    single_epiparameter = TRUE
)

# Extract the serial interval distribution
serial_interval_dist <- epiparameter_db(
    disease = "covid",
    epi_name = "serial",
    single_epiparameter = TRUE
)

# ============================================================================== #
# ESTIMATE INFECTIONS AND Rt WITH EPINOW2
# ============================================================================== #

# Extract parameters and maximum of the distribution for EpiNow2
incubation_params <- epiparameter::get_parameters(covid_incubation_dist)
incubation_max_days <- round(quantile(covid_incubation_dist, 0.999))  # Upper 99.9% range needed for EpiNow2

# Create a LogNormal object for the incubation period
incubation_lognormal <- EpiNow2::LogNormal(
    meanlog = incubation_params[["meanlog"]],
    sdlog = incubation_params[["sdlog"]],
    max = incubation_max_days
)

# Extract parameters and maximum of the distribution for EpiNow2
serial_interval_params <- epiparameter::get_parameters(serial_interval_dist)
serial_interval_max_days <- round(quantile(serial_interval_dist, 0.999)) # Upper 99.9% range needed for EpiNow2

# Create a LogNormal object for the serial interval
serial_interval_lognormal <- EpiNow2::LogNormal(
    meanlog = serial_interval_params[["meanlog"]],
    sdlog = serial_interval_params[["sdlog"]],
    max = serial_interval_max_days
)

# Create data with missing dates filled in for EpiNow2
reported_cases_weekly_epinow <- fill_missing(
    reported_cases_weekly,
    missing_dates = "accumulate",
    initial_accumulate = 1 # Don't model the first data point (to match EpiEstim method)
)

# Estimate infections using EpiNow2
estimates_epinow <- epinow(
    data = reported_cases_weekly_epinow,
    generation_time = generation_time_opts(serial_interval_lognormal),
    delays = delay_opts(incubation_lognormal),
    forecast = forecast_opts(horizon = 0, accumulate = 1), # Forecasting is turned off to match with EpiEstim
    rt = NULL, # Use the non-mechanistic model equivalent to EpiEstim
    verbose = FALSE
)

# Initial look at the output
plot(Rt_epinow$plots$R)

# ============================================================================== #
# ESTIMATE RT WITH EPIESTIM
# ============================================================================== #

# Prepare serial interval distribution. We'll reuse the serial interval distribution
# extracted earlier.
si_mean <- serial_interval_dist$summary_stats$mean
si_sd <- serial_interval_dist$summary_stats$sd

# Prepare the data by renaming the columns and converting the date to a Date object
epiestim_data <- reported_cases_weekly |>
    rename(I = confirm) |>
    mutate(dates = as.Date(date), .keep = "unused") |>
    as.data.frame()

# Estimate infections and Rt using weekly aggregated data
estimates_epiestim <- estimate_R(
    incid = epiestim_data,
    dt = 7L, # Aggregation window
    dt_out = 7L, # Estimation rolling window
    recon_opt = "naive",
    method = "parametric_si",
    config = make_config(list(mean_si = si_mean, std_si = si_sd))
)

# Initial look at the output
plot(estimates_epiestim$R, "R") # Rt estimates only

# ==============================================================================
# Comparing the results from EpiNow2 and EpiEstim
# ==============================================================================
# Extract and process the Rt estimates from EpiEstim output
epiestim_Rt <- mutate(estimates_epiestim$R, method = "EpiEstim")

# Add the dates in the complete data to the Rt estimates
Rt_ts_epiestim <- reported_cases_weekly_epinow |>
    mutate(lookup = 1:nrow(reported_cases_weekly_epinow)) |>
    left_join(
        epiestim_Rt,
        by = join_by(lookup == t_start)
    ) |>
    select(-c(lookup, confirm, accumulate))

# Extract and process the Rt estimates from EpiNoww output
Rt_ts_epinow <- estimates_epinow$estimates$summarised[variable == "R"][date >= min(Rt_ts_epiestim$date, na.rm = TRUE)]

