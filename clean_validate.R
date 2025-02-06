# HPRU workshop script

## Clean and validate and plot outbreak data

# Load required R packages ------------------------------------------------

# {simulist} needs to be installed from a feature branch for now
# remotes::install_github("epiverse-trace/simulist@messy")
library(simulist)
library(tibble)
library(english)
library(cleanepi)
library(numberize)
library(incidence2)
library(tidyr)
library(dplyr)

# Choose a seed that results in suitable and reproducible outbreak --------

set.seed(1)

# Simulate outbreak -------------------------------------------------------

linelist <- sim_linelist()

# optionally run commented out code below to can convert
# to tibble for tidier printing
# linelist <- tibble::as_tibble(linelist)

linelist

# Create messy line list data ---------------------------------------------

linelist <- messy(linelist)

# convert some of the ages into words
linelist$age <- english::words(as.numeric(linelist$age))

linelist

# Tag line list of data validation ----------------------------------------

# see what tags are available
linelist::tags_names()

# in this case the tags have the same name but line list columns can be
# named differently from the tag names
linelist <- linelist::make_linelist(
  x = linelist,
  date_onset = "date_onset",
  date_admission = "date_admission",
  date_outcome = "date_outcome"
)
linelist

# line list can be validated using tags
# this will error due to the line list being messy
linelist::validate_linelist(linelist)

# Scan line list data for issues ------------------------------------------

cleanepi::scan_data(linelist)

cleanepi::check_subject_ids(linelist, target_columns = "id", range = c(1, 350))

# Clean line list ---------------------------------------------------------

linelist$age <- numberize::numberize(linelist$age)
linelist$age

# zeros need to be changed to NA (in discussion)
linelist$age[linelist$age == 0L] <- NA

# routine cleaning steps to tidy column names and remove duplicated rows
clean_linelist <- linelist |>
  cleanepi::standardize_column_names() |>
  cleanepi::remove_constants() |>
  cleanepi::remove_duplicates()

date_columns <- colnames(linelist)[startsWith(colnames(linelist), "date_")]
linelist <- linelist |>
  cleanepi::standardize_dates(target_columns = date_columns)


# clean inconsistent sex using dictionary
linelist |> cleanepi::clean_using_dictionary(
  dictionary = data.frame())

# clean spelling mistakes using dictionary
linelist$case_type[agrep(pattern = "suspected", x = linelist$case_type)] <- "suspected"
linelist$case_type[agrep(pattern = "probable", x = linelist$case_type)] <- "probable"
linelist$case_type[agrep(pattern = "confirmed", x = linelist$case_type)] <- "confirmed"

linelist$outcome[agrep(pattern = "recovered", x = linelist$outcome)] <- "recovered"
linelist$outcome[agrep(pattern = "died", x = linelist$outcome)] <- "died"

# Validate clean line list ------------------------------------------------

# line list is now valid after cleaning
linelist::validate_linelist(linelist)


# Aggregate and visualise data --------------------------------------------

# aggregate to daily incidence data
daily <- incidence(x = linelist, date_index = "date_onset", interval = "daily", complete_dates = TRUE)
plot(daily)

# aggregate to epiweek incidence data
weekly <- incidence(x = linelist, date_index = "date_onset", interval = "epiweek", complete_dates = TRUE)
plot(weekly)

# transform line list to aggregate and plot onset, hospital admission and death
linelist <- linelist |>
  tidyr::pivot_wider(
    names_from = outcome,
    values_from = date_outcome
  )
linelist <- linelist |>
  dplyr::rename(
    date_death = died,
    date_recovery = recovered
  )

# aggregate onset, admission and death
weekly_chd <- incidence(
  linelist,
  date_index = c(
    onset = "date_onset",
    hospitalisation = "date_admission",
    death = "date_death"
  ),
  interval = "epiweek",
  complete_dates = TRUE
)
plot(weekly_chd)
