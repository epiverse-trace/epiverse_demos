# Simulating epidemic attack rates -------------------------------------------------
# with heterogeneous social contacts
# This script builds this vignette:
# https://epiverse-trace.github.io/finalsize/articles/varying_contacts.html
# now maintained in
# how-to guide: https://epiverse-trace.github.io/howto/analyses/simulate_transmission/finalsize-attack-rate-heterogeneity.html


# Load packages --------------------------------------------------------------------
library(finalsize)
library(socialmixr)
library(tidyverse)


# Simple quick calculation with homogenous mixing ----------------------------------
r0_input <- 2
finalsize::final_size(r0 = r0_input)

# Set up the transmission model -------------------------------------------
# Note this is the same structure as the `scenarios_with_uncertainty.R` example
# So contact_data object can be reused in both scripts

# load contact and population data from socialmixr::polymod
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 5, 18, 40, 65),
  symmetric = TRUE
)

# prepare contact matrix and demography vector for use in model
contact_matrix <- t(contact_data$matrix) # transpose so R0 calculated correctly inside model
demography_vector <- contact_data$demography$population
names(demography_vector) <- rownames(contact_matrix)

# scale the contact matrix so the largest eigenvalue is 1.0
# this is to ensure that the overall epidemic dynamics correctly reflect
# the assumed value of R0
contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))

# divide each row of the contact matrix by the corresponding demography
# this reflects the assumption that each individual in group {j} make contacts
# at random with individuals in group {i}
contact_matrix <- contact_matrix / demography_vector

n_demo_grps <- length(demography_vector)

# all individuals are equally and highly susceptible
n_susc_groups <- 1L
susc_guess <- 1.0

susc_uniform <- matrix(
  data = susc_guess,
  nrow = n_demo_grps,
  ncol = n_susc_groups
)

# Final size calculations also need to know the proportion of each demographic group {𝑖} 
# that falls into the susceptibility group {𝑗}. This distribution of age groups into 
# susceptibility groups can be represented by the demography-susceptibility distribution matrix.
p_susc_uniform <- matrix(
  data = 1.0,
  nrow = n_demo_grps,
  ncol = n_susc_groups
)


output <- finalsize::final_size(
  r0 = r0_input,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susc_uniform,
  p_susceptibility = p_susc_uniform
)

output

output %>% 
  mutate(demo_grp = as_factor(demo_grp)) %>% 
  ggplot(aes(x = demo_grp, y = p_infected)) +
  geom_col() +
  ylim(0,1) +
  labs(
    x = "Age group",
    y = "Proportion infected",
    title = "Final size of an SIR epidemic",
    subtitle = "Fully susceptible population"
  )
