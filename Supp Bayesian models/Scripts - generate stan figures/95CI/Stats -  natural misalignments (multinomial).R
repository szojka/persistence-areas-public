
#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 1d Bayes version
#-----------------------------------------------------------------------------.

###############################
# Load packages and functions
###############################

# library(readr)
# library(tidyverse)
# library(truncnorm)
# library(stringr)
# library(BiodiversityR)
# library(car)
# library(lme4)
# library(glmmTMB)
# library(visreg)
# library(Rmisc)
# library(patchwork)
# library(sf)
# library(mclogit)

source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Source - MAIN fitnessdata POST STAN.R"))

#--------------------------------------------------------------------------.
# Data cleaning
#--------------------------------------------------------------------------.

fig1B_dat <- plotlev %>%
  select(block, contingency, treatment, species, site, grid) %>%
  filter(treatment %in% "B") %>%
  distinct() %>%
  dplyr::group_by(block) %>%
  dplyr::mutate(replicates = n()) %>% 
  group_by(block, treatment, contingency) %>%
  dplyr::mutate(species_no = n()) %>%
  ungroup() %>%
  select(block, contingency, species_no, replicates, site, grid) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = species_no, values_fill = 0)

# variable types
str(fig1B_dat)
fig1B_dat$block <- as.factor(fig1B_dat$block )

#--------------------------------------------------------------------------.
## Fit a model (stan version)
#--------------------------------------------------------------------------.

mod_dat <- fig1B_dat

######################
# Define variables

# set up random effects
N <- as.integer(dim(mod_dat)[1])
J <- 4
S <- as.integer(6)
G <- as.integer(18)

# issue with stan wanting sequential blocks, while maintaining global ID of block across species:
mod_dat$block_stan_id <- match(mod_dat$block,
                             sort(unique(mod_dat$block))) # assigned integer values to first position
min(mod_dat$block_stan_id); max(mod_dat$block_stan_id)
block_number <- length(unique(mod_dat$block_stan_id)) # 423 entries
B <- as.integer(length(unique(mod_dat$block_stan_id))) # using block_number specific to each species data

Site <- as.integer(as.factor(mod_dat$site)) # random effects are incorporated into the likelihood loops

mod_dat$grid_id <- as.integer(mod_dat$grid)
Grid <- as.integer(mod_dat$grid_id) # grid_id is temporarily necessary for sequ., but when you append dataframe later it will not be an issue to keep using 'grid'
Block <- as.integer(mod_dat$block_stan_id)

y_mat <- as.matrix(mod_dat[, c(8,6,7,5)]) # order = DL SS_n SS_y ME

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms 95CI/nat_mis_fit.rdata"))

params <- c("alpha", "sigma_site", "sigma_grid")
stan_plot(nat_mis_fit, pars = params)
stan_trace(nat_mis_fit, pars = params) # Sigma block is why they did not converge. 

post_natmis <- rstan::extract(nat_mis_fit)
# need to extract different columns for each contingency
output_nat <- data.frame(contingency = c("DL", "SS_n", "SS_y", "ME"),
                         prediction = apply(post_natmis$probs, 3, median),
                         prediction_lwr = apply(post_natmis$probs, 3, quantile, probs = 0.025),
                         prediction_upr = apply(post_natmis$probs, 3, quantile, probs = 0.975))
