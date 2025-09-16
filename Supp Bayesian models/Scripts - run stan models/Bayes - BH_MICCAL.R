############################
# Miccal population model
############################

library(rstan)
library(tidyverse)
library(coda)

############################
# Data load & cleaning 
############################

# prepare data
load(here::here("Data/dat_final_doubled.Rdata"))
dat_final_doubled <- na.omit(dat_final_doubled)
#names(dat_final_doubled)

# case match so that site is a number
dat_final_doubled$site <- case_match(dat_final_doubled$site,
                                     'HUT' ~ 1,
                                     'KNOX' ~ 2,
                                     'LM' ~ 3,
                                     'RES' ~ 4,
                                     'TAIL' ~ 5,
                                     'WHIT' ~ 6)

# Would like my blocks and grids to run sequentially:
# will deal with this within species (as not all species have data within each block!)
dat_final_doubled$block <- droplevels(dat_final_doubled$block)

# grids:
dat_final_doubled$grid_id <- as.integer(dat_final_doubled$grid)

grid_lookup <- data.frame(
  original_grid = levels(dat_final_doubled$grid),
  grid_id = seq_along(levels(dat_final_doubled$grid))
)

# case_match so that abundance bins from 0-3 go to mid abundance of each bin:
dat_final_doubled$intraspecific_abundance <- case_match(dat_final_doubled$intraspecific_abundance,
                                                        0 ~ 0,
                                                        1 ~ 5,
                                                        2 ~ 50,
                                                        3 ~ 101) # choose 100 instead of 'half-way point of >101)

# miccal:
focal<-as.data.frame(subset(dat_final_doubled,species=="miccal"))
#head(focal)

# issue with stan wanting sequential blocks, while maintaining global ID of block across species:
focal$block_stan_id <- match(focal$block,
                             sort(unique(focal$block))) # assigned integer values to first position
min(focal$block_stan_id)
max(focal$block_stan_id)
block_number <- length(unique(focal$block_stan_id)) # 391 entries

# key for which block original ID corresponds to miccals' stan ID:
miccal_block_key <- focal %>%
  select(block, block_stan_id) %>%
  distinct()

##############################
# DEFINE DATA & VARIABLES:
##############################

# response is per_captia_seed_prod:
Fecundity <- as.integer(focal$per_capita_seed_prod) # response for likelihood (y_hat)

# set up random effects
N <- as.integer(dim(focal)[1])
P <- as.integer(6)
G <- as.integer(18)
B <- as.integer(block_number) # using block_number specific to each species data
Site <- as.integer(as.factor(focal$site)) # random effects are incorporated into the likelihood loops
Grid <- as.integer(focal$grid_id) # grid_id is temporarily necessary for sequ., but when you append dataframe later it will not be an issue to keep using 'grid'
Block <- as.integer(focal$block_stan_id)

# variables parameters are dependent on:
Abundance <-  as.integer(focal$intraspecific_abundance) # this and cover have same min and max (0,100)
Cover <-  as.integer(focal$cover)
Env <- as.vector(focal$green_index_scaled) # this green_index_scaled variable already has been z-transformed via the 'scale' function. I had also multiplied it by -1 so that harsh plots are negative numbers and productive plots are positive numbers

##############################
# STAN MODEL
##############################

initials <- list(
  site_effect = rep(0, P), # exp(0) = 1 means no effect
  grid_effect = rep(0, G),
  sigma_site = 0.01,
  sigma_grid = 0.01
)
initials1 <- rep(list(initials), 4)  # For 4 chains

PrelimFit <- stan(file = here::here("Supp Bayesian models/Stan models/BH - ZIH model nested.stan"), data = c("Fecundity", "N", "Abundance", "Cover", "P", "Site", "G","Grid","Env"),iter = 20000, chains = 4, thin = 2, init = initials1) # "B","Block",

Miccal_posteriors <- PrelimFit # every parameter converged

save(Miccal_posteriors, here::here("Supp Bayesian models/Stan rdata/Posteriors/Miccal_posteriors.rdata"))

##############################
# DIAGNOSTIC PLOTS
##############################

# First check the distribution of Rhats and effective sample sizes
# hist(summary(PrelimFit)$summary[,"Rhat"])
# hist(summary(PrelimFit)$summary[,"n_eff"])

# visualize parameters
# params <- c("lambda_0", "lambda_1", 'alpha_intra_0', 'alpha_intra_1', 'alpha_other_0', 'alpha_other_1')
# stan_plot(PrelimFit, pars = params)
# stan_trace(PrelimFit, pars = params)

# print parameters
#summary(PrelimFit)$summary[1:14,]

# check the correlation among key model parameters
# pairs(PrelimFit, pars = c("lambda_0", "lambda_1",
#                                 "alpha_intra_0", "alpha_intra_1",
#                                 "alpha_other_0", "alpha_other_1"))







