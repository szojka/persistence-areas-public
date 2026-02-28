
#-----------------------------------------------------------------------------.
# DESCRIPTION: Run stan model for Bayes version of Figure 3 at blocks (potential biases corrected)
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
# library(emmeans)

source("Supp Bayesian models/Scripts - run stan models/Run models for Bayes figures/95CI/Source - MAIN fitnessdata POST STAN.R")

#--------------------------------------------------------------------------.
# POPULATION SCALE (PLOT)
#--------------------------------------------------------------------------.

########################
# Data cleaning
########################

fig2B_dat <- plotlev %>%
  select(block, grid, site, contingency, treatment, species) %>%
  distinct() %>%
  dplyr::group_by(block, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(block, treatment, contingency) %>%
  dplyr::mutate(species_no = n()) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
# I do need zeros, so expand, fill in NA with 0
  pivot_wider(.,names_from = contingency, values_from = species_no) 
  # all NAs = 0 because I standardized by replicate
fig2B_dat[is.na(fig2B_dat)] <- 0
fig2B_dat <- pivot_longer(fig2B_dat, cols = 6:9, names_to = "contingency", values_to = "species_no")
# now create proportion as I have species_no as zero
fig2B_dat <- fig2B_dat %>%
  dplyr::mutate(prop = species_no/replicates)

fig2B_dat$contingency<-as.factor(fig2B_dat$contingency)
fig2B_dat$site<-as.factor(fig2B_dat$site)
fig2B_dat$grid<-as.factor(fig2B_dat$grid)
fig2B_dat$block<-as.factor(fig2B_dat$block)
fig2B_dat$treatment<-as.factor(fig2B_dat$treatment)

########################
## Fit models 
########################

# BRMS:
fig2B_dat <- fig2B_dat %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig2B_dat$prop); max(fig2B_dat$prop)

m.b.pts_brm <- brm(
   bf(prop ~ treatment*contingency  + (1 | block) + (1 | site:grid) + (1 | site),
      phi ~ 1),
   data = fig2B_dat,
   family = Beta(),
   chains = 4, cores = 4, iter = 2000
 )

 # save model
 save(m.b.pts_brm, file = here::here("Supp Bayesian models/Stan rdata/brms 95CI/m.b.pts_brm.rdata"))
