#-----------------------------------------------------------------------------.
# DESCRIPTION: Run stan model for Bayes version of Figure 3 at site level (potential biases corrected)
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

#source("Persistence thresholds/Supp Bayesian models/Scripts - generate stan figures/95CI/Source - MAIN fitnessdata POST STAN.R")

#--------------------------------------------------------------------------.
# REGIONAL SCALE (SITES)
#--------------------------------------------------------------------------.

########################
# Data cleaning
########################

fig2F_dat <- sitelev %>%
  select(site, contingency, treatment, species) %>%
  distinct() %>%
  dplyr::group_by(site, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(site, treatment, contingency) %>%
  dplyr::mutate(species_no = n()) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = species_no)
# all NAs = 0 because I standardized by replicate
fig2F_dat[is.na(fig2F_dat)] <- 0
fig2F_dat$SS_n <- rep(0, length.out = length(fig2F_dat$ME))
fig2F_dat$DL <- rep(0, length.out = length(fig2F_dat$ME))
fig2F_dat <- pivot_longer(fig2F_dat, cols = 4:7, names_to = "contingency", values_to = "species_no")
fig2F_dat <- fig2F_dat %>% dplyr::mutate(prop = species_no/replicates) 

fig2F_dat$contingency<-as.factor(fig2F_dat$contingency)
fig2F_dat$site<-as.factor(fig2F_dat$site)
fig2F_dat$treatment<-as.factor(fig2F_dat$treatment)

########################
## Fit models 
########################

# all contingencies at once

# zeros are causing issues in model:
fig2F_dat_nozero <- fig2F_dat %>% filter(!species_no == 0)
#this is ok to do bc zeros are only in SS_n and DL

# BRMS:
fig2F_dat_nozero <- fig2F_dat_nozero %>%
 # mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig2F_dat_nozero$prop); max(fig2F_dat_nozero$prop)

 m.s.pts_brm <- brm(
   bf(prop ~ treatment*contingency + (1 | site),
      phi ~ 1),
   data = fig2F_dat_nozero,
   family = Beta(),
   chains = 4, cores = 4, iter = 2000
 )

 # save model
 save(m.s.pts_brm, file = here::here("Persistence thresholds/Supp Bayesian models/Stan rdata/brms 95CI/m.s.pts_brm.rdata"))
