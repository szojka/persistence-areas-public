
#-----------------------------------------------------------------------------.
# DESCRIPTION: Run stan model for Bayes version of Figure 3 at grids (potential biases corrected)
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

#source("Scripts - stan persistence/95CI/Source - MAIN fitnessdata POST STAN.R")

#--------------------------------------------------------------------------.
# META-POP SCALE (GRIDS)
#--------------------------------------------------------------------------.

########################
# Data cleaning
########################

fig2D_dat <- gridlev %>%
  select(grid, site, contingency, treatment, species) %>%
  distinct() %>%
  dplyr::group_by(grid, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(grid, treatment, contingency) %>%
  dplyr::mutate(species_no = n()) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = species_no)
# all NAs = 0 because I standardized by replicate
fig2D_dat[is.na(fig2D_dat)] <- 0
fig2D_dat <- pivot_longer(fig2D_dat, cols = 5:8, names_to = "contingency", values_to = "species_no")
fig2D_dat <- fig2D_dat %>% dplyr::mutate(prop = species_no/replicates) 

fig2D_dat$contingency<-as.factor(fig2D_dat$contingency)
fig2D_dat$site<-as.factor(fig2D_dat$site)
fig2D_dat$grid<-as.factor(fig2D_dat$grid)
fig2D_dat$treatment<-as.factor(fig2D_dat$treatment)

########################
## Fit models 
########################

# all contingencies at once

# BRMS:
fig2D_dat <- fig2D_dat %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig2D_dat$prop); max(fig2D_dat$prop)

m.g.pts_brm <- brm(
   bf(prop ~ treatment*contingency + (1 | site:grid) + (1 | site),
      phi ~ 1),
   data = fig2D_dat,
   family = Beta(),
   chains = 4, cores = 4, iter = 2000
 )

 # save model
 save(m.g.pts_brm, file = here::here("Supp Bayesian models/Stan rdata/brms 95CI/m.g.pts_brm.rdata"))
