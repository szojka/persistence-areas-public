
#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 3 grid level Bayes version
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

#source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Source - MAIN fitnessdata POST STAN.R"))

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

fig2D_dat <- fig2D_dat %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig2D_dat$prop); max(fig2D_dat$prop)

########################
## Load models 
########################

load(here::here("Supp Bayesian models/Stan rdata/brms 95CI/m.g.pts_brm.rdata"))

# check model fit:
summary(m.g.pts_brm)
# plot(m.g.pts_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.g.pts_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_g <- obj$`treatment:contingency` # data housed in function
output_g <- output_g %>%
  select(treatment, contingency, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'grid')

# add pairwise comparison table
em <- emmeans(m.g.pts_brm, ~treatment | contingency, type = "response") 

g.contrast <- as.data.frame(pairs(em))

############################
# Create objects for table

# tab2. Pairs significance

grid.contrast <- as.data.frame(g.contrast)

grid.contrast$scale <- NA
grid.contrast$scale <- 'grid'

