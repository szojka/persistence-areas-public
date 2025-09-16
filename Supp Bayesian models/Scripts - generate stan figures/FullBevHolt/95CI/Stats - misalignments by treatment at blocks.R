
#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 3 blocks stats, Bayes version (predictions from full Beverton-Holt)
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

source(here::here("Supp Bayesian models/Scripts - generate stan figures/FullBevHolt/95CI/Source - MAIN fitnessdata POST STAN.R")) 

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

fig2B_dat <- fig2B_dat %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig2B_dat$prop); max(fig2B_dat$prop)

########################
## Load models 
########################

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.b.pts_brm.rdata"))

# check model fit:
summary(m.b.pts_brm)
# plot(m.b.pts_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.b.pts_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_b <- obj$`treatment:contingency` # data housed in function
output_b <- output_b %>%
  select(treatment, contingency, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'plot')

# add pairwise comparison table
em <- emmeans(m.b.pts_brm, ~treatment | contingency, type = "response") 

p.contrast <- as.data.frame(pairs(em))

############################
# Create objects for table

# tab2. Pairs significance

plot.contrast <- as.data.frame(p.contrast)

plot.contrast$scale <- NA
plot.contrast$scale <- 'plot'


