#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 3 site level Bayes version
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

# zeros are causing issues in model:
fig2F_dat_nozero <- fig2F_dat %>% filter(!species_no == 0)
#this is ok to do bc zeros are only in SS_n and DL

# BRMS:
fig2F_dat_nozero <- fig2F_dat_nozero %>%
  # mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig2F_dat_nozero$prop); max(fig2F_dat_nozero$prop)

########################
## Load models 
########################

load(here::here("Supp Bayesian models/Stan rdata/brms 95CI/m.s.pts_brm.rdata"))

# check model fit:
summary(m.s.pts_brm)
# plot(m.s.pts_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.s.pts_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_s <- obj$`treatment:contingency` # data housed in function
output_s <- output_s %>%
  select(treatment, contingency, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'site')


# DL and SS_n are both zero, fill that in:
zero_dat <- data.frame(treatment = rep(c('A','B'), times = 2), 
                       contingency = rep(c('DL', 'SS_n'), each = 2), 
                       prediction = rep(0, times = 4), 
                       prediction_lwr = rep(0, times = 4), 
                       prediction_upr = rep(0, times = 4),
                       scale = 'site')

output_s <- rbind(output_s,zero_dat)

# add pairwise comparison table
em <- emmeans(m.s.pts_brm, ~treatment | contingency, type = "response") 

s.contrast <- as.data.frame(pairs(em))
pairs(em)

############################
# Create objects for table

# tab2. Pairs significance

site.contrast <- as.data.frame(s.contrast)

site.contrast$scale <- NA
site.contrast$scale <- 'site'

