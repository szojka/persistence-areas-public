
#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 1d Bayes version (predictions from full Beverton-Holt)
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

source(here::here("Supp Bayesian models/Scripts - generate stan figures/FullBevHolt/95CI/Source - MAIN fitnessdata POST STAN.R")) 
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
  pivot_wider(.,names_from = contingency, values_from = species_no)
fig1B_dat[is.na(fig1B_dat)] <- 0
fig1B_dat <- pivot_longer(fig1B_dat, cols = 5:8, names_to = "contingency", values_to = "species_no")

# add colors for plotting into dataframe to make it easy in ggplot
fig1B_dat$colortreat <- NA
fig1B_dat$colortreat[fig1B_dat$contingency == "DL"]<- "slategrey"
fig1B_dat$colortreat[fig1B_dat$contingency == "ME"]<- "violet"
fig1B_dat$colortreat[fig1B_dat$contingency == "SS_n"]<- "#8c96c6"
fig1B_dat$colortreat[fig1B_dat$contingency == "SS_y"]<- "#88419d"

# variable types
fig1B_dat$block <- as.factor(fig1B_dat$block )
fig1B_dat$contingency <- as.factor(fig1B_dat$contingency)
fig1B_dat$colortreat <- as.factor(fig1B_dat$colortreat)

mod_dat <- fig1B_dat
mod_dat

mod_dat$contingency<-as.factor(mod_dat$contingency)
mod_dat$block<-as.factor(mod_dat$block)
mod_dat$site<-as.factor(mod_dat$site)
mod_dat$grid<-as.factor(mod_dat$grid)
mod_dat$prop <- mod_dat$species_no / mod_dat$replicates

# data frame I need:
head(mod_dat) 
mod_dat <- mod_dat %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(mod_dat$prop); max(mod_dat$prop)

#--------------------------------------------------------------------------.
## Load a model 
#--------------------------------------------------------------------------.

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.nat_brm.rdata"))

# check model fit:
summary(m.nat_brm)
# plot(m.nat_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.nat_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_nat <- obj$contingency # data housed in function
output_nat <- output_nat %>%
  select(contingency, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__')

# add pairwise comparison table
em <- emmeans(m.nat_brm, ~contingency, type = "response") 

p.contrast <- as.data.frame(pairs(em))
pairs(em)
