
#-----------------------------------------------------------------------------.
# DESCRIPTION: Run stan model for Bayes version of Figure 1d (full Beverton-Holt)
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

source(here::here("Supp Bayesian models/Scripts - run stan models/Run models for Bayes figures/FullBevHolt_95CI/Source - MAIN fitnessdata POST STAN.R"))

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

#--------------------------------------------------------------------------.
## Fit a model 
#--------------------------------------------------------------------------.

#accounting for potential correlation among observations with random effects

# BRMS:
mod_dat <- mod_dat %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(mod_dat$prop); max(mod_dat$prop)

m.nat_brm <- brm(
  bf(prop ~ contingency + (1 | block) + (1 | site:grid) + (1 | site),
     phi ~ 1),
  data = mod_dat,
  family = Beta(),
  chains = 4, cores = 4, iter = 2000
)
 
#save model
save(m.nat_brm, file = here::here("Supp Bayesian models/Stan rdata/FullBevHolt/95CI/m.nat_brm.rdata"))
