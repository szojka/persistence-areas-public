
#------------------------------------------------------------
## DESCRIPTION: Source code for running brms models
#------------------------------------------------------------

########################
# Load Packages 
########################

library(readr)
library(tidyverse)
library(truncnorm)
library(stringr)
library(car)
library(lme4)
library(glmmTMB)
library(visreg)
library(Rmisc)
library(patchwork)
library(ggeffects)
library(DHARMa)
library(gt)
library(emmeans)
library(mclogit)
library(rstan)
library(brms)

################################
# Set confidence intervals:

# for 95%
lwr <- 0.025
upr <- 0.975

# for finding estimates of Betas
odds_to_prob <- function(t) {1/(1+exp(-(t)))}

###########################
# # join tags to photos using:
tag_labels <- read_csv("Data/Tag labels.csv")

# Should retain order such that unique key (names) is maintained for each group

# make a key to connect 'names' to 'seeds_clean' (Step 2 & 3)
nkey <- dplyr::select(tag_labels, -plot)
nkey$names1<-paste(nkey$site, nkey$grid, sep="")
num <- rep(c(1:25), each = 8)
nkey$names <- paste(num, nkey$names1, sep="_")
nkey <- dplyr::select(nkey, -names1)

#---------------------------------------------
# Read in data with stan predictions
load(here::here("Supp Bayesian models/Stan rdata/pred_data_saved_95CI.rdata"))
plotlev <- pred_data_saved
  
# create y/n
plotlev$occurrence[plotlev$occurrence > 0] <- c("yes")
plotlev$occurrence[plotlev$occurrence == 0] <- c("no")

######################################################
# This changes depending on how we decide persistence:
######################################################

plotlev$persistence <- NA
plotlev$persistence[plotlev$persist_final >= 1] <- c("yes")
plotlev$persistence[plotlev$persist_final < 1] <- c("no")
plotlev$persistence[is.na(plotlev$persist_final)] <- c("no")

# categories ME, SS_n, SS_y, DL
plotlev$contingency <- with(plotlev, paste0(occurrence, persistence))

plotlev$contingency[plotlev$contingency == "yesyes"] <- c('SS_y')
plotlev$contingency[plotlev$contingency == "nono"] <- c('SS_n')
plotlev$contingency[plotlev$contingency == "yesno"] <- c('ME')
plotlev$contingency[plotlev$contingency == "noyes"] <- c("DL")

# need to add treatment back in
plotlev <- plotlev %>%
  dplyr::rename('treatment' = 'biotic_treatment')

#-------------------------------------------------------.
#### Repeat dataset for grid and site levels ####
#-------------------------------------------------------.

# instead of using seeds to sum, use predicted_fecundity predictions:
head(plotlev)
names(plotlev)

#### gridlev ####
# code datasheet for contingencies at grid level
gridlev1 <- plotlev %>%
  dplyr::select(-persistence, -contingency) %>%
  dplyr::group_by(treatment, species, grid) %>%
  dplyr::mutate(grid_n = n()) %>%
  dplyr::mutate(grid_occur = ifelse('yes' %in% occurrence, 'yes', 'no')) %>% 
  dplyr::ungroup() %>%
  dplyr::select(treatment, species, grid, grid_occur, grid_n, site) %>%
  distinct()
gridlev1 <- as.data.frame(gridlev1)

# join with the grid_pred data that has persistence:
load(here::here("Supp Bayesian models/Stan rdata/grid_pred_95CI.rdata"))

# need to add treatment back in
grid_pred <- grid_pred %>%
  dplyr::rename('treatment' = 'biotic_treatment')

gridlev <- left_join(grid_pred, gridlev1, by = c("grid", "treatment", "species"))

# check if absence at grid level is valid or NA problem:
gridlev %>%
  filter(grid_occur %in% 'no') 
# 9 entries
# BROHOR absent from all plots (A or B) in grids 27, 2, 16; VULMIC absent from all plots in grid 27, missing MICCAL from 28.

# allow if occurrence is yes in A then yes in B, and if yes in B then yes in A. Necessary bc NA in grid 5 B for BROHOR translated to absence even though in grid 5 A BROHOR was present.
gridlev$grid_occur[gridlev$grid == 5 & gridlev$treatment %in% "B"] <- 'yes'

gridlev <- gridlev %>%
  dplyr::rename('occurrence' = 'grid_occur') %>%
  dplyr::rename('grid_seed' = 'grid_predicted_fecundity') %>%
  dplyr::rename('grid_seed_lwr' = 'grid_predicted_fecundity_lwr')%>%
  dplyr::rename('grid_seed_upr' = 'grid_predicted_fecundity_upr') # for compatibility downstream 

# Persistence threshold is one seeds per individual, where seeds out >= seeds in. Use threshold that takes into account missing data

gridlev$persistence <- NA
# THRESHOLD AT LOWER CI
gridlev$persistence[gridlev$grid_seed_lwr >= gridlev$grid_n] <- c("yes") # threshold at grid_n
gridlev$persistence[gridlev$grid_seed_lwr < gridlev$grid_n] <- c("no")
# THRESHOLD AT MEDIAN:
# gridlev$persistence[gridlev$grid_seed >= gridlev$grid_n] <- c("yes") # threshold at grid_n
# gridlev$persistence[gridlev$grid_seed < gridlev$grid_n] <- c("no")

gridlev$contingency <- with(gridlev, paste0(occurrence, persistence))
gridlev$contingency[gridlev$contingency == "yesyes"] <- c('SS_y')
gridlev$contingency[gridlev$contingency == "nono"] <- c('SS_n')
gridlev$contingency[gridlev$contingency == "yesno"] <- c('ME')
gridlev$contingency[gridlev$contingency == "noyes"] <- c("DL")

rm(gridlev1, grid_pred)

#### sitelev ####
# contingency yes/no for occurrence and persistence at site level
sitelev1 <- plotlev %>%
  dplyr::select(-persistence, -contingency)  %>%
  dplyr::group_by(site, treatment, species) %>%
  dplyr::mutate(site_n = n()) %>%
 # dplyr::mutate(site_seed = sum(predicted_fecundity)) %>%
  dplyr::mutate(site_occur = ifelse('yes' %in% occurrence, 'yes', 'no')) %>% 
  dplyr::ungroup() %>%
  dplyr::select(treatment, species, site, site_occur, site_n) %>%
  distinct()
sitelev1 <- as.data.frame(sitelev1)

# join with the grid_pred data that has persistence:
load(here::here("Supp Bayesian models/Stan rdata/site_pred_95CI.rdata"))

# need to add treatment back in
site_pred <- site_pred %>%
  dplyr::rename('treatment' = 'biotic_treatment')

sitelev <- left_join(site_pred, sitelev1, by = c("site", "treatment", "species"))

# check site level occurrence is 100%
sitelev %>%
  filter(site_occur %in% 'no') # empty, good.

sitelev <- sitelev %>%
  dplyr::rename("occurrence" = "site_occur")  %>%
  dplyr::rename('site_seed' = 'site_predicted_fecundity') %>%
  dplyr::rename('site_seed_lwr' = 'site_predicted_fecundity_lwr')%>%
  dplyr::rename('site_seed_upr' = 'site_predicted_fecundity_upr')# for compatibility downstream 

sitelev$persistence <- NA
# THRESHOLD AT LOWER CI:
sitelev$persistence[sitelev$site_seed_lwr >= sitelev$site_n] <- c("yes")
sitelev$persistence[sitelev$site_seed_lwr < sitelev$site_n] <- c("no")
# THRESHOLD AT MEDIAN:
# sitelev$persistence[sitelev$site_seed >= sitelev$site_n] <- c("yes")
# sitelev$persistence[sitelev$site_seed < sitelev$site_n] <- c("no")

sitelev$contingency <- with(sitelev, paste0(occurrence, persistence))
sitelev$contingency[sitelev$contingency == "yesyes"] <- c('SS_y')
sitelev$contingency[sitelev$contingency == "nono"] <- c('SS_n')
sitelev$contingency[sitelev$contingency == "yesno"] <- c('ME')
sitelev$contingency[sitelev$contingency == "noyes"] <- c("DL")

rm(sitelev1, site_pred)

#---------------------------------------------------------------------------------.
#### Data for creating SAR and PAR matrices and fitting models: SAR dt_o dt_p ####
#---------------------------------------------------------------------------------.

# Note, in 'plotlev' some ab-cat = 0, occurrence = yes because A was cleared so paired B plot is used to assess occurrence

length(plotlev$tube)
length(unique(plotlev$tube)) # similar to tag

plotlev$block # similar to names

dt2 <- plotlev %>%
  dplyr::select(tube, block, species, treatment, grid, site, predicted_fecundity, persist_final, intraspecific_abundance, occurrence, persistence)
dt1 <- dt2
  
#------------------------------------------
# Making sure variables are correct type
#------------------------------------------

plotlev$persistence <- as.factor(plotlev$persistence)
plotlev$species <- as.factor(plotlev$species)
plotlev$treatment <- as.factor(plotlev$treatment)
plotlev$grid <- as.factor(plotlev$grid)
plotlev$site <- as.factor(plotlev$site)
plotlev$occurrence <- as.factor(plotlev$occurrence)
plotlev$contingency <- as.factor(plotlev$contingency)
plotlev$block <- as.factor(plotlev$block) # replaced names

gridlev$persistence <- as.factor(gridlev$persistence)
gridlev$species <- as.factor(gridlev$species)
gridlev$treatment <- as.factor(gridlev$treatment)
gridlev$grid <- as.factor(gridlev$grid)
gridlev$site <- as.factor(gridlev$site)
gridlev$occurrence <- as.factor(gridlev$occurrence)
gridlev$contingency <- as.factor(gridlev$contingency)

sitelev$persistence <- as.factor(sitelev$persistence)
sitelev$species <- as.factor(sitelev$species)
sitelev$treatment <- as.factor(sitelev$treatment)
sitelev$site <- as.factor(sitelev$site)
sitelev$occurrence <- as.factor(sitelev$occurrence)
sitelev$contingency <- as.factor(sitelev$contingency)
