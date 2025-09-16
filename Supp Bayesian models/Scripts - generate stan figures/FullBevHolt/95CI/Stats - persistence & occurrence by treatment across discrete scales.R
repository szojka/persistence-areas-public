
#------------------------------------------------------------------------------
# DESCRIPTION: Figure 2 stats, Bayes version (predictions from full Beverton-Holt)
#------------------------------------------------------------------------------

##############################
# Load packages and functions
##############################

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

source(here::here("Supp Bayesian models/Scripts - generate stan figures/FullBevHolt/95CI/Source - MAIN fitnessdata POST STAN.R")) 
#-------------------------------
# Block w/ & w/out neighbors
#-------------------------------

########################
# Data wrangling
########################
{
  # Plot level persistence
  
  fig2A_per <- plotlev %>%
    select(block, persistence, treatment, species, grid, site) %>%
    distinct() %>%
    dplyr::group_by(block, treatment) %>% # find number of species replicates in each plot
    dplyr::mutate(replicates = n()) %>%
    group_by(block, treatment, persistence) %>%
    dplyr::mutate(species_no = n()) %>%
    dplyr::mutate(prop = species_no/replicates) %>%
    select(-species) %>%
    distinct() %>%
    ungroup()
  
  fig2A_per_yes <- fig2A_per %>% 
    filter(persistence %in% "yes")%>%
    select(-persistence) # filter for yes, as 1 - proportion persisting = proportion not persisting
  fig2A_per_no <- fig2A_per %>% 
    filter(persistence %in% "no") %>% 
    filter(prop == 1) %>%
    select(-persistence)
  fig2A_per_no$prop[fig2A_per_no$prop == 1] <- 0    
  fig2A_per <- rbind(fig2A_per_yes,fig2A_per_no)  

  # Plot level occurrence
  
  fig2A_occ <- plotlev %>%
    filter(treatment %in% 'B') %>%
    select(block, occurrence, species, grid, site) %>%
    distinct() %>%
    dplyr::group_by(block) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(block, occurrence) %>%
    dplyr::mutate(species_no = n()) %>%
    dplyr::mutate(prop = species_no/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  
  fig2A_occ_yes <-fig2A_occ %>% # now filter for 'yes' 
    filter(occurrence %in% 'yes') %>%
    mutate(treatment = c("O")) %>%
    select(-occurrence)
  fig2A_occ_no <- fig2A_occ %>% 
    filter(occurrence %in% "no") %>% 
    mutate(treatment = c("O")) %>%
    filter(prop == 1) %>%
    select(-occurrence)
  fig2A_occ_no$prop[fig2A_occ_no$prop == 1] <- 0    
  fig2A_occ <- rbind(fig2A_occ_yes,fig2A_occ_no)  

  # put persistence and occurrence raw data together for figure and models
  fig_fig_2A <- rbind(fig2A_occ,fig2A_per)
  
  # correct variable types
  fig_fig_2A$block<- as.factor(fig_fig_2A$block)
  fig_fig_2A$treatment<- as.factor(fig_fig_2A$treatment)
  fig_fig_2A$site<- as.factor(fig_fig_2A$site)
  fig_fig_2A$grid<- as.factor(fig_fig_2A$grid)
  
  fig_fig_2A$treatment <- case_match(fig_fig_2A$treatment, 
                                     "A" ~ "suitability without neighbors", 
                                     "B" ~ "suitability with neighbors", 
                                     "O" ~ "occupancy", 
                                     .default = fig_fig_2A$treatment)
  fig_fig_2A$treatment <- as.factor(fig_fig_2A$treatment)
# 0s to 0.01, 1s to 0.99
fig_fig_2A <- fig_fig_2A %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig_fig_2A$prop); max(fig_fig_2A$prop)
}

########################
# Load models
########################

load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.plot_brm.rdata"))

# check model fit:
summary(m.plot_brm) 
# plot(m.plot_brm) # shows posteriors and chains

############################
# Set up viz and table data
############################

# Extract estimates:
obj <- conditional_effects(m.plot_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
vis.p <- obj$treatment # data housed in function
vis.p <- vis.p %>%
  select(treatment, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  mutate(scale = "plot")

# add pairwise comparison table, response should mirror predictions above ^
em <- emmeans(m.plot_brm, ~treatment, type = "response")

mplot.contrast <- as.data.frame(pairs(em))
pairs(em)

#-------------------------------
# Grid w/ & w/out neighbors
#-------------------------------

#################
# Data wrangling
#################

{
  # Grid level persistence
  
  fig2C_per <- gridlev %>%
    select(site, grid, persistence, treatment, species) %>%
    distinct() %>%
    dplyr::group_by(grid, treatment) %>% # find number of species replicates in each plot
    dplyr::mutate(replicates = n()) %>%
    group_by(grid, treatment, persistence) %>%
    dplyr::mutate(species_no = n()) %>%
    dplyr::mutate(prop = species_no/replicates) %>%
    select(-species) %>%
    distinct() %>%
    ungroup()
  
  fig2C_per_yes <- fig2C_per %>% 
    filter(persistence %in% "yes")%>%
    select(-persistence) # filter for yes, as 1 - proportion persisting = proportion not persisting
  fig2C_per_no <- fig2C_per %>% 
    filter(persistence %in% "no") %>% 
    filter(prop == 1) %>%
    select(-persistence)
  fig2C_per_no$prop[fig2C_per_no$prop == 1] <- 0    
  fig2C_per <- rbind(fig2C_per_yes,fig2C_per_no)  
  
  # Grid level occurrence
  
  fig2C_occ <- gridlev %>%
    filter(treatment %in% 'B') %>%
    select(occurrence, species, grid, site) %>%
    distinct() %>%
    dplyr::group_by(grid) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(grid, occurrence) %>%
    dplyr::mutate(species_no = n()) %>%
    dplyr::mutate(prop = species_no/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  
  fig2C_occ_yes <-fig2C_occ %>% # now filter for 'yes' 
    filter(occurrence %in% 'yes') %>%
    mutate(treatment = c("O")) %>%
    select(-occurrence)
  fig2C_occ_no <- fig2C_occ %>% 
    filter(occurrence %in% "no") %>% # idea is to keep zeros, so find where occurrence no = 1, and make the prop zero, then rbind to persistence yes dataframe
    mutate(treatment = c("O")) %>%
    filter(prop == 1) %>%
    select(-occurrence)
  fig2C_occ_no$prop[fig2C_occ_no$prop == 1] <- 0    
  fig2C_occ <- rbind(fig2C_occ_yes,fig2C_occ_no)  
 
  # together for figure and model
  fig_fig_2C <- rbind(fig2C_occ,fig2C_per)
  
  # correct variable types
  fig_fig_2C$treatment<- as.factor(fig_fig_2C$treatment)
  fig_fig_2C$site<- as.factor(fig_fig_2C$site)
  fig_fig_2C$grid<- as.factor(fig_fig_2C$grid)
  
  fig_fig_2C$treatment <- case_match(fig_fig_2C$treatment, 
                                     "A" ~ "suitability without neighbors", 
                                     "B" ~ "suitability with neighbors", 
                                     "O" ~ "occupancy", 
                                     .default = fig_fig_2C$treatment)
  fig_fig_2C$treatment <- as.factor(fig_fig_2C$treatment)
fig_fig_2C <- fig_fig_2C %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig_fig_2C$prop); max(fig_fig_2C$prop)
}

########################
# Load models
########################

load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.grid_brm.rdata"))

# check model fit:
summary(m.grid_brm) 
# plot(m.grid_brm) # shows posteriors and chains

############################
# Set up viz and table data
############################

# Extract estimates:
obj <- conditional_effects(m.grid_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
vis.g <- obj$treatment # data housed in function
vis.g <- vis.g %>%
  select(treatment, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  mutate(scale = "grid")

# add pairwise comparison table, response should mirror predictions above ^
em <- emmeans(m.grid_brm, ~treatment, type = "response") # specify green_index_scaled values

mgrid.contrast <- as.data.frame(pairs(em))
pairs(em)

#-------------------------------
# Site w/ & w/out neighbors
#-------------------------------

#################
# Data wrangling
#################

# Site level persistence
{
  fig2E_per <- sitelev %>%
    select(site,persistence, treatment, species) %>%
    distinct() %>%
    dplyr::group_by(site, treatment) %>% # find number of species replicates in each plot
    dplyr::mutate(replicates = n()) %>%
    group_by(site, treatment, persistence) %>%
    dplyr::mutate(species_no = n()) %>%
    dplyr::mutate(prop = species_no/replicates) %>%
    select(-species) %>%
    distinct() %>%
    ungroup()
  
  fig2E_per_yes <- fig2E_per %>% 
    filter(persistence %in% "yes")%>%
    select(-persistence) # filter for yes, as 1 - proportion persisting = proportion not persisting
  fig2E_per_no <- fig2E_per %>% 
    filter(persistence %in% "no") %>% 
    filter(prop == 1) %>%
    select(-persistence)
  fig2E_per_no$prop[fig2E_per_no$prop == 1] <- 0    
  fig2E_per <- rbind(fig2E_per_yes,fig2E_per_no)  
  
  # Site level occurrence
  
  fig2E_occ <- sitelev %>%
    filter(treatment %in% 'B') %>%
    select(site, occurrence, species) %>%
    distinct() %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(site, occurrence) %>%
    dplyr::mutate(species_no = n()) %>%
    dplyr::mutate(prop = species_no/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  
  fig2E_occ_yes <-fig2E_occ %>% # now filter for 'yes' 
    filter(occurrence %in% 'yes') %>%
    mutate(treatment = c("O")) %>%
    select(-occurrence)
  fig2E_occ_no <- fig2E_occ %>% 
    filter(occurrence %in% "no") %>% 
    mutate(treatment = c("O")) %>%
    filter(prop == 1) %>%
    select(-occurrence)
  fig2E_occ_no$prop[fig2E_occ_no$prop == 1] <- 0    
  fig2E_occ <- rbind(fig2E_occ_yes,fig2E_occ_no)
 
  # pull together for plotting and model
  fig_fig_2E <- rbind(fig2E_occ,fig2E_per)
  
  # correct variable types
  fig_fig_2E$treatment<- as.factor(fig_fig_2E$treatment)
  fig_fig_2E$site<- as.factor(fig_fig_2E$site)
  
  fig_fig_2E$treatment <- case_match(fig_fig_2E$treatment, 
                                     "A" ~ "suitability without neighbors", 
                                     "B" ~ "suitability with neighbors", 
                                     "O" ~ "occupancy", 
                                     .default = fig_fig_2E$treatment)
  fig_fig_2E$treatment <- as.factor(fig_fig_2E$treatment)
fig_fig_2E <- fig_fig_2E %>%
  mutate(prop = ifelse(prop == 0, 0.01, prop)) %>%
  mutate(prop = ifelse(prop == 1, 0.99, prop))
min(fig_fig_2E$prop); max(fig_fig_2E$prop)
}

########################
# Load models
########################

load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.site_brm.rdata"))

# check model fit:
summary(m.site_brm) # bulk ess > 2000  and Rhat = 1.00
# plot(m.site_brm) # shows posteriors and chains

############################
# Set up viz and table data
############################

obj <- conditional_effects(m.site_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
vis.s <- obj$treatment # data housed in function
vis.s <- vis.s %>%
  select(treatment, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  mutate(scale = "site")

# add pairwise comparison table, response should mirror predictions above ^
em <- emmeans(m.site_brm, ~treatment, type = "response")  # specify green_index_scaled values

msite.contrast <- as.data.frame(pairs(em))
pairs(em)

###############################################
# used for figure and tables
###############################################

summary <- rbind(vis.p,vis.g,vis.s)
names(summary)
summary <- as.data.frame(summary)
