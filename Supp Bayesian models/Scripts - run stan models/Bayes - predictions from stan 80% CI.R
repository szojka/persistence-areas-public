###################################
# Pulling out fits
###################################

# Predictions correct for potential biases, using 80% credible intervals 

# this script does not assess model fits at all, assuming all is well, it simply works the posteriors into dataframes used for figures and analyses in persistence-areas repo.

##############################################
# Making plot level predictions on fecundity
##############################################

# load libraries
library(rstan)
library(tidyverse)

load(here::here("Supp Bayesian models/Stan rdata/Posteriors/Brohor_posteriors.rdata")) # this model used a poisson likelihood
load(here::here("Supp Bayesian models/Stan rdata/Posteriors/Vulmic_posteriors.rdata")) # this model used a poisson likelihood
load(here::here("Supp Bayesian models/Stan rdata/Posteriors/Plaere_posteriors.rdata")) # this model used a zero-inflated poisson likelihood
load(here::here("Supp Bayesian models/Stan rdata/Posteriors/Miccal_posteriors.rdata")) # this model used a zero-inflated poisson likelihood

# prepare data as I did prior to model fitting
load(here::here("Data/dat_final_doubled.Rdata"))

# code start:
dat_final_doubled <- na.omit(dat_final_doubled)

# case match so that site is a number
dat_final_doubled$site <- case_match(dat_final_doubled$site,
                                     'HUT' ~ 1,
                                     'KNOX' ~ 2,
                                     'LM' ~ 3,
                                     'RES' ~ 4,
                                     'TAIL' ~ 5,
                                     'WHIT' ~ 6)

# case_match so that abundance bins from 0-3 go to mid abundance of each bin:
dat_final_doubled$intraspecific_abundance <- case_match(dat_final_doubled$intraspecific_abundance,
                                                        0 ~ 0,
                                                        1 ~ 5,
                                                        2 ~ 50,
                                                        3 ~ 101) # chose 100 instead of 'half-way point of >101' which isn't clear)

# load focal dataframes and name them so that we can join with species specific predictions
brohor_dat <- as.data.frame(subset(dat_final_doubled,species=="brohor"))
vulmic_dat <- as.data.frame(subset(dat_final_doubled,species=="vulmic"))
plaere_dat <- as.data.frame(subset(dat_final_doubled,species=="plaere"))
miccal_dat <- as.data.frame(subset(dat_final_doubled,species=="miccal"))

#-----------
# BROHOR

# include relevent estimates and credible intervals into each dat
post <- rstan::extract(Brohor_posteriors)

# lambda predictions
brohor_dat$lambda_ei <- apply(post$lambda_ei_out, 2, median) # the two indicates applied over columns
brohor_dat$lambda_ei_lower <- apply(post$lambda_ei_out, 2, quantile, probs = 0.1)
brohor_dat$lambda_ei_upper <- apply(post$lambda_ei_out, 2, quantile, probs = 0.9)

# intra predictions 
brohor_dat$alpha_intra_ei <- apply(post$alpha_intra_ei_out, 2, median)
brohor_dat$alpha_intra_ei_lower <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.1)
brohor_dat$alpha_intra_ei_upper <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.9)

# inter predictions
brohor_dat$alpha_other_ei <- apply(post$alpha_other_ei_out, 2, median)
brohor_dat$alpha_other_ei_lower <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.1)
brohor_dat$alpha_other_ei_upper <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.9)

# fecundity predictions (all accounting for site-level epsilons):

# predictions calculated with all parameters in BH
brohor_dat$y_hat <- apply(post$y_hat, 2, median)
brohor_dat$y_hat_lower <- apply(post$y_hat, 2, quantile, probs = 0.1)
brohor_dat$y_hat_upper <- apply(post$y_hat, 2, quantile, probs = 0.9)

# predictions calculated where intraspecific density accounted for
brohor_dat$y_hat_no_intra <- apply(post$y_hat_no_intra, 2, median)
brohor_dat$y_hat_no_intra_lower <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.1)
brohor_dat$y_hat_no_intra_upper <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.9)

# predictions calculated where so that no interactions; low density growth rate
brohor_dat$y_hat_baseline <- apply(post$y_hat_baseline, 2, median)
brohor_dat$y_hat_baseline_lower <- apply(post$y_hat_baseline, 2, quantile, probs = 0.1)
brohor_dat$y_hat_baseline_upper <- apply(post$y_hat_baseline, 2, quantile, probs = 0.9)

# we now need to save posts for calculating iteration specific probabilities 
# to scale from plot to grid persistence
brohor_post <- post

#------------
# VULMIC

# include relevent estimates and credible intervals into each dat
post <- rstan::extract(Vulmic_posteriors)

# lambda predictions
vulmic_dat$lambda_ei <- apply(post$lambda_ei_out, 2, median) # the two indicates applied over columns
vulmic_dat$lambda_ei_lower <- apply(post$lambda_ei_out, 2, quantile, probs = 0.1)
vulmic_dat$lambda_ei_upper <- apply(post$lambda_ei_out, 2, quantile, probs = 0.9)

# intra predictions 
vulmic_dat$alpha_intra_ei <- apply(post$alpha_intra_ei_out, 2, median)
vulmic_dat$alpha_intra_ei_lower <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.1)
vulmic_dat$alpha_intra_ei_upper <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.9)

# inter predictions
vulmic_dat$alpha_other_ei <- apply(post$alpha_other_ei_out, 2, median)
vulmic_dat$alpha_other_ei_lower <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.1)
vulmic_dat$alpha_other_ei_upper <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.9)

# fecundity predictions (all accounting for site-level epsilons):

# predictions calculated with all parameters in BH
vulmic_dat$y_hat <- apply(post$y_hat, 2, median)
vulmic_dat$y_hat_lower <- apply(post$y_hat, 2, quantile, probs = 0.1)
vulmic_dat$y_hat_upper <- apply(post$y_hat, 2, quantile, probs = 0.9)

# predictions calculated where intraspecific density accounted for
vulmic_dat$y_hat_no_intra <- apply(post$y_hat_no_intra, 2, median)
vulmic_dat$y_hat_no_intra_lower <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.1)
vulmic_dat$y_hat_no_intra_upper <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.9)

# predictions calculated where so that no interactions; low density growth rate
vulmic_dat$y_hat_baseline <- apply(post$y_hat_baseline, 2, median)
vulmic_dat$y_hat_baseline_lower <- apply(post$y_hat_baseline, 2, quantile, probs = 0.1)
vulmic_dat$y_hat_baseline_upper <- apply(post$y_hat_baseline, 2, quantile, probs = 0.9)

# we now need to save posts for calculating iteration specific probabilities 
# to scale from plot to grid persistence
vulmic_post <- post

#------------
# PLAERE

# include relevent estimates and credible intervals into each dat
post <- rstan::extract(Plaere_posteriors)

# lambda predictions
plaere_dat$lambda_ei <- apply(post$lambda_ei_out, 2, median) # the two indicates applied over columns
plaere_dat$lambda_ei_lower <- apply(post$lambda_ei_out, 2, quantile, probs = 0.1)
plaere_dat$lambda_ei_upper <- apply(post$lambda_ei_out, 2, quantile, probs = 0.9)

# intra predictions 
plaere_dat$alpha_intra_ei <- apply(post$alpha_intra_ei_out, 2, median)
plaere_dat$alpha_intra_ei_lower <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.1)
plaere_dat$alpha_intra_ei_upper <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.9)

# inter predictions
plaere_dat$alpha_other_ei <- apply(post$alpha_other_ei_out, 2, median)
plaere_dat$alpha_other_ei_lower <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.1)
plaere_dat$alpha_other_ei_upper <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.9)

# fecundity predictions (all accounting for site-level epsilons):

# predictions calculated with all parameters in BH
plaere_dat$y_hat <- apply(post$y_hat, 2, median)
plaere_dat$y_hat_lower <- apply(post$y_hat, 2, quantile, probs = 0.1)
plaere_dat$y_hat_upper <- apply(post$y_hat, 2, quantile, probs = 0.9)

# predictions calculated where intraspecific density accounted for
plaere_dat$y_hat_no_intra <- apply(post$y_hat_no_intra, 2, median)
plaere_dat$y_hat_no_intra_lower <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.1)
plaere_dat$y_hat_no_intra_upper <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.9)

# predictions calculated where so that no interactions; low density growth rate
plaere_dat$y_hat_baseline <- apply(post$y_hat_baseline, 2, median)
plaere_dat$y_hat_baseline_lower <- apply(post$y_hat_baseline, 2, quantile, probs = 0.1)
plaere_dat$y_hat_baseline_upper <- apply(post$y_hat_baseline, 2, quantile, probs = 0.9)

# we now need to save posts for calculating iteration specific probabilities 
# to scale from plot to grid persistence
plaere_post <- post

#--------------
# MICCAL

# include relevent estimates and credible intervals into each dat
post <- rstan::extract(Miccal_posteriors)

# lambda predictions
miccal_dat$lambda_ei <- apply(post$lambda_ei_out, 2, median) # the two indicates applied over columns
miccal_dat$lambda_ei_lower <- apply(post$lambda_ei_out, 2, quantile, probs = 0.1)
miccal_dat$lambda_ei_upper <- apply(post$lambda_ei_out, 2, quantile, probs = 0.9)

# intra predictions 
miccal_dat$alpha_intra_ei <- apply(post$alpha_intra_ei_out, 2, median)
miccal_dat$alpha_intra_ei_lower <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.1)
miccal_dat$alpha_intra_ei_upper <- apply(post$alpha_intra_ei_out, 2, quantile, probs = 0.9)

# inter predictions
miccal_dat$alpha_other_ei <- apply(post$alpha_other_ei_out, 2, median)
miccal_dat$alpha_other_ei_lower <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.1)
miccal_dat$alpha_other_ei_upper <- apply(post$alpha_other_ei_out, 2, quantile, probs = 0.9)

# fecundity predictions (all accounting for site-level epsilons):

# predictions calculated with all parameters in BH
miccal_dat$y_hat <- apply(post$y_hat, 2, median)
miccal_dat$y_hat_lower <- apply(post$y_hat, 2, quantile, probs = 0.1)
miccal_dat$y_hat_upper <- apply(post$y_hat, 2, quantile, probs = 0.9)

# predictions calculated where intraspecific density accounted for
miccal_dat$y_hat_no_intra <- apply(post$y_hat_no_intra, 2, median)
miccal_dat$y_hat_no_intra_lower <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.1)
miccal_dat$y_hat_no_intra_upper <- apply(post$y_hat_no_intra, 2, quantile, probs = 0.9)

# predictions calculated where so that no interactions; low density growth rate
miccal_dat$y_hat_baseline <- apply(post$y_hat_baseline, 2, median)
miccal_dat$y_hat_baseline_lower <- apply(post$y_hat_baseline, 2, quantile, probs = 0.1)
miccal_dat$y_hat_baseline_upper <- apply(post$y_hat_baseline, 2, quantile, probs = 0.9)

# we now need to save posts for calculating iteration specific probabilities 
# to scale from plot to grid persistence
miccal_post <- post

# ALL SPECIES TOGETHER PLOT SCALE:
pred_data <- rbind(brohor_dat,vulmic_dat,plaere_dat,miccal_dat)

# point estimate determines persistence
# already exists: y
pred_data <- pred_data %>%
  
  # Scenario 1 (y_hat_baseline)
  mutate(persist_1_point = ifelse(pred_data$y_hat_baseline >= 1, 1, 0)) %>%
  mutate(persist_1_cred = ifelse(pred_data$y_hat_baseline_lower >= 1, 1, 0)) %>%# lwr credible interval determines persistence
  
  # Scenario 2 (y_hat_no_intra)
  mutate(persist_2_point = ifelse(pred_data$y_hat_no_intra >= 1, 1, 0)) %>%
  mutate(persist_2_cred = ifelse(pred_data$y_hat_no_intra_lower >= 1, 1, 0)) %>%# lwr credible interval determines persistence
  
  # Scenario 3 (y_hat)
  mutate(persist_3_point = ifelse(pred_data$y_hat >= 1, 1, 0)) %>%
  mutate(persist_3_cred = ifelse(pred_data$y_hat_lower >= 1, 1, 0)) # lwr credible interval determines persistence

# Final predictions I will use that are specific to biotic treatment
pred_data$persist_final <- NA
pred_data$persist_final <- ifelse(pred_data$biotic_treatment %in% 'B', pred_data$persist_2_cred, pred_data$persist_1_cred) # if yes = control intra, if no = 

# now you need the predicted fecundity and lower and upper CI's to be able to plot plot level results:
pred_data$predicted_fecundity <- NA
pred_data$predicted_fecundity_lwr <- NA 
pred_data$predicted_fecundity_upr <- NA 
pred_data$predicted_fecundity <-  ifelse(pred_data$biotic_treatment %in% 'B', pred_data$y_hat_no_intra, pred_data$y_hat_baseline)
pred_data$predicted_fecundity_lwr <-  ifelse(pred_data$biotic_treatment %in% 'B', pred_data$y_hat_no_intra_lower, pred_data$y_hat_baseline_lower)
pred_data$predicted_fecundity_upr <-  ifelse(pred_data$biotic_treatment %in% 'B', pred_data$y_hat_no_intra_upper, pred_data$y_hat_baseline_upper)

# summarize to check that this is doing what I want
pred_data %>%
  select(biotic_treatment, persist_final, persist_1_cred, persist_2_cred) %>%
  head()

# save version with less columns
pred_data_saved <- select(pred_data, tube, block, grid, site, species, biotic_treatment, cover, intraspecific_abundance, replicate_spp_per_plot, green_index_scaled ,occurrence, per_capita_seed_prod, predicted_fecundity, predicted_fecundity_lwr, predicted_fecundity_upr, persist_final)

##############################################
# Making grid level predictions on fecundity
##############################################

# BROHOR

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- brohor_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for brohor_dat
identical(length(brohor_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(brohor_dat, grid, biotic_treatment )
head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_B <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_B <- rbind(grid_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above.

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- brohor_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for brohor_dat
identical(length(brohor_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_A <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_A <- rbind(grid_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and use apply to extract median and quantiles for each row

dat_a <- grid_sum_dat_A %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_baseline = median(grid_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- grid_sum_dat_B %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(grid_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and go ahead and make a single summarized dataframe...

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -grid, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$grid_predicted_fecundity <- NA
dat_ab$grid_predicted_fecundity_lwr <- NA
dat_ab$grid_predicted_fecundity_upr <- NA

dat_ab$grid_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$grid_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$grid_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_brohor <- select(dat_ab, grid, biotic_treatment, grid_predicted_fecundity, grid_predicted_fecundity_lwr, grid_predicted_fecundity_upr)

# VULMIC

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- vulmic_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for vulmic_dat
identical(length(vulmic_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(vulmic_dat, grid, biotic_treatment )
head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_B <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_B <- rbind(grid_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above.

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- vulmic_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for vulmic_dat
identical(length(vulmic_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_A <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_A <- rbind(grid_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and use apply to extract median and quantiles for each row

dat_a <- grid_sum_dat_A %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_baseline = median(grid_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- grid_sum_dat_B %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(grid_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and go ahead and make a single summarized dataframe...

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -grid, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$grid_predicted_fecundity <- NA
dat_ab$grid_predicted_fecundity_lwr <- NA
dat_ab$grid_predicted_fecundity_upr <- NA

dat_ab$grid_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$grid_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$grid_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_vulmic <- select(dat_ab, grid, biotic_treatment, grid_predicted_fecundity, grid_predicted_fecundity_lwr, grid_predicted_fecundity_upr)


# PLAERE

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- plaere_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for plaere_dat
identical(length(plaere_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(plaere_dat, grid, biotic_treatment )
head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_B <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_B <- rbind(grid_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above. 

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- plaere_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for plaere_dat
identical(length(plaere_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_A <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_A <- rbind(grid_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and use apply to extract median and quantiles for each row

dat_a <- grid_sum_dat_A %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_baseline = median(grid_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- grid_sum_dat_B %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(grid_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and go ahead and make a single summarized dataframe...

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -grid, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$grid_predicted_fecundity <- NA
dat_ab$grid_predicted_fecundity_lwr <- NA
dat_ab$grid_predicted_fecundity_upr <- NA

dat_ab$grid_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$grid_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$grid_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_plaere <- select(dat_ab, grid, biotic_treatment, grid_predicted_fecundity, grid_predicted_fecundity_lwr, grid_predicted_fecundity_upr)

# MICCAL

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- miccal_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for miccal_dat
identical(length(miccal_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(miccal_dat, grid, biotic_treatment )
head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_B <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_B <- rbind(grid_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above.

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- miccal_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for miccal_dat
identical(length(miccal_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
grid_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
grid_sum_dat_A <- data.frame()
column_names <- colnames(grid_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- grid_dat %>%
    group_by(grid, biotic_treatment) %>%
    reframe(
      grid_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  grid_sum_dat_A <- rbind(grid_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and use apply to extract median and quantiles for each row

dat_a <- grid_sum_dat_A %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_baseline = median(grid_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- grid_sum_dat_B %>%
  group_by(grid, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(grid_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(grid_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(grid_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(grid, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and go ahead and make a single summarized dataframe...

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -grid, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$grid_predicted_fecundity <- NA
dat_ab$grid_predicted_fecundity_lwr <- NA
dat_ab$grid_predicted_fecundity_upr <- NA

dat_ab$grid_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$grid_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$grid_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_miccal <- select(dat_ab, grid, biotic_treatment, grid_predicted_fecundity, grid_predicted_fecundity_lwr, grid_predicted_fecundity_upr)

dat_ab_brohor$species <- 'brohor'
dat_ab_vulmic$species <- 'vulmic'
dat_ab_plaere$species <- 'plaere'
dat_ab_miccal$species <- 'miccal'

grid_pred <- rbind(dat_ab_brohor,dat_ab_vulmic,dat_ab_plaere,dat_ab_miccal)

##############################################
# Making site level predictions on fecundity
##############################################

# BROHOR

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- brohor_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for brohor_dat
identical(length(brohor_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(brohor_dat, site, biotic_treatment )
head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_B <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_B <- rbind(site_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above.

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- brohor_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for brohor_dat
identical(length(brohor_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_A <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_A <- rbind(site_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and extract median and quantiles for each site

dat_a <- site_sum_dat_A %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_baseline = median(site_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- site_sum_dat_B %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(site_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and make a single summarized dataframe

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -site, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$site_predicted_fecundity <- NA
dat_ab$site_predicted_fecundity_lwr <- NA
dat_ab$site_predicted_fecundity_upr <- NA

dat_ab$site_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$site_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$site_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_brohor <- select(dat_ab, site, biotic_treatment, site_predicted_fecundity, site_predicted_fecundity_lwr, site_predicted_fecundity_upr)

# VULMIC

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- vulmic_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for vulmic_dat
identical(length(vulmic_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(vulmic_dat, site, biotic_treatment )
head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_B <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_B <- rbind(site_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above.

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- vulmic_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for vulmic_dat
identical(length(vulmic_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_A <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_A <- rbind(site_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and extract median and quantiles for each site

dat_a <- site_sum_dat_A %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_baseline = median(site_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- site_sum_dat_B %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(site_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and make a single summarized dataframe

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -site, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$site_predicted_fecundity <- NA
dat_ab$site_predicted_fecundity_lwr <- NA
dat_ab$site_predicted_fecundity_upr <- NA

dat_ab$site_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$site_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$site_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_vulmic <- select(dat_ab, site, biotic_treatment, site_predicted_fecundity, site_predicted_fecundity_lwr, site_predicted_fecundity_upr)


# PLAERE

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- plaere_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for plaere_dat
identical(length(plaere_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(plaere_dat, site, biotic_treatment )
head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_B <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_B <- rbind(site_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above.

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- plaere_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for plaere_dat
identical(length(plaere_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_A <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_A <- rbind(site_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and extract median and quantiles for each site

dat_a <- site_sum_dat_A %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_baseline = median(site_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- site_sum_dat_B %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(site_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and make a single summarized dataframe

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -site, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$site_predicted_fecundity <- NA
dat_ab$site_predicted_fecundity_lwr <- NA
dat_ab$site_predicted_fecundity_upr <- NA

dat_ab$site_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$site_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$site_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_plaere <- select(dat_ab, site, biotic_treatment, site_predicted_fecundity, site_predicted_fecundity_lwr, site_predicted_fecundity_upr)

# MICCAL

#-----------------------
# Biotic trt: no intra
#-----------------------

# extract the type of predicted fecundity we need
new_dat <- miccal_post$y_hat_no_intra
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for miccal_dat
identical(length(miccal_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

key <- select(miccal_dat, site, biotic_treatment )
head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_B <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_B <- rbind(site_sum_dat_B, temp)
  print((i/20000)*100)
  
}

# repeat for the y_hat_baseline persistence, then you can pull A out of that one and B out of the one above.

#------------------------------
# Abiotic trt: basline lambda
#------------------------------

# extract the type of predicted fecundity we need
new_dat <- miccal_post$y_hat_baseline
dim(new_dat) # rows are iterations (20000) and columns are my data.

# rearrange this posterior data so that you can append the grid and biotic_trt identity!
trans_new_dat <- t(new_dat)
trans_new_dat <- as.data.frame(trans_new_dat)

# find grid and biotic_trt information for miccal_dat
identical(length(miccal_dat$grid) , length(trans_new_dat$V1)) # must be true to go on

head(key)

# column append these two dataframes
site_dat <- cbind(key, trans_new_dat)

# Loop and calculate sum for each grid for each iteration:
site_sum_dat_A <- data.frame()
column_names <- colnames(site_dat)[3:20002]

for(i in 1:length(column_names)) { # this indexes column number
  
  temp <- site_dat %>%
    group_by(site, biotic_treatment) %>%
    reframe(
      site_sum = sum(.data[[column_names[i]]], na.rm = TRUE),
      iteration = column_names[i]
    )
  
  site_sum_dat_A <- rbind(site_sum_dat_A, temp)
  print((i/20000)*100)
  
}

# then I will want to summarize by iteration and extract median and quantiles for each site

dat_a <- site_sum_dat_A %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_baseline = median(site_sum)) %>%
  mutate(y_hat_baseline_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_baseline_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_baseline, y_hat_baseline_lower, y_hat_baseline_upper) %>%
  distinct()

dat_b <- site_sum_dat_B %>%
  group_by(site, biotic_treatment) %>%
  mutate(y_hat_no_intra = median(site_sum)) %>%
  mutate(y_hat_no_intra_lower = quantile(site_sum, probs = 0.1)) %>%
  mutate(y_hat_no_intra_upper = quantile(site_sum, probs = 0.9)) %>%
  ungroup() %>%
  select(site, biotic_treatment, y_hat_no_intra, y_hat_no_intra_lower, y_hat_no_intra_upper) %>%
  distinct()

# select the a and b of the biotic treatment and make a single summarized dataframe

# rm grid and biotic_treatment from B
dat_b <- select(dat_b, -site, -biotic_treatment)
# bind data
dat_ab <- cbind(dat_a, dat_b)

dat_ab$site_predicted_fecundity <- NA
dat_ab$site_predicted_fecundity_lwr <- NA
dat_ab$site_predicted_fecundity_upr <- NA

dat_ab$site_predicted_fecundity <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra, dat_ab$y_hat_baseline)
dat_ab$site_predicted_fecundity_lwr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_lower, dat_ab$y_hat_baseline_lower)
dat_ab$site_predicted_fecundity_upr <-  ifelse(dat_ab$biotic_treatment %in% 'B', dat_ab$y_hat_no_intra_upper, dat_ab$y_hat_baseline_upper)

dat_ab_miccal <- select(dat_ab, site, biotic_treatment, site_predicted_fecundity, site_predicted_fecundity_lwr, site_predicted_fecundity_upr)

dat_ab_brohor$species <- 'brohor'
dat_ab_vulmic$species <- 'vulmic'
dat_ab_plaere$species <- 'plaere'
dat_ab_miccal$species <- 'miccal'

site_pred <- rbind(dat_ab_brohor,dat_ab_vulmic,dat_ab_plaere,dat_ab_miccal)

################################################
# Save dataframe for persistence-areas
################################################

# pred_dat_saved
# write csv
save(pred_data_saved, file = here::here("Supp Bayesian models/Stan rdata/pred_data_saved_80CI.rdata"))# cluster

# grid_pred
save(grid_pred, file = here::here("Supp Bayesian models/Stan rdata/grid_pred_80CI.rdata"))# cluster

# site_pred
save(site_pred, file = here::here("Supp Bayesian models/Stan rdata/site_pred_80CI.rdata"))# cluster
