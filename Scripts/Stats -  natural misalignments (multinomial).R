
#------------------------------------------------------------------------------------.
# DESCRIPTION: Figure 1d (mis)alignments in natural conditions (i.e., with neighbors)
#------------------------------------------------------------------------------------.

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

source("Scripts/Source - MAIN fitnessdata.R")

#--------------------------------------------------------------------------.
# Data cleaning
#--------------------------------------------------------------------------.

fig1B_dat <- plotlev %>%
  select(names, contingency, treatment, species, site, grid) %>%
  filter(treatment %in% "B") %>%
  distinct() %>%
  dplyr::group_by(names) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(names, treatment, contingency) %>%
  dplyr::mutate(prop = n()/replicates) %>%
  ungroup() %>%
  select(names, contingency, prop, replicates, site, grid) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = prop)
fig1B_dat[is.na(fig1B_dat)] <- 0
fig1B_dat <- pivot_longer(fig1B_dat, cols = 5:8, names_to = "contingency", values_to = "prop")

# add colors for plotting into dataframe to make it easy in ggplot
fig1B_dat$colortreat <- NA
fig1B_dat$colortreat[fig1B_dat$contingency == "DL"]<- "slategrey"
fig1B_dat$colortreat[fig1B_dat$contingency == "ME"]<- "violet"
fig1B_dat$colortreat[fig1B_dat$contingency == "SS_n"]<- "#8c96c6"
fig1B_dat$colortreat[fig1B_dat$contingency == "SS_y"]<- "#88419d"

# variable types
fig1B_dat$names <- as.factor(fig1B_dat$names )
fig1B_dat$contingency <- as.factor(fig1B_dat$contingency)
fig1B_dat$colortreat <- as.factor(fig1B_dat$colortreat)

# # A CHECK: using replicates to fill in which NA's are actually 0s
# 
# # First for the easy cases where replicates include all 4 species:
# temp1 <- filter(fig1B_dat, replicates == 4)
# # for all these sum of props should == 1 (check)
# temp1$sum <- NA
# # must change NAs to zero temporarily to perform row sum
# temp1[is.na(temp1)] <- 0
# temp1$sum <-rowSums(temp1[3:6])
# # if all sum = 1, then all NA's are zero
# print(temp1$sum)
# # keep all these NA's as zero. 
# 
# # now for replicates that do not have all 4 species:
# temp2 <- filter(fig1B_dat, !replicates == 4) 
# # must change NAs to zero temporarily to perform row sum
# temp2[is.na(temp2)] <- 0
# # these stay NA because we cannot assume which category is zero
# temp2$sum <- NA
# temp2$sum <-rowSums(temp2[3:6])
# # unless the sum is 1, then the other categories can be assumed to be 0.
# # all NAs are actually 0
# #----------------------------------------------------------------------.

# fig1B_dat # proportion of contingencies in biotic conditions
mod_dat <- fig1B_dat
mod_dat

mod_dat$contingency<-as.factor(mod_dat$contingency)
mod_dat$names<-as.factor(mod_dat$names)
mod_dat$site<-as.factor(mod_dat$site)
mod_dat$grid<-as.factor(mod_dat$grid)

# data frame for multinomial model:
head(mod_dat) 

mod_dat3 <- mod_dat %>% 
  select(-colortreat) %>% # makes each row repetitive
  pivot_wider(., names_from = 'contingency', values_from = 'prop') %>% # make wide
  distinct()
mod_dat3

y_mat <- as.matrix(mod_dat3[, 5:8]) # make response matrix, order is DL SS_n SS_y ME

#--------------------------------------------------------------------------.
## Fit a model 
#--------------------------------------------------------------------------.

#accounting for potential correlation among observations with random effects

mfit <- mblogit(y_mat ~ 1, random = ~ 1|site/grid, data = mod_dat3) # converged

summary(mfit)

#--------------------------------------------------------------------------.
# write probability transformation (to go from alpha to probability)
#--------------------------------------------------------------------------.

# PROBABILITIES

## passing though probability transformation gives the category that's not baseline, calculate baseline by finding all then p = (1-all)

a <- c(summary(mfit)$coef[1,1],summary(mfit)$coef[2,1],summary(mfit)$coef[3,1]) # intercept estimates

prob_trans <- function(j, J=3) { # have to have J = 3 because indexing of alpha vector above
  exp(a[j])/(1 + sum(sapply(1:J, function(i) exp(a[i])))) # correct when summation occurs from 1:J always :)
}

SS_n <- prob_trans(j=1) # SS_n
SS_y <- prob_trans(j=2) # SS_y
ME <- prob_trans(j=3) # ME
DL <- 1 - sum(SS_n,SS_y,ME) # sum of those 3 cats = 0.8984995

# CONFIDENCE INTERVALS
# can go through same prob transformation, but first find upper and lower ci from summary output of 'alpha' +- 1.96*SE 

# version using proportional function 
#(original veresion with transformations found at bottom of script)
ci_prop <- function(level = 0.975, n, p) qt(level,df=n-1)*sqrt(p*(1-p)/n)

upr.ssn <- SS_n + ci_prop(n= length(y_mat[,3]), p=SS_n )
lwr.ssn <- SS_n - ci_prop(n= length(y_mat[,3]), p=SS_n )
upr.ssy <- SS_y + ci_prop(n= length(y_mat[,4]), p=SS_y )
lwr.ssy <- SS_y - ci_prop(n= length(y_mat[,4]), p=SS_y )
upr.me  <- ME + ci_prop(n= length(y_mat[,1]), p=ME )
lwr.me  <- ME - ci_prop(n= length(y_mat[,1]), p=ME )
upr.dl  <- DL + ci_prop(n= length(y_mat[,2]), p=DL )
lwr.dl  <- DL - ci_prop(n= length(y_mat[,2]), p=DL )

ci.upr <- c(upr.dl,upr.ssn,upr.ssy,upr.me)
ci.lwr <- c(lwr.dl,lwr.ssn,lwr.ssy,lwr.me)

# pull into dataframe
contingency <- c("DL","SS_n","SS_y","ME")
prob <- c(DL, SS_n, SS_y, ME)
output1 <- data.frame(contingency= contingency, prob=prob, ci.upr=ci.upr, ci.lwr=ci.lwr)
output1[output1 < 0] <- 0
output1

#--------------------------------------------------------------------------.
# WALD TESTS
#--------------------------------------------------------------------------.

#https://www.econometrics.blog/post/the-wilson-confidence-interval-for-a-proportion/
## to compare each contingency 
## don't need to compute z stats or p-vals when comparing to baseline cat, as these are found in the summary

## extract covariance matrix
V <- vcov(mfit)

## define the contrast vector
# 1) ask whether the odds of being in category 2 vs 3 is different from 1.
c_vec <- c(1, -1, 0) # this is how we tell z stat which categories we are comparing
                     # bc all summary has is each category compared to baseline (cat = 1)

## compute the z-stat
z_23 <- c_vec %*% coefficients(mfit)/
  sqrt(c_vec %*% V %*% c_vec) # z score can be reported w/ p-val in paper, explains log-ratio of probabilities

## p-value
pval_23 <- pnorm(abs(z_23), lower.tail = F)

# 2) ask whether the odds of being in category 4 vs 2 is different from 1.
c_vec <- c(1, 0, -1) 

## compute the z-stat
z_42 <- c_vec %*% coefficients(mfit)/
  sqrt(c_vec %*% V %*% c_vec) # z score can be reported w/ p-val in paper, explains log-ratio of probabilities

## p-value
pval_42 <- pnorm(abs(z_42), lower.tail = F)

# 3) ask whether the odds of being in category 4 vs 3 is different from 1.
c_vec <- c(0, 1, -1) 

## compute the z-stat
z_43 <- c_vec %*% coefficients(mfit)/
  sqrt(c_vec %*% V %*% c_vec) # z score can be reported w/ p-val in paper, explains log-ratio of probabilities

## p-value
pval_43 <- pnorm(abs(z_43), lower.tail = F)

# 4) comparisons to baseline

z_21 <- summary(mfit)$coef[1,3] # SS_n~DL
pval_21 <- summary(mfit)$coef[1,4] 
z_31 <- summary(mfit)$coef[2,3] # SS_y ~ DL
pval_31 <- summary(mfit)$coef[2,4]
z_41 <- summary(mfit)$coef[3,3] # ME ~ DL
pval_41 <- summary(mfit)$coef[3,4]

# gather into dataframe 
comparing <- c("aligned absent / dispersal limitation", #2-1
               "aligned present / dispersal limtation", #3-1
               "sink / dispersal limitation", # 4-1
               "aligned absent / aligned present", # 2-3
               "sink / aligned absent", # 4-2
               "sink / aligned present") # 4-3
z_score <- c(z_21,z_31,z_41,z_23,z_42,z_43)
p_val <- c(round(pval_21,3),round(pval_31,3),round(pval_41,3),round(pval_23,3),round(pval_42,3),round(pval_43,3))

dat_compare <- data.frame(comparing = comparing, z_score = z_score, p_val = p_val)
dat_compare
# DL = 1
# SS_n = 2
# SS_y = 3
# ME = 4

