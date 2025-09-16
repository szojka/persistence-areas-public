#-----------------------------------------------------------------------------.
# DESCRIPTION: Fitting a glm to site level data (n = 6) for both neighbor 
# treatments, to find the estimated proportions of each (mis)alignment for 
# Figure 3. Then, calculate confidence intervals around these proportions, 
# and compare proportions to assess significance using emmeans. 
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

#source("Scripts/Source - MAIN fitnessdata.R")

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
  dplyr::mutate(prop = n()/replicates) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = prop)
# all NAs = 0 because I standardized by replicate
fig2F_dat[is.na(fig2F_dat)] <- 0
fig2F_dat$SS_n <- rep(0, length.out = length(fig2F_dat$ME))
fig2F_dat$DL <- rep(0, length.out = length(fig2F_dat$ME))
fig2F_dat <- pivot_longer(fig2F_dat, cols = 4:7, names_to = "contingency", values_to = "prop")

fig2F_dat$contingency<-as.factor(fig2F_dat$contingency)
fig2F_dat$site<-as.factor(fig2F_dat$site)
fig2F_dat$treatment<-as.factor(fig2F_dat$treatment)

########################
## Fit models 
########################

################
# sinks
m.me.r <- glmmTMB(prop ~ treatment + (1|site), 
                  data = filter(fig2F_dat, contingency %in% 'ME'), 
                  weights = replicates, 
                  family = binomial())

# Check model fit 
testZeroInflation(m.me.r)
testDispersion(m.me.r) 
plot(fitted(m.me.r), residuals(m.me.r))
hist(residuals(m.me.r)) 
loc_sim_ouput <- simulateResiduals(m.me.r)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.me.r)

a <- Anova(m.me.r, type = 3)

# save anova values for table:
site.me.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('site','site'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('sink','sink'))

em <- emmeans(m.me.r, ~treatment, type = "response") # specify green_index_scaled values
em

site.me.contrast <- pairs(em)

################
# d. limitation
# m.dl.r <- glmmTMB(prop ~ treatment, 
#                   data = filter(fig2F_dat, contingency %in% 'DL'), 
#                   weights = replicates, 
#                   family = binomial()) # convergence issue, because all regional sinks = 0


################
# a. present
m.ap.r <- glmmTMB(prop ~ treatment + (1|site), 
                  data = filter(fig2F_dat, contingency %in% 'SS_y'), 
                  weights = replicates, 
                  family = binomial())

# Check model fit 
testZeroInflation(m.ap.r)
testDispersion(m.ap.r) 
plot(fitted(m.ap.r), residuals(m.ap.r))
hist(residuals(m.ap.r)) 
loc_sim_ouput <- simulateResiduals(m.ap.r)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.ap.r)

a <- Anova(m.ap.r, type = 3)

# save anova values for table:
site.ap.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('site','site'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned present','aligned present'))

em <- emmeans(m.ap.r, ~treatment, type = "response") # specify green_index_scaled values
em

site.ap.contrast <- pairs(em)

################
# a. absent
# m.aa.r <- glmmTMB(prop ~ treatment, 
#                   data = filter(fig2F_dat, contingency %in% 'SS_n'), 
#                   weights = replicates, 
#                   family = binomial()) # all zeros!

##########################
# Create objects for viz
##########################

vis.me.r <- ggpredict(m.me.r, 
                      terms = c("treatment"), 
                      type = "fe")
vis.me.r$contingency <- c("ME")
vis.me.r$scale <- c("site")

vis.ap.r <- ggpredict(m.ap.r, 
                      terms = c("treatment"), 
                      type = "fe")
vis.ap.r$contingency <- c("SS_y")
vis.ap.r$scale <- c("site")

vis.dl.r <- data.frame(x = c("A", "B"), 
                       predicted = c(0,0), 
                       std.error = c(0,0), 
                       conf.low = c(0,0), 
                       conf.high = c(0,0), 
                       group = c(1, 1), 
                       contingency = "DL", 
                       scale = "site")

vis.aa.r <- data.frame(x = c("A", "B"), 
                       predicted = c(0,0), 
                       std.error = c(0,0), 
                       conf.low = c(0,0), 
                       conf.high = c(0,0), 
                       group = c(1, 1), 
                       contingency = "SS_n", 
                       scale = "site")

############################
# Create objects for table

# tab1. ANOVA

# fill in last two contingencies that were zeros and thus couldn't be fit as models
site.dl.anova <-  data.frame(Chi.squared = c(NA, NA),
                             scale = c('site','site'),
                             Df =c(NA, NA),
                             P_value = c(NA, NA),
                             predictor = c("Intercept","Data type"),
                             contingency = c('dispersal limitation','dispersal limitation'))
site.aa.anova <-  data.frame(Chi.squared = c(NA, NA),
                             scale = c('site','site'),
                             Df =c(NA, NA),
                             P_value =  c(NA, NA),
                             predictor = c("Intercept","Data type"),
                             contingency = c('aligned absent','aligned absent'))
site.anova <- rbind(site.me.anova,site.dl.anova,site.ap.anova,site.aa.anova)

# tab2. Pairs significance

site.me.contrast <- as.data.frame(site.me.contrast)
site.dl.contrast <- data.frame(contrast = "A / B", odds.ratio = NA, SE = NA, df = NA, null = NA, z.ratio = NA, p.value = NA)
site.ap.contrast <- as.data.frame(site.ap.contrast)
site.aa.contrast <- data.frame(contrast = "A / B", odds.ratio = NA, SE = NA, df = NA, null = NA, z.ratio = NA, p.value = NA)

site.contrast <- rbind(site.me.contrast,site.dl.contrast,site.ap.contrast,site.aa.contrast)

site.contrast$contingency <- NA
site.contrast$contingency <- c("sink","dispersal limitation", "aligned present", "aligned absent")

site.contrast$scale <- NA
site.contrast$scale <- 'site'
