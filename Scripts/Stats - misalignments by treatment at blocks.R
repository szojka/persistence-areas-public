
#-----------------------------------------------------------------------------.
# DESCRIPTION: Fitting a glm to block level data (n = 450) for both neighbor 
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

source("Scripts/Source - MAIN fitnessdata.R")

#--------------------------------------------------------------------------.
# POPULATION SCALE (PLOT)
#--------------------------------------------------------------------------.

########################
# Data cleaning
########################

fig2B_dat <- plotlev %>%
  select(names, grid, site, contingency, treatment, species) %>%
  distinct() %>%
  dplyr::group_by(names, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(names, treatment, contingency) %>%
  dplyr::mutate(prop = n()/replicates) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = prop)
# all NAs = 0 because I standardized by replicate
fig2B_dat[is.na(fig2B_dat)] <- 0
fig2B_dat <- pivot_longer(fig2B_dat, cols = 6:9, names_to = "contingency", values_to = "prop")

fig2B_dat$contingency<-as.factor(fig2B_dat$contingency)
fig2B_dat$site<-as.factor(fig2B_dat$site)
fig2B_dat$grid<-as.factor(fig2B_dat$grid)
fig2B_dat$names<-as.factor(fig2B_dat$names)
fig2B_dat$treatment<-as.factor(fig2B_dat$treatment)

########################
## Fit models 
########################

################
# sinks
m.me.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1| grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'ME'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.me.p)
testDispersion(m.me.p) 
plot(fitted(m.me.p), residuals(m.me.p))
hist(residuals(m.me.p)) 
loc_sim_ouput <- simulateResiduals(m.me.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.me.p)

a <- Anova(m.me.p, type = 3)

# save anova values for table:
plot.me.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
           scale = c('plot','plot'),
           Df =c(a$Df[1],a$Df[2]),
           P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
           predictor = c("Intercept","Data type"),
           contingency = c('sinks','sinks'))


em <- emmeans(m.me.p, ~treatment, type = "response") # specify green_index_scaled values
em

plot.me.contrast <- pairs(em)
pairs(em)

################
# d. limitation
m.dl.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1| grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'DL'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.dl.p)
testDispersion(m.dl.p) 
plot(fitted(m.dl.p), residuals(m.dl.p))
hist(residuals(m.dl.p)) 
loc_sim_ouput <- simulateResiduals(m.dl.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.dl.p)

a <- Anova(m.dl.p, type = 3)

# save anova values for table:
plot.dl.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                         scale = c('plot','plot'),
                         Df =c(a$Df[1],a$Df[2]),
                         P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                         predictor = c("Intercept","Data type"),
                         contingency = c('dispersal limitation','dispersal limitation'))

em <- emmeans(m.dl.p, ~treatment, type = "response") # specify green_index_scaled values
em

plot.dl.contrast <- pairs(em)

################
# a. present
m.ap.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1| grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'SS_y'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.ap.p)
testDispersion(m.ap.p) 
plot(fitted(m.ap.p), residuals(m.ap.p))
hist(residuals(m.ap.p)) 
loc_sim_ouput <- simulateResiduals(m.ap.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.ap.p)

a <- Anova(m.ap.p, type = 3)

# save anova values for table:
plot.ap.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('plot','plot'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned present','aligned present'))

em <- emmeans(m.ap.p, ~treatment, type = "response") # specify green_index_scaled values
em

plot.ap.contrast <- pairs(em)

################
# a. absent
m.aa.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1 | grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'SS_n'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.aa.p)
testDispersion(m.aa.p) 
plot(fitted(m.aa.p), residuals(m.aa.p))
hist(residuals(m.aa.p)) 
loc_sim_ouput <- simulateResiduals(m.aa.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

summary(m.aa.p)

a <- Anova(m.aa.p, type = 3)

# save anova values for table:
plot.aa.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('plot','plot'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned absent','aligned absent'))


em <- emmeans(m.aa.p, ~treatment, type = "response") # specify green_index_scaled values
em

plot.aa.contrast <- pairs(em)

##########################
# Create objects for viz
##########################

vis.me.p <- ggpredict(m.me.p, 
                      terms = c("treatment"), 
                      type = "fe")
vis.me.p$contingency <- c("ME")
vis.me.p$scale <- c("plot")

vis.ap.p <- ggpredict(m.ap.p, 
                      terms = c("treatment"), 
                      type = "fe")
vis.ap.p$contingency <- c("SS_y")
vis.ap.p$scale <- c("plot")

vis.dl.p <- ggpredict(m.dl.p, 
                      terms = c("treatment"), 
                      type = "fe")
vis.dl.p$contingency <- c("DL")
vis.dl.p$scale <- c("plot")

vis.aa.p <- ggpredict(m.aa.p, 
                      terms = c("treatment"), 
                      type = "fe")
vis.aa.p$contingency <- c("SS_n")
vis.aa.p$scale <- c("plot")


############################
# Create objects for table

# tab1. ANOVA

plot.anova <- rbind(plot.me.anova,plot.dl.anova,plot.ap.anova,plot.aa.anova)

# tab2. Pairs significance

plot.me.contrast <- as.data.frame(plot.me.contrast)
plot.dl.contrast <- as.data.frame(plot.dl.contrast)
plot.ap.contrast <- as.data.frame(plot.ap.contrast)
plot.aa.contrast <- as.data.frame(plot.aa.contrast)

plot.contrast <- rbind(plot.me.contrast,plot.dl.contrast,plot.ap.contrast,plot.aa.contrast)

plot.contrast$contingency <- NA
plot.contrast$contingency <- c("sink","dispersal limitation", "aligned present", "aligned absent")

plot.contrast$scale <- NA
plot.contrast$scale <- 'plot'


