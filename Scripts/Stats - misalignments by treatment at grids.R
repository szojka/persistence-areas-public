
#-----------------------------------------------------------------------------.
# DESCRIPTION: Fitting a glm to grid level data (n = 18) for both neighbor 
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
# META-POP SCALE (GRIDS)
#--------------------------------------------------------------------------.

########################
# Data cleaning
########################

fig2D_dat <- gridlev %>%
  select(grid, site, contingency, treatment, species) %>%
  distinct() %>%
  dplyr::group_by(grid, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(grid, treatment, contingency) %>%
  dplyr::mutate(prop = n()/replicates) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = prop)
# all NAs = 0 because I standardized by replicate
fig2D_dat[is.na(fig2D_dat)] <- 0
fig2D_dat <- pivot_longer(fig2D_dat, cols = 5:8, names_to = "contingency", values_to = "prop")

fig2D_dat$contingency<-as.factor(fig2D_dat$contingency)
fig2D_dat$site<-as.factor(fig2D_dat$site)
fig2D_dat$grid<-as.factor(fig2D_dat$grid)
fig2D_dat$treatment<-as.factor(fig2D_dat$treatment)

########################
## Fit models 
########################

################
# sinks
m.me.g <- glmmTMB(prop ~ treatment + (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'ME'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.me.g)
testDispersion(m.me.g) 
plot(fitted(m.me.g), residuals(m.me.g))
hist(residuals(m.me.g)) 
loc_sim_ouput <- simulateResiduals(m.me.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.me.g)

a <- Anova(m.me.g, type = 3)

# save anova values for table:
grid.me.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('sink','sink'))

em <- emmeans(m.me.g, ~treatment, type = "response") # specify green_index_scaled values
em

grid.me.contrast <- pairs(em)

################
# d. limitation
m.dl.g <- glmmTMB(prop ~ treatment +  (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'DL'), 
                  weights = replicates, 
                  family = binomial(),
                  ziformula=~1, # adding a zero inflation term fixes convergence
                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))) 

# Check model fit 
testZeroInflation(m.dl.g)
testDispersion(m.dl.g) 
plot(fitted(m.dl.g), residuals(m.dl.g))
hist(residuals(m.dl.g)) 
loc_sim_ouput <- simulateResiduals(m.dl.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.dl.g)

a <- Anova(m.dl.g, type = 3)

# save anova values for table:
grid.dl.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('dispersal limitation','dispersal limitation'))


em <- emmeans(m.dl.g, ~treatment, type = "response") # specify green_index_scaled values
em

grid.dl.contrast <- pairs(em)

################
# a. present
m.ap.g <- glmmTMB(prop ~ treatment +  (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'SS_y'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.ap.g)
testDispersion(m.ap.g) 
plot(fitted(m.ap.g), residuals(m.ap.g))
hist(residuals(m.ap.g)) 
loc_sim_ouput <- simulateResiduals(m.ap.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.ap.g)


# save anova values for table:
grid.ap.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned present','aligned present'))

em <- emmeans(m.ap.g, ~treatment, type = "response") # specify green_index_scaled values
em

grid.ap.contrast <- pairs(em)

################
# a. absent
m.aa.g <- glmmTMB(prop ~ treatment +  (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'SS_n'), 
                  weights = replicates, 
                  family = binomial())

# Check model fit 
testZeroInflation(m.aa.g)
testDispersion(m.aa.g) 
plot(fitted(m.aa.g), residuals(m.aa.g))
hist(residuals(m.aa.g)) 
loc_sim_ouput <- simulateResiduals(m.aa.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.aa.g)

a <- Anova(m.aa.g, type = 3)

# save anova values for table:
grid.aa.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned absent','aligned absent'))


em <- emmeans(m.aa.g, ~treatment, type = "response") # specify green_index_scaled values
em

grid.aa.contrast <- pairs(em)

##########################
# Create objects for viz
##########################

vis.me.g <- ggpredict(m.me.g, 
                      terms = c("treatment"), 
                      type = "fe")
vis.me.g$contingency <- c("ME")
vis.me.g$scale <- c("grid")

vis.ap.g <- ggpredict(m.ap.g, 
                      terms = c("treatment"), 
                      type = "fe")
vis.ap.g$contingency <- c("SS_y")
vis.ap.g$scale <- c("grid")

vis.dl.g <- ggpredict(m.dl.g, 
                      terms = c("treatment"), 
                      type = "fe")
vis.dl.g$contingency <- c("DL")
vis.dl.g$scale <- c("grid")

vis.aa.g <- ggpredict(m.aa.g, 
                      terms = c("treatment"), 
                      type = "fe")
vis.aa.g$contingency <- c("SS_n")
vis.aa.g$scale <- c("grid")


############################
# Create objects for table

# tab1. ANOVA

grid.anova <- rbind(grid.me.anova,grid.dl.anova,grid.ap.anova,grid.aa.anova)

# tab2. Pairs significance

grid.me.contrast <- as.data.frame(grid.me.contrast)
grid.dl.contrast <- as.data.frame(grid.dl.contrast)
grid.ap.contrast <- as.data.frame(grid.ap.contrast)
grid.aa.contrast <- as.data.frame(grid.aa.contrast)

grid.contrast <- rbind(grid.me.contrast,grid.dl.contrast,grid.ap.contrast,grid.aa.contrast)

grid.contrast$contingency <- NA
grid.contrast$contingency <- c("sink","dispersal limitation", "aligned present", "aligned absent")

grid.contrast$scale <- NA
grid.contrast$scale <- 'grid'
