
#-----------------------------------------------------------------------------.
# SPECIES SPECIFIC PATTERNS OF MISALIGNMENTS
#-----------------------------------------------------------------------------.

################################
# Load packages and functions
################################

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
#library(mclogit)

source("Scripts/Source - MAIN fitnessdata.R")

{
  # Functions
prob_trans <- function(j, J=3) { # have to have J = 3 because indexing of alpha vector above
  exp(a[j])/(1 + sum(sapply(1:J, function(i) exp(a[i])))) # correct when summation occurs from 1:J 
}
ci_prop <- function(level = 0.975, n, p) qt(level,df=n-1)*sqrt(p*(1-p)/n)

# Colors
contingency_labs <- c( "i. aligned \n present", "ii. sinks","iii. dispersal \n limitation", "iv. aligned \n absent")
contingency_cols1 <- c("slategrey", "#8c96c6" , "#88419d","violet" )

  #-----------------------------------------------------------------------------
  ## plaere ####
  #-----------------------------------------------------------------------------
  
  ##################
  # DATA WRANGLING

  # proportion of each contingency in natural conditions (biotic plots)
  
  # Raw data component:
  fig1P_dat <- plotlev %>%
    select(names, contingency, treatment, species) %>%
    filter(treatment %in% "B" & species %in% "plaere") %>%
    distinct() %>%
    group_by(names, treatment, contingency) %>%
    dplyr::mutate(prop = n()) %>%
    ungroup() %>%
    select(names, contingency, prop) %>%
    distinct() %>%
    pivot_wider(.,names_from = contingency, values_from = prop)
  
  # solve NA problem: basically proportions don't add up to 1 because of some missing data. To delineate missing data from actual 0s, I only inputed 0 in other categories if one category == 1
  fig1P_dat1 <- fig1P_dat
  fig1P_dat1[is.na(fig1P_dat1)] <- 0
  fig1P_dat1$sum <- rowSums(fig1P_dat1[2:5])
  fig1P_dat2 <- filter(fig1P_dat1, sum != 1)
  fig1P_dat2[fig1P_dat2 == 0] <- NA
  fig1P_dat1 <- filter(fig1P_dat1, !names %in% fig1P_dat2$names)
  fig1P_dat <- rbind(fig1P_dat1,fig1P_dat2)
  fig1P_dat <- select(fig1P_dat, -sum)
  fig1P_dat <- pivot_longer(fig1P_dat, cols = c(2:5) ,names_to = 'contingency', values_to = 'prop')
  fig1P_dat <- na.omit(fig1P_dat)
  
  fig1P_dat$colortreat <- NA
  fig1P_dat$colortreat[fig1P_dat$contingency == "DL"]<- "slategrey"
  fig1P_dat$colortreat[fig1P_dat$contingency == "ME"]<- "violet"
  fig1P_dat$colortreat[fig1P_dat$contingency == "SS_n"]<- "#8c96c6"
  fig1P_dat$colortreat[fig1P_dat$contingency == "SS_y"]<- "#88419d"
  col_treat_long <- c(fig1P_dat$colortreat)
  
  ##################
  ## Fit a model 

  temp_key <- name_key %>% 
    select(names, grid, site)  %>% # attach to site grid plot
    distinct()
  mod_dat_p <- left_join(fig1P_dat, temp_key, by = "names")
  mod_dat_p$contingency<-as.factor(mod_dat_p$contingency)
  mod_dat_p$site<-as.factor(mod_dat_p$site)
  mod_dat_p$grid<-as.factor(mod_dat_p$grid)
  
  mod_dat_p <- mod_dat_p %>% 
    select(-colortreat) %>% # makes each row repetitive
    pivot_wider(., names_from = 'contingency', values_from = 'prop') %>% # make wide
    distinct()
  mod_dat_p
  
  y_mat <- as.matrix(mod_dat_p[, 4:7]) # make response matrix
  #accounting for potential correlation among observations with random effects
  
  mfit.p <- mblogit(y_mat ~ 1, random = ~ 1|site/grid, data = mod_dat_p) 

  summary(mfit.p)

  ##################
  # PROBABILITIES
  
  ## passing though probability transformation gives the category that's not baseline, calculate baseline by finding all then p = (1-all)
  
  a <- c(summary(mfit.p)$coef[1,1],summary(mfit.p)$coef[2,1],summary(mfit.p)$coef[3,1]) # intercept estimates
  
  SS_n <- prob_trans(j=1) # SS_n
  SS_y <- prob_trans(j=2) # SS_y
  ME <- prob_trans(j=3) # ME
  DL <- 1 - sum(SS_n,SS_y,ME) 
  
  ######################
  # CONFIDENCE INTERVALS
  
  # can go through same prob transformation, but first find upper and lower ci from summary output of 'alpha' +- 1.96*SE 
  
  # version using proportional function 

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
  output.p <- data.frame(contingency= contingency, prob=prob, ci.upr=ci.upr, ci.lwr=ci.lwr)
  output.p[output.p < 0] <- 0
  output.p
  
  ##############
  ## Plotting
  
  # Order contingincies to mirror the table:
  output.p_ordered <- as.data.frame(output.p)
  output.p_ordered$contingency <- as.factor(output.p_ordered$contingency)
  output.p_ordered$contingency <- factor(output.p_ordered$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  output.p_ordered <- arrange(output.p_ordered, match(output.p_ordered$contingency, levels(output.p_ordered$contingency)))
  
  fig1P_dat$contingency <- as.factor(fig1P_dat$contingency)
  fig1P_dat$contingency <- factor(fig1P_dat$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  fig1P_dat <- arrange(fig1P_dat, match(fig1P_dat$contingency, levels(fig1P_dat$contingency)))
  
  fig1p_all <- left_join(fig1P_dat, output.p_ordered, by = c("contingency"))
  
  col_treat_long <- c(fig1p_all$colortreat)
  man_list <- unique(fig1p_all$contingency)
  
  p <- ggplot(fig1p_all) +
    ylim(0,1) +
    labs(x="", y = "") +
    geom_jitter(aes(x = contingency, y = prop), 
                alpha = 0.1, 
                size = 1.5, 
                width = 0.2, 
                height = 0.03,
                color = "mediumpurple1") +
    geom_point( aes(x = contingency, y = prob), 
                color = "mediumpurple1",
                size = 4) +
    geom_linerange(aes(x = contingency, y = prob, ymin = ci.lwr, ymax = ci.upr), 
                  color = "mediumpurple1",
                  linewidth = 1) +
    geom_hline(yintercept = 0.25, linetype='dashed', col = 'grey') +
    theme_bw() +
    theme(text = element_text(size = 16),
          legend.position = 'none',
          axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Proportion of populations") +
    ggtitle("Plantago") +
    scale_x_discrete(labels = contingency_labs) +
    annotate(xmin = which(man_list=="SS_y")-0.5, xmax = which(man_list=="SS_y")+0.5,
                                                          ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate(xmin = which(man_list=="SS_n")-0.5, xmax = which(man_list=="SS_n")+0.5,
             ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate("text", x = 0.5, y = 0.99, label = "A", size = 6, fontface = "bold")
  
  p
  
  #---------------------------------------------------------------------------
  ## miccal ####
  #---------------------------------------------------------------------------
  
  # proportion of each contingency in natural conditions (biotic plots)
  
  # Raw data component:
  fig1M_dat <- plotlev %>%
    select(names, contingency, treatment, species) %>%
    filter(treatment %in% "B" & species %in% "miccal") %>%
    distinct() %>%
    group_by(names, treatment, contingency) %>%
    dplyr::mutate(prop = n()) %>%
    ungroup() %>%
    select(names, contingency, prop) %>%
    distinct() %>%
    pivot_wider(.,names_from = contingency, values_from = prop)
  
  # solve NA problem: basically proportions don't add up to 1 because of some missing data. To delineate missing data from actual 0s, I only inputed 0 in other categories if one category == 1
  fig1M_dat1 <- fig1M_dat
  fig1M_dat1[is.na(fig1M_dat1)] <- 0
  fig1M_dat1$sum <- rowSums(fig1M_dat1[2:5])
  fig1M_dat2 <- filter(fig1M_dat1, sum != 1)
  fig1M_dat2[fig1M_dat2 == 0] <- NA
  fig1M_dat1 <- filter(fig1M_dat1, !names %in% fig1M_dat2$names)
  fig1M_dat <- rbind(fig1M_dat1,fig1M_dat2)
  fig1M_dat <- select(fig1M_dat, -sum)
  fig1M_dat <- pivot_longer(fig1M_dat, cols = c(2:5) ,names_to = 'contingency', values_to = 'prop')
  fig1M_dat <- na.omit(fig1M_dat)
  
  fig1M_dat$colortreat <- NA
  fig1M_dat$colortreat[fig1M_dat$contingency == "DL"]<- "slategrey"
  fig1M_dat$colortreat[fig1M_dat$contingency == "ME"]<- "violet"
  fig1M_dat$colortreat[fig1M_dat$contingency == "SS_n"]<- "#8c96c6"
  fig1M_dat$colortreat[fig1M_dat$contingency == "SS_y"]<- "#88419d"
  col_treat_long <- c(fig1M_dat$colortreat)
  
  ##################
  ## Fit a model 
  
  temp_key <- name_key %>% 
    select(names, grid, site)  %>% # attach to site grid plot
    distinct()
  mod_dat_m <- left_join(fig1M_dat, temp_key, by = "names")
  mod_dat_m$contingency<-as.factor(mod_dat_m$contingency)
  mod_dat_m$site<-as.factor(mod_dat_m$site)
  mod_dat_m$grid<-as.factor(mod_dat_m$grid)
  
  mod_dat_m <- mod_dat_m %>% 
    select(-colortreat) %>% # makes each row repetitive
    pivot_wider(., names_from = 'contingency', values_from = 'prop') %>% # make wide
    distinct()
  mod_dat_m
  
  y_mat <- as.matrix(mod_dat_m[, c(6,4,7,5)]) # order as DL, SS_n, SS_y, ME
  #accounting for potential correlation among observations with random effects
  
  mfit.m <- mblogit(y_mat ~ 1, random = ~ 1|site/grid, data = mod_dat_m) 
  
  summary(mfit.m)
  
  ##################
  # PROBABILITIES
  
  ## passing though probability transformation gives the category that's not baseline, calculate baseline by finding all then p = (1-all)
  
  a <- c(summary(mfit.m)$coef[1,1],summary(mfit.m)$coef[2,1],summary(mfit.m)$coef[3,1]) # intercept estimates
  
  SS_n <- prob_trans(j=1) # SS_n
  SS_y <- prob_trans(j=2) # SS_y
  ME <- prob_trans(j=3) # ME
  DL <- 1 - sum(SS_n,SS_y,ME) 
  
  ######################
  # CONFIDENCE INTERVALS
  
  # version using proportional function 
  
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
  output.m <- data.frame(contingency= contingency, prob=prob, ci.upr=ci.upr, ci.lwr=ci.lwr)
  output.m[output.m < 0] <- 0
  output.m
  
  ##############
  ## Plotting
  
  # Order contingincies to mirror the table:
  output.m_ordered <- as.data.frame(output.m)
  output.m_ordered$contingency <- as.factor(output.m_ordered$contingency)
  output.m_ordered$contingency <- factor(output.m_ordered$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  output.m_ordered <- arrange(output.m_ordered, match(output.m_ordered$contingency, levels(output.m_ordered$contingency)))
  
  fig1M_dat$contingency <- as.factor(fig1M_dat$contingency)
  fig1M_dat$contingency <- factor(fig1M_dat$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  fig1M_dat <- arrange(fig1M_dat, match(fig1M_dat$contingency, levels(fig1M_dat$contingency)))
  
  fig1m_all <- left_join(fig1M_dat, output.m_ordered, by = c("contingency"))
  
  col_treat_long <- c(fig1m_all$colortreat)
  man_list <- unique(fig1m_all$contingency)
  
  m <- ggplot(fig1m_all) +
    ylim(0,1) +
    labs(x="", y = "") +
    geom_jitter(aes(x = contingency, y = prop), 
                alpha = 0.1, 
                size = 1.5, 
                width = 0.2, 
                height = 0.03,
                color = "mediumpurple1") +
    geom_point( aes(x = contingency, y = prob), 
                color = "mediumpurple1",
                size = 4) +
    geom_linerange(aes(x = contingency, y = prob, ymin = ci.lwr, ymax = ci.upr), 
                   color = "mediumpurple1", linewidth = 1) +
    geom_hline(yintercept = 0.25, linetype='dashed', col = 'grey') +
    theme_bw() +
    theme(text = element_text(size = 16),
          legend.position = 'none',
          axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) +
    labs(y = "") +
    ggtitle("Micropus") +
    scale_x_discrete(labels = contingency_labs) +
    annotate(xmin = which(man_list=="SS_y")-0.5, xmax = which(man_list=="SS_y")+0.5,
             ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate(xmin = which(man_list=="SS_n")-0.5, xmax = which(man_list=="SS_n")+0.5,
             ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate("text", x = 0.5, y = 0.99, label = "B", size = 6, fontface = "bold")
  
  m
  
  #-----------------------------------------------------------------------------
  ## brohor ####
  #-----------------------------------------------------------------------------
  
  # proportion of each contingency in natural conditions (biotic plots)
  
  # Raw data component:
  fig1Br_dat <- plotlev %>%
    select(names, contingency, treatment, species) %>%
    filter(treatment %in% "B" & species %in% "brohor") %>%
    distinct() %>%
    group_by(names, treatment, contingency) %>%
    dplyr::mutate(prop = n()) %>%
    ungroup() %>%
    select(names, contingency, prop) %>%
    distinct() %>%
    pivot_wider(.,names_from = contingency, values_from = prop)
  
  # solve NA problem: basically proportions don't add up to 1 because of some missing data. To delineate missing data from actual 0s, I only inputed 0 in other categories if one category == 1
  fig1Br_dat1 <- fig1Br_dat
  fig1Br_dat1[is.na(fig1Br_dat1)] <- 0
  fig1Br_dat1$sum <- rowSums(fig1Br_dat1[2:5])
  fig1Br_dat2 <- filter(fig1Br_dat1, sum != 1)
  fig1Br_dat2[fig1Br_dat2 == 0] <- NA
  fig1Br_dat1 <- filter(fig1Br_dat1, !names %in% fig1Br_dat2$names)
  fig1Br_dat <- rbind(fig1Br_dat1,fig1Br_dat2)
  fig1Br_dat <- select(fig1Br_dat, -sum)
  fig1Br_dat <- pivot_longer(fig1Br_dat, cols = c(2:5) ,names_to = 'contingency', values_to = 'prop')
  fig1Br_dat <- na.omit(fig1Br_dat)
  
  fig1Br_dat$colortreat <- NA
  fig1Br_dat$colortreat[fig1Br_dat$contingency == "DL"]<- "slategrey"
  fig1Br_dat$colortreat[fig1Br_dat$contingency == "ME"]<- "violet"
  fig1Br_dat$colortreat[fig1Br_dat$contingency == "SS_n"]<- "#8c96c6"
  fig1Br_dat$colortreat[fig1Br_dat$contingency == "SS_y"]<- "#88419d"
  col_treat_long <- c(fig1Br_dat$colortreat)
  
  ##################
  ## Fit a model 
  
  temp_key <- name_key %>% 
    select(names, grid, site)  %>% # attach to site grid plot
    distinct()
  mod_dat_b <- left_join(fig1Br_dat, temp_key, by = "names")
  mod_dat_b$contingency<-as.factor(mod_dat_b$contingency)
  mod_dat_b$site<-as.factor(mod_dat_b$site)
  mod_dat_b$grid<-as.factor(mod_dat_b$grid)
  
  mod_dat_b <- mod_dat_b %>% 
    select(-colortreat) %>% # makes each row repetitive
    pivot_wider(., names_from = 'contingency', values_from = 'prop') %>% # make wide
    distinct()
  mod_dat_b
  
  y_mat <- as.matrix(mod_dat_b[, c(5,4,6,7)]) # order as DL, SS_n, SS_y, ME
  #accounting for potential correlation among observations with random effects
  
  mfit.b <- mblogit(y_mat ~ 1, random = ~ 1|site/grid, data = mod_dat_b) 
  
  summary(mfit.b)
  
  ##################
  # PROBABILITIES
  
  ## passing though probability transformation gives the category that's not baseline, calculate baseline by finding all then p = (1-all)
  
  a <- c(summary(mfit.b)$coef[1,1],summary(mfit.b)$coef[2,1],summary(mfit.b)$coef[3,1]) # intercept estimates
  
  SS_n <- prob_trans(j=1) # SS_n
  SS_y <- prob_trans(j=2) # SS_y
  ME <- prob_trans(j=3) # ME
  DL <- 1 - sum(SS_n,SS_y,ME) 
  
  ######################
  # CONFIDENCE INTERVALS
  
  # version using proportional function 
  
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
  output.b <- data.frame(contingency= contingency, prob=prob, ci.upr=ci.upr, ci.lwr=ci.lwr)
  output.b[output.b < 0] <- 0
  output.b
  
  ##############
  ## Plotting
  
  # Order contingincies to mirror the table:
  output.b_ordered <- as.data.frame(output.b)
  output.b_ordered$contingency <- as.factor(output.b_ordered$contingency)
  output.b_ordered$contingency <- factor(output.b_ordered$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  output.b_ordered <- arrange(output.b_ordered, match(output.b_ordered$contingency, levels(output.b_ordered$contingency)))
  
  fig1Br_dat$contingency <- as.factor(fig1Br_dat$contingency)
  fig1Br_dat$contingency <- factor(fig1Br_dat$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  
  
fig1Br_dat <- arrange(fig1Br_dat, match(fig1Br_dat$contingency, levels(fig1Br_dat$contingency)))

  fig1b_all <- left_join(fig1Br_dat, output.b_ordered, by = c("contingency"))
  
  col_treat_long <- c(fig1b_all$colortreat)
  man_list <- unique(fig1b_all$contingency)
  
  b <- ggplot(fig1b_all) +
    ylim(0,1) +
    labs(x="", y = "") +
    geom_jitter(aes(x = contingency, y = prop), 
                alpha = 0.1, 
                size = 1.5, 
                width = 0.2, 
                height = 0.03,
                color = "mediumpurple1") +
    geom_point( aes(x = contingency, y = prob), 
                color = "mediumpurple1",
                size = 4) +
    geom_linerange(aes(x = contingency, y = prob, ymin = ci.lwr, ymax = ci.upr), 
                   color = "mediumpurple1", linewidth = 1) +
    geom_hline(yintercept = 0.25, linetype='dashed', col = 'grey') +
    theme_bw() +
    theme(text = element_text(size = 16),
          legend.position = 'none',
          axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) +
    labs(y = "Proportion of populations") +
    ggtitle("Bromus") +
    scale_x_discrete(labels = contingency_labs) +
    annotate(xmin = which(man_list=="SS_y")-0.5, xmax = which(man_list=="SS_y")+0.5,
             ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate(xmin = which(man_list=="SS_n")-0.5, xmax = which(man_list=="SS_n")+0.5,
             ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate("text", x = 0.5, y = 0.99, label = "C", size = 6, fontface = "bold")
  
  b
  
  #-----------------------------------------------------------------------------
  ## vulmic ####
  #-----------------------------------------------------------------------------
  
  # proportion of each contingency in natural conditions (biotic plots)
  
  # Raw data component:
  fig1V_dat <- plotlev %>%
    select(names, contingency, treatment, species) %>%
    filter(treatment %in% "B" & species %in% "vulmic") %>%
    distinct() %>%
    group_by(names, treatment, contingency) %>%
    dplyr::mutate(prop = n()) %>%
    ungroup() %>%
    select(names, contingency, prop) %>%
    distinct() %>%
    pivot_wider(.,names_from = contingency, values_from = prop)
  
  # solve NA problem: basically proportions don't add up to 1 because of some missing data. To delineate missing data from actual 0s, I only inputed 0 in other categories if one category == 1
  fig1V_dat1 <- fig1V_dat
  fig1V_dat1[is.na(fig1V_dat1)] <- 0
  fig1V_dat1$sum <- rowSums(fig1V_dat1[2:5])
  fig1V_dat2 <- filter(fig1V_dat1, sum != 1)
  fig1V_dat2[fig1V_dat2 == 0] <- NA
  fig1V_dat1 <- filter(fig1V_dat1, !names %in% fig1V_dat2$names)
  fig1V_dat <- rbind(fig1V_dat1,fig1V_dat2)
  fig1V_dat <- select(fig1V_dat, -sum)
  fig1V_dat <- pivot_longer(fig1V_dat, cols = c(2:5) ,names_to = 'contingency', values_to = 'prop')
  fig1V_dat <- na.omit(fig1V_dat)
  
  fig1V_dat$colortreat <- NA
  fig1V_dat$colortreat[fig1V_dat$contingency == "DL"]<- "slategrey"
  fig1V_dat$colortreat[fig1V_dat$contingency == "ME"]<- "violet"
  fig1V_dat$colortreat[fig1V_dat$contingency == "SS_n"]<- "#8c96c6"
  fig1V_dat$colortreat[fig1V_dat$contingency == "SS_y"]<- "#88419d"
  col_treat_long <- c(fig1V_dat$colortreat)

  
  ##################
  ## Fit a model 
  
  temp_key <- name_key %>% 
    select(names, grid, site)  %>% # attach to site grid plot
    distinct()
  mod_dat_v <- left_join(fig1V_dat, temp_key, by = "names")
  mod_dat_v$contingency<-as.factor(mod_dat_v$contingency)
  mod_dat_v$site<-as.factor(mod_dat_v$site)
  mod_dat_v$grid<-as.factor(mod_dat_v$grid)
  
  mod_dat_v <- mod_dat_v %>% 
    select(-colortreat) %>% # makes each row repetitive
    pivot_wider(., names_from = 'contingency', values_from = 'prop') %>% # make wide
    distinct()
  mod_dat_v
  
  y_mat <- as.matrix(mod_dat_v[, 4:7]) # order as DL, SS_n, SS_y, ME
  #accounting for potential correlation among observations with random effects
  
  mfit.v <- mblogit(y_mat ~ 1, random = ~ 1|site/grid, data = mod_dat_v) 
  
  summary(mfit.v)
  
  ##################
  # PROBABILITIES
  
  ## passing though probability transformation gives the category that's not baseline, calculate baseline by finding all then p = (1-all)
  
  a <- c(summary(mfit.v)$coef[1,1],summary(mfit.v)$coef[2,1],summary(mfit.v)$coef[3,1]) # intercept estimates
  
  SS_n <- prob_trans(j=1) # SS_n
  SS_y <- prob_trans(j=2) # SS_y
  ME <- prob_trans(j=3) # ME
  DL <- 1 - sum(SS_n,SS_y,ME) 
  
  ######################
  # CONFIDENCE INTERVALS

  # version using proportional function 
  
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
  output.v <- data.frame(contingency= contingency, prob=prob, ci.upr=ci.upr, ci.lwr=ci.lwr)
  output.v[output.v < 0] <- 0
  output.v
  
  ##############
  ## Plotting
  
  # Order contingincies to mirror the table:
  output.v_ordered <- as.data.frame(output.v)
  output.v_ordered$contingency <- as.factor(output.v_ordered$contingency)
  output.v_ordered$contingency <- factor(output.v_ordered$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  output.v_ordered <- arrange(output.v_ordered, match(output.v_ordered$contingency, levels(output.v_ordered$contingency)))
  
  fig1V_dat$contingency <- as.factor(fig1V_dat$contingency)
  fig1V_dat$contingency <- factor(fig1V_dat$contingency, levels = c(
    "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
  fig1V_dat <- arrange(fig1V_dat, match(fig1V_dat$contingency, levels(fig1V_dat$contingency)))
  
  fig1v_all <- left_join(fig1V_dat, output.v_ordered, by = c("contingency"))
  
  col_treat_long <- c(fig1v_all$colortreat)
  man_list <- unique(fig1v_all$contingency)
  
  v <- ggplot(fig1v_all) +
    ylim(0,1) +
    labs(x="", y = "") +
    geom_jitter(aes(x = contingency, y = prop), 
                alpha = 0.1, 
                size = 1.5, 
                width = 0.2, 
                height = 0.03,
                color = "mediumpurple1") +
    geom_point( aes(x = contingency, y = prob), 
                color = "mediumpurple1",
                size = 4) +
    geom_linerange(aes(x = contingency, y = prob, ymin = ci.lwr, ymax = ci.upr), 
                   linewidth = 1, color = "mediumpurple1") +
    geom_hline(yintercept = 0.25, linetype='dashed', col = 'grey') +
    theme_bw() +
    theme(text = element_text(size = 16),
          legend.position = 'none',
          axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) +
    labs(y = "") +
    ggtitle("Festuca") +
    scale_x_discrete(labels = contingency_labs) +
    annotate(xmin = which(man_list=="SS_y")-0.5, xmax = which(man_list=="SS_y")+0.5,
             ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate(xmin = which(man_list=="SS_n")-0.5, xmax = which(man_list=="SS_n")+0.5,
             ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
    annotate("text", x = 0.5, y = 0.99, label = "D", size = 6, fontface = "bold")
v
}

##########################################################
library(patchwork)

pdf('Figures/fig_supp_species.pdf', width = 9, height = 9)
(p + m) / (b + v)
dev.off()
##########################################################

# RESULTS: species follow similar patterns. miccal and vulmic show dominant 
# sink populations, whereas plaere and brohor show dominant aligned absent. 
# aligned present and dispersal limitation are least common for all four species.

##############################
# TABLE
##############################

# summarize all model estimates
output.p_ordered$species <- NA
output.p_ordered$species <- "Plantago erecta"
output.m_ordered$species <- NA
output.m_ordered$species <- "Micropus californicus"
output.b_ordered$species <- NA
output.b_ordered$species <- "Bromus hordeaceus"
output.v_ordered$species <- NA
output.v_ordered$species <- "Festuca microstachys"
output_ordered <- rbind(output.p_ordered,output.m_ordered,output.b_ordered,output.v_ordered)

# get the order and names of misalignments correct:
output_ordered$contingency <- as.factor(output_ordered$contingency)
output_ordered$contingency <- case_match(output_ordered$contingency,
                                  'SS_y' ~ 'i. aligned present',
                                  'ME' ~ 'ii. sink',
                                  'DL' ~ 'iii. dispersal limitation',
                                  'SS_n' ~ 'iv. aligned absent')
output_ordered <- output_ordered %>%
  arrange(species, contingency)

tab_mulinomial_estimates <- output_ordered |>
  dplyr::mutate(conf.int =paste(round(ci.lwr,2), round(ci.upr,2), sep = ", ")) |>
  dplyr::mutate(prob = round(prob,2)) |>
  dplyr::select(species, contingency, prob, conf.int) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    species = "Species",
    contingency = "(Mis)alignment",
    prob = "Estimated proportion",
    conf.int =  "95% CI") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))
tab_mulinomial_estimates |>
  gtsave(paste0(here::here(),"/Tables/15tab_multinomial_estimates.pdf")) 
