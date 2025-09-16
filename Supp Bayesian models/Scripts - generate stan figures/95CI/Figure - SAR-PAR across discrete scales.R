
#---------------------------------------------------------------------
# DESCRIPTION: Figures S14 = Figure 4 re-made using Bayesian estimates
# that control for biases and use 95% credible intervals
#----------------------------------------------------------------------

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
# library(ggeffects)

#-----------------------------------
# i. SHMEAR METHOD OF PERSISTENCE

source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats - SAR-PAR categories across discrete scales - i. averaged.R")) # loads main source too

# ii. ACCUMULATED METHOD OF PERSISTENCE
source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats - SAR-PAR categories across discrete scales - ii. accumulated.R")) # loads main source too


##################################################
# BRING RAW DATA TOGETHER FROM SCALING TYPES

# using means and raw data from par_all_discrete

par_all_discrete$treatment <- NA
par_all_discrete$treatment[par_all_discrete$type %in% 
                             c("Potential \n with neighbors","Diversity (SAR)","Realized \n with neighbors")] <- "With neighbors"
par_all_discrete$treatment[par_all_discrete$type %in% 
                             c("Potential \n without neighbors","Realized \n without neighbors")] <- "Without neighbors"

par_all_discrete$kind <- NA
par_all_discrete$kind[par_all_discrete$type %in% 
                        c("Potential \n with neighbors","Potential \n without neighbors")] <- "Potential (PAR)"

par_all_discrete$kind[par_all_discrete$type %in% 
                        c("Realized \n with neighbors","Realized \n without neighbors")] <- "Realized (PAR)"

par_all_discrete$kind[par_all_discrete$type %in% 
                        c("Diversity (SAR)")] <- "Diversity (SAR)"

par_all_discrete$scenario <- NA
par_all_discrete$scenario <- "Data"

par_all_discrete$scale_type <- NA
par_all_discrete$scale_type <- "Averaged"

# using par_all_accum

levels(par_all_accum$type)
par_all_accum$treatment <- NA
par_all_accum$treatment[par_all_accum$type %in% 
                             c("Potential \n with neighbors","Diversity (SAR)","Realized \n with neighbors")] <- "With neighbors"
par_all_accum$treatment[par_all_accum$type %in% 
                             c("Potential \n without neighbors","Realized \n without neighbors")] <- "Without neighbors"

par_all_accum$kind <- NA
par_all_accum$kind[par_all_accum$type %in% 
                        c("Potential \n with neighbors","Potential \n without neighbors")] <- "Potential (PAR)"

par_all_accum$kind[par_all_accum$type %in% 
                        c("Realized \n with neighbors","Realized \n without neighbors")] <- "Realized (PAR)"

par_all_accum$kind[par_all_accum$type %in% 
                        c("Diversity (SAR)")] <- "Diversity (SAR)"

par_all_accum$scenario <- NA
par_all_accum$scenario <- "Data"

par_all_accum$scale_type <- NA
par_all_accum$scale_type <- "Accumulated"

par_all <- rbind(par_all_accum, par_all_discrete)

###############################################
# BRING PREDICTIONS TOGETHER FROM SCALING TYPES

vis <- rbind(vis1,vis2)
names(vis)

#####################################
# DATA FIGURE

dodge <- position_dodge(width=0.5) 
par_dat <- ggplot(vis) + 
  facet_grid(cols = vars(scale_type), rows = vars(scenario)) +
  geom_point(mapping = aes(x = scale, y = prediction, color = kind, shape = treatment, group = type),
             position = dodge,inherit.aes = TRUE, size = 3) +
  geom_line(mapping = aes(x = scale, y = prediction, color = kind, linetype = treatment, group = type), position = dodge, linewidth = 1) +
  geom_pointrange(mapping = aes(x = scale, y = prediction, ymin = prediction_lwr, ymax = prediction_upr, color = kind, shape = treatment, group = type), position = dodge, inherit.aes = TRUE) +
  theme_bw() +
  labs(x = "Spatial scale", y = "Proportion of species", color = "", linetype = "") + # I want y to be proportion of species
  scale_linetype_manual(values = c("Without neighbors" = "dashed",
                                   "With neighbors" = "solid"), 
                        guide = "none") +
  scale_shape_manual(values = c("Without neighbors" = 2,
                                "With neighbors" = 19)) +
    color = guide_legend(title = "", order = 1, nrow = 3),
    shape = guide_legend(title = "", order = 2, nrow = 2),
    linetype = guide_legend(title = "", order = 2, nrow = 2) # shape + line share row 2
  ) +
  scale_color_manual(values = c("Diversity (SAR)" = "blue",
                                "Potential (PAR)" = "orchid",
                                "Realized (PAR)" = "orange",
                                guide = "none"
  )) +
  theme(legend.position = "top",
        text = element_text(size = 16)) +
  ylim(0, 1) 
par_dat

#---------------------------------------------------------------------------
library(patchwork)
png("Figures/Bayes/95CI/discrete_SAR-PAR_full.png", height = 5, width = 7, units = "in", res=800)
par_dat 
dev.off()
#---------------------------------------------------------------------------
