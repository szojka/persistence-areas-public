
#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 4 simulated conceptual figure demonstrating scaling type 
# (averaged suitable habitat vs. accumulated suitable habitat) for 
# persistence-area relationships (realized and potential) and 
# species-area relationships (diversity) for both with neighbors and without neighbors. 
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
# library(ggeffects)

#-----------------------------------
# i. AVERAGED METHOD OF PERSISTENCE

source("Scripts/Stats - SAR-PAR categories across discrete scales - i. averaged.R") # loads main source too

# ii. ACCUMULATED METHOD OF PERSISTENCE
source("Scripts/Stats - SAR-PAR categories across discrete scales - ii. accumulation.R") # loads main source too

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
  geom_point(mapping = aes(x = scale, y = predicted, color = kind, shape = treatment, group = type),
             position = dodge,inherit.aes = TRUE, size = 3) +
  geom_line(mapping = aes(x = scale, y = predicted, color = kind, linetype = treatment, group = type), position = dodge, linewidth = 1) +
  geom_pointrange(mapping = aes(x = scale, y = predicted, ymin = conf.low, ymax = conf.high, color = kind, shape = treatment, group = type), position = dodge, inherit.aes = TRUE) +
  theme_bw() +
  labs(x = "Spatial scale", y = "Proportion of species", color = "", linetype = "") + # I want y to be proportion of species
  scale_linetype_manual(values = c("Without neighbors" = "dashed",
                                   "With neighbors" = "solid"), 
                        guide = "none") +
  scale_shape_manual(values = c("Without neighbors" = 2,
                                "With neighbors" = 19)) +
  guides(color = "none") + 
  guides(shape = guide_legend(title = "", title.position = "left", direction = "horizontal"))+ 
  scale_color_manual(values = c("Diversity (SAR)" = "blue",
                                "Potential (PAR)" = "orchid",
                                "Realized (PAR)" = "orange", guide = "none"
  )) +
  theme(legend.position = "bottom",
        text = element_text(size = 16)) +
  ylim(0, 1) 
par_dat


################################################
# SIMULATE CONCEPTUAL COMPONENT

t <- 18*2 # 2 scaling types * 3 div cats * 3 scales * 2 scenarios

# In real data scales = x, predicted = predicted, group = div_type
scenario <- c(rep("Scenario i", times = t/2), rep("Scenario ii", times = t/2))
scale_type <- c(rep(c(rep("Accumulated", times = t/4), rep("Averaged", times = t/4)), times = 2))
scales <- c(rep(c(rep("plot", times = t/12),rep("grid", times = t/12),rep("site", times = t/12)),times = 4))
div_type <- c(rep(c("Diversity (SAR)","Realized (PAR)","Potential (PAR)"), times = t/3))
predicted <- c(12/16, 6/16, 7/16, 1, 1, 1, 1, 1, 1,
               12/16, 6/16, 7/16, 1, 1, 1, 1, 1, 1,
               12/16, 6/16, 7/16, 1, 1, 1, 1, 1, 1,
               12/16, 6/16, 7/16, 1, 0, 0, 1, 0, 0)

# Notes on structure of simulated data: 

# a accumulated plot diverity = 12/16
# a accumulated plot realized = 6/16
# a accumulated plot potential = 7/16

# a accumulated grid diverity = 1/1
# a accumulated grid realized = 1/1
# a accumulated grid potential = 1/1

# a accumulated site diverity = 1/1
# a accumulated site realized = 1/1
# a accumulated site potential = 1/1

# a averaged plot diverity = 12/16
# a averaged plot realized =6/16
# a averaged plot potential = 7/16

# a averaged grid diverity = 1/1
# a averaged grid realized = 1/1
# a averaged grid potential = 1/1

# a averaged site diverity = 1/1
# a averaged site realized = 1/1
# a averaged site potential = 1/1

#-----------

# b accumulated plot diverity = 12/16
# b accumulated plot realized = 6/16
# b accumulated plot potential = 7/16

# b accumulated grid diverity = 1
# b accumulated grid realized = 1
# b accumulated grid potential = 1

# b accumulated site diverity = 1
# b accumulated site realized = 1
# b accumulated site potential = 1

# b averaged plot diverity = 12/16
# b averaged plot realized =6/16
# b averaged plot potential = 7/16

# b averaged grid diverity = 0
# b averaged grid realized =0
# b averaged grid potential =0

# b averaged site diverity = 0
# b averaged site realized =0
# b averaged site potential =0

dat <- data.frame(scenario, scale_type,scales,div_type,predicted)

dat$scales <- factor(dat$scales, c("plot","grid","site"))

par_sim <- ggplot(dat, aes(x = scales, y = predicted, color = div_type, group = div_type)) + 
  geom_point(position = dodge,inherit.aes = TRUE, size = 3) +
  geom_line(aes(color = div_type, group = div_type), position = dodge, linewidth = 1) +
  ylim(0, 1) +
  scale_color_manual(values = c("Diversity (SAR)" = "blue", 
                                "Potential (PAR)" = "orchid",
                                "Realized (PAR)" = "orange"
  )) +
  guides(color=guide_legend(override.aes=list(shape=1), title = "", title.position = "left", direction = "horizontal"))+ 
  guides(theme(legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
               legend.text=element_text(color = "grey20",size = 16, angle = 0, hjust = 0, face = "plain"))) +
  labs(x = "", y = "Proportion of populations", color = "", linetype = "", title = "") +
  facet_grid(cols = vars(scale_type), rows = vars(scenario)) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size = 16))
par_sim

#---------------------------------------------------------------------------
library(patchwork)
png("Figures/discrete_SAR-PAR_full.png", height = 10, width = 7, units = "in", res=800)
par_sim / par_dat + plot_layout (heights = c (1.7, 1))
dev.off()
#---------------------------------------------------------------------------

