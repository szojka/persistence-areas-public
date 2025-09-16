
#-----------------------------------------------------------------------------
# DESCRIPTION:
# Supplemental material to show histogram oftransplant seed production 
# for each species, in conditions with and without neighbors
#-----------------------------------------------------------------------------

##################################
# RAW SEED PRODUCTION PER SPECIES
##################################

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

source(paste0(here::here(),"/Scripts/Source - MAIN fitnessdata.R"))


####################
# Figure S6
####################

seed_dat <- plotlev

# filter out germination == 0
seed_dat$status <- as.factor(seed_dat$status) # "no fruit"
levels(seed_dat$status) 

# filter for normal > 0 and no fruit and herbivory (i.e., situations where plant germinated)
seed_dat_germ1 <- seed_dat %>%
  filter(status %in% c("normal") & seed > 0)
seed_dat_germ2 <- seed_dat %>%
  filter(status %in% c("no fruit", "herbivory"))
seed_dat_germ <- rbind(seed_dat_germ1, seed_dat_germ2)   

seed_dat_germ$treatment <- case_match(seed_dat_germ$treatment, 
                               'A' ~ 'without neighbors',
                               'B' ~ 'with neighbors')

seed_dat_germ$species <- case_match(seed_dat_germ$species, 
                                 'plaere' ~ 'Plantago',
                                 'brohor' ~ 'Bromus',
                               'vulmic' ~ 'Festuca',
                               'miccal' ~ 'Micropus')

fig_seed <- ggplot(seed_dat_germ, aes(x = seed, fill = treatment, color = treatment)) +
  geom_histogram(alpha = 0.75, bins = 60) +
  scale_fill_manual(values = c( "mediumpurple1","honeydew4")) +
  scale_color_manual(values = c("mediumpurple1","honeydew4"), guide = FALSE) +
  labs(fill = "", x = "Seed production of germinated transplants", y = 'Count') +
  theme_bw() +
  facet_wrap(~species) +
  theme(legend.position = 'top',
        text = element_text(size = 16)) +
  geom_vline(xintercept = 2, color = 'red')
fig_seed

#######################
# determine how much of the data is cut off after 20 seeds on the x axis:
tot_n <- seed_dat_germ %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::select(species, freq) %>%
  distinct()
filt_n <- seed_dat_germ %>%
  filter(seed <= 20) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::select(species, freq) %>%
  distinct()
perc_p <- filt_n[,2]/tot_n[,2]

1 - perc_p  
# 1 0.13074205 Plantago
#2 0.02325581 Micropus
#3 0.02422145 Bromus
#4 0.11869436 Festuca

dat_text <- data.frame(
  label = c("2.4% of data omitted", "11.9% of data omitted", "2.3% of data omitted", "13.1% of data omitted"),
  species   = c('Bromus', 'Festuca', 'Micropus', 'Plantago'),
  treatment = 'without neighbors'
)

# labels for each facet
facet_labels <- data.frame(
  species = c("Bromus", "Festuca", "Micropus", "Plantago"),
  label = c("A", "B", "C", "D") 
)

# zoomed in:
fig_seed_zoomed <- ggplot(seed_dat_germ, aes(x = seed, fill = treatment, color = treatment)) +
  geom_bar(alpha = 0.75, just = 0) + # just 0 means 1 is to the right
  scale_fill_manual(values = c( "mediumpurple1","honeydew4")) +
  scale_color_manual(values = c("mediumpurple1","honeydew4"), guide = FALSE) +
  labs(fill = "", x = "Seed production of germinated transplants", y = 'Count') +
  theme_bw() +
  facet_wrap(~species) +
  theme(legend.position = 'top',
        text = element_text(size = 16),
        plot.title = element_text(size = 16)) +
  geom_vline(xintercept = 2, color = 'red') +
  xlim(-1, 20) + 
  geom_text(
    data    = dat_text,
    mapping = aes(x = 13, y = 80, label = label)
  ) +
  geom_text(data = facet_labels, aes(x = 1.5, y = Inf, label = label), 
            hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold", 
            inherit.aes = FALSE)
fig_seed_zoomed

#####################################################################
pdf('Figures/supp_fig_seedprod_zoomed.pdf', width = 7, height = 6)
fig_seed_zoomed
dev.off()
#####################################################################

#---------------------------------
# supporting stats
#---------------------------------

# What is the max seed production observed per species per treatment?
seed_dat_germ %>%
  dplyr::group_by(species, treatment) %>%
  dplyr::mutate(max_seed = max(seed)) %>%
  dplyr::mutate(mean_seed = mean(seed)) %>%
  dplyr::mutate(sd_seed = sd(seed)) %>%
  dplyr::mutate(n_seed = n()) %>%
  dplyr::mutate(lower_ci_seed = mean_seed-1.96*(sd_seed/sqrt(n_seed))) %>%
  dplyr::mutate(upper_ci_seed = mean_seed+1.96*(sd_seed/sqrt(n_seed))) %>%
  dplyr::select(species, treatment, max_seed, mean_seed, lower_ci_seed, upper_ci_seed) %>%
  distinct()

################################################
# overall mean seed production (not by species)
seed_dat_germ %>%
  dplyr::group_by(treatment) %>%
  dplyr::mutate(max_seed = max(seed)) %>%
  dplyr::mutate(mean_seed = mean(seed)) %>%
  dplyr::select(treatment, max_seed, mean_seed) %>%
  distinct()

#######################################################################
# Figure S5 - histogram of seed production per sown seed
#######################################################################
 
load(here::here("Data/dat_final_doubled.Rdata"))

dat_final_doubled$species <- case_match(dat_final_doubled$species, 
                                    'plaere' ~ 'Plantago',
                                    'brohor' ~ 'Bromus',
                                    'vulmic' ~ 'Festuca',
                                    'miccal' ~ 'Micropus')

dat_final_doubled$biotic_treatment <- case_match(dat_final_doubled$biotic_treatment, 
                                      'A' ~ 'without neighbors',
                                      'B' ~ 'with neighbors')

names(dat_final_doubled)

# labels for each facet
facet_labels <- data.frame(
  species = c("Bromus", "Festuca", "Micropus", "Plantago"),
  label = c("A", "B", "C", "D")  # Adjust labels as necessary
)

#######################
# how much of the data do we cut off after 20 seeds?
tot_n <- dat_final_doubled %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::select(species, freq) %>%
  distinct()
filt_n <- dat_final_doubled %>%
  filter(per_capita_seed_prod <= 20) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::select(species, freq) %>%
  distinct()
perc_p <- filt_n[,2]/tot_n[,2]

(1 - perc_p) *100 
  
# 1 3.1077348 Plantago
# 2 0.1351351 Micropus
# 3 0.4943503 Bromus
# 4 2.7009223 Festuca

dat_text <- data.frame(
  label = c("3.1% of data omitted", "0.14% of data omitted", "0.49% of data omitted", "2.7% of data omitted"),
  species   = c('Plantago', 'Micropus', 'Bromus', 'Festuca'),
  biotic_treatment = 'without neighbors'
)

dat_final_doubled %>%
  filter(per_capita_seed_prod == 0) %>%
  group_by(species) %>%
  mutate(n = n()) %>%
  select(species,n) %>%
  distinct

fig_seed_doubled <- ggplot(dat_final_doubled, aes(x = per_capita_seed_prod, fill = biotic_treatment, color = biotic_treatment)) +
  geom_bar(alpha = 0.75, just = 0) + # just 0 means 1 is to the right
  scale_fill_manual(values = c( "mediumpurple1","honeydew4")) +
  scale_color_manual(values = c("mediumpurple1","honeydew4"), guide = FALSE) +
  labs(fill = "", x = "Seed production of each sown seed", y = 'Count') +
  theme_bw() +
  facet_wrap(~species) +
  theme(legend.position = 'top',
        text = element_text(size = 16)) +
  geom_vline(xintercept = 1, color = 'red') +
  geom_text(
    data    = dat_text,
    mapping = aes(x = 13, y = 1500, label = label)
  ) +
  geom_text(data = facet_labels, aes(x = 1.5, y = Inf, label = label), 
            hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold", 
            inherit.aes = FALSE) +
   xlim(-1, 20) 
fig_seed_doubled

pdf(here::here('Figures/supp_fig_seedprod_doubled.pdf'), width = 7, height = 6)
fig_seed_doubled
dev.off()
