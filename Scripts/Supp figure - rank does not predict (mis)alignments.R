
#------------------------------------------------------
# DESCRIPTION:
# Supplemental material to see whether the rank abundance 
# of focal species predicts the proportion of (mis)alignments 
# for each species with neighbors
#------------------------------------------------------

##################################
# SUPPLEMENTAL MATERIAL
##################################

# library(tidyverse)
# library(stringr)
# library(car)
# library(lme4)
# library(glmmTMB)
# library(visreg)
# library(patchwork)
# library(ggeffects)
# library(DHARMa)

source("Scripts/Source - MAIN fitnessdata.R")
source("Scripts/Supp figure - species specific misalignments.R")

###########################
# TOTAL ABUNDANCE DATA
###########################

# prepare plotlev abundances by name by species
plotlev$ab_cat <- as.factor(plotlev$ab_cat)
fig1_ab <- plotlev %>%
  filter(treatment %in% 'B') %>%
  select(names, species, occurrence, plot, grid, site) %>%
  mutate(presence = ifelse(occurrence == 'no', 0, 1))

# prepare species proportion of contingencies
fig1P_dat$species <- 'plaere'
fig1V_dat$species <- 'vulmic'
fig1M_dat$species <- 'miccal'
fig1Br_dat$species <- 'brohor'
#fig1Br_dat <- select(fig1Br_dat, -replicates)
spp_prop_dat <- rbind(fig1P_dat,fig1V_dat,fig1M_dat,fig1Br_dat)
spp_prop_dat$species <- as.factor(spp_prop_dat$species)

# link together plotlev and proportion of contingencies for each species
review_dat <- full_join(fig1_ab, spp_prop_dat, by = c('species','names'), relationship = 'many-to-many')

# rank abundance step
# to organize x axis in species from most abundant to least:
review_dat <- na.omit(review_dat)
review_dat %>%
  dplyr::select(presence, species) %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(sum_presence = sum(presence))

rank <- c('Festuca', 'Micropus', 'Plantago', 'Bromus') # most to least presences

# pull out each contingency and model:
mr1 <- glmmTMB(prop ~ species*contingency + (1|site/grid/plot), 
               data = review_dat,
               family = binomial()
               )

vis_mr1 <- ggpredict(mr1,
                     terms = c("species", "contingency"),
                     type = 'fixed'); plot(vis_mr1)

vis_mr1$x <- case_match(vis_mr1$x, #species
                     "brohor" ~ "Bromus",
                     "miccal" ~ "Micropus",
                     "plaere" ~ "Plantago",
                     "vulmic" ~ "Festuca")
vis_mr1$x <- factor(vis_mr1$x, levels = c(rank)) # by rank
vis_mr1$group <- case_match(vis_mr1$group, #contingency
                            "SS_y" ~ "Aligned present",
                            "ME" ~ "Sink",
                            "DL" ~ "Dispersal limitation",
                            "SS_n" ~ "Aligned absent")
vis_mr1$group <- factor(vis_mr1$group, levels = c("Aligned present", "Sink","Dispersal limitation", "Aligned absent"))

library(viridis)

facet_labels <- data.frame(
  group = unique(vis_mr1$group),
  label = c("A", "B", "C", "D")  # Adjust labels as necessary
)

review_fig <- ggplot(data = vis_mr1, mapping = aes(x = x, y = predicted, color = x)) +
  geom_line(color = "lightgrey", aes(group = group)) +
  geom_point(size = 3) + 
  geom_linerange(aes(y = predicted, ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~group) +
  labs(x = "Species by rank abundance (high to low)", y = "Proportion of each (mis)alignment", color = "") +
  scale_color_viridis(option = "H", direction = 1, discrete = TRUE,
                      labels = c("Bromus","Plantago", "Micropus","Festuca")) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16)) +
  geom_text(data = facet_labels, aes(x = 1, y = Inf, label = label), 
            hjust = 1.5, vjust = 1.5, size = 6, fontface = "bold", 
            inherit.aes = FALSE)
review_fig

###################################################################
pdf('Figures/fig_supp_review.pdf', width = 6, height = 5)
review_fig
dev.off()
###################################################################
