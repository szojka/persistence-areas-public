
#-----------------------------------------------------------------------------
# DESCRIPTION:
# Supplementary material to see how effectively we 
# reduced the cover in cleared plots in comparison 
# to the intact paired plot. 
#-----------------------------------------------------------------------------


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

############################
# download read in csv version of data

cover1 <- read_csv(paste0(here::here(),"/Data/Plot photo cover estimates 2023.csv"))

# rename some things
names(cover1)[1] <- "plot"
names(cover1)[2] <- 'image'
names(cover1)[3] <- "treatment"
cover1 <- cover1[,-8]

cover1$plot <- as.factor(cover1$plot)
cover1$treatment <- as.factor(cover1$treatment)
cover1$image <- as.factor(cover1$image)

# summary stats:

cover <- cover1 %>%
  pivot_longer(., cols = 4:7, names_to = "category", values_to = 'cover')
cover$category <- as.factor(cover$category)

# magnitudes of difference in cover
cover %>%
  group_by(treatment, category) %>%
  dplyr::summarize(min = min(cover), max = max(cover))

# total cover includes litter this year, last year, and live biomass
total_cover <- cover1
total_cover$total <- NA
total_cover$total <- rowSums(total_cover[5:7])

# Link total_cover with another key to get grid using image and plot.
key_for_cover <- plotlev %>%
  dplyr::select(image, grid, site, plot, names) %>%
  distinct()
key_for_cover$image <- as.factor(key_for_cover$image) 

# join
cover_green_all <- left_join(total_cover, key_for_cover, by = c("plot", "image"),
                             relationship = "many-to-many")
cover_green_all  

# to use beta family need 0 < y < 1
cover_green_all$prop <- NA
cover_green_all$prop <- cover_green_all$total/100

cover_green_all$prop[cover_green_all$prop == 0] <- 0.01
cover_green_all$prop[cover_green_all$prop >= 1] <- 0.99

cover_green_all <- na.omit(cover_green_all) 

# match a to without neighbors and b to with neighbors
cover_green_all$treatment <- case_match(cover_green_all$treatment,
                                        "A" ~ "Without neighbors",
                                        "B" ~ "With neighbors")

# calculate reduction in cover at each plot:
cover_reduction <- cover_green_all %>%
  group_by(names) %>%
  select(plot, total, grid, site, names, treatment) %>%
  pivot_wider(names_from = 'treatment', values_from = 'total') %>%
  mutate(reduction = `With neighbors`-`Without neighbors`/ `With neighbors`)
  
######################################
# plot raw cover reduction in boxplot 

fig_r <- ggplot(data = cover_reduction, aes(x = grid, y = reduction)) +
  geom_boxplot(fill = "midnightblue", alpha = 0.6) +
  theme_classic() + 
  labs(x = "Grid number", y = "Cover reduced \n by clearing (%)", fill = "") + 
  ylim(0,100) +
  theme(text = element_text(size = 16),
        legend.position = 'top') +
  annotate("text", x = 1, y = 99, label = "B", size = 6, fontface = "bold")
fig_r

fig_c <- ggplot(data = cover_green_all, aes(x = grid, y = total, fill = treatment, color = treatment)) +
  geom_boxplot(alpha = 0.5) +
  theme_classic() + 
  labs(x = "", y = "Total cover (%)", fill = "") + 
  ylim(0,100) +
  scale_fill_manual(values = c( "mediumpurple1","honeydew4")) +
  scale_color_manual(values = c( "mediumpurple1","honeydew4"), guide = 'none') +
  theme(text = element_text(size = 16),
        legend.position = 'top')+
  annotate("text", x = 1, y = 99, label = "A", size = 6, fontface = "bold")
fig_c

##############################################################################
library(patchwork)

pdf('Figures/fig_supp_cover.pdf', width = 7, height = 9)
fig_c / fig_r
dev.off()
##############################################################################
