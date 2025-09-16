#------------------------------------------------
# Supplemental figure for annual climate averages
# across 20 years at Knoxville weather statison, CA
#------------------------------------------------

library(tidyverse)

clim_dat <- read_csv("Data/climate_data.csv")
class(clim_dat$Year)

fig_temp <- ggplot(clim_dat, aes(x = Year, y = Temperature, group = 1)) +
  geom_point(color = 'darkred', size = 3) +
  geom_line(color = 'darkred', linewidth = 1)+
  theme_bw() +
  labs(x = "Year",y = "Annual temperature (C)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
  text = element_text(size = 16)) +
  ylim(10,20) +
  annotate("text", x = 1, y = 20, label = "A", size = 6, fontface = "bold") +
  scale_x_discrete(breaks = clim_dat$Year[seq(1, length(clim_dat$Year), by = 3)])
fig_temp

fig_precip <- ggplot(clim_dat, aes(x = Year, y = Precipitation, group = 1)) +
  geom_point(color = 'lightblue4', size = 3) +
  geom_line(color = 'lightblue4', linewidth = 1)+
  theme_bw() +
  labs(x = "Year",y = "Annual precipitation (cm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16)) +
  ylim(0,180) + 
  annotate("text", x = 1, y = 180, label = "B", size = 6, fontface = "bold") +
  scale_x_discrete(breaks = clim_dat$Year[seq(1, length(clim_dat$Year), by = 3)])
fig_precip

library(patchwork)

###################################################################################
pdf('Figures/fig_supp_enviro.pdf', width = 10, height = 5)
fig_temp + fig_precip
dev.off()
####################################################################################