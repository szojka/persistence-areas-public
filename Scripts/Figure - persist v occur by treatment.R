
#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 2: proportion of occurrence and persistence 
# with and without neighbors across scales
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

# source("Scripts/Source - MAIN fitnessdata.R") # loaded in following script:
source("Scripts/Stats - persistence & occurrence by treatment across discrete scales.R")

#################################
# TABLES
#################################

# Estimates
tab_occ.v.per_estimates <- summary |>
  dplyr::mutate(conf.int =paste(round(conf.low,2), round(conf.high,1), sep = ", ")) |>
  dplyr::mutate(predicted = round(predicted,2)) |>
  dplyr::select(scale, x, predicted, conf.int) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    x = "Data type",
    predicted = 'Estimated proportion',
    conf.int = "95% CI",
    scale = 'Scale level') |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_occ.v.per_estimates |>
  gtsave(paste0(here::here(),"/Tables/3tab_occ.v.per_estimates.pdf")) 

# Anovas

all.anova <- rbind(p.anova,g.anova,s.anova)

all_anovas <- all.anova |>
  dplyr::select(scale, predictor, Chi.squared, Df, P_value) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    predictor = 'Predictor',
    scale = 'Scale level',
    P_value = "P-value"
  ) |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

all_anovas |>
  gtsave(paste0(here::here(),"/Tables/4tab_occ.v.per_anova.pdf")) 

# pairwise contrasts

contrastz <- rbind(mplot.contrast, gplot.contrast, splot.contrast)

contrastz$scale <- NA
contrastz$scale <- c(rep('plot',times = 3),rep('grid', times = 3), rep('site',times = 3))

tab_occ.v.per_contrasts <- contrastz |>
  dplyr::select(scale, contrast, z.ratio, p.value) |>
  # dplyr::mutate(odds.ratio = round(odds.ratio,3)) |>
  dplyr::mutate(z.ratio = round(z.ratio,3)) |>
  dplyr::mutate(p.value = round(p.value,3)) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    p.value = "P-value",
    contrast = 'Pairwise comparison',
    scale = "Scale level",
    z.ratio = "z-ratio") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_occ.v.per_contrasts |>
  gtsave(paste0(here::here(),"/Tables/5tab_occ.v.per_contrasts.pdf")) 

####################
# Figure
####################

###################
# combine datasets
fig_fig_2A$scale <- "plot"

fig_fig_2C$scale <- "grid"
fig_fig_2C$names <- NA
  
fig_fig_2E$scale <- "site"
fig_fig_2E$names <- NA
fig_fig_2E$grid <- NA
  
fig_fig <- rbind(fig_fig_2A,fig_fig_2C,fig_fig_2E)
fig_fig$scale <- as.factor(fig_fig$scale)
fig_fig$scale <- factor(fig_fig$scale,
                        levels = c("plot","grid","site"))

# use summary from above
names(summary) <- c("treatment","prop","std.error", "conf.low", "conf.high", "group", "scale" )
summary$conf.low[summary$conf.low == 0] <- 1
summary$conf.low[round(summary$conf.low,4) == 0.0000] <- 1 # where there is no variation remove confidence intervals

###############
# All scales

#col_treat_long <- c(fig_fig_2A$colortreat)
figa <- ggplot(fig_fig, aes(x = scale, y = prop, fill = treatment, color = treatment)) +
  ylim(0,1) +
  geom_point(col = "deepskyblue4", alpha = 0.05, size = 1.5, position = position_jitterdodge(jitter.width = 0.8)) +
  geom_point(data = summary, aes(x = scale, y = prop), col = "deepskyblue4", size = 3) +
  geom_linerange(data = summary, aes(x = scale, y = prop, ymin = conf.low, ymax = conf.high), col = "deepskyblue4", linewidth = 1)+
  geom_line(data = summary, aes(x = scale, y = prop, group = treatment), col = "deepskyblue4", linewidth = 1) +
  theme_bw() +
  theme(plot.title = element_text(size = 16),
        text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs( x = "", y = "Proportion of species") +
  facet_wrap(~treatment)
figa

#---------------------------------------------------------------------------
jpeg('Figures/fig_occ_persist_by_treatment.jpeg', width = 10, height = 5, units = 'in', res = 600)
figa
dev.off()
#---------------------------------------------------------------------------


