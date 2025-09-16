
#--------------------------------------------------
# DESCRIPTION: Figure 3 comparing the effect of neighbors 
# on the proportion of species misaligned across scales
#--------------------------------------------------

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

source(paste0(here::here(),"/Scripts/Stats - misalignments by treatment at blocks.R"))
source(paste0(here::here(),"/Scripts/Stats - misalignments by treatment at grids.R"))
source(paste0(here::here(),"/Scripts/Stats - misalignments by treatment at sites.R")) # warnings here are ok

# combine data so that scale is new column, and other variables are treatment, replicates, contingency, prop, colortreat

fig2B_dat_new <-  fig2B_dat %>%
  select(treatment, replicates, contingency, prop) %>%
  mutate(scale = 'plot')
fig2D_dat_new <- fig2D_dat %>%
  select(treatment, replicates, contingency, prop) %>%
  mutate(scale = 'grid')
fig2F_dat_new <- fig2F_dat %>%
  select(treatment, replicates, contingency, prop) %>%
  mutate(scale = 'site')

scales_figdat <- rbind(fig2B_dat_new,fig2D_dat_new,fig2F_dat_new)
scales_figdat$scale <- as.factor(scales_figdat$scale)
scales_figdat$scale <- factor(scales_figdat$scale,
                              levels = c('plot','grid','site'))
levels(scales_figdat$scale)

# new objects that replace output = 
output_all <- rbind(vis.me.r,vis.dl.r,vis.ap.r,vis.aa.r,
      vis.me.g,vis.dl.g,vis.ap.g,vis.aa.g,
      vis.me.p,vis.dl.p,vis.ap.p,vis.aa.p)

output_all <- as.data.frame(output_all)

# order factor levels
output_all$scale <- as.factor(output_all$scale)
output_all$scale <- factor(output_all$scale,
                              levels = c('plot','grid','site'))
levels(output_all$scale)
names(output_all)[1] <- "treatment"

# Make labels what I want
output_all$treatment<-mapvalues(output_all$treatment, from=c("A","B"),
                                to=c("Without neighbors","With neighbors"))
scales_figdat$treatment<-mapvalues(scales_figdat$treatment, from=c("A","B"),
                                   to=c("Without neighbors","With neighbors"))
output_all$contingency<-mapvalues(output_all$contingency, from=c("SS_y","ME","DL","SS_n"),
                                  to=c("i. aligned present","ii. sink","iii. dispersal limitation","iv. aligned absent"))
scales_figdat$contingency<-mapvalues(scales_figdat$contingency, from=c("SS_y","ME","DL","SS_n"),
                                     to=c("i. aligned present","ii. sink","iii. dispersal limitation","iv. aligned absent"))

# control order of contingencies
output_all$contingency <- factor(output_all$contingency,
                           levels = c("i. aligned present","ii. sink","iii. dispersal limitation","iv. aligned absent"))
scales_figdat$contingency <- factor(scales_figdat$contingency,
                           levels = c("i. aligned present","ii. sink","iii. dispersal limitation","iv. aligned absent"))

output_all$conf.high[output_all$conf.high > 1] <- 1
output_all$predicted <- round(output_all$predicted,2)
output_all$conf.low <- round(output_all$conf.low,2)
output_all$conf.high <- round(output_all$conf.high,2)

# change conf.high from 1 to 0 for dispersal limitation grid, as 0 occurrences
output_all$conf.high[output_all$scale %in% "grid" & output_all$contingency %in% "iii. dispersal limitation" & output_all$treatment %in% "With neighbors"] <- 0 
  
# FIGURE

# COMBINE output_all and scales_figdat
scales_figdat_all <- left_join(scales_figdat, output_all, by = c("contingency", "treatment", "scale"))

fig2_points <- ggplot(scales_figdat_all) + 
  ylim(0,1) +
  labs(x="Spatial scale", y = "Proportion of species") +
  facet_wrap(~contingency) +
  geom_jitter(aes(x = scale, y = prop, color = treatment, shape = treatment), 
              position = position_jitterdodge(jitter.width = 0.6), 
              #col = col_treat_long,
              alpha = 0.1, 
              size = 1.5) +
  geom_jitter(aes(x = scale, y = predicted, group = treatment, color = treatment),
              position = position_dodge(width = 0.8), 
              size = 4,#col = col_treat_long
              ) +
  geom_line(aes(x = scale, y = predicted, group = treatment, color = treatment),
            position = position_dodge(width = 0.8), 
            linewidth = 1, 
            alpha = 0.7#,col = col_treat_mod
            ) +
  geom_linerange(aes(x = scale, y = predicted, ymin = conf.high, ymax = conf.low, color = treatment),
            position = position_dodge(width = 0.8),#col = col_treat_mod,
            linewidth = 1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'top',
        text = element_text(size = 16)) +
  scale_color_manual(values = c( "honeydew4","mediumpurple1"))
fig2_points  

jpeg('Figures/fig_2_scaleonx.jpeg', width = 8.5, height = 7, units = 'in', res = 600)
fig2_points
dev.off()

pdf('Figures/fig_2_scaleonx.pdf', width = 8.5, height = 7)
fig2_points
dev.off()

##############################
# TABLE OUTPUT
###############################

# Estimates will be easy to table using output_all
# Anovas found in each associated Stats script 

# make W lower case for consistency within tables 
output_all$treatment<-mapvalues(output_all$treatment, from=c("Without neighbors","With neighbors"),
                                to=c("without neighbors","with neighbors"))

paste(round(output_all$conf.low,2), round(output_all$conf.high,1), sep = ", ")

tab_neighbor_estimates <- output_all |>
  dplyr::arrange(contingency,scale) |>
  dplyr::mutate(conf.int =paste(round(conf.low,2), round(conf.high,1), sep = ", ")) |>
  dplyr::mutate(predicted = round(predicted,2)) |>
  dplyr::select(contingency, scale, treatment, predicted, conf.int) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    contingency = "(Mis)alignment",
   treatment = "Neighbor treatment",
    predicted = 'Estimated proportion',
    conf.int = "95% CI",
    scale = 'Scale level') |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_neighbor_estimates |>
  gtsave(paste0(here::here(),"/Tables/6tab_neighbor_estimates.pdf")) 

###############################
# ANOVA and CONTRASTS

# tab1. ANOVA

anovas <- rbind(plot.anova, grid.anova, site.anova)

# get the order and names of misalignments correct:
anovas$contingency <- as.factor(anovas$contingency)
levels(anovas$contingency)

anovas$contingency <- case_match(anovas$contingency,
                                     'aligned present' ~ 'i. aligned present',
                                     'sink' ~ 'ii. sink',
                                     'sinks' ~ 'ii. sink',
                                     'dispersal limitation' ~ 'iii. dispersal limitation',
                                     'aligned absent' ~ 'iv. aligned absent')

# order plot -> site
anovas$scale <- factor(anovas$scale, levels = c(
  "plot", "grid", "site"), ordered = TRUE) 
levels(anovas$scale)

anovas <- anovas %>%
 arrange(scale, contingency)

tab_anovas <- anovas |>
  dplyr::select(scale, contingency, predictor, Chi.squared, Df, P_value) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    predictor = 'Predictor',
    scale = 'Scale level',
    P_value = "P-value",
    contingency = '(Mis)alignment') |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_anovas |>
  gtsave(paste0(here::here(),"/Tables/7tab_neighbor_anova.pdf")) 

# tab2. Pairs significance

contrastz <- rbind(plot.contrast, grid.contrast, site.contrast)

# get the order and names of misalignments correct:
contrastz$contingency <- as.factor(contrastz$contingency)
levels(contrastz$contingency)

contrastz$contingency <- case_match(contrastz$contingency,
                                 'aligned present' ~ 'i. aligned present',
                                 'sink' ~ 'ii. sink',
                                 'dispersal limitation' ~ 'iii. dispersal limitation',
                                 'aligned absent' ~ 'iv. aligned absent')

# order plot -> site
contrastz$scale <- factor(contrastz$scale, levels = c(
  "plot", "grid", "site"), ordered = TRUE) 
levels(contrastz$scale)

contrastz <- contrastz %>%
  arrange(scale, contingency)

tab_contrasts <- contrastz |>
  dplyr::select(scale, contingency, z.ratio, p.value) |>
 # dplyr::mutate(odds.ratio = round(odds.ratio,3)) |>
  dplyr::mutate(z.ratio = round(z.ratio,3)) |>
  dplyr::mutate(p.value = round(p.value,3)) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    z.ratio = "z-ratio",
    p.value = "P-value",
    contingency = '(Mis)alignment',
    scale = "Scale level") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_contrasts |>
  gtsave(paste0(here::here(),"/Tables/8tab_neighbor_contrasts.pdf")) 
