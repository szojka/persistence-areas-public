
#---------------------------------------------------------------------
# DESCRIPTION: Figures S13 = Figure 3 re-made using Bayesian estimates
# that control for biases and use 95% credible intervales
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

source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats - misalignments by treatment at blocks.R"))
source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats - misalignments by treatment at grids.R"))
source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats - misalignments by treatment at sites.R")) # warnings here are ok

# new figure with scale on the x-axis, and proportion still on the y, facet by a and b

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
output_all <- rbind(output_b, output_g, output_s)

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

output_all$prediction_upr[output_all$prediction_upr  > 1] <- 1
output_all$prediction  <- round(output_all$prediction ,2)
output_all$prediction_lwr  <- round(output_all$prediction_lwr ,2)
output_all$prediction_upr  <- round(output_all$prediction_upr ,2)

# FIGURE FOR ALL

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
  geom_jitter(aes(x = scale, y = prediction, group = treatment, color = treatment),
              position = position_dodge(width = 0.8), 
              size = 4,#col = col_treat_long
              ) +
  geom_line(aes(x = scale, y = prediction, group = treatment, color = treatment),
            position = position_dodge(width = 0.8), 
            linewidth = 1, 
            alpha = 0.7#,col = col_treat_mod
            ) +
  geom_linerange(aes(x = scale, y = prediction, ymin = prediction_lwr, ymax = prediction_upr, color = treatment),
            position = position_dodge(width = 0.8),#col = col_treat_mod,
            linewidth = 1) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'top',
        text = element_text(size = 16)) +
  scale_color_manual(values = c( "honeydew4","mediumpurple1"))
fig2_points  

jpeg('Figures/Bayes/95CI/fig_2_scaleonx.jpeg', width = 8.5, height = 7, units = 'in', res = 600)
fig2_points
dev.off()

##############################
# TABLE OUTPUT
###############################

# Estimates will be easy to table using output_all
# Anovas found in each associated Stats script 

# make W lower case for consistency within tables 
# output_all$treatment<-mapvalues(output_all$treatment, from=c("A","B"),
#                                 to=c("without neighbors","with neighbors"))
# output_all$contingency<-mapvalues(output_all$contingency, from=c("ME","SS_n","SS_y", "DL"),
#                                 to=c("sinks","aligned absent","aligned present", "dispersal limitation"))

paste(round(output_all$prediction_lwr,2), round(output_all$prediction_upr,1), sep = ", ")

tab_neighbor_estimates <- output_all |>
  dplyr::arrange(contingency,scale) |>
  dplyr::mutate(conf.int =paste(round(prediction_lwr,2), round(prediction_upr,1), sep = ", ")) |>
  dplyr::mutate(predicted = round(prediction,2)) |>
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

# chrome update means I now save tables as pdfs using pagedown:
gt_file <- tempfile(fileext = ".html")
gt::gtsave(tab_neighbor_estimates, filename = gt_file) 
pagedown::chrome_print(input = gt_file,
                       output = here::here("Tables/Bayes/95CI/6_STAN_tab_neighbor_estimates.pdf"))
