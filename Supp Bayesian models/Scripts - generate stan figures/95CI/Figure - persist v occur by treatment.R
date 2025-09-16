
#---------------------------------------------------------------------
# DESCRIPTION: Figures S12 = Figure 2 re-made using Bayesian estimates
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

# source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Source - MAIN fitnessdata POST STAN.R"))  # loaded in following script:
source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats - persistence & occurrence by treatment across discrete scales.R"))

#################################
# TABLES
#################################

# Estimates
tab_occ.v.per_estimates <- summary |>
  dplyr::mutate(conf.int =paste(round(prediction_lwr,2), round(prediction_upr,1), sep = ", ")) |>
  dplyr::mutate(predicted = round(prediction,2)) |>
  dplyr::select(scale, treatment, predicted, conf.int) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    treatment = "Data type",
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
gt::gtsave(tab_occ.v.per_estimates, filename = gt_file) 
pagedown::chrome_print(input = gt_file,
                       output = here::here("Tables/y_hat/95CI/3_STAN_tab_occ.v.per_estimates.pdf"))

####################
# Figure
####################

###################
# combine datasets
fig_fig_2A$scale <- "plot"

fig_fig_2C$scale <- "grid"
fig_fig_2C$block <- NA
  
fig_fig_2E$scale <- "site"
fig_fig_2E$block <- NA
fig_fig_2E$grid <- NA
  
fig_fig <- rbind(fig_fig_2A,fig_fig_2C,fig_fig_2E)
fig_fig$scale <- as.factor(fig_fig$scale)
fig_fig$scale <- factor(fig_fig$scale,
                        levels = c("plot","grid","site"))

# use summary from above
names(summary) <- c("treatment","prop", "conf.low", "conf.high", "scale" )
summary$conf.low[summary$conf.low == 0] <- 1
#summary$conf.low[round(summary$conf.low,4) == 0.0000] <- 1 # where there is no variation remove confidence intervals

###############
# All scales

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
jpeg('Figures/Bayes/95CI/fig_occ_persist_by_treatment.jpeg', width = 10, height = 5, units = 'in', res = 600)
figa
dev.off()
#---------------------------------------------------------------------------


