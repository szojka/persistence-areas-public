
#---------------------------------------------------------------------
# DESCRIPTION: Figures S11B = Figure 1d re-made using Bayesian estimates
# that control for biases and use 95% credible intervals
#----------------------------------------------------------------------

{
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
  # library(here)
  
source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats -  natural misalignments (multinomial).R"))
source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Stats - persistence & occurrence by treatment across discrete scales.R"))
  

contingency_labs <- c( "i. aligned \npresent", "ii. sink","iii. dispersal \nlimitation", "iv. aligned \nabsent")

#-----------------------------------------------------------------------------.
# Figure 1C ####
# proportion of populations occurring, vs proportion of populations persisting (note: won't add to 1)

fig_fig_2A$treatment <- as.factor(fig_fig_2A$treatment )
levels(fig_fig_2A$treatment)
dat_fig <- fig_fig_2A %>%
  dplyr::filter(treatment %in% c("occupancy","suitability with neighbors" ))
dat_fig$treatment <- case_match(dat_fig$treatment, # make the name fit
                                "suitability with neighbors" ~ "suitability \n with neighbors",
                                "occupancy"~ "occupancy")
summary$scale <- as.factor(summary$scale)
dat_mod <- summary %>%
  dplyr::filter(treatment %in% c("occupancy","suitability with neighbors" ) &
                  scale %in% "plot")
dat_mod$treatment <- case_match(dat_mod$treatment, # make the name fit
                                "suitability with neighbors" ~ "suitability \n with neighbors",
                                "occupancy"~ "occupancy")

figa <- ggplot() + # should be prop for raw data
  geom_jitter(data = dat_fig, aes(x = treatment, y = prop), color = "deepskyblue4", alpha = 0.1, size = 1, width = 0.2, height = 0.03) + 
  theme_bw() +
  geom_point(data = dat_mod, aes(x = treatment, y = prediction), color = "deepskyblue4", size = 4) +
  geom_linerange(data = dat_mod, aes(x = treatment, y = prediction, ymin = prediction_lwr, ymax = prediction_upr), color = "deepskyblue4", linewidth = 1)+
  theme(text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) +
  labs( x = "", y = "Proportion of species")+
  ylim(0,1)
figa

#-----------------------------------------------------------------------------.
# Figure 1D ####
# contingency proportions in natural conditions (biotic plots)

# Order contingincies to mirror the table:
plot_dat <- fig1B_dat %>% # need to make contingency longer
  pivot_longer(., cols = 5:8, names_to = 'contingency', values_to = 'prop') %>%
  mutate(prop = prop/replicates)
output1_ordered <- as.data.frame(output_nat)
output1_ordered$contingency <- as.factor(output1_ordered$contingency)
output1_ordered$contingency <- factor(output1_ordered$contingency, levels = c(
  "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
levels(output1_ordered$contingency)
output1_ordered <- arrange(output1_ordered, match(output1_ordered$contingency, levels(output1_ordered$contingency)))

plot_dat$contingency <- as.factor(plot_dat$contingency)
plot_dat$contingency <- factor(plot_dat$contingency, levels = c(
  "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
levels(plot_dat$contingency)

plot_dat <- arrange(plot_dat, match(plot_dat$contingency, levels(plot_dat$contingency)))

man_list <- unique(plot_dat$contingency)

figb <- ggplot() +
  ylim(0,1) +
  labs(x="", y = "") +
  geom_jitter(data = plot_dat, aes(x = contingency, y = prop), 
              alpha = 0.1, 
              size = 1, 
              width = 0.2, 
              height = 0.03, 
              color = "mediumpurple1") +
  geom_point(data = output1_ordered, aes(x = contingency, y = prediction),size = 4, color = "mediumpurple3") +
  geom_linerange(data = output1_ordered, aes(x = contingency, y = prediction, ymin = prediction_lwr, ymax = prediction_upr), 
                linewidth = 1,
                color = "mediumpurple3") +
  geom_hline(yintercept = 0.25, linetype='dashed', col = 'grey') +
  theme_bw() +
  theme(text = element_text(size = 16),
        legend.position = 'none',
        axis.text.x = element_text (angle = 45, vjust = 1, hjust=1)) +
  annotate(xmin = which(man_list=="SS_y")-0.5, xmax = which(man_list=="SS_y")+0.5,
           ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
  annotate(xmin = which(man_list=="SS_n")-0.5, xmax = which(man_list=="SS_n")+0.5,
           ymin = -Inf, ymax = Inf, geom = 'rect', alpha = 0.2) +
  scale_x_discrete(labels = contingency_labs)
figb
}

#-----------------------------------------------------------------------------.
# graph
jpeg('Figures/Bayes/95CI/fig_1.jpeg', width = 8.5, height = 4.5, units = 'in', res = 300)
figa + figb  + plot_layout(widths = c(1,2))
dev.off()
#-----------------------------------------------------------------------------.


##############################
# TABLE
##############################

# get the order and names of misalignments correct:
output_nat$contingency <- as.factor(output_nat$contingency)
output_nat$contingency <- case_match(output_nat$contingency,
                                  'SS_y' ~ 'i. aligned present',
                                  'ME' ~ 'ii. sink',
                                  'DL' ~ 'iii. dispersal limitation',
                                  'SS_n' ~ 'iv. aligned absent')
output_nat <- output_nat %>%
  arrange(contingency)

tab_mulinomial_estimates <- output_nat |>
  dplyr::mutate(conf.int =paste(round(prediction_lwr,2), round(prediction_upr,2), sep = ", ")) |>
  dplyr::mutate(prob = round(prediction,2)) |>
  dplyr::select(contingency, prob, conf.int) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    contingency = "(Mis)alignment",
    prob = "Estimated proportion",
    conf.int =  "95% CI") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

# chrome update means I now save tables as pdfs using pagedown:
gt_file <- tempfile(fileext = ".html")
gt::gtsave(tab_mulinomial_estimates, filename = gt_file) 
pagedown::chrome_print(input = gt_file,
                       output = here::here("Tables/Bayes/95CI/1_STAN_tab_multinomial_estimates.pdf"))

