
#----------------------------------------------------------------------------
# FIGURE 1d - MISALIGNMENTS IN NATURAL CONDITIONS (WITH BIOTIC INTERACTIONS)
#----------------------------------------------------------------------------
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
  
 # source(paste0(here::here(), "/Scripts/Stats - natural persistence v occurrence.R")) defunct
  source(paste0(here::here(),"/Scripts/Stats -  natural misalignments (multinomial).R"))
  source("Scripts/Stats - persistence & occurrence by treatment across discrete scales.R")
  
  
# aesthetics for Figure 1
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
names(summary)[1] <- "treatment"
dat_mod <- summary %>%
  dplyr::filter(treatment %in% c("occupancy","suitability with neighbors" ) &
                  scale %in% "plot")
dat_mod$treatment <- case_match(dat_mod$treatment, # make the name fit
                                "suitability with neighbors" ~ "suitability \n with neighbors",
                                "occupancy"~ "occupancy")

figa <- ggplot() + # should be prop for raw data
  geom_jitter(data = dat_fig, aes(x = treatment, y = prop), color = "deepskyblue4", alpha = 0.1, size = 1, width = 0.2, height = 0.03) + 
  theme_bw() +
  geom_point(data = dat_mod, aes(x = treatment, y = predicted), color = "deepskyblue4", size = 4) +
  geom_linerange(data = dat_mod, aes(x = treatment, y = predicted, ymin = conf.low, ymax = conf.high), color = "deepskyblue4", linewidth = 1)+
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
fig1B_dat
output1_ordered <- as.data.frame(output1)
output1_ordered$contingency <- as.factor(output1_ordered$contingency)
output1_ordered$contingency <- factor(output1_ordered$contingency, levels = c(
  "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
levels(output1_ordered$contingency)
output1_ordered <- arrange(output1_ordered, match(output1_ordered$contingency, levels(output1_ordered$contingency)))

fig1B_dat$contingency <- as.factor(fig1B_dat$contingency)
fig1B_dat$contingency <- factor(fig1B_dat$contingency, levels = c(
  "SS_y", "ME", "DL","SS_n"), ordered = TRUE) 
levels(fig1B_dat$contingency)

fig1B_dat <- arrange(fig1B_dat, match(fig1B_dat$contingency, levels(fig1B_dat$contingency)))

fig1b_all <- left_join(fig1B_dat, output1_ordered, by = 'contingency')
levels(fig1b_all$contingency) 

col_treat_long <- c(fig1b_all$colortreat)
man_list <- unique(fig1b_all$contingency)

figb <- ggplot(fig1b_all) +
  ylim(0,1) +
  labs(x="", y = "") +
  geom_jitter(aes(x = contingency, y = prop), 
              alpha = 0.1, 
              size = 1, 
              width = 0.2, 
              height = 0.03, 
              color = "mediumpurple1") +
  geom_point(aes(x = contingency, y = prob),size = 4, color = "mediumpurple3") +
  geom_linerange(data = output1_ordered, aes(x = contingency, y = prob, ymin = ci.lwr, ymax = ci.upr), 
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
# save graph
jpeg('Figures/fig_1.jpeg', width = 8.5, height = 4.5, units = 'in', res = 300)
figa + figb  + plot_layout(widths = c(1,2))
dev.off()
#-----------------------------------------------------------------------------.


##############################
# TABLE
##############################

# get the order and names of misalignments correct:
output1$contingency <- as.factor(output1$contingency)
output1$contingency <- case_match(output1$contingency,
                                  'SS_y' ~ 'i. aligned present',
                                  'ME' ~ 'ii. sink',
                                  'DL' ~ 'iii. dispersal limitation',
                                  'SS_n' ~ 'iv. aligned absent')
output1 <- output1 %>%
  arrange(contingency)

tab_mulinomial_estimates <- output1 |>
  dplyr::mutate(conf.int =paste(round(ci.lwr,2), round(ci.upr,2), sep = ", ")) |>
  dplyr::mutate(prob = round(prob,2)) |>
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
tab_mulinomial_estimates |>
 gtsave(paste0(here::here(),"/Tables/1tab_multinomial_estimates.pdf")) 

tab_wald <- dat_compare |>
  dplyr::mutate(z_score = round(z_score,4)) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    p_val = "P-value",
    comparing = 'Pairwise comparison',
    z_score = "z-score") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_wald |>
 gtsave(paste0(here::here(),"/Tables/2tab_wald.pdf")) 






