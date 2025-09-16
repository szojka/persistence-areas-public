#####################################################
# Supplemental figures from Bayesian analyses
#####################################################

# Includes:

# 1 - difference graph between Beverton-Holt vs pop model (=unwanted species interactions), Bev-Holt vs empirical (=model choice and/or stochasiticity). Corresponds to Figure S15

# 2 - plotted estimated fecundity at three scales with 95% CI, corresponds to Figures S8-S10

###################
# Libraries
library(rstan)
library(tidyverse)

##############################################################
# load model objects to extract posteriors and compare them

source(here::here("Supp bayesian models/Scripts - generate stan figures/95CI/Stats - persistence & occurrence by treatment across discrete scales.R"))
head(summary)
post_popmodel <- summary
post_popmodel$model <- 'popmodel'

source(here::here("Supp bayesian models/Scripts - generate stan figures/FullBevHolt/95CI/Stats - persistence & occurrence by treatment across discrete scales.R"))
head(summary)
post_bevholt <- summary
post_bevholt$model <- 'bevholt'

# load original empirical data
source(here::here("Scripts/Stats - persistence & occurrence by treatment across discrete scales.R"))
original <- summary
original$model <- 'original'
# different titles here
original <- select(original, -std.error, -group)
names(original) <- c('treatment', "prediction", "prediction_lwr", "prediction_upr", "scale", "model")

######################################
# Create difference dataframe

dat_models <- rbind(post_popmodel, post_bevholt, original)

# expand predictions into their 'model name'

dat_models <- dat_models %>%
  select(-prediction_lwr, -prediction_upr) %>%
  filter(!treatment %in% ("occupancy")) %>%
  pivot_wider(., names_from = 'model', values_from = 'prediction') %>%
  mutate(colonizer = popmodel - bevholt) %>%
  mutate(stochasticity = bevholt - original) %>%
  mutate(`total difference` = original - popmodel) %>%
  select(treatment, scale, colonizer, stochasticity, `total difference`) %>%
  pivot_longer(cols = c(3:5), values_to = 'difference', names_to = 'process')
  
min(dat_models$difference)
max(dat_models$difference)

dat_models$scale <- factor(dat_models$scale, levels = c("plot", "grid", "site"))

# visualize differences

# Row one, facets A and B
fig_diff1a <- ggplot(filter(dat_models, process %in% 'total difference' & treatment %in% 'suitability with neighbors'), aes(x = scale, y = difference)) +
  geom_segment( aes(x=scale,y=0, yend=difference), color = "mediumpurple") +
  geom_point(size=2, alpha=0.7, shape=21, stroke=2, color = "mediumpurple") +
  geom_hline(yintercept = 0, linetype = "dashed", "darkgrey") +
  theme_bw()+
  ylab(bquote('Empirical' - 'Model'[' no conspecifics'])) +
  labs(x = "", color = "") +
=  ylim(-1,1) +
  theme(legend.position = 'top',
        text = element_text(size = 16)) +
  annotate("text", x = 2, y = .75, label = "Empirical analyses estimate a higher proportion of \n suitable habitat compared to population models.", color = "slategrey") +
  annotate("text", x = 2, y = -.75, label = "Empirical analyses estimate a lower proportion of \n suitable habitat compared to population models.", color = "slategrey")+
  annotate("text", x = 0.5, y = 1, label = "A", color = 'black', fontface = 'bold', size = 7)
fig_diff1a

fig_diff1b <- ggplot(filter(dat_models, process %in% 'total difference' & treatment %in% 'suitability without neighbors'), aes(x = scale, y = difference)) +
  geom_segment( aes(x=scale,y=0, yend=difference), color = "honeydew4") +
  geom_point(size=2, alpha=0.7, shape=21, stroke=2, color = "honeydew4") +
  geom_hline(yintercept = 0, linetype = "dashed", "darkgrey") +
  theme_bw()+
  ylab(bquote('Empirical' - 'Model'[' no neighbors'])) +
  labs(x = "", color = "") +
  ylim(-1,1) +
  theme(legend.position = 'top',
        text = element_text(size = 16)) +
  annotate("text", x = 2, y = .75, label = "Empirical analyses estimate a higher proportion of \n suitable habitat compared to population models.", color = "slategrey") +
  annotate("text", x = 2, y = -.75, label = "Empirical analyses estimate a lower proportion of \n suitable habitat compared to population models.", color = "slategrey")+
  annotate("text", x = 0.5, y = 1, label = "B", color = 'black', fontface = 'bold', size = 7)
fig_diff1b

# Row 2, Facets C and D
fig_diff2a <- ggplot(filter(dat_models, process %in% 'colonizer' & treatment %in% 'suitability with neighbors'), aes(x = scale, y = difference, color = treatment)) +
  geom_segment( aes(x=scale,y=0, yend=difference), color = 'mediumpurple') +
  geom_point(size=2, alpha=0.7, shape=21, stroke=2,color = 'mediumpurple') +
  geom_hline(yintercept = 0, linetype = "dashed", "darkgrey") +
  theme_bw()+
  ylab(bquote('Model'[' no conspecifics'] - 'Model'['full Beverton-Holt'])) +
  labs(x = "", color = "") +
  ylim(-1,1) +
  theme(legend.position = 'none',
        text = element_text(size = 16)) +
  annotate("text", x = 2, y = .75, label = "Conspecifics of focal species in 'with neighbors' plots \n competitively impact LDGR of transplants.", color = "slategrey") +
  annotate("text", x = 2, y = -.75, label = "Conspecifics of focal species in 'with neighbors' plots \n facilitatively impact LDGR of transplants.", color = "slategrey") +
  annotate("text", x = 0.5, y = 1, label = "C", color = 'black', fontface = 'bold', size = 7)
fig_diff2a

fig_diff2b <- ggplot(filter(dat_models, process %in% 'colonizer' & treatment %in% 'suitability without neighbors'), aes(x = scale, y = difference, color = treatment)) +
  geom_segment( aes(x=scale,y=0, yend=difference), color = 'honeydew4') +
  geom_point(size=2, alpha=0.7, shape=21, stroke=2,color = 'honeydew4') +
  geom_hline(yintercept = 0, linetype = "dashed", "darkgrey") +
  theme_bw()+
  # 'Model predictions'['no neighbors']
  ylab(bquote('Model'[' no neighbors'] - 'Model'['full Beverton-Holt'])) +
  labs(x = "",  color = "") +
  ylim(-1,1) +
  theme(legend.position = 'none',
        text = element_text(size = 16)) +
  annotate("text", x = 2, y = .75, label = "Colonizers of 'without neighbors' plots are net competitive, \n resulting in detection of fewer competitive effects \n on habitat suitability.", color = "slategrey") +
  annotate("text", x = 2, y = -.75, label = "Colonizers of 'without neighbors' plots are net facilitative, \n resulting in detection of fewer facilitative effects \n on habitat suitability.", color = "slategrey") +
  annotate("text", x = 0.5, y = 1, label = "D", color = 'black', fontface = 'bold', size = 7)
fig_diff2b

# Row 3, Facets E and F
fig_diff3a <- ggplot(filter(dat_models, process %in% 'stochasticity'  & treatment %in% 'suitability with neighbors'), aes(x = scale, y = difference, color = treatment)) +
  geom_segment( aes(x=scale,y=0, yend=difference), color = "mediumpurple") +
  geom_point(size=2, alpha=0.7, shape=21, stroke=2, color = "mediumpurple") +
  geom_hline(yintercept = 0, linetype = "dashed", "darkgrey") +
  theme_bw()+
  ylab(bquote('Model'[' full Beverton-Holt'] - 'Empirical')) +
  labs(x = "", color = "") +
  ylim(-1,1) +
  theme(legend.position = 'none',
        text = element_text(size = 16)) +
  annotate("text", x = 2, y = .75, label = "Modelling population dynamics (stochasticity & \n choice of Beverton-Holt) estimates a higher proportion \n of habitat suitability compared to empirical analyses.", color = "slategrey") +
  annotate("text", x = 2, y = -.75, label = "Modelling population dynamics estimates a lower proportion \n of habitat suitability compared to empirical analyses.", color = "slategrey") +
  annotate("text", x = 0.5, y = 1, label = "E", color = 'black', fontface = 'bold', size = 7)
fig_diff3a

fig_diff3b <- ggplot(filter(dat_models, process %in% 'stochasticity'  & treatment %in% 'suitability without neighbors'), aes(x = scale, y = difference, color = treatment)) +
  geom_segment( aes(x=scale,y=0, yend=difference), color = "honeydew4") +
  geom_point(size=2, alpha=0.7, shape=21, stroke=2, color = "honeydew4") +
  geom_hline(yintercept = 0, linetype = "dashed", "darkgrey") +
  theme_bw()+
  ylab(bquote('Model'[' full Beverton-Holt'] - 'Empirical')) +
  labs(x = "", color = "") +
  ylim(-1,1) +
  theme(legend.position = 'none',
        text = element_text(size = 16)) +
  annotate("text", x = 2, y = .75, label = "Modelling population dynamics (stochasticity & \n choice of Beverton-Holt) estimates a higher proportion \n of habitat suitability compared to empirical analyses.", color = "slategrey") +
  annotate("text", x = 2, y = -.75, label = "Modelling population dynamics estimates a lower proportion \n of habitat suitability compared to empirical analyses.", color = "slategrey") +
  annotate("text", x = 0.5, y = 1, label = "F", color = 'black', fontface = 'bold', size = 7)
fig_diff3b

#---------------------------------------------------------------------------
library(patchwork)
library(cowplot)
{
stacked_plots <- (fig_diff1a + fig_diff1b) / (fig_diff2a + fig_diff2b) / (fig_diff3a + fig_diff3b) + 
  theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))  # ensure some margin

# Add the shared y-axis label on the left using `plot_grid`
with_y_label <- plot_grid(
  ggdraw() + draw_label("Difference in proportion suitable habitat", 
                        angle = 90, size = 18, hjust = 0.5),
  stacked_plots,
  ncol = 2,
  rel_widths = c(0.05, 1)  # space for y label
)

# Add shared x-axis label underneath
combined_plot <- plot_grid(
  with_y_label,
  ggdraw() + draw_label("Spatial scale", size = 18, hjust = 0.5, x = 0.55),
  ncol = 1,
  rel_heights = c(1, 0.03)  # adjust height of label
)

legend_plot <- ggplot(data.frame(x=1:2, y=1:2, group=c("With neighbors", "Without neighbors")),
                      aes(x, y, color=group)) +
  geom_point() +
  scale_color_manual(values=c("mediumpurple", "honeydew4")) +
  theme(legend.position = "top",
        text = element_text(size = 16)) +
  labs(color = "Biotic treatment")
legend <- get_legend(legend_plot)

final_plot <- plot_grid(
  legend,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)
final_plot
}

#----------------------
pdf('Figures/supp_diff_predictions.pdf', width = 12, height = 14)
final_plot
dev.off()
#---------------------------------------------------------------------------

########################################################################
# Predictions from the pop model where we control for neighbors:

# in the 96CI folder, need to predictions from the objects plotlev, grilev, sitelev:
# predictions and predictions_lwr, predictions_upr

source(here::here("Supp bayesian models/Scripts - generate stan figures/95CI/Source - MAIN fitnessdata POST STAN.R"))

names(plotlev)

# check that these estimates are all > or equal to 0 
filter(plotlev, predicted_fecundity_lwr < 0)

# change species codes to Genus names
estimates_dat <- plotlev
estimates_dat$species <- case_match(estimates_dat$species, 
                                    'brohor' ~ "Bromus", 
                                      'vulmic' ~ "Festuca", 
                                      'miccal' ~ "Micropus", 
                                      'plaere' ~ "Plantago")
estimates_dat$species <- as.factor(estimates_dat$species)

# prepare facet tags
facet_labels <- data.frame(
  species = c("Bromus", "Festuca", "Micropus", "Plantago"),
  label = c("A", "B", "C", "D")  # Adjust labels as necessary
)

# shrink y axis with log 10
a <- ggplot(estimates_dat) +
  geom_linerange(aes(x = block, y = predicted_fecundity, ymin = predicted_fecundity_lwr, ymax = predicted_fecundity_upr, color = treatment)) +
geom_jitter(aes(x = block, y = predicted_fecundity, color = treatment)) +
  theme_bw() +
  facet_wrap(~species) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  scale_color_manual(values = c("honeydew4","mediumpurple"),
                     labels = c("without neighbors", 'with neighbors')) + 
  theme(legend.position = 'top',
        text = element_text(size = 16),
        axis.text.x=element_blank()) +
  labs(x = "Replicate (block)", y = "Predicted seed production", color = "") +
  geom_text(data = facet_labels, aes(x = -Inf, y = Inf, label = label), 
            size = 6, fontface = "bold", hjust = 0, vjust = 1.2,
            inherit.aes = FALSE) +
  scale_y_log10() 
a 

# grid 
filter(gridlev, grid_seed_lwr < 0)

estimates_dat <- gridlev
estimates_dat$species <- case_match(estimates_dat$species, 
                                    'brohor' ~ "Bromus", 
                                    'vulmic' ~ "Festuca", 
                                    'miccal' ~ "Micropus", 
                                    'plaere' ~ "Plantago")
estimates_dat$species <- as.factor(estimates_dat$species)

threshold_grid <- median(estimates_dat$grid_n)
threshold_grid

b <- ggplot(estimates_dat) +
  geom_linerange(aes(x = grid, y = grid_seed, ymin = grid_seed_lwr, ymax = grid_seed_upr, color = treatment)) +
  geom_point(aes(x = grid, y = grid_seed, color = treatment)) +
  theme_bw() +
  facet_wrap(~species) +
  geom_hline(yintercept = threshold_grid, linetype = "dashed", color = "darkred") +
  scale_color_manual(values = c("honeydew4","mediumpurple"),
                     labels = c("without neighbors", 'with neighbors')) + 
  theme(legend.position = 'top',
        text = element_text(size = 16),
        axis.text.x=element_blank()) +
  labs(x = "Replicate (grid)", y = "Predicted seed production", color = "") +  
  scale_y_log10() +
  geom_text(data = facet_labels, aes(x = -Inf, y = Inf, label = label), 
            size = 6, fontface = "bold", hjust = 0, vjust = 1.2,
            inherit.aes = FALSE) 
b

# site
filter(sitelev, site_seed_lwr < 0)

estimates_dat <- sitelev
estimates_dat$species <- case_match(estimates_dat$species, 
                                    'brohor' ~ "Bromus", 
                                    'vulmic' ~ "Festuca", 
                                    'miccal' ~ "Micropus", 
                                    'plaere' ~ "Plantago")
estimates_dat$species <- as.factor(estimates_dat$species)

threshold_site <- median(estimates_dat$site_n)
threshold_site

c <- ggplot(estimates_dat) +
  geom_linerange(aes(x = site, y = site_seed, ymin = site_seed_lwr, ymax = site_seed_upr, color = treatment)) +
  geom_point(aes(x = site, y = site_seed, color = treatment)) +
  theme_bw() +
  facet_wrap(~species) +
  geom_hline(yintercept =threshold_site, linetype = "dashed", color = "darkred") +
  scale_color_manual(values = c("honeydew4","mediumpurple"),
                     labels = c("without neighbors", 'with neighbors')) + 
  theme(legend.position = 'top',
        text = element_text(size = 16),
        axis.text.x=element_blank()) +
  labs(x = "Replicate (site)", y = "Predicted seed production", color = "") +  
  scale_y_log10() +
  geom_text(data = facet_labels, aes(x = -Inf, y = Inf, label = label), 
            size = 6, fontface = "bold", hjust = 0, vjust = 1.2,
            inherit.aes = FALSE)
c

#--------------------------------------
pdf('Figures/supp_pop_model_lwr_CI_plot.pdf', width = 12, height = 12)
a
dev.off()

pdf('Figures/supp_pop_model_lwr_CI_grid.pdf', width = 8, height = 8)
b
dev.off()

pdf('Figures/supp_pop_model_lwr_CI_site.pdf',width = 8, height = 8)
c
dev.off()
#--------------------------------------
