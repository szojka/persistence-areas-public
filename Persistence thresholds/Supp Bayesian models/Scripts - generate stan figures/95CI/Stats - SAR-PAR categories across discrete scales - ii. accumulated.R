
#-------------------------------------------------------------------------------
# DESCRIPTION: Figure 4A Bayes version
#-------------------------------------------------------------------------------

##############################
# Load packages and functions
##############################

# library(readr)
# library(tidyverse)
# library(truncnorm)
# library(stringr)
# library(BiodiversityR)
# library(car)
# library(lme4)
# library(glmmTMB)
# library(rstan)
# library(brms)
# library(visreg)
# library(Rmisc)
# library(patchwork)
# library(sf)

source(here::here("Supp Bayesian models/Scripts - generate stan figures/95CI/Source - MAIN fitnessdata POST STAN.R"))

#----------------
#Plot
#----------------

# For poa and pob 
dt_p <- dt1 %>% # dt1 is plotlev (length 3060) at 
  dplyr::select(-occurrence, -intraspecific_abundance ) 
dt_o <- dt1 %>%
  dplyr::select(-persistence, -predicted_fecundity, -persist_final) #%>% dplyr::filter(!occurrence%in%'no')
dt_op <- left_join(dt_o, dt_p, by = c('species','treatment','block','tube', 'grid', 'site'), 
                   relationship = "many-to-many")
dt_op <- dt_op %>% #dplyr::filter(occurrence == 'yes') %>% # removed this because we want the zeros that come from occ = no
  dplyr::select(-intraspecific_abundance, -tube, -predicted_fecundity, -persistence  )
dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
dt_op$occurrence <- as.numeric(dt_op$occurrence)

# if statements: if both occ and persistence = 0, then 0, if occ = 0, then 0, if persist = 0 , then 0, if pers AND occ have 1, then 1
dt_op$realized <- NA
dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persist_final  == 0, 0, 
                         ifelse(dt_op$occurrence == 1 & dt_op$persist_final  == 1, 1, NA))

# rename dataframe so that I can keep this dataframe for grids and scripts
op <- dt_op %>% 
  select(-occurrence, -persist_final ) %>%
  distinct() %>%
  group_by(block,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
  dplyr::mutate(replicates = n()) %>%
  tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species

# NA to zeros
op[is.na(op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
op$species_no <- rowSums(op[c(6:9)])

op$grid <- as.factor(op$grid)
op$site <- as.factor(op$site)
op$treatment <- as.factor(op$treatment)

# good data for POA or POB

# for OB PA and PB  
temp <- plotlev # easier coding than yes/no
temp$occurrence <- case_match( temp$occurrence, 
                               'yes' ~ "1",
                               'no' ~ "0",
                               .default =  temp$occurrence)
#temp$persistence <- as.numeric(temp$persistence)
temp$occurrence <- as.numeric(temp$occurrence)

temp_p <- temp %>% # temporary data frame to manipulate persistence (yes or no)
  dplyr::select(block,  species, treatment, persist_final , grid, site) %>%
  dplyr::distinct() %>%
  dplyr::group_by(block, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  distinct()
p <- pivot_wider(temp_p, names_from = "species", values_from = "persist_final")
p[is.na(p)] <- 0
p$type <- "P"
p$species_no <- rowSums(p[c(6:9)])

temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
  dplyr::select(block, species, treatment, occurrence, grid, site) %>%
  dplyr::distinct() %>%
  dplyr::group_by(block, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  distinct()
o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
o[is.na(o)] <- 0
o$type <- "O"
o$species_no <- rowSums(o[c(6:9)])

o$type <- as.factor(o$type)
p$type <- as.factor(p$type)

# combine all types together, using persistence and occurrence but getting more specific.
pob <- op %>% 
  dplyr::filter(treatment %in% c('B')) %>%
  dplyr::mutate(type = "POB") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment) %>%
  dplyr::distinct()

poa <- op %>% 
  dplyr::filter(treatment %in% c('A')) %>%
  dplyr::mutate(type = "POA") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()

ob <- o %>% 
  dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
  dplyr::mutate(type = "OB") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()

pb <- p %>% 
  dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
  dplyr::mutate(type = "PB") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()

pa <- p %>% 
  dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
  dplyr::mutate(type = "PA") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()


#type as predictor, with treatment implicit.
plot_par <- rbind(pob,ob,pb,poa,pa) # all these are types.

# make response variable proportions and remove species
plot_par <- plot_par %>%
  dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
  mutate(species_prop = species_no/replicates)

# make type a factor
plot_par$type <- as.factor(plot_par$type)
plot_par$scale <- c("plot")

#-------------------------------
# GRID & SITE
#-------------------------------

# Changing the way I do this. Summing persistence grouped by grid, and if > 0, persistence = 1

dt_op %>%
  dplyr::ungroup() %>%
  dplyr::select(-occurrence, -persist_final) %>%
  dplyr::group_by(treatment, grid, species) %>%
  reframe(sum_grid = sum(realized), treatment=treatment, grid=grid, species=species) %>%
  distinct() # there are fair amounts of zeros here so I trust this.

dt_op %>%
  dplyr::ungroup() %>%
  dplyr::select(-occurrence, -persist_final, -grid) %>%
  dplyr::group_by(treatment, site, species) %>%
  reframe(sum_site = sum(realized), treatment=treatment, site=site, species=species) %>%
  distinct() # also zeros

# if there is a realized within a grid and species, make grid_realized a 1
dt_op_scale <- dt_op %>% 
  dplyr::ungroup() %>%
  dplyr::select(-occurrence, -persist_final) %>%
  dplyr::group_by(treatment, grid, species) %>%
  dplyr::mutate(sum_grid_realized = sum(realized)) %>% # sum by grid
  dplyr::mutate(grid_realized = ifelse(sum_grid_realized > 0, 1, 0)) %>%
  # now for site
  dplyr::ungroup() %>%
  dplyr::group_by(treatment, site, species) %>%
  dplyr::mutate(sum_site_realized = sum(realized)) %>% # sum by grid
  dplyr::mutate(site_realized = ifelse(sum_site_realized > 0, 1, 0))
  

# break up by grid, site, calculate replicates and species rowsums, then proportions
dt_op_grid <- dt_op_scale %>%
  ungroup() %>%
  dplyr::select(treatment, grid, site, species, grid_realized) %>%
  dplyr::distinct() %>%
  dplyr::group_by(treatment, grid) %>%
  dplyr::mutate(replicates = n()) %>%
  pivot_wider(., names_from = "species", values_from = "grid_realized")
dt_op_grid[is.na(dt_op_grid)] <- 0
dt_op_grid$type <- "OP"
dt_op_grid$species_no <- rowSums(dt_op_grid[c(5:8)])

dt_op_grid$type <- as.factor(dt_op_grid$type)
dt_op_grid$site <- as.factor(dt_op_grid$site)
dt_op_grid$treatment <- as.factor(dt_op_grid$treatment)
dt_op_grid$grid <- as.factor(dt_op_grid$grid)
dt_op_grid$scale <- "grid"
dt_op_grid$block <- NA

dt_op_site <- dt_op_scale %>%
  dplyr::ungroup() %>%
  dplyr::select(treatment, site, species, site_realized) %>%
  dplyr::distinct() %>%
  dplyr::group_by(treatment, site) %>%
  dplyr::mutate(replicates = n()) %>%
  pivot_wider(., names_from = "species", values_from = "site_realized")
dt_op_site[is.na(dt_op_site)] <- 0
dt_op_site$type <- "OP"
dt_op_site$species_no <- rowSums(dt_op_site[c(4:7)])

dt_op_site$type <- as.factor(dt_op_site$type)
dt_op_site$site <- as.factor(dt_op_site$site)
dt_op_site$treatment <- as.factor(dt_op_site$treatment)
dt_op_site$scale <- "site"
dt_op_site$block <- NA
dt_op_site$grid <- NA

dt_op_site_grid <- rbind(dt_op_site,dt_op_grid)

# if there is a potential within a grid and species, make grid_potential a 1
temp_p %>%
  dplyr::ungroup() %>%
  dplyr::group_by(treatment, grid, species) %>%
  reframe(sum_grid = sum(persist_final), treatment=treatment, grid=grid, species=species) %>%
  distinct() %>%
  filter(sum_grid == 0) # 63 zeros
temp_p %>%
  dplyr::ungroup() %>%
  dplyr::group_by(treatment, site, species) %>%
  reframe(sum_site = sum(persist_final), treatment=treatment, site=site, species=species) %>%
  distinct() %>%
  filter(sum_site == 0) # 20 zeros

temp_p_scale <- temp_p %>%
  dplyr::group_by(treatment, grid, species) %>%
  dplyr::mutate(sum_grid_p = sum(persist_final)) %>%
  #dplyr::mutate(grid_persist = ifelse(1 %in% persist_final, 1, 0)) %>%
  dplyr::mutate(grid_persist = ifelse(sum_grid_p > 0, 1, 0)) %>%
  # now site:
  dplyr::ungroup() %>%
  dplyr::group_by(treatment, site, species) %>%
  dplyr::mutate(sum_site_p = sum(persist_final)) %>%
  dplyr::mutate(site_persist = ifelse(sum_site_p > 0, 1, 0))
  #dplyr::mutate(site_persist = ifelse(1 %in% persist_final, 1, 0)) # changed from grid_persist to persistence

# break up by grid, site, calculate replicates and species rowsums, then proportions
p_grid <- temp_p_scale %>%
  dplyr::ungroup() %>%
  dplyr::select(treatment, grid, site, species, grid_persist) %>%
  dplyr::distinct() %>%
  dplyr::group_by(treatment, grid) %>%
  dplyr::mutate(replicates = n()) %>%
  pivot_wider(., names_from = "species", values_from = "grid_persist")
p_grid[is.na(p_grid)] <- 0
p_grid$type <- "P"
p_grid$species_no <- rowSums(p_grid[c(5:8)])

p_grid$type <- as.factor(p_grid$type)
p_grid$site <- as.factor(p_grid$site)
p_grid$treatment <- as.factor(p_grid$treatment)
p_grid$grid <- as.factor(p_grid$grid)
p_grid$scale <- "grid"
p_grid$block <- NA

p_site <- temp_p_scale %>%
  ungroup() %>%
  dplyr::select(treatment, site, species, site_persist) %>%
  dplyr::distinct() %>%
  dplyr::group_by(treatment, site) %>%
  dplyr::mutate(replicates = n()) %>%
  pivot_wider(., names_from = "species", values_from = "site_persist")
p_site[is.na(p_site)] <- 0
p_site$type <- "P"
p_site$species_no <- rowSums(p_site[c(4:7)])

p_site$type <- as.factor(p_site$type)
p_site$site <- as.factor(p_site$site)
p_site$treatment <- as.factor(p_site$treatment)
p_site$scale <- "site"
p_site$grid <- NA
p_site$block <- NA

p_site_grid <- rbind(p_site,p_grid)


# use same method for calculating occurrence as averaged:

temp <- gridlev # easier coding than yes/no
temp$occurrence <- case_match( temp$occurrence, 
                               'yes' ~ "1",
                               'no' ~ "0",
                               .default =  temp$occurrence)
#temp$persistence <- as.numeric(temp$persistence)
temp$occurrence <- as.numeric(temp$occurrence)

o_grid <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
  dplyr::select(species, treatment, occurrence, grid, site) %>%
  dplyr::distinct() %>%
  dplyr::group_by(grid, treatment) %>%
  dplyr::mutate(replicates = n())
o_grid <- pivot_wider(o_grid, names_from = "species", values_from = "occurrence")
o_grid[is.na(o_grid)] <- 0 
o_grid$type <- "O"
o_grid$species_no <- rowSums(o_grid[c(5:8)]) # should be identical to averaged method diversity (SAR)

o_grid$type <- as.factor(o_grid$type)
o_grid$site <- as.factor(o_grid$site)
o_grid$treatment <- as.factor(o_grid$treatment)
o_grid$grid <- as.factor(o_grid$grid)
o_grid$scale <- "grid"
o_grid$block <- NA

# site level
temp <- sitelev %>% # easier coding than yes/no
  select(-contingency, -site_n, -site_seed)
temp$occurrence <- case_match( temp$occurrence, 
                               'yes' ~ "1",
                               #'no' ~ "0",
                               .default =  temp$occurrence)
temp$occurrence <- as.numeric(temp$occurrence)

o_site <- temp %>%
  dplyr::select(treatment, site, species, occurrence) %>%
  dplyr::distinct() %>%
  dplyr::group_by(treatment, site) %>%
  dplyr::mutate(replicates = n()) %>%
  pivot_wider(., names_from = "species", values_from = "occurrence")
o_site[is.na(o_site)] <- 0
o_site$type <- "O"
o_site$species_no <- rowSums(o_site[c(4:7)])

o_site$type <- as.factor(o_site$type)
o_site$site <- as.factor(o_site$site)
o_site$treatment <- as.factor(o_site$treatment)
o_site$scale <- "site"
o_site$grid <- NA
o_site$block <- NA

o_site_grid <- rbind(o_site,o_grid)

# ORGANIZE INTO PB, PA, OB, POA, POB

# using dt_op_site_grid, p_site_grid, o_site_grid
# then combine with plot_par

# combine all types together, using persistence and occurrence but getting more specific.
pob <- dt_op_site_grid %>%  
  dplyr::filter(treatment %in% c('B')) %>%
  dplyr::mutate(type = "POB") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment) %>%
  dplyr::distinct()

poa <- dt_op_site_grid %>% 
  dplyr::filter(treatment %in% c('A')) %>%
  dplyr::mutate(type = "POA") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()

ob <- o_site_grid %>% 
  dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
  dplyr::mutate(type = "OB") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()

pb <- p_site_grid %>% 
  dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
  dplyr::mutate(type = "PB") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()

pa <- p_site_grid %>% 
  dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
  dplyr::mutate(type = "PA") %>%
  dplyr::ungroup() %>%
  dplyr::select(-treatment)%>%
  dplyr::distinct()


# type as predictor, with treatment implicit.
grid_site_par <- rbind(pob,ob,pb,poa,pa) # all these are types.

# make response variable proportions and remove species
grid_site_par <- grid_site_par %>%
  dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
  mutate(species_prop = species_no/replicates)

plot_par$scale <- "plot"

par_all_accum <- rbind(plot_par,grid_site_par)

# order scales plot grid site
par_all_accum$scale <- factor(par_all_accum$scale, 
                                 levels = c('plot', 'grid', 'site'))

# rename ob pa poa pob
par_all_accum$type <- case_match(par_all_accum$type, 
                                    'OB' ~ "Diversity (SAR)",
                                    'PA' ~ "Potential \n without neighbors",
                                    'POA' ~ "Realized \n without neighbors",
                                    'PB' ~ "Potential \n with neighbors",
                                    'POB' ~ "Realized \n with neighbors",
                                    .default =  par_all_accum$type)

par_all_accum$type<- as.factor(par_all_accum$type)
par_all_accum$block<- as.factor(par_all_accum$block)

############################
# summary stats and checks
############################

par_all_accum <- par_all_accum %>%
  dplyr::group_by(type,scale) %>%
  dplyr::mutate(mean = mean(species_prop)) 
ggplot(par_all_accum) + 
  geom_histogram(aes(x=species_prop)) +
  facet_grid(rows = vars(scale), cols = vars(type)) + # realized are both zero inflated or more poisson like
  geom_vline(aes(xintercept = mean), color = "red")

par_all_accum %>%
  dplyr::group_by(type, scale) %>%
  dplyr::mutate(min = min(species_prop)) %>%
  dplyr::mutate(max = max(species_prop)) %>%
  dplyr::mutate(std = sd(species_prop)) %>%
  dplyr::mutate(N = n()) %>%
  dplyr::select(type, scale, min, max, std, mean, N) %>%
  distinct()

###############
# load models 
###############

# Loading .rdata models

############################
# ABIOTIC & BIOTIC TOGETHER

#------------
# Plot

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms 95CI/m.p_brm.rdata"))

# check model fit:
summary(m.p_brm)
# plot(m.p_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.p_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_p <- obj$type # data housed in function
output_p <- output_p %>%
  select(type, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'plot')

# add pairwise comparison table
em <- emmeans(m.p_brm, ~type, type = "response") # specify green_index_scaled values

p.contrast <- as.data.frame(pairs(em))

#------------
# Grid

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms 95CI/m.g_brm.rdata"))

# check model fit:
summary(m.g_brm)
# plot(m.g_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.g_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_g <- obj$type # data housed in function
output_g <- output_g %>%
  select(type, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'grid')

# add pairwise comparison table
em <- emmeans(m.g_brm, ~type, type = "response") 

g.contrast <- as.data.frame(pairs(em))

#------------
# Site

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms 95CI/m.s_brm.rdata"))

# check model fit:
summary(m.s_brm)
# posterior_summary(m.s_brm)
# plot(m.s_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.s_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_s <- obj$type # data housed in function
output_s <- output_s %>%
  select(type, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'site')

# add pairwise comparison table
em <- emmeans(m.s_brm, ~type, type = "response") 

s.contrast <- as.data.frame(pairs(em))

##############################################
# PULL PREDICTIONS TOGETHER WITH GGPREDICT
##############################################

vis1 <- rbind(output_p,output_g,output_s)
vis1$scale <- factor(vis1$scale, c("plot","grid","site")) # order group into increasing scales

names(vis1) # all correct classes.
vis1 <- as.data.frame(vis1)
# add scenario, kind, treatment, and scale_type

vis1$treatment <- NA
vis1$treatment[vis1$type %in% 
                          c("Potential \n with neighbors","Diversity (SAR)","Realized \n with neighbors")] <- "With neighbors"
vis1$treatment[vis1$type %in% 
                          c("Potential \n without neighbors","Realized \n without neighbors")] <- "Without neighbors"

vis1$kind <- NA
vis1$kind[vis1$type %in% 
                     c("Potential \n with neighbors","Potential \n without neighbors")] <- "Potential (PAR)"

vis1$kind[vis1$type %in% 
                     c("Realized \n with neighbors","Realized \n without neighbors")] <- "Realized (PAR)"

vis1$kind[vis1$type %in% 
                     c("Diversity (SAR)")] <- "Diversity (SAR)"

vis1$scenario <- NA
vis1$scenario <- "Data"

vis1$scale_type <- NA
vis1$scale_type <- "Accumulated"

#############################
# TABLE OUTPUT
#############################

#####################
# 1. Estimates

# vis1, organize by scale, type, predicted, confidence intervals "[{conf.low}, {conf.high}]]
# discard 'treatment' and 'scenario' and 'kind'

paste(round(vis1$prediction_lwr,2), round(vis1$prediction_upr,1), sep = ", ")

tab_acc_PAR_estimates <- vis1 |>
  dplyr::mutate(conf.int =paste(round(prediction_lwr,2), round(prediction_upr,1), sep = ", ")) |>
  dplyr::mutate(predicted = round(prediction,2)) |>
  # dplyr::select(-kind, -scenario, -treatment, -std.error, -scale_type, -conf.low, -conf.high) |>
  dplyr::select(type, scale, predicted, conf.int) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    type = 'Data type',
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
gt::gtsave(tab_acc_PAR_estimates, filename = gt_file) 
pagedown::chrome_print(input = gt_file,
                      output = here::here("Tables/Bayes/95CI/9_STAN_tab_acc_PAR_estimates.pdf"))

