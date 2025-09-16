
#-------------------------------------------------------------------------------
# DESCRIPTION: Figure 4F stats, Bayes version (predictions from full Beverton-Holt)
#-------------------------------------------------------------------------------

# using plotlev, gridlev, sitelev

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
# library(visreg)
# library(Rmisc)
# library(patchwork)

source(here::here("Supp Bayesian models/Scripts - generate stan figures/FullBevHolt/95CI/Source - MAIN fitnessdata POST STAN.R")) 
#--------------------------------------------
# Plot level
#--------------------------------------------

########################
# Data wrangling
########################

# For poa and pob 
dt_p <- dt1 %>% 
  dplyr::select(-occurrence, -intraspecific_abundance) 
dt_o <- dt1 %>%
  dplyr::select(-persistence, -persist_final, -predicted_fecundity ) 
dt_op <- left_join(dt_o, dt_p, by = c('species','treatment','block','tube', 'grid', 'site'), relationship = 'many-to-many')
dt_op <- dt_op %>% 
  dplyr::select(-intraspecific_abundance, -tube, -predicted_fecundity )
dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
dt_op$occurrence <- as.numeric(dt_op$occurrence)
dt_op$persist_final<- as.numeric(dt_op$persist_final)

dt_op$realized <- NA
dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persist_final == 0, 0, 
                        ifelse(dt_op$occurrence == 1 & dt_op$persist_final == 1, 1, NA))

dt_op <- dt_op %>% 
  select(-occurrence, -persist_final, -persistence) %>%
  distinct() %>%
  group_by(block,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
  dplyr::mutate(replicates = n()) %>%
  tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species

# NA to zeros
dt_op[is.na(dt_op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
dt_op$species_no <- rowSums(dt_op[c(6:9)])

dt_op$grid <- as.factor(dt_op$grid)
dt_op$site <- as.factor(dt_op$site)
dt_op$treatment <- as.factor(dt_op$treatment)

# good data for POA or POB
  
# for OB PA and PB  
temp <- plotlev # easier coding than yes/no
  temp$occurrence <- case_match( temp$occurrence, 
                                  'yes' ~ "1",
                                  'no' ~ "0",
                                  .default =  temp$occurrence)
 # temp$persistence <- as.numeric(temp$persistence)
  temp$occurrence <- as.numeric(temp$occurrence)
  
  temp_p <- temp %>% # temporary data frame to manipulate persistence (yes or no)
    dplyr::select(block,  species, treatment, persist_final, grid, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(block, treatment) %>%
    dplyr::mutate(replicates = n()) %>%
    distinct()
  temp_p <- pivot_wider(temp_p, names_from = "species", values_from = "persist_final")
  temp_p[is.na(temp_p)] <- 0
  temp_p$type <- "P"
  temp_p$species_no <- rowSums(temp_p[c(6:9)])
  
  temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
    dplyr::select(block, species, treatment, occurrence, grid, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(block, treatment) %>%
    dplyr::mutate(replicates = n()) %>%
    distinct()
  temp_o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
  temp_o[is.na(temp_o)] <- 0
  temp_o$type <- "O"
  temp_o$species_no <- rowSums(temp_o[c(6:9)])
 
  temp_o$type <- as.factor(temp_o$type)
  temp_p$type <- as.factor(temp_p$type)

  # combine all types together, using persistence and occurrence but getting more specific.
  pob <- dt_op %>% 
    dplyr::filter(treatment %in% c('B')) %>%
    dplyr::mutate(type = "POB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment) %>%
    dplyr::distinct()
  
  poa <- dt_op %>% 
    dplyr::filter(treatment %in% c('A')) %>%
    dplyr::mutate(type = "POA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  
  ob <- temp_o %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
    dplyr::mutate(type = "OB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()

  pb <- temp_p %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
    dplyr::mutate(type = "PB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()

  pa <- temp_p %>% 
    dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
    dplyr::mutate(type = "PA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()

  
  # type as predictor, with treatment implicit.
  plot_par <- rbind(pob,ob,pb,poa,pa) # all these are types.
  
  # make response variable proportions and remove species
  plot_par <- plot_par %>%
    dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
    mutate(species_prop = species_no/replicates)
  
  # make type a factor
  plot_par$type <- as.factor(plot_par$type)
  plot_par$scale <- c("plot")

#--------------------------------------------
# Grid level
#--------------------------------------------

########################
# Data wrangling
########################

  # Create 'realized' data for poa and pob 
  
  dt_op <-  gridlev %>% dplyr::select(-grid_n, -grid_seed, -contingency)
  
  dt_op$occurrence <- as.character(dt_op$occurrence)
  dt_op$persistence <- as.character(dt_op$persistence)
  
  dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
  dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
  dt_op$occurrence <- as.numeric(dt_op$occurrence)
  dt_op$persistence[dt_op$persistence=='yes'] <- 1
  dt_op$persistence[dt_op$persistence == 'no'] <- 0
  dt_op$persistence<- as.numeric(dt_op$persistence)
  
  dt_op$realized <- NA
  dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persistence == 0, 0, 
                           ifelse(dt_op$occurrence == 1 & dt_op$persistence == 1, 1, NA))
  
  dt_op <- dt_op %>% 
    select(-occurrence, -persistence, -grid_seed_lwr,-grid_seed_upr) %>%
    distinct() %>%
    group_by(grid,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
    dplyr::mutate(replicates = n()) %>%
    tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species
  
  # NA to zeros
  dt_op[is.na(dt_op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
  dt_op$species_no <- rowSums(dt_op[c(5:8)])
  
  dt_op$grid <- as.factor(dt_op$grid)
  dt_op$site <- as.factor(dt_op$site)
  dt_op$treatment <- as.factor(dt_op$treatment)

 # DATA "POTENTIAL" FOR PA PB
  
  temp <- gridlev # easier coding than yes/no
  temp$persistence <- case_match( temp$persistence, 
                                  'yes' ~ "1",
                                  'no' ~ "0",
                                  .default =  temp$persistence)
  temp$occurrence <- case_match( temp$occurrence, 
                                 'yes' ~ "1",
                                 'no' ~ "0",
                                 .default =  temp$occurrence)
  temp$persistence <- as.numeric(temp$persistence)
  temp$occurrence <- as.numeric(temp$occurrence)
  
  temp_p <- temp %>% # temporary data frame to manipulate persistence (yes or no)
    dplyr::select(species, treatment, persistence, grid, site) %>%
    dplyr::distinct()  %>%
    dplyr::group_by(grid, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_p <- pivot_wider(temp_p, names_from = "species", values_from = "persistence")
  temp_p[is.na(temp_p)] <- 0 # this is fine, replicates have already been calculated, turn to 0 just to allow rowSum function to work.
  temp_p$type <- "P"
  temp_p$species_no <- rowSums(temp_p[c(5:8)]) # only places where there were NAs in resevoir grid 26.
  
  # Check that things are correct is that potential should always be equal or higher than realized
  dt_op %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no)) 

  temp_p %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no))

  # DATA "DIVERSITY" FOR OB
  
  temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
    dplyr::select(species, treatment, occurrence, grid, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(grid, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
  temp_o[is.na(temp_o)] <- 0 
  temp_o$type <- "O"
  temp_o$species_no <- rowSums(temp_o[c(5:8)])

  # combine all types together, using persistence and occurrence but getting more specific.
  pob <- dt_op %>% 
    dplyr::filter(treatment %in% c('B')) %>%
    dplyr::mutate(type = "POB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  poa <- dt_op %>% 
    dplyr::filter(treatment %in% c('A')) %>%
    dplyr::mutate(type = "POA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  ob <- temp_o %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
    dplyr::mutate(type = "OB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pb <- temp_p %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
    dplyr::mutate(type = "PB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pa <- temp_p %>% 
    dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
    dplyr::mutate(type = "PA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  
  # type as predictor, with treatment implicit.
  grid_par <- rbind(pob,ob,pb,poa,pa) # all these are types.
  
  # make response variable proportions and remove species
  grid_par <- grid_par %>%
    dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
    mutate(species_prop = species_no/replicates) 
  
  # make type a factor
  grid_par$type <- as.factor(grid_par$type)
  grid_par$scale <- c("grid")
  
#--------------------------------------------
# Site level
#--------------------------------------------

########################
# Data wrangling
########################

  # For poa and pob 
  dt_op <- sitelev %>% dplyr::select(-site_n, -site_seed, -contingency) 

  dt_op$occurrence <- as.character(dt_op$occurrence)
  dt_op$persistence <- as.character(dt_op$persistence)
  
  dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
  dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
  dt_op$occurrence <- as.numeric(dt_op$occurrence)
  dt_op$persistence[dt_op$persistence=='yes'] <- 1
  dt_op$persistence[dt_op$persistence == 'no'] <- 0
  dt_op$persistence<- as.numeric(dt_op$persistence)
  
  dt_op$realized <- NA
  dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persistence == 0, 0, 
                           ifelse(dt_op$occurrence == 1 & dt_op$persistence == 1, 1, NA))
  
  dt_op <- dt_op %>% 
    dplyr::select(-occurrence, -persistence,-site_seed_lwr,-site_seed_upr) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
    dplyr::mutate(replicates = n()) %>%
    tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species
   
  # NA to zeros
  dt_op[is.na(dt_op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
  dt_op$species_no <- rowSums(dt_op[c(4:7)])
  
  dt_op$site <- as.factor(dt_op$site)
  dt_op$treatment <- as.factor(dt_op$treatment)
  # good data for POA or POB
  
  # for PA PB OB
  temp <- sitelev %>% # easier coding than yes/no
    dplyr::select(-contingency, -site_n, -site_seed)
  temp$persistence <- case_match( temp$persistence, 
                                  'yes' ~ "1",
                                  'no' ~ "0",
                                  .default =  temp$persistence)
  temp$occurrence <- case_match( temp$occurrence, 
                                 'yes' ~ "1",
                                 #'no' ~ "0",
                                 .default =  temp$occurrence)
  temp$persistence <- as.numeric(temp$persistence)
  temp$occurrence <- as.numeric(temp$occurrence)
  
  temp_p <- temp %>% # temporary data frame to manipulate persistence (yes or no)
    dplyr::select(species, treatment, persistence, site) %>%
    dplyr::distinct()  %>%
    dplyr::group_by(site, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_p <- pivot_wider(temp_p, names_from = "species", values_from = "persistence")
  temp_p[is.na(temp_p)] <- 0
  temp_p$type <- "P"
  temp_p$species_no <- rowSums(temp_p[c(4:7)])
  
  temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
    dplyr::select(species, treatment, occurrence, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
  temp_o[is.na(temp_o)] <- 0
  temp_o$type <- "O"
  temp_o$species_no <- rowSums(temp_o[c(4:7)])
  
  # combine all types together, using persistence and occurrence but getting more specific.
  pob <- dt_op %>% 
    dplyr::filter(treatment %in% c('B')) %>%
    dplyr::mutate(type = "POB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  poa <- dt_op %>% 
    dplyr::filter(treatment %in% c('A')) %>%
    dplyr::mutate(type = "POA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  ob <- temp_o %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
    dplyr::mutate(type = "OB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pb <- temp_p %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
    dplyr::mutate(type = "PB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pa <- temp_p %>% 
    dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
    dplyr::mutate(type = "PA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  
  # type as predictor, with treatment implicit.
  site_par <- rbind(pob,ob,pb,poa,pa) # all these are types.
  
  # make response variable proportions and remove species
  site_par <- site_par %>%
  dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
  dplyr::mutate(species_prop = species_no/replicates) 
  
  # make type a factor
  site_par$type <- as.factor(site_par$type)
  site_par$scale <- c("site")

  # Check that things are correct is that potential should always be equal or higher than realized
  dt_op %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no)) 

  temp_p %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no))

  
################################################
# BRING ALL RAW DATAFRAMES TOGETHER FOR PLOTTING

grid_par$block <- NA
site_par$block <- NA
site_par$grid <- NA

names(plot_par)
names(grid_par)
names(site_par)
par_all_discrete <- rbind(plot_par, grid_par, site_par)

# Change variable types / names for all dataframes
par_all_discrete$scale <- as.factor(par_all_discrete$scale)

# order scales plot grid site
par_all_discrete$scale <- factor(par_all_discrete$scale, 
                                 levels = c('plot', 'grid', 'site'))

# rename ob pa poa pob
par_all_discrete$type <- case_match(par_all_discrete$type, 
                                'OB' ~ "Diversity (SAR)",
                                'PA' ~ "Potential \n without neighbors",
                                'POA' ~ "Realized \n without neighbors",
                                'PB' ~ "Potential \n with neighbors",
                                'POB' ~ "Realized \n with neighbors",
                                .default =  par_all_discrete$type)

par_all_discrete$type<- as.factor(par_all_discrete$type)
par_all_discrete$block<- as.factor(par_all_discrete$block)

############################
# summary stats and checks
############################
par_all_discrete <- par_all_discrete %>%
  dplyr::group_by(type,scale) %>%
  dplyr::mutate(mean = mean(species_prop)) 
ggplot(par_all_discrete) + 
 geom_histogram(aes(x=species_prop)) +
  facet_grid(rows = vars(scale), cols = vars(type)) + # realized are both zero inflated 
  geom_vline(aes(xintercept = mean), color = "red")

par_all_discrete %>%
  dplyr::group_by(type, scale) %>%
  dplyr::mutate(min = min(species_prop)) %>%
  dplyr::mutate(max = max(species_prop)) %>%
  dplyr::mutate(std = sd(species_prop)) %>%
  dplyr::mutate(N = n()) %>%
  dplyr::select(type, scale, min, max, std, mean, N) %>%
distinct()

#----------------------------------------------------------------------------
# MODELLING TYPE AS FACTOR RESPONSE.

###############
# Fit models 
###############

############################
# ABIOTIC & BIOTIC TOGETHER

#------------
# Plot

dat <- par_all_discrete %>%
  dplyr::filter(scale %in% 'plot') %>%
  mutate(species_prop = ifelse(species_prop == 0, 0.01, species_prop)) %>%
  mutate(species_prop = ifelse(species_prop == 1, 0.99, species_prop))
min(dat$species_prop); max(dat$species_prop)

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.p_avg_brm.rdata"))

# check model fit:
summary(m.p_avg_brm)
# plot(m.p_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.p_avg_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_p <- obj$type # data housed in function
output_p <- output_p %>%
  select(type, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'plot')

# add pairwise comparison table
em <- emmeans(m.p_avg_brm, ~type, type = "response") 

p.contrast <- as.data.frame(pairs(em))

#------------
# Grid

dat <- par_all_discrete %>%
  dplyr::filter(scale %in% 'grid') %>%
  mutate(species_prop = ifelse(species_prop == 0, 0.01, species_prop)) %>%
  mutate(species_prop = ifelse(species_prop == 1, 0.99, species_prop))
min(dat$species_prop); max(dat$species_prop)

# load model
load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.g_avg_brm.rdata"))

# check model fit:
summary(m.g_avg_brm)
# plot(m.g_avg_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.g_avg_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_g <- obj$type # data housed in function
output_g <- output_g %>%
  select(type, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'grid')

# add pairwise comparison table
em <- emmeans(m.g_avg_brm, ~type, type = "response") 

g.contrast <- as.data.frame(pairs(em))

#------------
# Site

dat <- par_all_discrete %>%
  dplyr::filter(scale %in% 'site') %>%
  mutate(species_prop = ifelse(species_prop == 0, 0.01, species_prop)) %>%
  mutate(species_prop = ifelse(species_prop == 1, 0.99, species_prop))
min(dat$species_prop); max(dat$species_prop)

# load 
load(here::here("Supp Bayesian models/Stan rdata/brms FullBevHolt/95CI/m.s_avg_brm.rdata"))

# check model fit:
summary(m.s_avg_brm)
# plot(m.s_brm) # shows posteriors and chains

# Extract estimates:
obj <- conditional_effects(m.s_avg_brm, prob = 0.95, robust = FALSE) # last arg when FALSE = mean when TRUE = median
output_s <- obj$type # data housed in function
output_s <- output_s %>%
  select(type, estimate__, lower__, upper__) %>%
  dplyr::rename('prediction' = 'estimate__') %>%
  dplyr::rename('prediction_lwr' = 'lower__') %>%
  dplyr::rename('prediction_upr' = 'upper__') %>%
  dplyr::mutate(scale = 'site')

# add pairwise comparison table
em <- emmeans(m.s_avg_brm, ~type, type = "response") 

s.contrast <- as.data.frame(pairs(em))

##############################################
# PULL PREDICTIONS TOGETHER WITH GGPREDICT
##############################################

vis2 <- rbind(output_p,output_g,output_s)
vis2$scale <- factor(vis2$scale, c("plot","grid","site")) # order group into increasing scales

names(vis2) # all correct classes.

vis2 <- as.data.frame(vis2)
# add scenario, kind, treatment, and scale_type

vis2$treatment <- NA
vis2$treatment[vis2$type %in% 
                 c("Potential \n with neighbors","Diversity (SAR)","Realized \n with neighbors")] <- "With neighbors"
vis2$treatment[vis2$type %in% 
                 c("Potential \n without neighbors","Realized \n without neighbors")] <- "Without neighbors"

vis2$kind <- NA
vis2$kind[vis2$type %in% 
            c("Potential \n with neighbors","Potential \n without neighbors")] <- "Potential (PAR)"

vis2$kind[vis2$type %in% 
            c("Realized \n with neighbors","Realized \n without neighbors")] <- "Realized (PAR)"

vis2$kind[vis2$type %in% 
            c("Diversity (SAR)")] <- "Diversity (SAR)"

vis2$scenario <- NA
vis2$scenario <- "Data"

vis2$scale_type <- NA
vis2$scale_type <- "Averaged"

#############################
# TABLE OUTPUT
#############################

# vis2, organize by scale, type, predicted, confidence intervals "[{conf.low}, {conf.high}]]
# discard 'treatment' and 'scenario' and 'kind'

paste(round(vis2$prediction_lwr,2), round(vis2$prediction_upr,1), sep = ", ")

tab_avg_PAR_estimates <- vis2 |>
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
 gt::gtsave(tab_avg_PAR_estimates, filename = gt_file) 
 pagedown::chrome_print(input = gt_file,
                        output = here::here("Tables/Bayes/FullBevHolt/95CI/12_STAN_tab_avg_PAR_estimates.pdf"))
 
