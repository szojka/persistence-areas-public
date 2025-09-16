
#-------------------------------------------------------------------------------
# DESCRIPTION: create SAR PAR categories that accumulate across scales: diversity, 
# realized-PAR, and potential-PAR, with and without neighbors. 
# accumulated persistence calculated as once it persists or both, stays that way 
# through higher scales
# Corresponds to figure 4E.
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
# library(visreg)
# library(Rmisc)
# library(patchwork)
# library(sf)

source("Scripts/Source - MAIN fitnessdata.R") # must run to reset data I used in 'averaged' script.

#----------------
#Plot
#----------------

# For poa and pob 
# need to have separate dataframes for just Occurrence info and just Persistence info. 
dt_p <- dt1 %>% 
  dplyr::select(-occurrence, -intraspecific_abundance ) 
dt_o <- dt1 %>%
  dplyr::select(-persistence, -y_hat ) 
dt_op <- left_join(dt_o, dt_p, by = c('species','treatment','names','tag', 'grid', 'site'), 
                   relationship = "many-to-many")
dt_op <- dt_op %>% 
  dplyr::select(-intraspecific_abundance, -tag, -y_hat )
dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
dt_op$occurrence <- as.numeric(dt_op$occurrence)
dt_op$persistence[dt_op$persistence=='yes'] <- 1
dt_op$persistence[dt_op$persistence == 'no'] <- 0
dt_op$persistence<- as.numeric(dt_op$persistence)

# realized PAR = aligned present.
# if statements: if both occ and persistence = 0, then 0, if occ = 0, then 0, if persist = 0 , then 0, if pers AND occ have 1, then 1
dt_op$realized <- NA
dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persistence == 0, 0, 
                         ifelse(dt_op$occurrence == 1 & dt_op$persistence == 1, 1, NA))

# rename dataframe so that I can keep this dataframe for grids and scripts
op <- dt_op %>% 
  select(-occurrence, -persistence) %>%
  distinct() %>%
  group_by(names,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
  dplyr::mutate(replicates = n()) %>%
  tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species

# NA to zeros
op[is.na(op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
op$species_no <- rowSums(op[c(6:9)])

op$grid <- as.factor(op$grid)
op$site <- as.factor(op$site)
op$treatment <- as.factor(op$treatment)

# for OB PA and PB  
temp <- plotlev # easier coding than yes/no
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
  dplyr::select(names,  species, treatment, persistence, grid, site) %>%
  dplyr::distinct() %>%
  dplyr::group_by(names, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  distinct()
p <- pivot_wider(temp_p, names_from = "species", values_from = "persistence")
p[is.na(p)] <- 0
p$type <- "P"
p$species_no <- rowSums(p[c(6:9)])

temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
  dplyr::select(names, species, treatment, occurrence, grid, site) %>%
  dplyr::distinct() %>%
  dplyr::group_by(names, treatment) %>%
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


# type as predictor, with treatment implicit.
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

# Differences from averaged method start here.

# if there is a realized within a grid and species, make grid_realized a 1
dt_op <- dt_op %>% 
  dplyr::ungroup() %>%
  dplyr::select(-persistence, -occurrence) %>%
  dplyr::group_by(treatment, grid, species) %>%
  dplyr::mutate(grid_realized = ifelse(1 %in% realized, 1, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(treatment, site, species) %>%
  dplyr::mutate(site_realized = ifelse(1 %in% realized, 1, 0)) # changed from grid_realized to realized

# break up by grid, site, calculate replicates and species rowsums, then calculate proportions
dt_op_grid <- dt_op %>%
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
dt_op_grid$names <- NA

dt_op_site <- dt_op %>%
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
dt_op_site$names <- NA
dt_op_site$grid <- NA

dt_op_site_grid <- rbind(dt_op_site,dt_op_grid)


# if there is a potential within a grid and species, make grid_potential a 1
temp_p <- temp_p %>%
  dplyr::group_by(treatment, grid, species) %>%
  dplyr::mutate(grid_persist = ifelse(1 %in% persistence , 1, 0)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(treatment, site, species) %>%
  dplyr::mutate(site_persist = ifelse(1 %in% persistence, 1, 0)) # changed from grid_persist to persistence

# break up by grid, site, calculate replicates and species rowsums, then proportions
p_grid <- temp_p %>%
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
p_grid$names <- NA

p_site <- temp_p %>%
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
p_site$names <- NA

p_site_grid <- rbind(p_site,p_grid)


# use same method for calculating occurrence as averaged:

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

o_grid <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
  dplyr::select(species, treatment, occurrence, grid, site) %>%
  dplyr::distinct() %>%
  dplyr::group_by(grid, treatment) %>%
  dplyr::mutate(replicates = n())
o_grid <- pivot_wider(o_grid, names_from = "species", values_from = "occurrence")
o_grid[is.na(o_grid)] <- 0 
o_grid$type <- "O"
o_grid$species_no <- rowSums(o_grid[c(5:8)]) # should be identical to shmear method diversity (SAR)

o_grid$type <- as.factor(o_grid$type)
o_grid$site <- as.factor(o_grid$site)
o_grid$treatment <- as.factor(o_grid$treatment)
o_grid$grid <- as.factor(o_grid$grid)
o_grid$scale <- "grid"
o_grid$names <- NA

# site level
temp <- sitelev %>% # easier coding than yes/no
  select(-contingency, -site_n, -site_seed)
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
o_site$names <- NA

o_site_grid <- rbind(o_site,o_grid)

# ORGANIZE INTO PB, PA, OB, POA, POB

# using dt_op_site_grid, p_site_grid, o_site_grid
# then combine with plot_par

# combine all types together, using persistence and occurrence but getting more specific.
pob <- dt_op_site_grid %>%  # note that some sites are missing brohor and miccal entirely from realized.
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
par_all_accum$names<- as.factor(par_all_accum$names)

############################
# summary stats and checks
############################

par_all_accum <- par_all_accum %>%
  dplyr::group_by(type,scale) %>%
  dplyr::mutate(mean = mean(species_prop)) 
ggplot(par_all_accum) + 
  geom_histogram(aes(x=species_prop)) +
  facet_grid(rows = vars(scale), cols = vars(type)) + # realized are both zero inflated
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
# Fit models 
###############

############################
# ABIOTIC & BIOTIC TOGETHER

#------------
# Plot
# 
# remove  + (1|site:grid) for now as won't work with higher scales
m.p <- glmmTMB::glmmTMB(species_prop ~ type + (1|names) + (1|site:grid) + (1|site),
                          data = dplyr::filter(par_all_accum, scale %in% 'plot'),
                          family = binomial(),
                          weights = replicates)

# Check model fit
testZeroInflation(m.p)
testDispersion(m.p)
plot(fitted(m.p), residuals(m.p))
hist(residuals(m.p))
loc_sim_ouput <- simulateResiduals(m.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.p)

Anova(m.p, type = 3)

v <- ggpredict(m.p, terms = c("type"), 
          type = "fe"); plot(v)

# add pairwise comparison table
em <- emmeans(m.p, ~type, type = "response") 
em

p.contrast <- as.data.frame(pairs(em))
pairs(em)

#------------
# Grid

m.g <- glmmTMB::glmmTMB(species_prop ~ type + (1|site:grid) + (1|site),
                          data = dplyr::filter(par_all_accum, scale %in% 'grid'),
                          family = binomial(),
                          weights = replicates)


# Check model fit
testZeroInflation(m.g)
testDispersion(m.g)
plot(fitted(m.g), residuals(m.g))
hist(residuals(m.g))
loc_sim_ouput <- simulateResiduals(m.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # looks great,


summary(m.g)

Anova(m.g, type = 3)

v <- ggpredict(m.g, terms = c("type"), 
               type = "fe"); plot(v)

# add pairwise comparison table
em <- emmeans(m.g, ~type, type = "response") # specify green_index_scaled values
em

g.contrast <- as.data.frame(pairs(em))
pairs(em)

#------------
# Site

m.s <- glmmTMB::glmmTMB(species_prop ~ type + (1|site),
                          data = dplyr::filter(par_all_accum, scale %in% 'site'),
                          family = binomial(),
                          weights = replicates)

# Check model fit
testZeroInflation(m.s)
testDispersion(m.s)
plot(fitted(m.s), residuals(m.s))
hist(residuals(m.s))
loc_sim_ouput <- simulateResiduals(m.s)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.s)

Anova(m.s, type = 3)

v <- ggpredict(m.s, terms = c("type"), 
               type = "fe"); plot(v) 

# add pairwise comparison table
em <- emmeans(m.s, ~type, type = "response") 
em

s.contrast <- as.data.frame(pairs(em))
pairs(em)

##############################################
# PULL PREDICTIONS TOGETHER WITH GGPREDICT
##############################################

vis1p <- ggpredict(m.p, 
                    terms = c("type"), 
                    type = "fe")
vis1p$group <- 'plot'

vis1g <- ggpredict(m.g,
                    terms = c("type"), 
                    type = "fe")
vis1g$group <- 'grid'

vis1s <- ggpredict(m.s, 
                    terms = c("type"), 
                    type = "fe")
vis1s$group <- 'site'
# take confs to 1, as no error can really be estimated when all values are 1.
vis1s$conf.low[vis1s$x %in% 'Diversity (SAR)'] <- 1 
vis1s$conf.low[vis1s$x %in% 'Potential \n with neighbors'] <- 1
vis1s$conf.low[vis1s$x %in% 'Potential \n without neighbors'] <- 1
vis1s$conf.low[vis1s$x %in% 'Realized \n without neighbors'] <- 1

vis1 <- rbind(vis1p,vis1g,vis1s)
vis1$group <- factor(vis1$group, c("plot","grid","site")) # order group into increasing scales

names(vis1) # all correct classes.
vis1 <- as.data.frame(vis1)
# add scenario, kind, treatment, and scale_type

vis1$treatment <- NA
vis1$treatment[vis1$x %in% 
                          c("Potential \n with neighbors","Diversity (SAR)","Realized \n with neighbors")] <- "With neighbors"
vis1$treatment[vis1$x %in% 
                          c("Potential \n without neighbors","Realized \n without neighbors")] <- "Without neighbors"

vis1$kind <- NA
vis1$kind[vis1$x %in% 
                     c("Potential \n with neighbors","Potential \n without neighbors")] <- "Potential (PAR)"

vis1$kind[vis1$x %in% 
                     c("Realized \n with neighbors","Realized \n without neighbors")] <- "Realized (PAR)"

vis1$kind[vis1$x %in% 
                     c("Diversity (SAR)")] <- "Diversity (SAR)"

vis1$scenario <- NA
vis1$scenario <- "Data"

vis1$scale_type <- NA
vis1$scale_type <- "Accumulated"

names(vis1)[1] <- "type"
names(vis1)[6] <- "scale"

#############################
# TABLE OUTPUT
#############################

#####################
# 1. Estimates

# vis1, organize by scale, type, predicted, confidence intervals "[{conf.low}, {conf.high}]]
# discard 'treatment' and 'scenario' and 'kind'

paste(round(vis1$conf.low,2), round(vis1$conf.high,1), sep = ", ")

tab_acc_PAR_estimates <- vis1 |>
  dplyr::mutate(conf.int =paste(round(conf.low,2), round(conf.high,1), sep = ", ")) |>
  dplyr::mutate(predicted = round(predicted,2)) |>
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

tab_acc_PAR_estimates |>
  gtsave(paste0(here::here(),"/Tables/9tab_acc_PAR_estimates.pdf")) 

######################
# 2. ANOVA 

# three models one for each scale
#m.p, m.g, m.s

a1 <- Anova(m.p, type = 3) 
a2 <- Anova(m.g, type = 3)
a3 <- Anova(m.s, type = 3)

anova_acc <- data.frame(Intercept = c(a1$Chisq[1], a2$Chisq[1], a3$Chisq[1], 
                                      a1$Df[1], a2$Df[1],a3$Df[1],
                                      a1$`Pr(>Chisq)`[1],a2$`Pr(>Chisq)`[1],a3$`Pr(>Chisq)`[1]),
                        type = c(a1$Chisq[2], a2$Chisq[2], a3$Chisq[2], 
                                 a1$Df[2], a2$Df[2],a3$Df[2],
                                 a1$`Pr(>Chisq)`[2],a2$`Pr(>Chisq)`[2],a3$`Pr(>Chisq)`[2]),
                        scale = c(rep(c('plot','grid','site'),times = 3)),
                        stat = c(rep('Chi squared', times = 3), rep('df', times = 3), rep('P-value', times = 3))
                        )

anova_acc <- data.frame(Chi.squared = round(c(a1$Chisq[1],a1$Chisq[2],a2$Chisq[1],a2$Chisq[2],a3$Chisq[1],a3$Chisq[2]),3),
           scale = c('plot','plot','grid','grid','site','site'),
           Df =c(a1$Df[1],a1$Df[2],a2$Df[1],a2$Df[2],a3$Df[1],a3$Df[2]),
           P_value =  round(c(a1$`Pr(>Chisq)`[1],a1$`Pr(>Chisq)`[2],a2$`Pr(>Chisq)`[1],a2$`Pr(>Chisq)`[2],a3$`Pr(>Chisq)`[1],a3$`Pr(>Chisq)`[2]),3),
           predictor = c(rep(c("Intercept","Data type"), times = 3)))

tab_acc_anova <- anova_acc |>
  dplyr::select(scale, predictor, Chi.squared, Df, P_value) |>
  gt() |>
  tab_header( title = "",
              subtitle = "")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    predictor = 'Predictor',
    scale = 'Scale level',
    P_value = "P-value") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_acc_anova |>
  gtsave(paste0(here::here(),"/Tables/10tab_acc_anova.pdf")) 

############################
# 3. Contrasts

contrastz <- rbind(p.contrast, g.contrast, s.contrast)
contrastz$scale <- NA
contrastz$scale <- c(rep('plot',times = 10),rep('grid',times = 10), rep('site',times = 10))

levels(contrastz$contrast)
# remove \n from the level contrasts
contrastz$contrast <- case_match(contrastz$contrast,  
            "Diversity (SAR) / Potential \n with neighbors"  ~ "SAR / Potential-PAR with neighbors",            
            "Diversity (SAR) / Potential \n without neighbors"  ~ "SAR / Potential-PAR without neighbors",       
            "Diversity (SAR) / Realized \n with neighbors"  ~  "SAR / Realized-PAR with neighbors",             
            "Diversity (SAR) / Realized \n without neighbors" ~ "SAR / Realized-PAR without neighbors",
            "Potential \n with neighbors / Potential \n without neighbors"  ~ "Potential-PAR with neighbors / Potential-PAR without neighbors",
            "Potential \n with neighbors / Realized \n with neighbors"  ~ "Potential-PAR with neighbors / Realized-PAR with neighbors",     
            "Potential \n with neighbors / Realized \n without neighbors" ~ "Potential-PAR with neighbors / Realized-PAR without neighbors" ,  
            "Potential \n without neighbors / Realized \n with neighbors"  ~  "Potential-PAR without neighbors / Realized-PAR with neighbors", 
            "Potential \n without neighbors / Realized \n without neighbors" ~ "Potential-PAR without neighbors / Realized-PAR without neighbors",
            "Realized \n with neighbors / Realized \n without neighbors"  ~"Realized-PAR with neighbors / Realized-PAR without neighbors" ) 

tab_accPAR_contrasts <- contrastz |>
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

tab_accPAR_contrasts |>
  gtsave(paste0(here::here(),"/Tables/11tab_acc_PAR_contrasts.pdf")) 

