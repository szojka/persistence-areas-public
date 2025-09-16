
#------------------------------------------------------------
## DESCRIPTION: Source code for clean data files for scales
#------------------------------------------------------------

########################
# Load Packages 
########################

library(readr)
library(tidyverse)
library(truncnorm)
library(stringr)
library(car)
library(lme4)
library(glmmTMB)
library(visreg)
library(Rmisc)
library(patchwork)
library(sf)
library(ggeffects)
library(DHARMa)
library(gt)
library(emmeans)
library(mclogit)
library(here)

###########################

# https://gt.albert-rapp.de/getting_started.html
# necessary for saving but shouldn't have to use more than once.
#Sys.setenv(CHROMOTE_CHROME = "C:/Users/Megan Szojka/AppData/Local/Google/Chrome/Application/chrome.exe")
#chromote::find_chrome()

# set chrome for gt tables
Sys.setenv(CHROMOTE_CHROME = "C:/Users/Megan Szojka/AppData/Local/Google/Chrome/Application/chrome.exe") # replace with local chrome path

########################
#### Clean dataset ####

seeds_unclean <- read_csv("Data/Transplant_master data_Oct2020.csv")

seeds_unclean$status <- as.factor(as.character(seeds_unclean$status))
levels(seeds_unclean$status) # four levels

seeds_unclean$status <- as.character(as.factor(seeds_unclean$status)) # make status = NA into factor
seeds_unclean <- filter(seeds_unclean, !status%in%c('missing')) # keep herbivory - meaningful zero like 'zero' category for rocks
seeds_unclean$status[is.na(seeds_unclean$status)] <- c("normal")
seeds_unclean$seed[is.na(seeds_unclean$seed)] <- 0
seeds_clean <- ungroup(seeds_unclean)

#--------------------------------------.
#### Abundance and occurrence data ####

# abundance categories = 0,1,2,3
# key = 0, 1-10, 11-100, 101-1000

# preliminary miccal and plaere only
photos1 <- read_csv("Data/Abundance & Occurance 2019_final.csv")
photos1 <- select(photos1, -notes)
# gather plaere and miccal to species
photos <- gather(photos1, species, ab_cat, 6:9) # abundance categories
# join tags to photos using:
tag_labels <- read_csv("Data/Tag labels.csv")
key <- left_join(tag_labels, photos, by = c("species", "treatment", "grid", "site", "plot"), 
                 relationship = "many-to-many")

#### UTM & dat ####

utm <- read_csv("Data/UTM_szojkalocs_id.csv")
utm2 <- dplyr::select(utm, -name) 
utm2 <- as.matrix(utm2)

# Starting point Sampling

max(utm$utm_x) #557021.5
min(utm$utm_x) #547533.3

max(utm$utm_y) #4302411
min(utm$utm_y) #4297447

# find distances between my plot points in "utm" and each other (so that I always start a run on a data point)
nrun <- 450
utm2 <- sp::SpatialPoints(utm2)

# convert data to an sf object
utm2_sf <- st_as_sf(utm2, coords = c("UTM_x", "UTM_y"), crs = 32633) # Replace with the correct column names and CRS

# calculate distances between each point
x <- st_distance(utm2_sf)

# create a distance dataset
names <- c(utm2_sf$name) # Assuming 'name' is a column in your sf object

dist.dat <- as.data.frame(x)

# Should retain order such that unique key (names) is maintained for each group

# make a key to connect 'names' to 'seeds_clean' (Step 2 & 3)
nkey <- dplyr::select(tag_labels, -plot)
nkey$names1<-paste(nkey$site, nkey$grid, sep="")
num <- rep(c(1:25), each = 8)
nkey$names <- paste(num, nkey$names1, sep="_")
nkey <- dplyr::select(nkey, -names1)

# change seed column to persistence
pers <- dplyr::select(seeds_clean, -status)
#pers <- pivot_wider(pers, names_from = "species", values_from = "seed")
name_key <- left_join(pers, nkey, by = c("tag","species","treatment","grid","site"), 
                      relationship = "many-to-many")
name_key <- na.omit(name_key)

# temporary filter to make persistence column
# for each species, in each 'unique' value, sumarize into new column
name_key$unique <- paste(name_key$names, name_key$treatment, sep = "_")

p_sum <- name_key %>%
  dplyr::select(species, seed, unique)

p_sum$seed[p_sum$seed > 0] <- 1

plot_list <- p_sum %>% 
  group_by(unique) %>% 
  add_count() %>% 
  filter(n == 4)

uni <- p_sum %>% 
  filter(unique %in% plot_list$unique) %>%
  distinct() %>%
  pivot_wider(names_from = "species", values_from = "seed") %>%
  mutate(persist = rowSums(.[2:5]))

uni_key1 <- dplyr::select(uni, unique, persist)
uni_key <- left_join(uni_key1, name_key, by = "unique", 
                     relationship = "many-to-many")
uni_key <- uni_key %>%
  dplyr::select(-unique, -seed, -tag, -species)
rm(uni, uni_key1, p_sum, plot_list, pers)

# collapse dataframe (dist.dat) into two columns where column names are third column indicating value 
dat1 <- pivot_longer(dist.dat, cols = c(1:nrun), names_to = "run")
names <- unique(nkey$names)
dat1$names <- rep(names, each=nrun)
names(dat1)[2] <- c("area")
dat <- left_join(dat1, uni_key, by = "names", relationship = "many-to-many")
dat <- distinct(dat)

# checking all is well with this dataset
str(dat)
dat$run <- as.factor(as.character(dat$run))
levels(dat$run) # good: 450 exist
dat$treatment <- as.factor(as.character(dat$treatment))
levels(dat$treatment) # good
dat <- na.omit(dat)

# area and percent of persisting species are on way difference scales
# change area to km (/1000)
dat$area <- dat$area/1000

# scale category for joy divisions
dat$scale <- cut(dat$area, breaks = 15, labels = c(0:14), include.lowest = T)
dat$scale <- as.numeric(as.factor(dat$scale)) # do this so my loop below will run
table(dat$scale) # frequency of obs w/in each level is well distributed

dat$run <- as.character(as.factor(dat$run))

#-------------------------------------------------------.
#### Dataset for persistence & occurrence (plotlev) ####

# abundance categories = 0,1,2,3
# ab_cat = 0, 1-10, 11-100, 101-1000

# join seeds_clean to photo abundances
plotlev1 <- left_join(key, seeds_clean, by = c("tag", "species", "treatment", "grid", "site"), 
                      relationship = "many-to-many")
plotlev1 <- dplyr::filter(plotlev1, !status%in%NA) # NAs represent observations I do not have

# make sure occurrence data for abiotic plot reflects biotic since I cleared abiotic
# trying to say: if occurrence in A is < occurrence for B per plot, then A occurrence = B occurrence
plotlev1$ab_cat <- as.numeric(as.character(plotlev1$ab_cat))

# assign occurence values in B treatment to A (cleared)

plots <- levels(factor(plotlev1$plot))
grids <- as.numeric(levels(factor(plotlev1$grid)))
spp <- levels(factor(plotlev1$species))

plotlev1$occ <- plotlev1$ab_cat

plotlev1 %>% filter(species == "plaere") %>% na.omit() -> plotlev2

for(j in grids) {
  plots <- levels(factor(plotlev2$plot[plotlev2$grid == j]))
  for(i in plots) {
  if(length(plotlev2$ab_cat[plotlev2$treatment %in% c("A","B") & plotlev2$plot == i & plotlev2$grid == j]) == 2) {
  if(plotlev2$ab_cat[plotlev2$treatment == "A" & plotlev2$plot == i & plotlev2$grid == j] <= plotlev2$ab_cat[plotlev2$treatment == "B" & plotlev2$plot == i & plotlev2$grid == j]){
    
    plotlev2$occ[plotlev2$treatment == "A" & plotlev2$plot == i & plotlev2$grid == j] <- plotlev2$occ[plotlev2$treatment == "B" & plotlev2$plot == i & plotlev2$grid == j]
  }
    if(plotlev2$ab_cat[plotlev2$treatment == "B" & plotlev2$plot == i & plotlev2$grid == j] <= plotlev2$ab_cat[plotlev2$treatment == "A" & plotlev2$plot == i & plotlev2$grid == j]){
      
      plotlev2$occ[plotlev2$treatment == "B" & plotlev2$plot == i & plotlev2$grid == j] <- plotlev2$occ[plotlev2$treatment == "A" & plotlev2$plot == i & plotlev2$grid == j] 
      
    }
    }
  }
}

# run loop again for each species then rbind
plotlev1 %>% filter(species == "miccal") %>% na.omit() -> plotlev3

for(j in grids) {
  plots <- levels(factor(plotlev3$plot[plotlev3$grid == j]))
  for(i in plots) {
    if(length(plotlev3$ab_cat[plotlev3$treatment %in% c("A","B") & plotlev3$plot == i & plotlev3$grid == j]) == 2) {
      if(plotlev3$ab_cat[plotlev3$treatment == "A" & plotlev3$plot == i & plotlev3$grid == j] <= plotlev3$ab_cat[plotlev3$treatment == "B" & plotlev3$plot == i & plotlev3$grid == j]){
        
        plotlev3$occ[plotlev3$treatment == "A" & plotlev3$plot == i & plotlev3$grid == j] <- plotlev3$occ[plotlev3$treatment == "B" & plotlev3$plot == i & plotlev3$grid == j]
      }
      if(plotlev3$ab_cat[plotlev3$treatment == "B" & plotlev3$plot == i & plotlev3$grid == j] <= plotlev3$ab_cat[plotlev3$treatment == "A" & plotlev3$plot == i & plotlev3$grid == j]){
        
        plotlev3$occ[plotlev3$treatment == "B" & plotlev3$plot == i & plotlev3$grid == j] <- plotlev3$occ[plotlev3$treatment == "A" & plotlev3$plot == i & plotlev3$grid == j] 
        
      }
    }
  }
}

plotlev1 %>% filter(species == "brohor") %>% na.omit() -> plotlev4

for(j in grids) {
  plots <- levels(factor(plotlev4$plot[plotlev4$grid == j]))
  for(i in plots) {
    if(length(plotlev4$ab_cat[plotlev4$treatment %in% c("A","B") & plotlev4$plot == i & plotlev4$grid == j]) == 2) {
      if(plotlev4$ab_cat[plotlev4$treatment == "A" & plotlev4$plot == i & plotlev4$grid == j] <= plotlev4$ab_cat[plotlev4$treatment == "B" & plotlev4$plot == i & plotlev4$grid == j]){
        
        plotlev4$occ[plotlev4$treatment == "A" & plotlev4$plot == i & plotlev4$grid == j] <- plotlev4$occ[plotlev4$treatment == "B" & plotlev4$plot == i & plotlev4$grid == j]
      }
      if(plotlev4$ab_cat[plotlev4$treatment == "B" & plotlev4$plot == i & plotlev4$grid == j] <= plotlev4$ab_cat[plotlev4$treatment == "A" & plotlev4$plot == i & plotlev4$grid == j]){
        
        plotlev4$occ[plotlev4$treatment == "B" & plotlev4$plot == i & plotlev4$grid == j] <- plotlev4$occ[plotlev4$treatment == "A" & plotlev4$plot == i & plotlev4$grid == j] 
        
      }
    }
  }
}

plotlev1 %>% filter(species == "vulmic") %>% na.omit() -> plotlev5

for(j in grids) {
  plots <- levels(factor(plotlev5$plot[plotlev5$grid == j]))
  for(i in plots) {
    if(length(plotlev5$ab_cat[plotlev5$treatment %in% c("A","B") & plotlev5$plot == i & plotlev5$grid == j]) == 2) {
      if(plotlev5$ab_cat[plotlev5$treatment == "A" & plotlev5$plot == i & plotlev5$grid == j] <= plotlev5$ab_cat[plotlev5$treatment == "B" & plotlev5$plot == i & plotlev5$grid == j]){
        
        plotlev5$occ[plotlev5$treatment == "A" & plotlev5$plot == i & plotlev5$grid == j] <- plotlev5$occ[plotlev5$treatment == "B" & plotlev5$plot == i & plotlev5$grid == j]
      }
      if(plotlev5$ab_cat[plotlev5$treatment == "B" & plotlev5$plot == i & plotlev5$grid == j] <= plotlev5$ab_cat[plotlev5$treatment == "A" & plotlev5$plot == i & plotlev5$grid == j]){
        
        plotlev5$occ[plotlev5$treatment == "B" & plotlev5$plot == i & plotlev5$grid == j] <- plotlev5$occ[plotlev5$treatment == "A" & plotlev5$plot == i & plotlev5$grid == j] 
        
      }
    }
  }
}

plotlev <- rbind(plotlev2, plotlev3,plotlev4, plotlev5)
rm(plotlev1,plotlev2,plotlev3,plotlev4, plotlev5)
names(plotlev)[11] <- "occurrence" 

# create y/n
plotlev$occurrence[plotlev$occurrence > 0] <- c("yes")
plotlev$occurrence[plotlev$occurrence == 0] <- c("no")

######################################################
# This changes depending on persistence threshold:
######################################################

plotlev$persistence <- NA
plotlev$persistence[plotlev$seed >= 2] <- c("yes")
plotlev$persistence[plotlev$seed < 2] <- c("no")
plotlev$persistence[is.na(plotlev$seed)] <- c("no")

# categories ME, SS_n, SS_y, DL
plotlev$contingency <- with(plotlev, paste0(occurrence, persistence))

plotlev$contingency[plotlev$contingency == "yesyes"] <- c('SS_y')
plotlev$contingency[plotlev$contingency == "nono"] <- c('SS_n')
plotlev$contingency[plotlev$contingency == "yesno"] <- c('ME')
plotlev$contingency[plotlev$contingency == "noyes"] <- c("DL")

# there is a single date in the ab_cat -> fix
plotlev$ab_cat <- as.numeric(plotlev$ab_cat)
# observation 1262 is culprit -> NA = 2
plotlev$ab_cat[is.na(plotlev$ab_cat)] <- 2
plotlev <- left_join(plotlev, nkey, by = c("tag","species","treatment","grid","site"), 
                     relationship = "many-to-many")

# Important note:
# occurrence differing from ab_cat is normal and fine. ab_cat is raw, occurrance considers that abiotic treatment was cleared of naturally occuring species

# calculate how many situations exist where A does not equal B
plotlev %>%
  select(occurrence, treatment, species, grid, site, plot) %>%
  group_by(species, grid, site, plot)

  #if occurrence in A = occurrence in B say y or else say no
# tt <- plotlev  %>%
#   select(occurrence, treatment, species, grid, site, plot) %>%
#   pivot_wider(names_from = "treatment", values_from = "occurrence") %>%
#   mutate(same = c(B == A)) # all either true or NA due to missing tag

#-------------------------------------------------------.
#### Repeat dataset for grid and site levels ####

#### gridlev ####
# code datasheet for contingencies at grid level
gridlev1 <- plotlev %>%
  dplyr::select(-ab_cat, -persistence, -contingency) %>%
  dplyr::group_by(treatment, species, grid) %>%
  dplyr::mutate(grid_n = n()) %>%
  dplyr::mutate(grid_seed = sum(seed)) %>%
  dplyr::mutate(grid_occur = ifelse('yes' %in% occurrence, 'yes', 'no')) %>% 
  dplyr::ungroup() %>%
  dplyr::select(-tag,-plot,-image,-occurrence,-status,-seed,-names) %>%
  distinct()
gridlev <- as.data.frame(gridlev1)

# check if absence at grid level is valid or NA problem:
gridlev %>%
  filter(grid_occur %in% 'no') 
# 9 entries
# BROHOR absent from all plots (A or B) in grids 27, 2, 16; VULMIC absent from all plots in grid 27.

# allow if occurrence is yes in A then yes in B, and if yes in B then yes in A. Necessary bc NA in grid 5 B for BROHOR translated to absence even though in grid 5 A BROHOR was present.
gridlev$grid_occur[gridlev$grid == 5 & gridlev$treatment %in% "B"] <- 'yes'

names(gridlev)[7] <- 'occurrence'

# Persistence threshold is two seeds per individual, where seeds out >= seeds in. Use threshold that takes into account missing data

gridlev$persistence <- NA
gridlev$persistence[gridlev$grid_seed >= gridlev$grid_n*2] <- c("yes") # threshold at grid_n*2 bc 2 seeds per tube
gridlev$persistence[gridlev$grid_seed < gridlev$grid_n*2] <- c("no")

gridlev$contingency <- with(gridlev, paste0(occurrence, persistence))
gridlev$contingency[gridlev$contingency == "yesyes"] <- c('SS_y')
gridlev$contingency[gridlev$contingency == "nono"] <- c('SS_n')
gridlev$contingency[gridlev$contingency == "yesno"] <- c('ME')
gridlev$contingency[gridlev$contingency == "noyes"] <- c("DL")

rm(gridlev1)

#### sitelev ####
# contingency yes/no for occurrence and persistence at site level
sitelev1 <- plotlev %>%
  dplyr::select(-ab_cat, -persistence, -contingency)  %>%
  dplyr::group_by(site, treatment, species) %>%
  dplyr::mutate(site_n = n()) %>%
  dplyr::mutate(site_seed = sum(seed)) %>%
  dplyr::mutate(site_occur = ifelse('yes' %in% occurrence, 'yes', 'no')) %>% 
  dplyr::ungroup() %>%
  dplyr::select(-tag,-plot,-grid,-image,-occurrence,-status,-seed,-names) %>%
  distinct()
sitelev <- as.data.frame(sitelev1)

# check site level occurrence is 100%
sitelev %>%
  filter(site_occur %in% 'no') # empty, good.

names(sitelev)[6] <- "occurrence"

sitelev$persistence <- NA
sitelev$persistence[sitelev$site_seed >= sitelev$site_n*2] <- c("yes")
sitelev$persistence[sitelev$site_seed < sitelev$site_n*2] <- c("no")

sitelev$contingency <- with(sitelev, paste0(occurrence, persistence))
sitelev$contingency[sitelev$contingency == "yesyes"] <- c('SS_y')
sitelev$contingency[sitelev$contingency == "nono"] <- c('SS_n')
sitelev$contingency[sitelev$contingency == "yesno"] <- c('ME')
sitelev$contingency[sitelev$contingency == "noyes"] <- c("DL")

rm(sitelev1)

#---------------------------------------------------------------------------------.
#### Data for creating SAR and PAR matrices and fitting models: SAR dt_o dt_p ####

# Note, in 'plotlev' some ab-cat = 0, occurrence = yes because A was cleared so paired B plot is used to assess occurrence
dt2 <- plotlev %>%
  dplyr::select(tag, species, treatment, grid, site, seed, ab_cat, occurrence, persistence)
dt1 <- left_join(dt2, nkey, by = c("tag", "species", "treatment", "grid", "site"), 
                 relationship = "many-to-many")
dt_p <- dt1 %>%
  dplyr::select(-occurrence, -ab_cat) #%>% dplyr::filter(!persistence%in% 'no')
dt_o <- dt1 %>%
  dplyr::select(-persistence, -seed) #%>% dplyr::filter(!occurrence%in%'no')

dt_p <- left_join(dt_p, dat1, by = "names", 
                  relationship = "many-to-many") # merge with locations data
dt_p <- dplyr::select(dt_p, -tag,-grid,-site)
dt_p$area <- round(dt_p$area, 0)
dt_p <- as.data.frame(dt_p)
dt_p <- distinct(dt_p)

dt_o <- left_join(dt_o, dat1, by = "names", 
                  relationship = "many-to-many") # merge with locations data
dt_o <- dplyr::select(dt_o,-grid,-site)
dt_o$area <- round(dt_o$area, 0)
dt_o <- as.data.frame(dt_o)
dt_o <- distinct(dt_o)

#---------------------------.
#### Bootstrapping data ####
# P: OLD ####
# For contingencies:

## Local
dl <- dplyr::select(plotlev, species, grid, treatment, contingency)
# Abiotic
dla <- dplyr::filter(dl, treatment%in%"A")
dla <- dplyr:: select(dla, -treatment)
# Biotic
# specify treatment
dlb <- dplyr::filter(dl, treatment%in%"B")
dlb <- dplyr:: select(dlb, -treatment)

## Grid
dg <- dplyr::select(gridlev, grid, species, treatment, contingency)
# Abiotic
dga <- dplyr::filter(dg, treatment%in%"A")
dga <- dplyr:: select(dga, -treatment)
# Biotic
dgb <- dplyr::filter(dg, treatment%in%"B")
dgb <- dplyr:: select(dgb, -treatment)

## Site
ds <- dplyr::select(sitelev, site, species, treatment, contingency)
# Abiotic
dsa <- dplyr::filter(ds, treatment%in%"A")
dsa <- dplyr:: select(dsa, -treatment)
dsa <- distinct(dsa)
# Biotic
dsb <- dplyr::filter(ds, treatment%in%"B")
dsb <- dplyr:: select(dsb, -treatment)

# For persistence & occurrence:

## Local
bl <- dplyr::select(plotlev, species, grid, treatment, occurrence, persistence)
bl <- pivot_longer(bl, cols = c(4,5), names_to = 'type',values_to = 'value')
bl$value[bl$value== 'yes'] <- 1
bl$value[bl$value== 'no'] <- 0
# Abiotic
bla <- bl %>% dplyr::filter(treatment%in%"A")%>%
  dplyr:: select(-treatment)
# Biotic
blb <- bl %>% dplyr::filter(treatment%in%"B")%>%
  dplyr:: select(-treatment)

## Grid
bg <- dplyr::select(gridlev, species, grid, treatment, occurrence, persistence)
bg <- pivot_longer(bg, cols = c(4,5), names_to = 'type',values_to = 'value')
bg$value[bg$value== 'yes'] <- 1
bg$value[bg$value== 'no'] <- 0
# Abiotic
bga <- bg %>% dplyr::filter(treatment%in%"A")%>%
  dplyr:: select(-treatment)
# Biotic
bgb <- bg %>% dplyr::filter(treatment%in%"B")%>%
  dplyr:: select(-treatment)

## Site
bs <- dplyr::select(sitelev, species, site, treatment, occurrence, persistence)
bs <- pivot_longer(bs, cols = c(4,5), names_to = 'type',values_to = 'value')
bs$value[bs$value== 'yes'] <- 1
bs$value[bs$value== 'no'] <- 0
# Abiotic
bsa <- bs %>% dplyr::filter(treatment%in%"A")%>%
  dplyr:: select(-treatment)
# Biotic
bsb <- bs %>% dplyr::filter(treatment%in%"B")%>%
  dplyr:: select(-treatment)


# Making sure variables are correct type
plotlev$persistence <- as.factor(plotlev$persistence)
plotlev$species <- as.factor(plotlev$species)
plotlev$treatment <- as.factor(plotlev$treatment)
plotlev$grid <- as.factor(plotlev$grid)
plotlev$site <- as.factor(plotlev$site)
plotlev$plot <- as.factor(plotlev$plot)
plotlev$occurrence <- as.factor(plotlev$occurrence)
plotlev$contingency <- as.factor(plotlev$contingency)
plotlev$names <- as.factor(plotlev$names)

gridlev$persistence <- as.factor(gridlev$persistence)
gridlev$species <- as.factor(gridlev$species)
gridlev$treatment <- as.factor(gridlev$treatment)
gridlev$grid <- as.factor(gridlev$grid)
gridlev$site <- as.factor(gridlev$site)
gridlev$occurrence <- as.factor(gridlev$occurrence)
gridlev$contingency <- as.factor(gridlev$contingency)

sitelev$persistence <- as.factor(sitelev$persistence)
sitelev$species <- as.factor(sitelev$species)
sitelev$treatment <- as.factor(sitelev$treatment)
sitelev$site <- as.factor(sitelev$site)
sitelev$occurrence <- as.factor(sitelev$occurrence)
sitelev$contingency <- as.factor(sitelev$contingency)


#######################################
# How many plants produces > 20 seeds?
# note must divide by 2 for the 2 seeds planted

hist(plotlev$seed)
plotlev %>%
 # ungroup()
  dplyr::mutate(total_ind = n()) %>% # total number of observations
  dplyr::filter(seed >= 40) %>% # over 20 seeds per captia = 20*2seed
  dplyr::mutate(no = n()) %>%
  dplyr::mutate(perc = no/total_ind) %>%# number of observations with seed > 20 divided by total
  dplyr::select(no, total_ind, perc) %>%
  distinct()
  
# What is the average seed production per ind?
plotlev %>%
  dplyr::mutate(avg_seed = mean(seed)/2) %>%
  dplyr::select(avg_seed) %>%
  distinct()

# checks that all loaded:
# for posting reproducible code start with these files
#View(plotlev)
#View(gridlev) 
#View(sitelev)
#View(df.scales)
