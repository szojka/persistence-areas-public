#########################################################################################################################
# CODE ASSOCIATED WITH MANUSCRIPT: "Species occupancy is a poor predictor of habitat suitability across spatial scales"

## AUTHORS: Megan Szojka, Lauren G. Shoemaker, and Rachel M. Germain
#########################################################################################################################

# OVERVIEW:

This project addressed how persistence and occurrence of four species change across scales in Californian serpentine grasslands. We planted two seeds of each of our focal species in two plots within a block. On plot was cleared of biotic material (without neighbors treatment) and one plot was left intact (with neighbors treatment). There were 25 blocks in a 5 block by 5 block grid, spaced one meter apart. There were three grids per site, and six sites spanning 45 km2 across the reserve. We assessed how the proportion of persisting species within each plot (Fig. 1a,c), and how this changed with and without neighbors across scales (Fig. 2). Our main goal was to combine persistence and occurrence data to assess when these variables were aligned or misaligned (Fig. 1b,d). We compared these (mis)aligned outcomes with and without neighbors (Fig. 3). Finally, building off of species-area relationships, we parameterized novel persistence area relationships: realized-PAR composed of when species occur and persist, and the potential-PAR composed of when species can persist (Fig. 4). There is a script for each figures, and a script for each analysis. As such, the below scripts are names based on whether they are the statistics and data wrangling scripts (Stats -) associated with a figure, or are a figure (Figure - ).

# HOW TO RUN PROJECT:

All scripts depend on 'Source - MAIN fitness.R'. Each figure script sources all necessary source and stats scripts within it. If one wants to run the source or a stats script independently of the figure scripts, there are hashed out packages at the top that must be loaded for that specific script. Specific descriptions of each script's function are found below. Figures are saved to the 'Figures' folder in the main repository, Tables are saved to the 'Tables' folder in the main repository, and all data used in the scripts are housed within the 'Data' folder in the main repository.

In regards to analyses described in the "Test of persistence using Bayesian fits from a Beverton-Holt model" section of the supplemental materials, we have a separate README.md that can be found within the Folder "/Supp Bayesian models/".

Similarly, the supplemental materials that testing the robustness of our results to more permissive persistence thresholds can be found in the sub folder "/Persistence thresholds/", and also has it's own README.md.

# DESCRIPTION OF SCRIPTS:

## Source - MAIN fitness.R

Reading in and wrangling the data. This is sourced at the top of each script.

## Figure - natural misalignments.R

Script to make figure 1, representing the proportion of species that are (mis)aligned with neighbors at the plot scale.

## Figure - persist v occur by treatment.R

Script to make Figure 2, representing the proportion of species occurring, or persisting with and without neighbors across scales.

## Figure - misalignments w processes point.R

Script to make Figure 3, representing the the proportion of species in each (mis)alignment category, with and without neighbors, across scales.

## Figure - SAR-PAR across discrete scales.R

Script to make Figure 4, representing the accumulation or averaging scaling methods, for persistence-area relationships with and without neighbors.

## Supp figure - abundance dependence.R

Script for supplementary figure on how per-capita seed production of each species changes or does not change across conspecific abundance.

## Supp figure - species specific misalignments.R

Script for supplementary figure on how (mis)alignments in figure 1 d break down for each focal species individually.

## Supp figure - temp precip

Using average annual temperature and precipitation data from the Western Regional Climate Center 2022 reflective of the annual plant growing season from October to September from 2001 to 2020. Script corresponds to supplementary figure on environmental conditions. 

## Supp figure - reduction in cover

Comparing the percent vegetative cover between paired plots with neighbors, and plots without neighbors. script corresponds to supplementary figure on total cover and cover reduced.

## Supp figure - histogram seed production.R
Transplant seed production for each species, in conditions with and without neighbors.

## Supp figure - rank does not predict (mis)alignments.R
Compares how well species by rank abundance predict the proportion of (mis)alignments for each species with neighbors.

## Stats - natural persistence v occurrence.R

Assessing the proportion of species that persistence with neighbors or occurr at the plot level. Models are fit to find the estimated proportion and 95% CIs for figure 1c.  

## Stats -  natural misalignments (multinomial).R

Assessing the proportion of species that fall into each (mis)alignment category at the plot level. Multinomial models are fit to find estimated proportions for each category and 95% CIs in figure 1d. 

## Stats - persistence & occurrence by treatment across discrete scales.R

Assessing the proportion of species that persist with and without neigbors, or occur at plot scales, grid scales, and site scales. Models were fit separately for occurrence data, persistence with neighbors, and persistence without neighbors, and were fit separately at plot, grid, and site scales for a total of 9 models. Estimated proportions and 95% CI from the models were used in figure 2.

## Stats - misalignments by treatment at blocks.R

Assessing the proportion of species that are (mis)aligned with and without neighbors at the plot scale (n = 450). Models were fit separately for each (mis)alignment, where proportion of species was predicted by neighbor treatment. Estimates and 95% CI shown in figure 3 at the plot scale were from these models. 

## Stats - misalignments by treatment at grids.R

Assessing the proportion of species that are (mis)aligned with and without neighbors at the grid scale (n = 18). Models were fit separately for each (mis)alignment, where proportion of species was predicted by neighbor treatment. Estimates and 95% CI shown in figure 3 at the grid scale were from these models. 

## Stats - misalignments by treatment at sites.R

Assessing the proportion of species that are (mis)aligned with and without neighbors at the site scale (n = 6). Models were fit separately for each (mis)alignment, where proportion of species was predicted by neighbor treatment. Estimates and 95% CI shown in figure 3 at the site scale were from these models. 

## Stats - SAR-PAR categories across discrete scales - ii. accumulated.R

This script uses the accumulating method of allowing suitable habitat to scale to assess the proportion of species that could persist (potetial-PAR), do persist and occur (realized-PAR), and occur (diversity-SAR) across scales and with and without neighbors. Three models were fit: the proportion of species predicted by each PAR or SAR consisting of five levels (diversity-SAR, potential-PAR with neighbors, potential-PAR without neighbors, realized-PAR with neighbors, realized-PAR without neighbors) at 1. plot scales, 2. grid scales, 3. site scales. Model estimates and 95% CIs are used in figure 4e.

## Stats - SAR-PAR categories across discrete scales - i. averaged.R

This script uses the averaged method of allowing suitable habitat to scale to assess the proportion of species that could persist (potetial-PAR), do persist and occur (realized-PAR), and occur (diversity-SAR) across scales and with and without neighbors. Three models were fit: the proportion of species predicted by each PAR or SAR consisting of five levels (diversity-SAR, potential-PAR with neighbors, potential-PAR without neighbors, realized-PAR with neighbors, realized-PAR without neighbors) at 1. plot scales, 2. grid scales, 3. site scales. Model estimates and 95% CIs are used in figure 4f.

# DESCRIPTION OF DATA FILES:

Here we summarize the contents of each data file and the meaning of each variable name within the data files. First we summarize column names that are found in more than one data file to reduce reducdency. 

tag : a factor with levels as transplant number. This number is not entirely unique to each transplant because of duplicate tags, but each combination or tag within a grid is unique. 

plot	: a factor with levels as plot number. Plot number is not entirely unique to each plot due to naming issues with transplant tags, but each combination of tag, plot and grid is unique.

grid	: a factor with 18 levels, where each grid has a unique number.

site	: a factor with 6 levels named after the area they are found in the reserve: "HUT", "WHIT", "LM", "RES", "TAIL", "KNOX".

treatment	: a factor with 2 levels representing plot neighbor status, either "A" for abiotic or without neighbors, and "B" for biotic or with neighbors.

species : a factor with four levels: Micropus californicus "miccal", Bromus hordeaceus "brohor", Plantago erecta "plaere", or Festuca microstachys "vulmic".

## Transplant_master data_Oct2020.csv

Seed production data for each transplant (n=3600). Used in script "Source - MAIN fitness.R".
	
status  : a factor used to categorize transplant state when we returned to quantify seed production. There are six factor levels which translate to the following: "missing" meaning the pipe each transplant was housed within was disturbed or otherwise lost and no transplant was found, "zero" meaning the transplant was never planted in the first place do to a rock or another obstruction, "no fruit" meaning the transplant germinated by naturally did not produce seed, "no pipe" meaning the pipe each transplant was housed within was disturbed or otherwise lost but the transplant was found, "herbivory" meaning the seed head was obviously eaten, "blank" indicating everything happened as expected.

seed  : the number of seeds that each transplant produced. The threshold for persist was 2 or more seeds for each transplant.

## Abundance & occurrence 2019_final.csv

Observed abundance of each block falling into four categories (0 = 0 individuals, 1 = between 1-10 individuals, 2 = between 11-100 individuals, 3 = above 101 individuals) for each species in each plot. Used in script "Source - MAIN fitness.R".

image	: a factor with each image number representing a block in our experiment (plot without neighbors + paired plot with neighbors).

miccal	: six letter code for Micropus californicus. Values under this column represent what abundance category characterized each plot (0 through 3 described above).

brohor	: six letter code for Bromus hordeaceus. Values under this column represent what abundance category characterized each plot (0 through 3 described above).

vulmic  : six letter code for Festuca (previously Vulpia) microstachys. Values under this column represent what abundance category characterized each plot (0 through 3 described above).

plaere  : six letter code for Plantago erecta. Values under this column represent what abundance category characterized each plot (0 through 3 described above).

## Tag labels.csv

Data that links the transplant unique code to its plot, grid, site, and neighbor treatment.

## UTM_szojkalocs_id.csv

name : a factor with unique levels for each combination of block (1-25), grid (1-18), and site ("HUT", "WHIT", "LM", "RES", "TAIL", or "KNOX") separated by an underscore.

utm.x : numerical location in UTM x coordinates.

utm.y : numerical location in UTM y coordinates.

## climate_data.csv

Year : growing season period (October to September)

Temperature : in Celsius

Precipitation : in cm

# SUPPORTING INFORMATION AND CONTACT:

1. code DOI: 10.5281/zenodo.17186785

2. although not required, the authors would appreciate if you would notify author Megan Szojka mszojka@uoguelph.ca if you are planning on using this data in a project for publication, as other projects are ongoing.

##################################################################


