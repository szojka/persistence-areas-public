#########################################################################################################################
# CODE ASSOCIATED WITH SUPPLEMENTARY MATERIALS of the manuscript: 
# "Species occupancy is a poor predictor of habitat suitability across spatial scales"
## AUTHORS: Megan Szojka, Lauren G. Shoemaker, and Rachel M. Germain
#########################################################################################################################

# OVERVIEW:

To evaluate the robustness of our persistence ($\lambda \geq 1$) threshold calculated using raw empirical fecundity, we fit a dynamical population model for each of our species and used the subsequent parameter estimates to quantify population growth of each species at each spatial scale. We did so using Bayesian methods, implemented in Stan. We then used these model predictions to test for sources of bias in our original analyses by controlling for the presence of conspecifics or neighbors that had colonized the cleared plots.

# HOW TO RUN PROJECT:

All scripts depend on 'Source - MAIN fitness POST STAN.R'. Each figure script sources all necessary source and stats scripts within it. If one wants to run the source or a stats script independently of the figure scripts, there are hashed out packages at the top that must be loaded for that specific script. Specific descriptions of each script's function are found below. Figures are saved to the 'Figures/Bayes/' folder in the main repository, Tables are saved to the 'Tables/Bayes/' folder in the main repository.

The order to recreate our supplementary analysis is as follows:

1) Fit Bayesian population models to our focal species' seed production data using the following four scripts found within the folder /Supp Baysian models/Scripts - run stan models/ (i) Bayes - BH_BROHOR.R, (ii) Bayes - BH_MICCAL.R, (iii) Bayes - BH_MICCAL.R, (iv) Bayes - BH_VULMIC.R.

2) Extract the posteriors of those models, and create versions of predictions with 95% Credible Intervals following the supplementary text (version 1 and 2, or version 3). Specifically, under the folder '/Supp Baysian models/Scripts - run stan models/', script 'Bayes - FullBevHolt predictions from stan 95% CI.R' generates .RData objects following version 3 of the model at the three spatial scales: fullBH_pred_data_saved_95CI.rdata, fullBH_grid_pred_95CI.rdata, fullBH_site_pred_95CI.rdata. Likewise, within the folder '/Supp Baysian models/Scripts - run stan models/', the script 'Bayes - predictions from stan 95% CI.R' generates .RData objects following combined versions 1 and 2 of the model at the three spatial scales: pred_data_saved_95CI.rdata, grid_pred_95CI.rdata, site_pred_95CI.rdata.

3) Use these predictions to recreate figures and analyses in the main text, using brms instead of glmmTMB. These scripts mirror those found in the main analyses (see README.md in main repository), and are found within the folder '/Supp Bayesian models/Scripts - generate stan figures/'

# FOLDER DESCRIPTION:

There are four folders we use to organize the workflow of these scripts. The below folders contain the same sub-structure. Any subfolder named '95CI/' uses predictions from model versions 1 and 2 that control for unwanted neighbors (see text for details), whereas the subfolders named 'FullBevHolt' uses predictions from model version 3 (see text for details).

## Supp Baysian models/Stan models/

Contains .stan scripts of the two population models we use to fit to our empricial data ('BH - nested.stan', 'BH - ZIH model nested.stan'), and one multinomial model.

## Supp Baysian models/Scripts - run stan models/

Contains R scripts that prepare data to run the stan or brms models, and then saves these model objects as .rdata files.

From the supplemental materials, model version 1, 2 and 3 are built using scripts in the main folder. The subfolder named '/Scripts - run stan models/Run models for Bayes figures/' house R scripts used to run brms models, in place of the original GLMM models, and use these outputs to recreate our Figures. 

## Supp Baysian models/Stan rdata/

Contains the .rdata files saved by the scripts in the 'Scripts - run stan models' folder.

## Supp Baysian models/Stan - generate stan figures/

These R scripts load .rdata brms models to recreate figures. There is an additional script named 'Supp figures - Bayesian posteriors.R' which is used to make supplementary figures S8, S9, S10, and S15. 

##################################################################


