
// To allow for nested random effects easily, I've switched to an additive version of accounting for random intercepts (see Statistical Rethinking), whereas previously I was using a multiplicative way, where epsilon was a site level scalar.

data {
  int<lower = 1> N;             // Number of plots (NEEDS TO BE LISTED FIRST)
  int Site[N];                  // Site ID for each observation
  int Grid[N];                  // Grid index per observation

  int Fecundity[N];             // Focal species fecundity
  int Abundance[N];             // Intraspecific abundance
  int Cover[N];                 // Heterospecific cover
  int<lower = 1> P;             // Number of sites
  vector[N] Env;                // Environmental covariate
  int<lower=1> G;               // Number of grids
}

parameters {
  real lambda_0;
  real lambda_1;
  real alpha_intra_0;
  real alpha_intra_1;
  real alpha_other_0;
  real alpha_other_1;

  vector[P] site_effect;    // Random site intercept estimate
  real<lower=0> sigma_site;   // Site-level SD
  vector[G] grid_effect;  // Grid
  real<lower=0> sigma_grid;       
  
}

model {
  // create a vectors for parameters
  vector[N] lambda_ei;
  vector[N] alpha_intra_ei;
  vector[N] alpha_other_ei;
  vector[N] F_hat2;

  // Priors
  lambda_0 ~ normal(0, 1);
  lambda_1 ~ normal(0, 1);
  alpha_intra_0 ~ cauchy(0, 1); // alpha intercept to cauchy
  alpha_intra_1 ~ normal(0, 1);
  alpha_other_0 ~ cauchy(0, 1); // alpha intercept to cauchy
  alpha_other_1 ~ normal(0, 1);

  site_effect ~ normal(0, sigma_site);
  grid_effect ~ normal(0, sigma_grid);

  // below necessary either formation of rnd effects
  sigma_site ~ normal(0, 1); // no longer exponential (which is weakly informed), now this is half-normal (normal but greater than zero)
  sigma_grid ~ normal(0, 1);

  // Likelihood loop
  // Model
  for (i in 1:N) {
    lambda_ei[i] = exp(lambda_0 + lambda_1 * Env[i]);
    alpha_intra_ei[i] = exp(alpha_intra_0 + alpha_intra_1 * Env[i]);
    alpha_other_ei[i] = exp(alpha_other_0 + alpha_other_1 * Env[i]);
    
    // Beverton Holt
    real F_hat = lambda_ei[i] / (1 + alpha_intra_ei[i] * Abundance[i] + alpha_other_ei[i] * Cover[i]); // when F-hat is here it means that its not saved (intermediate parameter)
    
    // Additive random effect on raw fecundity scale
    real total_rnd_effects = exp(site_effect[Site[i]] +
               grid_effect[Grid[i]]); // already exponentiated here
               
    
    F_hat2[i] = F_hat * total_rnd_effects;  
  }

  // Full likelihood
  Fecundity ~ poisson(F_hat2); // natural scale (no log link)
}

generated quantities {
  vector[N] y_hat; // full Bev-Holt
  vector[N] y_hat_no_intra; // Bev-Holt no intra
  vector[N] y_hat_baseline; // no intra no inter

  vector[N] lambda_ei_out;
  vector[N] alpha_intra_ei_out;
  vector[N] alpha_other_ei_out;

  for (i in 1:N) {
    lambda_ei_out[i] = exp(lambda_0 + lambda_1 * Env[i]);
    alpha_intra_ei_out[i] = exp(alpha_intra_0 + alpha_intra_1 * Env[i]);
    alpha_other_ei_out[i] = exp(alpha_other_0 + alpha_other_1 * Env[i]);

    real F_hat = lambda_ei_out[i] / (1 + alpha_intra_ei_out[i] * Abundance[i] + alpha_other_ei_out[i] * Cover[i]);
    real F_hat_no_intra = lambda_ei_out[i] / (1 + alpha_other_ei_out[i] * Cover[i]);
    real F_hat_baseline = lambda_ei_out[i];
    real total_rnd_effects = exp(site_effect[Site[i]] +
               grid_effect[Grid[i]]); 
               
  // the following predicted fecundities account for random effects
    y_hat[i] = F_hat * total_rnd_effects;  
    y_hat_no_intra[i] = F_hat_no_intra * total_rnd_effects;
    y_hat_baseline[i] = F_hat_baseline * total_rnd_effects;
  }
}

