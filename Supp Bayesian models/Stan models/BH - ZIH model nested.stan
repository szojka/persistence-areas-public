// Beverton-Holt with additive random effects and zero-inflation (non-centered version)

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
  real<lower=0, upper=1> psi;         // probability of structural zero

  vector[P] site_effect;    // Random site intercept estimate
  real<lower=0> sigma_site;   // Site-level SD
  vector[G] grid_effect;  // Grid
  real<lower=0> sigma_grid;       

  real lambda_0;
  real lambda_1;
  real alpha_intra_0;
  real alpha_intra_1;
  real alpha_other_0;
  real alpha_other_1;
}

model {
  // create a vector of predictions
  vector[N] lambda_ei;
  vector[N] alpha_intra_ei;
  vector[N] alpha_other_ei;
  vector[N] F_hat2;
 
  // set priors
  psi ~ beta(1, 1);  // flat prior on zero inflation parameter

  site_effect ~ normal(0, sigma_site);
  grid_effect ~ normal(0, sigma_grid);

  // below necessary either formation of rnd effects
  sigma_site ~ normal(0, 1); // no longer exponential (which is weakly informed), now this is half-normal (normal but greater than zero)
  sigma_grid ~ normal(0, 1);

  lambda_0 ~ normal(0, 1); // changed from gamma for interpretability
  lambda_1 ~ normal(0, 1); 
  alpha_intra_0 ~ cauchy(0, 1); // alpha intercept to cauchy
  alpha_intra_1 ~ normal(0, 1);
  alpha_other_0 ~ cauchy(0, 1); // alpha intercept to cauchy
  alpha_other_1 ~ normal(0, 1);

  // Likelihood loop
  // implement the biological model
  for (i in 1:N) {
    // environmental dependencies for lambda and alphas (equation for a line)
    lambda_ei[i] = exp(lambda_0 + lambda_1 * Env[i]);
    alpha_intra_ei[i] = exp(alpha_intra_0 + alpha_intra_1 * Env[i]);
    alpha_other_ei[i] = exp(alpha_other_0 + alpha_other_1 * Env[i]);

    // Beverton Holt
    real F_hat = lambda_ei[i] / (1 + alpha_intra_ei[i] * Abundance[i] + alpha_other_ei[i] * Cover[i]);

    // Additive random effect on raw fecundity scale
    real total_rnd_effects = exp(site_effect[Site[i]] +
                                 grid_effect[Grid[i]]); // already exponentiated here

    // these mu's follow the same structure of F_hats 
    real mu_rnd = F_hat * total_rnd_effects; //this part accounts for random effect of block, site, grid

    // add for zero-inflation Poisson likelihood component:
    if (Fecundity[i] == 0) {
      target += log_sum_exp(
        bernoulli_lpmf(1 | psi),                              // structural zero (zero from zero-inflation process)
        bernoulli_lpmf(0 | psi) + poisson_lpmf(0 | mu_rnd)    // sampling zero (true zero from poisson process)
      );
    } else {
      target += bernoulli_lpmf(0 | psi) + poisson_lpmf(Fecundity[i] | mu_rnd); // if fecundity > 0 it must come from poisson process
    }
  }
}

// save predictions
generated quantities {
  vector[N] y_hat; // ZIP adjusted full Bev-Holt
  vector[N] y_hat_no_intra; // ZIP adjusted Bev-Holt no intra
  vector[N] y_hat_baseline; // ZIP adjusted lambda (no intra no inter)

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

    // Multiply by random site effect and by (1 - psi) for ZIP expected value
    // Note that (1 - psi) is the log probability of not being a zero-inflated zero (aka structural zero)
    y_hat[i] = (1 - psi) * F_hat * total_rnd_effects;
    y_hat_no_intra[i] = (1 - psi) * F_hat_no_intra * total_rnd_effects;
    y_hat_baseline[i] = (1 - psi) * F_hat_baseline * total_rnd_effects;
  }
}
