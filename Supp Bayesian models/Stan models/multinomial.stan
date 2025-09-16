
// This stan script is to code a multinomial model, as BRMS was failing me in this regard (proportions need to = 1)


data {
  int<lower=1> N;             // number of observations (rows)
  int<lower=1> J;             // number of possible categories (4 (mis)alignments)
  int y_mat[N, J];                // counts of species per category (each row sums to total species sampled)
  
  int Site[N];                  // Site ID for each observation
  int Grid[N];                  // Grid index per observation
  //int Block[N];                 // Block index per observation
  
  int<lower = 1> S;             // Number of sites
  int<lower=1> G;               // Number of grids
  //int<lower=1> B;               // Number of blocks
}

parameters {
  vector[J-1] alpha;          // estimates for each categories J-1 = 3
  
  vector[S] site_effect;    // Random site intercept estimate
  real<lower=0> sigma_site;   // Site-level SD
  vector[G] grid_effect;  // Grid
  real<lower=0> sigma_grid;       
}

model {
  // priors
  alpha ~ normal(0, 2);
  // priors for random effects
  site_effect ~ normal(0, sigma_site);
  grid_effect ~ normal(0, sigma_grid);
   // priors for sigmas of random effects
  sigma_site ~ normal(0, 1); 
  sigma_grid ~ normal(0, 1);

  // likelihood
  for (i in 1:N) {
    
    // estimtes
    vector[J] eta;
    eta[1] = 0; // reference category
    
    // Additive random effects
    real total_rnd_effects = site_effect[Site[i]] +
               grid_effect[Grid[i]]; 

    // find eta estimates for a row i
    for (j in 2:J){
      eta[j] = alpha[j-1] + total_rnd_effects; 
    }

    y_mat[i] ~ multinomial(softmax(eta)); // caculates estimated proportion for each row and Category with random effects
    
  }
  
}

generated quantities {
  matrix[N, J] probs;    // predicted probabilities for each observation

  for (i in 1:N) {
    vector[J] eta;
    real total_rnd_effects = site_effect[Site[i]] +
                grid_effect[Grid[i]]; // block_effect[Block[i]]

    eta[1] = 0;
    for (j in 2:J) {
      eta[j] = alpha[j-1] + total_rnd_effects;
    }
    probs[i] = softmax(eta)';  // ' is column vector to row vector
  }
}
