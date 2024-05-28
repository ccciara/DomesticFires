data {
  int<lower=0> N;                                                        // number of spatial units or neighbourhoods in Campina Grande
  int<lower=0> N_edges;                                                  // number of edges connecting adjacent areas using Queens contiguity
  array[N_edges] int<lower=1, upper=N> node1;                            // list of index areas showing which spatial units are neighbours
  array[N_edges] int<lower=1, upper=N> node2;                            // list of neighbouring areas showing the connection to index spatial unit
  array[N] int<lower=0> Y;                                               // dependent variable i.e., number of households infested
  vector<lower=1>[N] offset;                                             // total number of houses in neoighbourhoods are used as an offset variable
  int<lower=3> K;                                                        // number of independent variables i.e., K = 4
  matrix[N, K] X;                                                        // independent variables in matrix form i.e., temp, prec, ndvi, and urban
}

transformed data {
    vector[N] log_offset = log(offset);
}

parameters {
  real alpha;                                                            // intercept
  vector[K] beta;                                                        // covariates
  real<lower=0> sigma;                                                   // overall standard deviation
  real<lower=0, upper=1> rho;                                            // proportion unstructured vs. spatially structured variance
  vector[N] theta;                                                       // unstructure random effects (heterogeneous)
  vector[N] phi;                                                         // spatial random effects
}

transformed parameters {
  vector[N] combined;                                                    // combined random effect i.e., unstructure and structured

  combined = sqrt(1 - rho) * theta + sqrt(rho) * phi;                    // formulation for the combined random effect i.e., unstructure and structured
}

model {
  Y ~ poisson_log(log_offset + alpha + X * beta + combined * sigma);     // likelihood function: multivariable Poisson ICAR regression model
  alpha ~ normal(0.0, 1.0);                                              // prior for alpha: weakly informative
  beta ~ normal(0.0, 1.0);                                               // prior for betas: weakly informative
  theta ~ normal(0.0, 1.0);                                              // prior for theta: weakly informative
  sigma ~ normal(0.0, 1.0);                                              // prior for sigma: weakly informative
  rho ~ beta(0.5, 0.5);                                                  // prior for rho: weakly informative

  target += -0.5 * dot_self(phi[node1] - phi[node2]);                    // prior for phi: using a conditional autoregressive distribution
  sum(phi) ~ normal(0, 0.001 * N);                                       // soft sum-to-zero constraint on phi, and its the equivalent to mean(phi) ~ normal(0,0.001)
}

generated quantities {
  vector[N] eta = alpha + X * beta + combined * sigma;                  // compute eta and exponentiate into mu                   
  vector[N] rr_mu = exp(eta);                                           // output the neighbourhood-specific relative risks in mu
  vector[K] rr_beta = exp(beta);                                        // output the risk ratios for each coefficient for the independent variables
  real rr_alpha = exp(alpha);                                           // output the risk ratios for the intercept
}
