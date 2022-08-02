// Schaefer surplus production model
// Abdullah Habibi 2022-06-29

data{
int<lower = 0> n_years; // the number of years
vector[n_years] years; // years
vector[n_years] cpue; // ideally abundance index from observation
vector[n_years] harvest; // vector of harvest not using catch since cates since apparantly catch is a reserved word
}

transformed data {
  vector[n_years] log_CPUE;
  log_CPUE = log(cpue);
}

parameters{
  real<lower = -2, upper = 0> log_r; // intrinsic growth rate
  real log_K; //carrying capacity
  real <lower = 1e-9, upper = 1> q; // catchability
  real<lower = 0> sigma_obs; // observation error
}

transformed parameters{
real r;
real K;
vector[n_years] EstBt; // vector of biomass
vector[n_years] EstCPUE; // vector of estimated cpue
vector[n_years] log_EstCPUE; // vector of log estimated cpue
r = exp(log_r);
K = exp(log_K);
EstBt[1] = K;

  for (i in 2:n_years){ // put the population model in transformed params to speed up computations
    EstBt[i] = fmax(1e-6,(EstBt[i - 1] + r * EstBt[i - 1] * (1 - EstBt[i - 1]/K) - harvest[i - 1]));
  } // close loop

  EstCPUE = EstBt * q; //
  log_EstCPUE = log(EstCPUE); // log transform
}

model{ // put only sampling priors statements in this section to increase efficiency
  log_CPUE ~ normal(log_EstCPUE, sigma_obs);
  log_K ~ normal(8,2);
  sigma_obs ~ cauchy(0,2.5);
}

generated quantities{
vector[n_years] pp_log_CPUE;
vector[n_years] pp_CPUE;
real MSY;
real EMSY;

 for (i in 1:n_years){ // posterior predictions for CPUE
    pp_log_CPUE[i] = normal_rng(log_EstCPUE[i], sigma_obs);
 }

pp_CPUE = exp(pp_log_CPUE);

   // Other ouputs
   MSY = r*K/4;
   EMSY = r/(2*q);
}
