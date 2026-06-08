data {

  int<lower=1> N;
  int<lower=1> S;

  // Specimen-level response

  vector[N] wing;

  // Species-level response

  vector[S] log_range;

  // Species index

  array[N] int<lower=1, upper=S> species;

  // Sex

  int<lower=1> n_sex;
  array[N] int<lower=1, upper=n_sex> sex;

  // Source

  int<lower=1> n_source;
  array[N] int<lower=1, upper=n_source> source;

  // Phylogenetic Cholesky matrix

  matrix[S,S] L_A;
}

parameters {

  // Global wing mean

  real alpha;

  vector[n_sex] beta_sex;

  vector[n_source] beta_source;

  // Species means
 
  vector[S] species_mean_raw;

  real mu_mean;

  real<lower=0> sigma_mean;

  // Latent phylogenetic log-CV

  real mu_logCV;

  real<lower=0> sigma_logCV;

  vector[S] z_logCV;

  // Range-size model

  real alpha_R;

  real beta_CV;

  real<lower=0> sigma_R;
}

transformed parameters {

  // Species mean wing lengths

  vector[S] log_species_mean;

  vector[S] species_mean;

  // Species log-CVs

  vector[S] logCV;

  // Species SDs

  vector[S] sigma_species;

  // Mean model

  log_species_mean =
    mu_mean +
    sigma_mean * species_mean_raw;

  species_mean =
    exp(log_species_mean);

  // Phylogenetic log-CV model

  logCV =
    mu_logCV +
    sigma_logCV * (L_A * z_logCV);

  // Convert CV -> SD

  for(i in 1:S){

    sigma_species[i] =
      species_mean[i] *
      exp(logCV[i]);

  }
}

model {

  // Priors

  alpha ~ normal(0, 5);

  beta_sex ~ normal(0, 2);

  beta_source ~ normal(0, 2);

  // Species mean priors

  mu_mean ~ normal(log(100), 1);

  sigma_mean ~ exponential(1);

  species_mean_raw ~ normal(0,1);

  // log-CV priors

  mu_logCV ~ normal(log(0.05), 1);

  sigma_logCV ~ exponential(1);

  z_logCV ~ normal(0,1);

  // Range-size priors

  alpha_R ~ normal(0,5);

  beta_CV ~ normal(0,1);

  sigma_R ~ exponential(1);

  // Specimen-level model

  for(n in 1:N){

    real mu_n;

    mu_n =
      species_mean[species[n]] +
      beta_sex[sex[n]] +
      beta_source[source[n]];

    wing[n] ~ normal(
      mu_n,
      sigma_species[species[n]]
    );
  }

  // Species-level range-size model

  log_range ~ normal(
    alpha_R +
    beta_CV * logCV,
    sigma_R
  );
}

generated quantities {

  // Back-transformed CVs

  vector[S] CV_species;

  // Predicted range size

  vector[S] predicted_range;

  // Mean CV

  real mean_CV;

  for(i in 1:S){

    CV_species[i] =
      exp(logCV[i]);

    predicted_range[i] =
      alpha_R +
      beta_CV * logCV[i];
  }

  mean_CV =
    mean(CV_species);
}
