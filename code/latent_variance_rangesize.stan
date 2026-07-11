data {

  int<lower=1> N;
  int<lower=1> S;

  vector[N] wing;

  vector[S] log_range;
  vector[S] density;
  vector[S] ESI;

  array[S] int<lower=0,upper=1> has_range;
  array[S] int<lower=0,upper=1> has_density;
  array[S] int<lower=0,upper=1> has_ESI;

  array[N] int<lower=1, upper=S> species;

  int<lower=1> n_source;
  array[N] int<lower=1, upper=n_source> source;

  real<lower=0> trait_mean;

  matrix[S,S] L_A;
}

parameters {

  // Random source effects

  vector[n_source] source_raw;

  real<lower=0> sigma_source;

  // Species means

  vector[S] species_mean_raw;

  real mu_mean;

  real<lower=0> sigma_mean;

  // Latent phylogenetic log-CV

  real mu_logCV;

  real<lower=0> sigma_logCV;

  vector[S] z_logCV;

  // Range-size model

  real alpha_range;

  real beta_CV_range;

  real<lower=0> sigma_range;

  // Density model

  real alpha_density;

  real beta_CV_density;

  real<lower=0> sigma_density;

  // ESI model

  real alpha_ESI;

  real beta_CV_ESI;

  real<lower=0> sigma_ESI;

}

transformed parameters {

  vector[n_source] source_effect;

  vector[S] log_species_mean;

  vector[S] species_mean;

  vector[S] logCV;

  vector[S] sigma_species;

  source_effect =
    sigma_source *
    source_raw;

  log_species_mean =
    mu_mean +
    sigma_mean *
    species_mean_raw;

  species_mean =
    exp(log_species_mean);

  logCV =
    mu_logCV +
    sigma_logCV *
    (L_A * z_logCV);

  for(i in 1:S){

    sigma_species[i] =
      species_mean[i] *
      exp(logCV[i]);

  }

}

model {

  //----------------------------------------------------------
  // Priors
  //----------------------------------------------------------

  source_raw ~ normal(0,1);

  sigma_source ~ exponential(1);

  mu_mean ~ normal(log(trait_mean),1);

  sigma_mean ~ exponential(1);

  species_mean_raw ~ normal(0,1);

  mu_logCV ~ normal(log(0.05),1);

  sigma_logCV ~ exponential(1);

  z_logCV ~ normal(0,1);

  alpha_range ~ normal(0,5);

  beta_CV_range ~ normal(0,1);

  sigma_range ~ exponential(1);

  alpha_density ~ normal(0,5);

  beta_CV_density ~ normal(0,1);

  sigma_density ~ exponential(1);

  alpha_ESI ~ normal(0,5);

  beta_CV_ESI ~ normal(0,1);

  sigma_ESI ~ exponential(1);

  //----------------------------------------------------------
  // Specimen-level model
  //----------------------------------------------------------

  for(n in 1:N){

    real mu_n;

    mu_n =

      species_mean[species[n]] +

      source_effect[source[n]];

    wing[n] ~

      normal(

        mu_n,

        sigma_species[species[n]]

      );

  }

  //----------------------------------------------------------
  // Species-level models
  //----------------------------------------------------------

  for(i in 1:S){

    if(has_range[i] == 1){

      log_range[i] ~

        normal(

          alpha_range +

          beta_CV_range *

          logCV[i],

          sigma_range

        );

    }

    if(has_density[i] == 1){

      density[i] ~

        normal(

          alpha_density +

          beta_CV_density *

          logCV[i],

          sigma_density

        );

    }

    if(has_ESI[i] == 1){

      ESI[i] ~

        normal(

          alpha_ESI +

          beta_CV_ESI *

          logCV[i],

          sigma_ESI

        );

    }

  }

}

generated quantities {

  vector[S] CV_species;

  vector[S] predicted_range;

  vector[S] predicted_density;

  vector[S] predicted_ESI;

  real mean_CV;

  for(i in 1:S){

    CV_species[i] =
      exp(logCV[i]);

    predicted_range[i] =

      alpha_range +

      beta_CV_range *

      logCV[i];

    predicted_density[i] =

      alpha_density +

      beta_CV_density *

      logCV[i];

    predicted_ESI[i] =

      alpha_ESI +

      beta_CV_ESI *

      logCV[i];

  }

  mean_CV =
    mean(CV_species);

}
