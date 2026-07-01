library(tidyverse)

set.seed(123)

####Number of random species
S <- 20

####Assigning random mean and CV values to each species
species_df <- tibble(
  species = paste0("sp_", 1:S),

  true_mean =
    runif(S, 50, 150),

  true_CV =
    runif(S, 0.02, 0.10)
)

####Assigning expected values of the parameters of the range-size model

#alpha_R is the intercept of the range-size model
alpha_R_true <- 10
#beta_CV is the slope, and the main estimate of interest
beta_CV_true <- 1.5
#sigma_R is the standard deviation of the residuals in the range-size model
sigma_R_true <- 0.3

####Now we generate the true range-size values for each species based on the formula of the range-size model
species_df <- species_df %>%
  mutate(
    log_range =
      alpha_R_true +
      beta_CV_true * log(true_CV) +
      rnorm(
        S,
        0,
        sigma_R_true
      )
  )

####We then assign random sample sizes to each species, which will be used to generate the specimen-level data

species_df$n_specimens <-
  sample(
    c(
      10,
      20,
      30,
      50,
      100
    ),
    S,
    replace = TRUE
  )

####Assigning effects for sex and source

beta_sex_true <- 3

beta_source_true <- c(
  field = 0,
  museum = -2
)

####Now we generate the specimen-level data for each species, based on the true mean and CV values, as well as accounting for the effects of Sex and source

sim_data <- map_dfr(
  1:S,
  function(i){

    n_i <- species_df$n_specimens[i]

    sex_i <- sample(
      c("F","M"),
      n_i,
      replace = TRUE
    )

    source_i <- sample(
      c("field","museum"),
      n_i,
      replace = TRUE
    )

    mean_i <- species_df$true_mean[i]

    cv_i <- species_df$true_CV[i]

    sd_i <- mean_i * cv_i

    mu_i <-
      mean_i +
      ifelse(
        sex_i == "M",
        beta_sex_true,
        0
      ) +
      ifelse(
        source_i == "museum",
        beta_source_true["museum"],
        0
      )

    tibble(

      species =
        species_df$species[i],

      sex =
        sex_i,

      source =
        source_i,

      wing =
        rnorm(
          n_i,
          mu_i,
          sd_i
        )
    )

  }
)

####Add the previously generated log_range values to the specimen-level data

sim_data <- left_join(
  sim_data,
  species_df %>%
    select(
      species,
      log_range
    ),
  by = "species"
)

####Make everything a factor or integer for Stan

sim_data$species <- factor(sim_data$species)
sim_data$sex <- factor(sim_data$sex)
sim_data$source <- factor(sim_data$source)

species_id <- as.integer(sim_data$species)
sex_id <- as.integer(sim_data$sex)
source_id <- as.integer(sim_data$source)

#Indexing species_df to match the order of species in sim_data
species_levels <- levels(sim_data$species)

species_df <- species_df %>%
  arrange(
    match(
      species,
      species_levels
    )
  )

#Checking if the species numbers match
all(
  species_df$species ==
  species_levels
)


####Now to put all the data together into a list for Stan

stan_data_sim <- list(

  N = nrow(sim_data),

  S = nlevels(sim_data$species),

  wing = sim_data$wing,

  log_range =
    species_df$log_range,

  species = species_id,

  n_sex =
    nlevels(sim_data$sex),

  sex =
    sex_id,

  n_source =
    nlevels(sim_data$source),

  source =
    source_id,

  #This is a basic identity matrix to act as a placeholder for the phylo-covariance matrix
  L_A =
    diag(nlevels(sim_data$species))
)

####Compiling the model
library(cmdstanr)

mod <- cmdstan_model(
  "latent_variance_rangesize.stan"
)

####Fitting the model

fit <- mod$sample(
  data = stan_data_sim,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.95,
  max_treedepth = 15,
  refresh = 100
)

####Post-hoc checks and summaries

posterior_cv <- fit$summary(
  "CV_species"
)

posterior_cv$species <- levels(sim_data$species)

cv_compare <- left_join(
  posterior_cv %>%
    dplyr::select(
      species,
      estimated_CV = mean
    ),
  species_df %>%
    dplyr::select(
      species,
      true_CV
    ),
  by = "species"
)

#To test true CV vs estimated CV correlation
cor.test(
  cv_compare$estimated_CV,
  cv_compare$true_CV
)

#To check the beta_CV estimate and see whether it captures the true value of 1.5
fit$summary("beta_CV")

#TO see whether the sample size relates to the uncertainty in the CV estimates as expected
posterior_cv$width <-
  posterior_cv$q95 -
  posterior_cv$q5

plot(
  species_df$n_specimens,
  posterior_cv$width
)
