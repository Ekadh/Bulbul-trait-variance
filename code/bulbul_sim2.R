library(tidyverse)
library(ape)
library(MASS)

set.seed(123)

#=================================================
# NUMBER OF SPECIES
#=================================================

S <- 100

#=================================================
# SIMULATED PHYLOGENY
#=================================================

phy <- rtree(S)

A <- ape::vcv(
  phy,
  corr = TRUE
)

L_A <- t(
  chol(A)
)

#=================================================
# PHYLOGENETIC CV VALUES
#=================================================

phy_effect <- MASS::mvrnorm(
  n = 1,
  mu = rep(0,S),
  Sigma = A
)

true_CV <-
  exp(
    log(0.05) +
    0.4 * phy_effect
  )

#=================================================
# SPECIES MEANS
#=================================================

true_mean <- runif(
  S,
  50,
  150
)

#=================================================
# TRUE RANGE-SIZE RELATIONSHIP
#=================================================

alpha_R_true <- 10

beta_CV_true <- 1.5

sigma_R_true <- 0.3

species_df <- tibble(

  species =
    paste0(
      "sp_",
      1:S
    ),

  true_mean =
    true_mean,

  true_CV =
    true_CV,

  log_range =
    alpha_R_true +
    beta_CV_true *
    log(true_CV) +
    rnorm(
      S,
      0,
      sigma_R_true
    )
)

#=================================================
# SAMPLE SIZES
#=================================================

species_df$n_specimens <-
  sample(
    c(
      20,
      50,
      100,
      200
    ),
    S,
    replace = TRUE
  )

#--------------------------------------------------
# Sex and source effects
#--------------------------------------------------

beta_sex_true <- 3

beta_source_true <- c(
  field = 0,
  museum = -2
)

#--------------------------------------------------
# Generate specimen data
#--------------------------------------------------

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

#--------------------------------------------------
# Add range size
#--------------------------------------------------

sim_data <- left_join(
  sim_data,
  species_df %>%
    select(
      species,
      log_range
    ),
  by = "species"
)

#===========================================================
# FACTORS
#===========================================================

sim_data$species <- factor(sim_data$species)
sim_data$sex <- factor(sim_data$sex)
sim_data$source <- factor(sim_data$source)

#===========================================================
# INDICES
#===========================================================

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

all(
  species_df$species ==
  species_levels
)


#===========================================================
# STAN DATA
#===========================================================

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

  # Identity phylogeny
  L_A =
    diag(nlevels(sim_data$species))
)

#===========================================================
# COMPILE MODEL
#===========================================================
library(cmdstanr)

mod <- cmdstan_model(
  "latent_variance_rangesize.stan"
)

#===========================================================
# FIT MODEL
#===========================================================

fit <- mod$sample(

  data = stan_data_sim,

  chains = 4,
  parallel_chains = 4,

  iter_warmup = 1000,
  iter_sampling = 2000,

  adapt_delta = 0.95,
  adapt_tree_depth = 15,
  refresh = 100
)

#===========================================================
# POST-HOC CHECKS
#===========================================================

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

cor.test(
  cv_compare$estimated_CV,
  cv_compare$true_CV
)

fit$summary("beta_CV")

posterior_cv$width <-
  posterior_cv$q95 -
  posterior_cv$q5

plot(
  species_df$n_specimens,
  posterior_cv$width
)
